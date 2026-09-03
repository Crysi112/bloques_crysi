/*
 * serial_win.h — Comunicación serie con protocolo HIL binario (Windows).
 *
 * Se compila como parte de block_core.c (incluido inline).
 * Solo se activa si _WIN32 está definido.
 *
 * Protocolo:
 *   PC → MCU:  [0xAA] [0x55] [duty_hi] [duty_lo] [XOR]  (5 bytes)
 *   MCU → PC:  [0x41|0x55] [angle_hi] [angle_lo] [rpm_b0..b3] [XOR]  (8 bytes)
 *   El primer byte es el sync de datos: 0x41 = firmware v2 (estimador
 *   nuevo + duty directo, el actual), 0x55 = firmware v1. Se aceptan ambos.
 */
#ifndef SERIAL_WIN_H
#define SERIAL_WIN_H
#ifdef _WIN32
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#define HIL_SYNC_CMD   0xAA
#define HIL_SYNC_CMD2  0x55
#define HIL_SYNC_DATA_V1 0x55   /* firmware v1 (estimador viejo) */
#define HIL_SYNC_DATA_V2 0x41   /* firmware v2 (duty directo) — el actual */
#define HIL_CMD_LEN    5
#define HIL_DATA_LEN   8
typedef struct {
  HANDLE hSerial;
  int    connected;
  float  angle_deg;
  float  rpm;
  uint16_t raw_angle;
  unsigned char rx[16];   /* FIFO de bytes recibidos (hasta 2 tramas) */
  int    rx_len;          /* bytes pendientes en rx */
} SerialHIL;
static int serial_open(SerialHIL *s, int com_port, int baudrate) {
  char portName[32];
  DCB dcbSerialParams = {0};
  COMMTIMEOUTS timeouts = {0};
  sprintf(portName, "\\\\.\\COM%d", com_port);
  /* FILE_SHARE_READ|WRITE: permite reabrir el puerto desde el mismo
     proceso (el 2º modelo de una prueba abría COMx mientras el 1er
     handle seguía vivo y fallaba → RPM 0 en la fase coast-down). */
  s->hSerial = CreateFileA(portName, GENERIC_READ | GENERIC_WRITE,
                           FILE_SHARE_READ | FILE_SHARE_WRITE,
                           NULL, OPEN_EXISTING, 0, NULL);
  if (s->hSerial == INVALID_HANDLE_VALUE) {
    s->connected = 0;
    return -1;
  }
  dcbSerialParams.DCBlength = sizeof(dcbSerialParams);
  if (!GetCommState(s->hSerial, &dcbSerialParams)) {
    CloseHandle(s->hSerial);
    s->connected = 0;
    return -2;
  }
  dcbSerialParams.BaudRate = (DWORD)baudrate;
  dcbSerialParams.ByteSize = 8;
  dcbSerialParams.StopBits = ONESTOPBIT;
  dcbSerialParams.Parity   = NOPARITY;
  if (!SetCommState(s->hSerial, &dcbSerialParams)) {
    CloseHandle(s->hSerial);
    s->connected = 0;
    return -3;
  }
  /* Timeouts: ReadFile retorna en max ~5ms si no llegan datos */
  timeouts.ReadIntervalTimeout         = MAXDWORD; /* retorna si hay algún dato */
  timeouts.ReadTotalTimeoutConstant    = 5;         /* 5 ms máximo total */
  timeouts.ReadTotalTimeoutMultiplier  = 0;
  timeouts.WriteTotalTimeoutConstant   = 50;
  timeouts.WriteTotalTimeoutMultiplier = 0;
  SetCommTimeouts(s->hSerial, &timeouts);
  /* Limpiar buffers */
  PurgeComm(s->hSerial, PURGE_RXCLEAR | PURGE_TXCLEAR);
  s->connected  = 1;
  s->angle_deg  = 0.0f;
  s->rpm        = 0.0f;
  s->raw_angle  = 0;
  s->rx_len     = 0;
  return 0;
}
static void serial_close(SerialHIL *s) {
  if (s->connected) {
    /* Enviar duty = 0 (parar motor): [AA 55 00 00 xor] */
    unsigned char cmd[HIL_CMD_LEN];
    DWORD written;
    cmd[0] = HIL_SYNC_CMD;
    cmd[1] = HIL_SYNC_CMD2;
    cmd[2] = 0; cmd[3] = 0;
    cmd[4] = cmd[0] ^ cmd[1] ^ cmd[2] ^ cmd[3];
    WriteFile(s->hSerial, cmd, HIL_CMD_LEN, &written, NULL);
    CloseHandle(s->hSerial);
    s->connected = 0;
  }
}
/* Envía duty (±14399 cuentas PWM crudas) al MCU: [AA 55 hi lo xor] */
static int serial_send_duty(SerialHIL *s, int16_t duty) {
  unsigned char cmd[HIL_CMD_LEN];
  DWORD written;
  if (!s->connected) return -1;
  cmd[0] = HIL_SYNC_CMD;
  cmd[1] = HIL_SYNC_CMD2;
  cmd[2] = (unsigned char)((uint16_t)duty >> 8);
  cmd[3] = (unsigned char)((uint16_t)duty & 0xFF);
  cmd[4] = cmd[0] ^ cmd[1] ^ cmd[2] ^ cmd[3];
  if (!WriteFile(s->hSerial, cmd, HIL_CMD_LEN, &written, NULL))
    return -1;
  return (written == HIL_CMD_LEN) ? 0 : -1;
}
/* Lee una trama de datos del MCU. Retorna 0 si OK, -1 si fallo/timeout.
 *
 * La recepción es asíncrona: el MCU transmite una trama de 8 bytes cada
 * 2 ms mientras que la simulación avanza a su propio ritmo. Para no
 * desincronizarnos ("sync hunting" que producía tramas rotas cuando se
 * leía 1 byte por paso), acumulamos en un FIFO y buscamos la trama
 * completa dentro de él. Si sobran bytes de una trama previa los
 * descartamos hasta encontrar 0x55.
 */
static int serial_recv_data(SerialHIL *s) {
  DWORD bytesRead = 0;
  unsigned char xor_chk;
  int i;
  float rpm_val;
  if (!s->connected) return -1;
  /* 1. Buscar una trama completa en lo que ya hay en el FIFO */
  for (i = 0; i < 2; i++) {
    while (s->rx_len >= HIL_DATA_LEN) {
      if (s->rx[0] == HIL_SYNC_DATA_V1 || s->rx[0] == HIL_SYNC_DATA_V2) {
        /* Verificar checksum */
        xor_chk = 0;
        for (int j = 0; j < 7; j++) xor_chk ^= s->rx[j];
        if (xor_chk == s->rx[7]) {
          /* Decodificar */
          s->raw_angle = ((uint16_t)s->rx[1] << 8) | s->rx[2];
          s->angle_deg = (float)s->raw_angle * (360.0f / 4096.0f);
          memcpy(&rpm_val, &s->rx[3], 4);
          s->rpm = rpm_val;
          /* Consumir la trama del FIFO */
          s->rx_len -= HIL_DATA_LEN;
          if (s->rx_len > 0) {
            memmove(s->rx, s->rx + HIL_DATA_LEN, s->rx_len);
          }
          return 0;
        }
      }
      /* No era sync o checksum malo: descartar 1 byte y resincronizar */
      s->rx_len -= 1;
      memmove(s->rx, s->rx + 1, s->rx_len);
    }
    if (i == 0) {
      int space = (int)sizeof(s->rx) - s->rx_len;
      if (space > 0) {
        if (ReadFile(s->hSerial, s->rx + s->rx_len, space, &bytesRead, NULL) && bytesRead > 0) {
          s->rx_len += (int)bytesRead;
        }
      }
    }
  }
  return -1;  /* sin trama completa disponible — guardamos lo leído */
}
#endif /* _WIN32 */
#endif /* SERIAL_WIN_H */