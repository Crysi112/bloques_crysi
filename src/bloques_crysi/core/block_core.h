#ifndef BLOQUES_CORE_H
#define BLOQUES_CORE_H
#ifdef _WIN32
  #define API __declspec(dllexport)
#else
  #define API
#endif
/* Opcodes */
enum {
  OP_SRC_CONST = 1,
  OP_SRC_STEP = 2,
  OP_SRC_RAMP = 3,
  OP_SRC_SIN = 4,
  OP_SRC_TRIF = 5,
  OP_GAIN = 6,
  OP_SUM = 7,
  OP_CLARKE = 8,
  OP_INV_CLARKE = 9,
  OP_PARK = 10,
  OP_INV_PARK = 11,
  OP_INTEGRADOR = 12,
  OP_TF = 13,
  OP_PID = 14,
  OP_MAQ_INDUCCION = 15,
  OP_MAQ_SINCRONA = 16,
  OP_MAQ_PMAC = 17,
  OP_MAQ_CC = 18,
  OP_POT_BUCK = 19,
  OP_POT_BOOST = 20,
  OP_POT_BUCKBOOST = 21,
  OP_POT_RECT_3F = 22,
  OP_POT_INV_3F = 23,
  OP_EJE_MECANICO = 24,
  OP_PWM_1F = 25,
  OP_PWM_SPWM = 26,
  OP_PWM_SVPWM = 27,
  OP_HW_SERIAL = 28,
  OP_BATERIA = 29,
  OP_SATURAR = 30,
  OP_RELAY = 31,
  OP_PULSO_RECT = 32,
  OP_PLL = 33,
  OP_POT_INV_1F = 34,
  OP_CARGA_RL_3F = 35,
  OP_TRANSFORMADOR = 36,
  OP_PANEL_SOLAR = 37,
  OP_MEDIDOR_POTENCIA = 38,
  OP_MAQ_DC_PM = 39,
  OP_INTERRUPTOR = 40,
  OP_DIODO = 41,
  OP_PUENTE_INV_3F = 42,
  OP_PUENTE_INV_1F = 43,
  OP_MUX = 44,
  OP_DEMUX = 45,
  OP_LUT1D = 46,
  OP_LUT2D = 47,
  OP_LUT3D = 48,
  OP_LOGICO = 49,
  OP_RELACIONAL = 50,
  OP_LIM_RAPIDEZ = 51,
  OP_RETENEDOR = 52,
  OP_MAQ_ESTADOS = 53,
  OP_MASA_TERMICA = 54,
  OP_RES_TERMICA = 55,
  OP_ENGRANAJE = 56,
  OP_EJE_FLEXIBLE = 57,
  OP_EMBRAGUE = 58,
  OP_SRC_CSV = 59,
  OP_FALLO_PROG = 60,
  OP_FALLO_EVENTO = 61,
  OP_MULTIPLICADOR = 62,
  OP_SAT_VECTORIAL = 63,
  OP_SRC_TABLE = 64,
  OP_BATERIA_ECM = 65,
  OP_CALCULO_IDC = 66,
  OP_VEHICULO = 67,
  OP_RESISTENCIA = 68,
  OP_INDUCTOR = 69,
  OP_CAPACITOR = 70,
  OP_QD = 71,
  OP_CARGA_PQ_3F = 72,
  OP_CARGA_PQ_1F = 73,
  OP_MNA = 74,
  OP_MUTUAL_INDUCTOR = 75,
  OP_VCVS = 76,
  OP_VCCS = 77,
};
typedef struct ModeloC ModeloC;
typedef struct BloqueC BloqueC;
struct BloqueC {
  int op;
  int n_in;
  long long *in_idx;
  int n_out;
  long long *out_idx;
  int n_param;
  double *param;
  int n_state;
  double *state;
  int n_ws;
  double *ws;
  double dt;
  /* Multitasa / Discreto */
  double Ts;
  double t_next_update;
  /* Punteros a funciones del bloque */
  void (*init)(BloqueC *bl, ModeloC *m);
  void (*eval_estatico)(BloqueC *bl, ModeloC *m, double *maxdelta);
  void (*deriv)(BloqueC *bl, ModeloC *m, const double *x, double *dx);
  void (*out)(BloqueC *bl, ModeloC *m, const double *x);
  void (*update)(BloqueC *bl, ModeloC *m);
};
struct ModeloC {
  int n_bloques;
  BloqueC *bloques;
  int n_sig;
  double *sig;
  int n_alg;
  long long *alg_list;
  int max_iter;
  double tol;
  double w_opt;
  int method;  /* 0 = Euler, 1 = RK4 */
  double t;
  double t_fin;
  double dt;
  int error_flag;  /* 1 si el lazo algebraico no convergio en algun paso */
};
/* ---------- HIL workspace size helper ----------
 * HIL_WS_DOUBLES: numero de doubles que ocupa SerialHIL en ws[] mas
 * 2 slots de control (ws[0]=init, ws[1]=last_t). Se calcula en tiempo
 * de compilacion para que Python no tenga que adivinar el tamano.
 * La formula redondea hacia arriba al double mas cercano.
 */
#ifdef _WIN32
#include "serial_win.h"
#define HIL_WS_DOUBLES \
    (2 + (int)((sizeof(SerialHIL) + sizeof(double) - 1) / sizeof(double)))
#else
#define HIL_WS_DOUBLES 16   /* fallback conservador en plataformas no-Win32 */
#endif

API void m_sim_run(ModeloC *m, int n_steps, int n_rec, long long *rec_idx,
                   double *rec_buf);
API int m_sim_iniciar(ModeloC *m);
API int m_sim_paso(ModeloC *m);
API int m_sim_guardar(ModeloC *m, double *buf);
API void m_sim_restaurar(ModeloC *m, const double *buf);
API int m_hil_ws_size(void);   /* devuelve HIL_WS_DOUBLES en tiempo de ejecucion */
#endif
