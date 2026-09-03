import ctypes
import ctypes.wintypes
import struct
import time
import matplotlib.pyplot as plt
from bloques_crysi import (
    Modelo, MotorHardwareCH32, FuenteEscalon, FuenteRampa, FuenteSeno,
    FalloProgramado,
)

try:
    COM_PORT
except NameError:
    COM_PORT = int(input("Puerto COM (ej. 3): ").strip())

def _enviar_stop():
    try:
        k32  = ctypes.windll.kernel32
        port = f"\\\\.\\COM{COM_PORT}".encode()
        h = k32.CreateFileA(port, 0xC0000000, 0, None, 3, 0, None)
        if h and h != ctypes.wintypes.HANDLE(-1).value:
            pkt = bytes([0xAA, 0x55, 0x00, 0x00, 0xFF])
            written = ctypes.wintypes.DWORD(0)
            k32.WriteFile(h, pkt, 5, ctypes.byref(written), None)
            time.sleep(0.15)
            k32.CloseHandle(h)
    except Exception:
        pass

def _graficar(res, nombre_ref, nombre_rpm, titulo):
    import numpy as np
    t    = res.t
    rpm  = np.asarray(res[nombre_rpm], dtype=float)
    duty = np.asarray(res[nombre_ref], dtype=float)

    rpm_max = max(np.percentile(np.abs(rpm), 99) * 1.3, 50.0)
    rpm_clip = np.clip(rpm, -rpm_max, rpm_max)

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax2 = ax1.twinx()

    l1, = ax1.plot(t, rpm_clip,  color='orange',    lw=2,   label='Velocidad (RPM)')
    l2, = ax2.plot(t, duty,      color='royalblue', lw=1.5, linestyle='--', label='Duty cycle')

    ax1.set_xlabel("Tiempo (s)")
    ax1.set_ylabel("Velocidad (RPM)", color='orange')
    ax2.set_ylabel("Duty cycle",      color='royalblue')
    ax1.set_ylim(-rpm_max * 0.1, rpm_max)
    ax2.set_ylim(-0.1, 1.1)
    ax1.grid(True)
    ax1.set_title(titulo)
    ax1.legend(handles=[l1, l2], loc='upper left')
    plt.tight_layout()
    plt.show()

def prueba_escalon():
    print("\n--- Prueba de ESCALÓN ---")
    amp  = float(input("  Amplitud (0.0–1.0)             [0.8]: ") or "0.8")
    t0   = float(input("  Tiempo del escalón [s]          [2.0]: ") or "2.0")
    tfin = float(input("  Duración total de la prueba [s] [6.0]: ") or "6.0")

    print("  Deteniendo motor previo...")
    _enviar_stop()
    time.sleep(1.5)

    m = Modelo(dt=0.002)
    ref   = m.add(FuenteEscalon("Ref", valor_final=amp, t_paso=t0, valor_inicial=0.0))
    motor = m.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m.conectar(ref.salida, motor.duty)

    print(f"\n  Corriendo {tfin}s — el motor arrancará en t = {t0}s...")
    try:
        res = m.run(t_fin=tfin, registrar=[ref, motor.sensorVelocidad()])
    finally:
        _enviar_stop()
        print("  Motor detenido.")

    _graficar(res, "Ref", "Motor_rpm", f"Escalón {amp*100:.0f}% en t={t0}s")

def prueba_rampa():
    print("\n--- Prueba de RAMPA ---")
    k    = float(input("  Pendiente (duty/s)              [0.15]: ") or "0.15")
    t0   = float(input("  Inicio de la rampa [s]           [1.0]: ") or "1.0")
    tfin = float(input("  Duración total de la prueba [s]  [8.0]: ") or "8.0")

    print("  Deteniendo motor previo...")
    _enviar_stop()
    time.sleep(1.5)

    m = Modelo(dt=0.002)
    ref   = m.add(FuenteRampa("Ref", pendiente=k, t_inicio=t0, offset=0.0))
    motor = m.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m.conectar(ref.salida, motor.duty)

    print(f"\n  Corriendo {tfin}s — rampa inicia en t = {t0}s...")
    try:
        res = m.run(t_fin=tfin, registrar=[ref, motor.sensorVelocidad()])
    finally:
        _enviar_stop()
        print("  Motor detenido.")

    _graficar(res, "Ref", "Motor_rpm", f"Rampa k={k} desde t={t0}s")

def prueba_seno():
    print("\n--- Prueba de SENO ---")
    amp  = float(input("  Amplitud (0.0–1.0)              [0.7]: ") or "0.7")
    frec = float(input("  Frecuencia [Hz]                 [0.3]: ") or "0.3")
    tfin = float(input("  Duración total de la prueba [s] [10.0]: ") or "10.0")

    print("  Deteniendo motor previo...")
    _enviar_stop()
    time.sleep(1.5)

    m = Modelo(dt=0.002)
    ref   = m.add(FuenteSeno("Ref", amplitud=amp, frecuencia=frec, fase=0.0, offset=0.0))
    motor = m.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m.conectar(ref.salida, motor.duty)

    print(f"\n  Corriendo {tfin}s — seno {amp}×sin(2π×{frec}t)...")
    try:
        res = m.run(t_fin=tfin, registrar=[ref, motor.sensorVelocidad()])
    finally:
        _enviar_stop()
        print("  Motor detenido.")

    _graficar(res, "Ref", "Motor_rpm", f"Seno {amp}, {frec} Hz")

def prueba_identificacion():
    import numpy as np

    print("\n--- Identificación de parámetros del motor ---")
    Vcc  = float(input("  Voltaje de la fuente [V]         [12.0]: ") or "12.0")
    R    = float(input("  Resistencia medida, rotor bloqueado [Ω] [21.3]: ") or "21.3")
    duty = float(input("  Duty del escalón (0–1)            [0.6]: ") or "0.6")
    t_paso  = float(input("  Tiempo del escalón [s]            [0.5]: ") or "0.5")
    t_fin_a = float(input("  Duración escalón [s]              [8.0]: ") or "8.0")
    t_fin_b = float(input("  Duración coast-down [s]           [8.0]: ") or "8.0")

    print("  Deteniendo motor previo...")
    _enviar_stop()
    time.sleep(1.5)

    print(f"\n  Fase A — escalón {duty*100:.0f}% en t={t_paso}s...")
    m = Modelo(dt=0.002)
    ref   = m.add(FuenteEscalon("Ref", valor_final=duty, t_paso=t_paso, valor_inicial=0.0))
    motor = m.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m.conectar(ref.salida, motor.duty)
    try:
        resA = m.run(t_fin=t_fin_a, registrar=[ref, motor.sensorVelocidad()])
    finally:
        m.cerrar_hw()
        _enviar_stop()
        time.sleep(1.5)

    ta   = np.asarray(resA.t, dtype=float)
    rpmA = np.asarray(resA["Motor_rpm"], dtype=float)

    t_b = 1.0
    print(f"  Fase B — coast-down (se corta el duty en t={t_b}s)...")
    m2 = Modelo(dt=0.002)
    ref2   = m2.add(FuenteEscalon("Ref2", valor_final=0.0, t_paso=t_b, valor_inicial=duty))
    motor2 = m2.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m2.conectar(ref2.salida, motor2.duty)
    try:
        resB = m2.run(t_fin=t_fin_b, registrar=[ref2, motor2.sensorVelocidad()])
    finally:
        m2.cerrar_hw()
        _enviar_stop()
        print("  Motor detenido.")

    tb   = np.asarray(resB.t, dtype=float)
    rpmB = np.asarray(resB["Motor_rpm"], dtype=float)

    def _mediana_movil(x, w=25):
        n = len(x)
        out = np.empty_like(x)
        half = w // 2
        for i in range(n):
            lo = max(0, i - half)
            hi = min(n, i + half + 1)
            out[i] = np.median(x[lo:hi])
        return out

    rpmA_f = _mediana_movil(np.asarray(rpmA, dtype=float))
    rpmB_f = _mediana_movil(np.asarray(rpmB, dtype=float))
    rpmA_f = np.nan_to_num(rpmA_f, nan=0.0)
    rpmB_f = np.nan_to_num(rpmB_f, nan=0.0)

    mask_ss  = ta > t_fin_a - 1.5
    rpm_ss   = float(np.median(rpmA_f[mask_ss]))
    if rpm_ss <= 0:
        print("  ERROR: RPM estable no válido (motor no llegó a girar?).")
        return

    rpm_max = float(np.max(rpmA_f))
    if rpm_max > 1.15 * rpm_ss:
        print(f"  (aviso: pico de arranque {rpm_max:.0f} RPM vs estable {rpm_ss:.0f} "
              f"— artefacto de medición, ignorado en el ajuste)")

    t0  = t_paso
    sel = ta >= t0
    t_tr   = ta[sel]
    y_tr   = np.clip(rpmA_f[sel] / rpm_ss, 0.0, 1.0)
    i_85   = np.argmax(y_tr >= 0.85)
    if i_85 == 0 and (y_tr[0] < 0.85).all():
        i_85 = len(t_tr)
    subida = t_tr[:i_85]
    y_sub  = y_tr[:i_85]
    lo, hi = y_sub >= 0.25, y_sub <= 0.75
    if np.sum(lo & hi) >= 5:
        p = np.polyfit(subida[lo & hi], np.log(1.0 - y_sub[lo & hi]), 1)
        tau_m = float(-1.0 / p[0])
    elif len(y_sub) > 0:
        i63 = np.argmax(y_sub >= 0.632)
        tau_m = float(subida[max(i63 - 1, 0)] - t0) if i63 > 0 else np.nan
    else:
        tau_m = np.nan

    pre = (tb >= t_b - 0.5) & (tb <= t_b + 0.02)
    rpm0 = float(np.median(rpmB_f[pre])) if pre.sum() > 0 else 0.0

    sel_b = tb >= t_b + 0.05
    t_cd  = tb[sel_b]
    r_cd  = rpmB_f[sel_b]
    ok    = r_cd > 30
    t_cd  = t_cd[ok]; r_cd = r_cd[ok]
    rpm_final = float(r_cd.min()) if len(r_cd) > 0 else 0.0

    tau_d = np.nan
    if len(t_cd) >= 15 and rpm0 > 50:
        lo2 = (r_cd <= 0.95 * rpm0) & (r_cd >= 0.15 * rpm0)
        if lo2.sum() < 5:
            lo2 = np.ones_like(r_cd, dtype=bool)
        p2    = np.polyfit(t_cd[lo2], np.log(np.clip(r_cd[lo2], 1e-9, None)), 1)
        tau_d = float(-1.0 / p2[0])
        if not np.isfinite(tau_d) or tau_d <= 0:
            sel2   = (r_cd >= 0.5 * rpm0) & (r_cd <= 0.9 * rpm0)
            if sel2.sum() >= 5:
                p3    = np.polyfit(t_cd[sel2], np.log(np.clip(r_cd[sel2], 1e-9, None)), 1)
                tau_d = float(-1.0 / p3[0])

    print(f"  (coast-down: rpm0={rpm0:.0f}, rpm_final={rpm_final:.0f}, "
          f"puntos={len(t_cd)}, τd = {tau_d if np.isfinite(tau_d) else float('nan'):.3f} s)")

    omega_ss = rpm_ss * 2.0 * np.pi / 60.0
    Ke       = Vcc * duty / omega_ss
    Kt       = Ke
    I_ss     = (Vcc * duty - Ke * omega_ss) / R
    b        = (Kt * I_ss) / omega_ss if omega_ss > 0 else np.nan
    J        = tau_m * (b * R + Kt * Ke) / R if (tau_m > 0 and R > 0) else np.nan
    b_cd = (J / tau_d) if (np.isfinite(tau_d) and tau_d > 0 and J > 0) else np.nan
    I_start  = Vcc / R
    T_stall  = Kt * I_start
    Km       = omega_ss / duty

    print("\n" + "=" * 55)
    print("  PARÁMETROS IDENTIFICADOS (Vcc = %.1f V, R = %.2f Ω)" % (Vcc, R))
    print("=" * 55)
    print(f"  RPM estable        : {rpm_ss:9.1f} RPM   ({omega_ss:8.2f} rad/s)")
    print(f"  Const. velocidad Ke: {Ke:9.4f} V·s/rad")
    print(f"  Const. torque Kt   : {Kt:9.4f} N·m/A   (≈ Ke en SI)")
    print(f"  Corriente régimen  : {I_ss:9.4f} A    (en vacío ≈0: b NO es fiable por balance eléctrico)")
    print(f"  Fricción viscosa b : {b:9.5f} N·m·s/rad   ← usar coast-down)")
    print(f"  τm (escalón)       : {tau_m:9.3f} s")
    print(f"  Inercia J (escalón): {J:9.2e} kg·m²")
    if np.isfinite(tau_d):
        print(f"  τd (coast-down)    : {tau_d:9.3f} s")
        print(f"  b (coast-down)     : {b_cd:9.5f} N·m·s/rad   ← fiabilidad alta)")
        print(f"  J_verif            : {J:9.2e} kg·m²   (consistencia J = b_cd·τd = {b_cd*tau_d:.2e})")
    else:
        print("  τd (coast-down)    :     N/D — no se pudo ajustar el decaimiento")
    print(f"  Corriente arranque : {I_start:9.3f} A   (rotor bloqueado)")
    print(f"  T_stall teórico    : {T_stall:9.4f} N·m")
    print(f"  Ganancia Km        : {Km:9.3f} rad/s por unidad de duty")
    print("=" * 55)
    print("  Modelo:  ω(s) = Km / (τm·s + 1) · duty(s)")
    print(f"          Km = {Km:9.3f} rad/s/duty,  τm = {tau_m:9.3f} s")
    print("=" * 55)

    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(ta, rpmA, color='lightgray', lw=1, label='RPM cruda')
    ax1.plot(ta, rpmA_f, color='orange', lw=1.5, label='RPM suavizada')
    ax1.axhline(rpm_ss, color='gray', ls='--', lw=1, label=f'RPM estable = {rpm_ss:.0f}')
    ax1.set_xlabel("Tiempo (s)"); ax1.set_ylabel("RPM"); ax1.set_title("Escalón")
    ax1.grid(True); ax1.legend(loc='lower right')
    ax2.plot(tb, rpmB, color='lightgray', lw=1, label='RPM cruda')
    ax2.plot(tb, rpmB_f, color='orange', lw=1.5, label='RPM suavizada')
    if np.isfinite(tau_d):
        t_fit  = np.linspace(t_b, t_fin_b, 200)
        rpm_fit = rpm0 * np.exp(-(t_fit - t_b) / tau_d)
        ax2.plot(t_fit, rpm_fit, 'b--', lw=1.5, label=f'τd = {tau_d:.2f} s')
    ax2.set_xlabel("Tiempo (s)"); ax2.set_ylabel("RPM"); ax2.set_title("Coast-down")
    ax2.grid(True); ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

def prueba_fallo():
    print("\n--- Prueba de FALLO PROGRAMADO (corte de duty) ---")
    amp  = float(input("  Amplitud (0.0–1.0)              [0.6]: ") or "0.6")
    t0   = float(input("  Tiempo del fallo [s]            [4.0]: ") or "4.0")
    tfin = float(input("  Duración total de la prueba [s] [8.0]: ") or "8.0")

    print("  Deteniendo motor previo...")
    _enviar_stop()
    time.sleep(1.5)

    m = Modelo(dt=0.002)
    ref  = m.add(FuenteEscalon("Ref", valor_final=amp, t_paso=0.5,
                               valor_inicial=0.0))
    fallo = m.add(FalloProgramado("Fallo", t_fallo=t0, valor=0.0, modo=0))
    motor = m.add(MotorHardwareCH32("Motor", puerto_com=COM_PORT, baudrate=115200))
    m.conectar(ref.salida, fallo.entrada)
    m.conectar(fallo.salida, motor.duty)

    print(f"\n  Corriendo {tfin}s — el duty se corta a 0 en t = {t0}s...")
    try:
        res = m.run(t_fin=tfin, registrar=[fallo, motor.sensorVelocidad()])
    finally:
        _enviar_stop()
        print("  Motor detenido.")

    _graficar(res, "Fallo", "Motor_rpm", f"Fallo programado en t={t0}s")

if __name__ == "__main__":
    print("=" * 45)
    print("  HIL Motor DC — Pruebas de entrada")
    print("=" * 45)
    print("  1. Escalón")
    print("  2. Rampa")
    print("  3. Seno")
    print("  4. Identificación de parámetros")
    print("  5. Fallo programado (corte de duty)")
    op = input("Selecciona (1/2/3/4/5): ").strip()
    {"1": prueba_escalon, "2": prueba_rampa, "3": prueba_seno,
     "4": prueba_identificacion, "5": prueba_fallo}.get(
        op, lambda: print("Opción inválida."))()
