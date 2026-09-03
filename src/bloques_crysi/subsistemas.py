from bloques_crysi import (
    Modelo, PID, TransformadaQD, SaturarVectorial,
    InvPark, InvClarke, PuenteInversorTrifasico,
    Suma, Ganancia, Multiplicador, Multiplexor,
    FuenteConstante, Puerto,
)
def _par(bloque, tipo, offset, n):
    return Puerto(bloque, tipo, offset, n)
def crear_control_foc(modelo: Modelo, nombre: str,
                      Kp_id: float = 0.0, Ki_id: float = 200.0,
                      Kp_iq: float = 0.0, Ki_iq: float = 200.0,
                      Vdc: float = 400.0, Ld: float = 1e-3, Lq: float = 1e-3,
                      lam_m: float = 0.1, P: int = 6) -> object:
    pref = f"{nombre}_"
    qd = modelo.add(TransformadaQD(f"{pref}qd"))
    the_blk = modelo.add(Ganancia(f"{pref}the_pass", 1.0))
    modelo.conectar(the_blk.salida, qd.th)
    cero_a = modelo.add(FuenteConstante(f"{pref}vabc_cero_a", 0.0))
    cero_b = modelo.add(FuenteConstante(f"{pref}vabc_cero_b", 0.0))
    cero_c = modelo.add(FuenteConstante(f"{pref}vabc_cero_c", 0.0))
    mux_vabc = modelo.add(Multiplexor(f"{pref}vabc_mux", n_canales=3))
    modelo.conectar(cero_a.salida, mux_vabc.entradas[0])
    modelo.conectar(cero_b.salida, mux_vabc.entradas[1])
    modelo.conectar(cero_c.salida, mux_vabc.entradas[2])
    modelo.conectar(mux_vabc.salida, qd.vabc)
    id_ref = modelo.add(FuenteConstante(f"{pref}id_ref", 0.0))
    iq_ref = modelo.add(FuenteConstante(f"{pref}iq_ref", 0.0))
    sum_d = modelo.add(Suma(f"{pref}sum_d", (1.0, -1.0)))
    sum_q = modelo.add(Suma(f"{pref}sum_q", (1.0, -1.0)))
    modelo.conectar(id_ref.salida, _par(sum_d, "ent", 0, 1))
    modelo.conectar(_par(qd, "sal", 3, 1), _par(sum_d, "ent", 1, 1))
    modelo.conectar(_par(qd, "sal", 2, 1), _par(sum_q, "ent", 1, 1))
    pid_d = modelo.add(PID(f"{pref}pid_d", Kp=Kp_id, Ki=Ki_id, u_min=-300, u_max=300))
    pid_q = modelo.add(PID(f"{pref}pid_q", Kp=Kp_iq, Ki=Ki_iq, u_min=-300, u_max=300))
    modelo.conectar(sum_d.salida, pid_d.entrada)
    modelo.conectar(sum_q.salida, pid_q.entrada)
    k_we = modelo.add(Ganancia(f"{pref}k_we", P / 2.0))
    lam_blk = modelo.add(FuenteConstante(f"{pref}lam", lam_m))
    mult_we_lq = modelo.add(Multiplicador(f"{pref}mult_we_lq"))
    mult_we_ld = modelo.add(Multiplicador(f"{pref}mult_we_ld"))
    mult_we_lm = modelo.add(Multiplicador(f"{pref}mult_we_lm"))
    modelo.conectar(k_we.salida, _par(mult_we_lq, "ent", 0, 1))
    modelo.conectar(_par(qd, "sal", 2, 1), _par(mult_we_lq, "ent", 1, 1))
    modelo.conectar(k_we.salida, _par(mult_we_ld, "ent", 0, 1))
    modelo.conectar(_par(qd, "sal", 3, 1), _par(mult_we_ld, "ent", 1, 1))
    modelo.conectar(k_we.salida, _par(mult_we_lm, "ent", 0, 1))
    modelo.conectar(lam_blk.salida, _par(mult_we_lm, "ent", 1, 1))
    lq_blk = modelo.add(Ganancia(f"{pref}lq", Lq))
    ld_blk = modelo.add(Ganancia(f"{pref}ld", Ld))
    modelo.conectar(mult_we_lq.salida, lq_blk.entrada)
    modelo.conectar(mult_we_ld.salida, ld_blk.entrada)
    sum_vd = modelo.add(Suma(f"{pref}sum_vd", (1.0, -1.0)))
    sum_vq = modelo.add(Suma(f"{pref}sum_vq", (1.0, 1.0, 1.0)))
    modelo.conectar(pid_d.salida, _par(sum_vd, "ent", 0, 1))
    modelo.conectar(lq_blk.salida, _par(sum_vd, "ent", 1, 1))
    modelo.conectar(pid_q.salida, _par(sum_vq, "ent", 0, 1))
    modelo.conectar(ld_blk.salida, _par(sum_vq, "ent", 1, 1))
    modelo.conectar(mult_we_lm.salida, _par(sum_vq, "ent", 2, 1))
    VMAX = Vdc / (3 ** 0.5)
    sat = modelo.add(SaturarVectorial(f"{pref}sat", Vmax=VMAX))
    modelo.conectar(sum_vd.salida, _par(sat, "ent", 0, 1))
    modelo.conectar(sum_vq.salida, _par(sat, "ent", 1, 1))
    ipark = modelo.add(InvPark(f"{pref}ipark"))
    iclarke = modelo.add(InvClarke(f"{pref}iclarke"))
    puente = modelo.add(PuenteInversorTrifasico(f"{pref}puente", promediado=True))
    vdc_src = modelo.add(FuenteConstante(f"{pref}vdc", Vdc))
    modelo.conectar(_par(sat, "sal", 0, 1), _par(ipark, "ent", 0, 1))
    modelo.conectar(_par(sat, "sal", 1, 1), _par(ipark, "ent", 1, 1))
    modelo.conectar(the_blk.salida, _par(ipark, "ent", 2, 1))
    modelo.conectar(ipark.salida, iclarke.entrada)
    modelo.conectar(vdc_src.salida, _par(puente, "ent", 0, 1))
    modelo.conectar(_par(iclarke, "sal", 0, 1), _par(puente, "ent", 1, 1))
    modelo.conectar(_par(iclarke, "sal", 1, 1), _par(puente, "ent", 2, 1))
    modelo.conectar(_par(iclarke, "sal", 2, 1), _par(puente, "ent", 3, 1))
    class _FOCInterface:
        def __init__(self, qd, sum_d, sum_q, k_we, puente, the_blk, Vdc):
            self.corrientes_abc = qd.iabc
            self.tension_abc = puente.salida
            self.th_e = the_blk.entrada
            self.omega_m = k_we.entrada
            self.id_ref = Puerto(sum_d, "ent", 0, 1)
            self.iq_ref = Puerto(sum_q, "ent", 0, 1)
            self.Vdc = Vdc
    return _FOCInterface(qd, sum_d, sum_q, k_we, puente, the_blk, Vdc)
def crear_driver_modelo(modelo: Modelo, nombre: str,
                        Kp: float = 150.0, Ki: float = 20.0, Kd: float = 0.0,
                        Tf: float = 0.01, T_min: float = -300.0, T_max: float = 300.0,
                        gear_ratio: float = 8.0, wheel_radius: float = 0.33) -> object:
    pref = f"{nombre}_"
    vref_blk = modelo.add(Ganancia(f"{pref}vref_pass", 1.0))
    sum_v = modelo.add(Suma(f"{pref}sum_v", (1.0, -1.0)))
    pid_v = modelo.add(PID(f"{pref}pid_v", Kp=Kp, Ki=Ki, Kd=Kd, Tf=Tf,
                           u_min=T_min, u_max=T_max))
    k_omega = modelo.add(Ganancia(f"{pref}k_omega", gear_ratio / wheel_radius))
    modelo.conectar(vref_blk.salida, _par(sum_v, "ent", 0, 1))
    modelo.conectar(vref_blk.salida, k_omega.entrada)
    modelo.conectar(sum_v.salida, pid_v.entrada)
    class _DriverInterface:
        def __init__(self, vref_blk, sum_v, pid_v, k_omega):
            self.v_ref = vref_blk.entrada
            self.v_veh = Puerto(sum_v, "ent", 1, 1)
            self.T_ref = pid_v.salida
            self.omega_ref = k_omega.salida
    return _DriverInterface(vref_blk, sum_v, pid_v, k_omega)
