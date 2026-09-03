def _ruta(rel):
    for p in (os.path.join(os.getcwd(), 'ejemplos', rel),
              os.path.join(os.getcwd(), rel)):
        if os.path.exists(p):
            return p
    return os.path.join(os.getcwd(), rel)

import os

import numpy as np

from bloques_crysi import (
    Modelo, FuenteConstante, FuenteCSV, FalloProgramado, FalloEvento,
    LimitadorRapidez, FiltroPasoBajo, MaquinaEstados, RetenedorDisparado,
    EjeFlexible, Embrague, Engranaje, MasaTermica, ResistenciaTermica,
    Multiplexor, Multiplicador, Tabla1D, Logico, Relacional,
    FuncionTransferencia, Suma, Ganancia, Scope, Puerto,
)

DT = 1e-3
T_END = 30.0
CSV_PERFIL = _ruta("data/perfil_velocidad.csv")
T_AMB = 25.0

m = Modelo(dt=DT)
perfil = m.add(FuenteCSV("perfil", CSV_PERFIL))
limit = m.add(LimitadorRapidez("limit", subida=60.0, bajada=60.0))
mot = m.add(FuncionTransferencia("mot", num=[1.0], den=[0.05, 1.0]))
eje = m.add(EjeFlexible("eje", K=30.0, B=1.0))
inercia = m.add(FuncionTransferencia("inercia", num=[1.0],
                                     den=[0.05, 0.0]))
m.conectar(perfil.salida, limit.entrada)
m.conectar(limit.salida, mot.entrada)
m.conectar(mot.salida, eje.w1)
m.conectar(inercia.salida, eje.w2)
m.conectar(eje.salida, inercia.entrada)

mux = m.add(Multiplexor("mux", n_canales=2))
m.conectar(perfil.salida, mux.entradas[0])
m.conectar(inercia.salida, mux.entradas[1])
sc = m.add(Scope("scope1", anchos=[2, 1],
                 guiones=["w_ref, w2 (rad/s)", "T_eje (N.m)"]))
m.conectar(mux.salida, sc.canales[0])
m.conectar(eje.salida, sc.canales[1])

res = m.run(t_fin=10.0, registrar=[sc, inercia, eje, perfil])
print("Paso 1 — cadena mecanica minima (perfil -> motor -> eje -> inercia)")
print("=" * 62)
print(f"w_ref final : {res['perfil'][-1]:6.2f} rad/s")
print(f"w2 final    : {res['inercia'][-1]:6.2f} rad/s  "
      f"(sigue al perfil con el filtro del motor)")
print(f"T_eje final : {res['eje'][-1]:6.2f} N.m  (resorte torsional "
      f"K = 30 N.m/rad)")

m = Modelo(dt=DT)
perfil = m.add(FuenteCSV("perfil", CSV_PERFIL))
uno = m.add(FuenteConstante("uno", 1.0))
cero = m.add(FuenteConstante("cero", 0.5))
ochenta5 = m.add(FuenteConstante("ochenta5", 85.0))
tc = m.add(FuenteConstante("tc", 4.0))
fallo_carga = m.add(FalloProgramado("fallo_carga", t_fallo=18.0,
                                    valor=2.0, modo=1))
limit = m.add(LimitadorRapidez("limit", subida=60.0, bajada=60.0))
filtro = m.add(FiltroPasoBajo("filtro", fc=8.0, orden=1))
mot = m.add(FuncionTransferencia("mot", num=[1.0], den=[0.05, 1.0]))
eje = m.add(EjeFlexible("eje", K=30.0, B=1.0))
emb = m.add(Embrague("emb", T_max=15.0, umbral=0.5))
inercia = m.add(FuncionTransferencia("inercia", num=[1.0],
                                     den=[0.05, 0.0]))
eng = m.add(Engranaje("eng", relacion=2.0))
and_emb = m.add(Logico("and_emb", opcion="AND", n_entradas=2))
ref_activa = m.add(Relacional("ref_activa", opcion=">"))
temp_ok = m.add(Relacional("temp_ok", opcion="<"))
suma_t = m.add(Suma("suma_t", signos=[1.0, -1.0]))

m.conectar(perfil.salida, limit.entrada)
m.conectar(limit.salida, mot.entrada)
m.conectar(mot.salida, eje.w1)
m.conectar(inercia.salida, eje.w2)
m.conectar(eje.salida, emb.entrada)
m.conectar(and_emb.salida, emb.control)
m.conectar(emb.salida, Puerto(suma_t, "ent", 0, 1))
m.conectar(tc.salida, fallo_carga.entrada)
m.conectar(fallo_carga.salida, Puerto(suma_t, "ent", 1, 1))
m.conectar(suma_t.salida, inercia.entrada)
m.conectar(inercia.salida, filtro.entrada)
m.conectar(filtro.salida, eng.w1)
m.conectar(eje.salida, eng.T1)
m.conectar(perfil.salida, ref_activa.a)
m.conectar(cero.salida, ref_activa.b)
m.conectar(uno.salida, temp_ok.a)
m.conectar(ochenta5.salida, temp_ok.b)
m.conectar(ref_activa.salida, and_emb.entradas[0])
m.conectar(temp_ok.salida, and_emb.entradas[1])

sc = m.add(Scope("scope2", anchos=[2, 1, 1],
                 guiones=["w_ref, w2 filtrada (rad/s)", "T_eje (N.m)",
                          "w_carga (rad/s)"]))
mux = m.add(Multiplexor("mux", n_canales=2))
m.conectar(perfil.salida, mux.entradas[0])
m.conectar(filtro.salida, mux.entradas[1])
m.conectar(mux.salida, sc.canales[0])
m.conectar(eje.salida, sc.canales[1])
m.conectar(Puerto(eng, "sal", 0, 1), sc.canales[2])

res = m.run(t_fin=20.0, registrar=[sc, filtro, eng, eje, emb])
print("Paso 2 — embrague con logica AND + carga + engranaje 2:1")
print("=" * 62)
print(f"T_eje final   : {res['eje'][-1]:6.2f} N.m  (el embrague "
      f"transmite el par del motor)")
print(f"w_carga final : {res['eng'][-1, 0]:6.2f} rad/s  "
      f"(= w2/2 por el engranaje)")
print(f"par de carga  : 4 N.m  (+2 N.m a los 18 s por el "
      f"FalloProgramado)")

m = Modelo(dt=DT)

perfil = m.add(FuenteCSV("perfil", CSV_PERFIL))
cero = m.add(FuenteConstante("cero", 0.5))
uno = m.add(FuenteConstante("uno", 1.0))
ochenta5 = m.add(FuenteConstante("ochenta5", 85.0))
amb = m.add(FuenteConstante("amb", T_AMB))
tc = m.add(FuenteConstante("tc", 4.0))
fallo_carga = m.add(FalloProgramado("fallo_carga", t_fallo=18.0,
                                    valor=2.0, modo=1))
fallo_t = m.add(FalloEvento("fallo_t", umbral=85.0, valor=0.0, modo=0))

limit = m.add(LimitadorRapidez("limit", subida=60.0, bajada=60.0))
filtro = m.add(FiltroPasoBajo("filtro", fc=8.0, orden=1))
sup = m.add(MaquinaEstados(
    "sup", n_estados=3, n_entradas=2,
    transiciones=[(0, 1, 0, ">", 0.5),
                  (1, 2, 1, ">", 85.0),
                  (2, 1, 1, "<", 80.0)],
    estado_inicial=0))
reten = m.add(RetenedorDisparado("reten", umbral=0.5, valor_inicial=T_AMB))

mot = m.add(FuncionTransferencia("mot", num=[1.0], den=[0.05, 1.0]))
eje = m.add(EjeFlexible("eje", K=30.0, B=10.0))
emb = m.add(Embrague("emb", T_max=15.0, umbral=0.5))
inercia = m.add(FuncionTransferencia("inercia", num=[1.0], den=[0.05, 0.0]))
eng = m.add(Engranaje("eng", relacion=2.0))
masa = m.add(MasaTermica("masa", C_th=40.0, T_inicial=T_AMB, n_entradas=2,
                         T_amb=T_AMB, R_amb=0.0))
rt = m.add(ResistenciaTermica("rt", R=0.45))
gr = m.add(Ganancia("gr", -1.0))

mux = m.add(Multiplexor("mux", n_canales=3))
mux_p = m.add(Multiplexor("mux_p", n_canales=2))
mult = m.add(Multiplicador("mult"))
mux_c = m.add(Multiplexor("mux_c", n_canales=2))
mult_c = m.add(Multiplicador("mult_c"))
kperd = m.add(Ganancia("kperd", 0.4))
gf = m.add(Ganancia("gf", 0.25))
kf = m.add(Suma("kf", signos=[1.0, 1.0]))
carga_tab = m.add(Tabla1D("carga_tab", [-160, -0.5, 0.5, 160],
                          [4.0, 0.0, 0.0, 4.0]))
tabla = m.add(Tabla1D("tabla", [0, 40, 80, 120, 160],
                      [0.55, 0.82, 0.90, 0.92, 0.88]))
and_emb = m.add(Logico("and_emb", opcion="AND", n_entradas=2))
ref_activa = m.add(Relacional("ref_activa", opcion=">"))
temp_ok = m.add(Relacional("temp_ok", opcion="<"))
en_marcha = m.add(Relacional("en_marcha", opcion=">"))
suma_t = m.add(Suma("suma_t", signos=[1.0, -1.0]))

m.conectar(perfil.salida, fallo_t.senal)
m.conectar(masa.salida, fallo_t.disparo)
m.conectar(fallo_t.salida, limit.entrada)
m.conectar(limit.salida, mot.entrada)

m.conectar(mot.salida, eje.w1)
m.conectar(inercia.salida, eje.w2)
m.conectar(eje.salida, emb.entrada)
m.conectar(and_emb.salida, emb.control)
m.conectar(emb.salida, Puerto(suma_t, "ent", 0, 1))
m.conectar(tc.salida, fallo_carga.entrada)
m.conectar(filtro.salida, carga_tab.entrada)
m.conectar(fallo_carga.salida, gf.entrada)
m.conectar(uno.salida, Puerto(kf, "ent", 0, 1))
m.conectar(gf.salida, Puerto(kf, "ent", 1, 1))
m.conectar(carga_tab.salida, mux_c.entradas[0])
m.conectar(kf.salida, mux_c.entradas[1])
m.conectar(mux_c.salida, mult_c.entrada)
m.conectar(mult_c.salida, Puerto(suma_t, "ent", 1, 1))
m.conectar(suma_t.salida, inercia.entrada)
m.conectar(inercia.salida, filtro.entrada)
m.conectar(filtro.salida, eng.w1)
m.conectar(eje.salida, eng.T1)
m.conectar(filtro.salida, mux.entradas[1])

m.conectar(mot.salida, mux_p.entradas[0])
m.conectar(emb.salida, mux_p.entradas[1])
m.conectar(mux_p.salida, mult.entrada)
m.conectar(mult.salida, kperd.entrada)
m.conectar(kperd.salida, masa.entradas[0])
m.conectar(rt.salida, gr.entrada)
m.conectar(gr.salida, masa.entradas[1])
m.conectar(masa.salida, rt.T1)
m.conectar(amb.salida, rt.T2)

m.conectar(perfil.salida, sup.entradas[0])
m.conectar(masa.salida, sup.entradas[1])
m.conectar(masa.salida, reten.senal)
m.conectar(ref_activa.salida, reten.trigger)
m.conectar(perfil.salida, ref_activa.a)
m.conectar(cero.salida, ref_activa.b)
m.conectar(masa.salida, temp_ok.a)
m.conectar(ochenta5.salida, temp_ok.b)
m.conectar(mot.salida, en_marcha.a)
m.conectar(cero.salida, en_marcha.b)
m.conectar(ref_activa.salida, and_emb.entradas[0])
m.conectar(temp_ok.salida, and_emb.entradas[1])
m.conectar(mot.salida, mux.entradas[0])
m.conectar(masa.salida, mux.entradas[2])
m.conectar(mot.salida, tabla.entrada)

sc = m.add(Scope("scope", anchos=[3, 1, 1, 1],
                 guiones=["w1, w2 filtrada, T (C)", "estado supervisor",
                          "T_eje (N.m)", "w_carga (rad/s)"]))
m.conectar(mux.salida, sc.canales[0])
m.conectar(sup.salida, sc.canales[1])
m.conectar(eje.salida, sc.canales[2])
m.conectar(Puerto(eng, "sal", 0, 1), sc.canales[3])

res = m.run(t_fin=T_END, registrar=[
    sc, sup, masa, reten, tabla, filtro, eng, eje, fallo_t, perfil, mot,
])

t = res.t
est = np.asarray(res["sup"])
T = np.asarray(res["masa"])
print("Paso 3 — sistema completo: termica + supervisor + fallos")
print("=" * 66)
nombres = {0: "DETENIDO", 1: "ARRANQUE", 2: "FALLO"}
cambio = np.where(np.diff(est) != 0)[0]
for i in cambio:
    print("  t = %6.2f s : %9s -> %s" % (t[i], nombres[int(est[i])],
          nombres[int(est[i + 1])]))
print("  T maxima       : %6.1f C   (proteccion a 85 C)" % T.max())
print("  T final        : %6.1f C" % T[-1])
print("  T al arrancar  : %6.1f C  (RetenedorDisparado en ref > 0.5)" %
      float(res["reten"][-1]))
print("  eficiencia media: %5.1f %%  (Tabla1D vs w1)" %
      (float(np.mean(res["tabla"])) * 100))
print("  w2 filtrada    : %7.2f rad/s" % float(res["filtro"][-1]))
print("  w_carga final  : %7.2f rad/s  (Engranaje 2:1)" % res["eng"][-1, 0])
print("  par en el eje  : %7.2f N.m" % float(res["eje"][-1]))
