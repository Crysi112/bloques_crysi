from bloques_crysi import Modelo, CargaPQTrifasica
from bloques_crysi.red import RedOpenDSS

red = RedOpenDSS(nombre="IEEE13", v_slack_kv_ll=4.16, f_hz=60.0)
red.dss.text("Set VoltageBases=[4.16, 0.48]")

red.definir_linecode("601", 3,
    r_mat=[[0.3465], [0.1560, 0.3375], [0.1580, 0.1535, 0.3414]],
    x_mat=[[1.0179], [0.5017, 1.0478], [0.4236, 0.3849, 1.0348]],
    b_mat_us=[[6.2998], [-1.9958, 5.9597], [-1.2595, -0.7417, 5.6386]]
)
red.definir_linecode("602", 3,
    r_mat=[[0.7526], [0.1580, 0.7475], [0.1560, 0.1535, 0.7436]],
    x_mat=[[1.1814], [0.4236, 1.1983], [0.5017, 0.3849, 1.2112]],
    b_mat_us=[[5.6990], [-1.0817, 5.1795], [-1.6905, -0.6588, 5.4246]]
)
red.definir_linecode("603", 2,
    r_mat=[[1.3294], [0.2066, 1.3238]],
    x_mat=[[1.3471], [0.4591, 1.3569]],
    b_mat_us=[[4.7097], [-0.8999, 4.6658]]
)
red.definir_linecode("604", 2,
    r_mat=[[1.3238], [0.2066, 1.3294]],
    x_mat=[[1.3569], [0.4591, 1.3471]],
    b_mat_us=[[4.6658], [-0.8999, 4.7097]]
)
red.definir_linecode("605", 1,
    r_mat=[[1.3292]], x_mat=[[1.3475]], b_mat_us=[[4.5193]]
)
red.definir_linecode("606", 3,
    r_mat=[[0.7982], [0.3192, 0.7891], [0.2849, 0.3192, 0.7982]],
    x_mat=[[0.4463], [0.0328, 0.4041], [-0.0143, 0.0328, 0.4463]],
    b_mat_us=[[96.8897], [0.0, 96.8897], [0.0, 0.0, 96.8897]]
)
red.definir_linecode("607", 1,
    r_mat=[[1.3425]], x_mat=[[0.5124]], b_mat_us=[[88.9912]]
)

red.agregar_banco_reguladores("RG60", bus_hv="650", bus_lv="RG60",
                               kva_fase=5000/3, kv_ll=4.16,
                               v_set=122.0, pt_ratio=20.0, ct_rating=700.0,
                               band=2.0, r_set=3.0, x_set=9.0,
                               xhl=0.01, r_perc=0.01)

red.agregar_transformador("XFM1", bus_hv="633", bus_lv="634", kvs=[4.16, 0.48], kvas=[500, 500], conns=["delta", "wye"], xhl=2.0, r_perc=1.1)
red.agregar_capacitor("Cap675", bus="675", kvar=600, kv=4.16, fases=3)
red.agregar_capacitor("Cap611", bus="611.3", kvar=100, kv=2.4018, fases=1)
red.agregar_capacitor("Cap675", bus="675", kvar=600, kv=4.16, fases=3)
red.agregar_capacitor("Cap611", bus="611.3", kvar=100, kv=2.4018, fases=1)

red.agregar_linea("L_RG60_632", "RG60.2.1.3", "632.2.1.3", longitud_ft=2000, linecode="601", fases=3)
red.agregar_linea("L_632_645",  "632.3.2", "645.3.2", longitud_ft=500, linecode="603", fases=2)
red.agregar_linea("L_645_646",  "645.3.2", "646.3.2", longitud_ft=300, linecode="603", fases=2)
red.agregar_linea("L_632_633",  "632.3.1.2", "633.3.1.2", longitud_ft=500, linecode="602", fases=3)
red.agregar_linea("L_632_671",  "632.2.1.3", "671.2.1.3", longitud_ft=2000, linecode="601", fases=3)
red.agregar_linea("L_671_680",  "671.2.1.3", "680.2.1.3", longitud_ft=1000, linecode="601", fases=3)
red.agregar_linea("L_671_684",  "671.1.3", "684.1.3", longitud_ft=300, linecode="604", fases=2)
red.agregar_linea("L_684_611",  "684.3", "611.3", longitud_ft=300, linecode="605", fases=1)
red.agregar_linea("L_684_652",  "684.1", "652.1", longitud_ft=800, linecode="607", fases=1)
red.agregar_linea("Sw_671_692", "671.2.1.3", "692.2.1.3", longitud_ft=1,    linecode="601", fases=3)
red.agregar_linea("L_692_675",  "692", "675", longitud_ft=500, linecode="606", fases=3)

red.agregar_carga("L_634_A", "634.1", kw=160, kvar=110, kv=0.277, fases=1, model=1)
red.agregar_carga("L_634_B", "634.2", kw=120, kvar=90,  kv=0.277, fases=1, model=1)
red.agregar_carga("L_634_C", "634.3", kw=120, kvar=90,  kv=0.277, fases=1, model=1)
red.agregar_carga("L_645_B", "645.2", kw=170, kvar=125, kv=2.402, fases=1, model=1)
red.agregar_carga("L_646_BC", "646.2.3", kw=230, kvar=132, kv=4.16, fases=1, model=2)
red.agregar_carga("L_652_A", "652.1", kw=128, kvar=86, kv=2.402, fases=1, model=2)
red.agregar_carga("L_671_AB", "671.1.2", kw=385, kvar=220, kv=4.16, fases=1, model=1)
red.agregar_carga("L_671_BC", "671.2.3", kw=385, kvar=220, kv=4.16, fases=1, model=1)
red.agregar_carga("L_671_CA", "671.3.1", kw=385, kvar=220, kv=4.16, fases=1, model=1)
red.agregar_carga("L_675_A", "675.1", kw=485, kvar=190, kv=2.402, fases=1, model=1)
red.agregar_carga("L_675_B", "675.2", kw=68,  kvar=60,  kv=2.402, fases=1, model=1)
red.agregar_carga("L_675_C", "675.3", kw=290, kvar=212, kv=2.402, fases=1, model=1)
red.agregar_carga("L_692_CA", "692.3.1", kw=170, kvar=151, kv=4.16, fases=1, model=5)
red.agregar_carga("L_611_C",  "611.3",   kw=170, kvar=80,  kv=2.402, fases=1, model=5)
red.agregar_carga("L_dist_A", "671.1", kw=17,  kvar=10, kv=2.402, fases=1, model=1)
red.agregar_carga("L_dist_B", "671.2", kw=66,  kvar=38, kv=2.402, fases=1, model=1)
red.agregar_carga("L_dist_C", "671.3", kw=117, kvar=68, kv=2.402, fases=1, model=1)

backend = red.compilar()
backend.dss.text("RegControl.RG60_1.Enabled=No")
backend.dss.text("RegControl.RG60_2.Enabled=No")
backend.dss.text("RegControl.RG60_3.Enabled=No")
backend.dss.text("Transformer.RG60_1.Taps=[1.0, 1.0625]")
backend.dss.text("Transformer.RG60_2.Taps=[1.0, 1.05]")
backend.dss.text("Transformer.RG60_3.Taps=[1.0, 1.06875]")
backend.dss.text("Solve")
for f in (1, 2, 3):
    print(f"RG60_{f}:", backend.dss.text(f"? Transformer.RG60_{f}.Taps"))
print("XFM1 Conns:", backend.dss.text("? Transformer.XFM1.Conns"))
print("XFM1 %Rs:", backend.dss.text("? Transformer.XFM1.%Rs"))
print("XFM1 buses:", backend.dss.text("? Transformer.XFM1.Buses"))
from bloques_crysi import Modelo, CargaPQTrifasica

m = Modelo(dt=1e-5)
carga_634 = m.add(CargaPQTrifasica("Carga_634_EMT", p_w=400e3, q_var=290e3))

cosim = m.acoplar_red(
    backend=backend,
    bus_pcc=["634"],
    elemento=[carga_634],
    v_nominal_ll=[480.0],
    dt_red=0.03333,
    tol_convergencia_v=1.0,
    max_iter_ventana=50,
    relajacion=0.3,
    reemplazar_cargas=True
)

res = cosim.run(t_fin=0.5)
cosim.reporte()
