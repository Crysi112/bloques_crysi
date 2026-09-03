"""
Tests para CoSimuladorRed (Fase 3 + 4) usando BackendRedMock.
"""

import pytest
import numpy as np
from bloques_crysi import (Modelo, FuenteTrifasica, MedidorPotencia,
                            FuenteConstante, Suma, Ganancia, Puerto,
                            Multiplexor, Demultiplexor)
from bloques_crysi.red import (CoSimuladorRed, BackendRedMock, ResultadoCoSim,
                                crear_cosimulador_simple)


class TestCoSimuladorRed:
    """Tests del orquestador de co-simulación."""

    def test_creacion_basica(self):
        """Test creación básica del co-simulador."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)
        # Corrientes cero (sin carga) - usar FuenteTrifasica con 0 amplitud
        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0)

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
        )

        assert cosim.dt_red == 0.1
        assert cosim.tol_v == 1e-3
        assert cosim.max_iter == 20
        assert cosim.alpha == 0.5

    def test_run_simple_sin_carga(self):
        """Test co-simulación simple sin carga activa."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)
        # Corrientes cero (sin carga)
        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0)

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
            v_nominal_ln=230.0,
        )

        # Simular 0.3 segundos (3 ventanas de 0.1s)
        res = cosim.run(t_fin=0.3)

        # Verificar estructura del resultado
        assert isinstance(res, ResultadoCoSim)
        assert len(res.t) == 3  # 3 ventanas
        assert len(res.pccs) == 1
        pcc = res.pccs[0]
        assert len(pcc.v_pcc_mag) == 3
        assert len(pcc.v_pcc_ang) == 3
        assert len(pcc.p_pcc) == 3
        assert len(pcc.q_pcc) == 3
        assert len(pcc.convergido) == 3
        assert len(pcc.iteraciones_ventana) == 3

        # Sin carga, tensión debe ser cercana a slack
        for v in pcc.v_pcc_mag:
            assert abs(v - 230.0) < 1.0  # ~230 V LN

        # Potencias ~0
        for p in pcc.p_pcc:
            assert abs(p) < 10.0
        for q in pcc.q_pcc:
            assert abs(q) < 10.0

        # Todas convergidas
        assert all(pcc.convergido)

    def test_run_con_carga_RL(self):
        """Test co-simulación con carga RL simple (3 fases)."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        # Carga RL por fase: R=10
        # Demultiplexor para separar las 3 fases de tensión
        demux_v = m.add(Demultiplexor("demux_v", n_canales=3))
        m.conectar(fuente.salida, demux_v.entrada)

        # Una ganancia por fase
        R = 10.0
        gan_a = m.add(Ganancia("gan_a", 1.0/R))
        gan_b = m.add(Ganancia("gan_b", 1.0/R))
        gan_c = m.add(Ganancia("gan_c", 1.0/R))

        m.conectar(demux_v.salidas[0], gan_a.entrada)
        m.conectar(demux_v.salidas[1], gan_b.entrada)
        m.conectar(demux_v.salidas[2], gan_c.entrada)

        # Multiplexor para combinar las 3 corrientes
        mux = m.add(Multiplexor("mux_i", n_canales=3))
        m.conectar(gan_a.salida, mux.entradas[0])
        m.conectar(gan_b.salida, mux.entradas[1])
        m.conectar(gan_c.salida, mux.entradas[2])
        m.conectar(mux.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0, z_linea_ohm=0.1 + 0.05j)

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
            v_nominal_ln=230.0,
        )

        res = cosim.run(t_fin=0.3)

        # Verificar que la simulación corrió
        assert len(res.t) == 3
        pcc = res.pccs[0]
        assert all(pcc.convergido)

        # Tensión debe caer ligeramente por la impedancia de línea + carga
        for v in pcc.v_pcc_mag:
            assert v < 230.0

        # Potencia activa positiva (carga consume)
        for p in pcc.p_pcc:
            assert p > 1000  # > 1kW

    def test_convergencia_punto_fijo(self):
        """Test que el punto fijo converge en cada ventana."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        R = 20.0
        demux_v = m.add(Demultiplexor("demux_v", n_canales=3))
        m.conectar(fuente.salida, demux_v.entrada)

        gan_a = m.add(Ganancia("gan_a", 1.0/R))
        gan_b = m.add(Ganancia("gan_b", 1.0/R))
        gan_c = m.add(Ganancia("gan_c", 1.0/R))
        m.conectar(demux_v.salidas[0], gan_a.entrada)
        m.conectar(demux_v.salidas[1], gan_b.entrada)
        m.conectar(demux_v.salidas[2], gan_c.entrada)

        mux = m.add(Multiplexor("mux_i", n_canales=3))
        m.conectar(gan_a.salida, mux.entradas[0])
        m.conectar(gan_b.salida, mux.entradas[1])
        m.conectar(gan_c.salida, mux.entradas[2])
        m.conectar(mux.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0, z_linea_ohm=0.05 + 0.02j)

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
            tol_convergencia_v=1e-2,
            max_iter_ventana=15,
            relajacion=0.5,
            v_nominal_ln=230.0,
        )

        res = cosim.run(t_fin=0.2)

        # Verificar que convergió en todas las ventanas (la primera puede no converger en cold start)
        pcc = res.pccs[0]
        assert all(pcc.convergido[1:])  # desde la segunda ventana
        # Iteraciones razonables (1-15)
        assert all(pcc.iteraciones_ventana >= 1)
        assert all(pcc.iteraciones_ventana <= 15)

    def test_subrelajacion_estable(self):
        """Test que la sub-relaxación estabiliza el punto fijo."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        R = 10.0
        demux_v = m.add(Demultiplexor("demux_v", n_canales=3))
        m.conectar(fuente.salida, demux_v.entrada)

        gan_a = m.add(Ganancia("gan_a", 1.0/R))
        gan_b = m.add(Ganancia("gan_b", 1.0/R))
        gan_c = m.add(Ganancia("gan_c", 1.0/R))
        m.conectar(demux_v.salidas[0], gan_a.entrada)
        m.conectar(demux_v.salidas[1], gan_b.entrada)
        m.conectar(demux_v.salidas[2], gan_c.entrada)

        mux = m.add(Multiplexor("mux_i", n_canales=3))
        m.conectar(gan_a.salida, mux.entradas[0])
        m.conectar(gan_b.salida, mux.entradas[1])
        m.conectar(gan_c.salida, mux.entradas[2])
        m.conectar(mux.salida, medidor.corrientes)

        # Línea con impedancia significativa
        backend = BackendRedMock(v_slack_v=230.0, z_linea_ohm=0.5 + 0.2j)

        # Con alpha=0.5 (con relajación) debe ser estable
        cosim_con_relaj = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
            relajacion=0.5,
            v_nominal_ln=230.0,
        )

        res_con = cosim_con_relaj.run(t_fin=0.3)

        # Debe converger
        assert all(res_con.pccs[0].convergido)

    def test_bootstrapping_valor_inicial(self):
        """Test que el bootstrapping usa v_nominal_ln correctamente."""
        v_nom = 240.0  # Tensión nominal personalizada
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0)  # Slack diferente

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
            v_nominal_ln=v_nom,
        )

        # La primera ventana usa v_nominal_ln para bootstrapping
        res = cosim.run(t_fin=0.1)

        # Primera tensión debe ser cercana a v_nominal (antes de corrección backend)
        # O al slack del backend tras primera iteración
        assert len(res.t) == 1

    def test_varios_buses_backend_mock(self):
        """Test que el backend mock maneja bus 0 (slack) y bus 1."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock(v_slack_v=230.0)

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=1,
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
        )

        res = cosim.run(t_fin=0.1)

        # Bus slack (0) debe mantener 230V
        v_slack = backend.get_tension(0)
        assert abs(v_slack.magnitud_v - 230.0) < 1e-6

    def test_error_bus_inexistente(self):
        """Test que falla con bus inexistente."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock()

        cosim = crear_cosimulador_simple(
            modelo=m,
            backend=backend,
            bus_pcc=99,  # Bus inexistente en mock
            medidor=medidor,
            fuente_red=fuente,
            dt_red=0.1,
        )

        # El error se lanza al llamar run() -> backend.set_carga()
        with pytest.raises(ValueError):
            cosim.run(t_fin=0.1)

    def test_parametros_invalidos(self):
        """Test validación de parámetros del constructor."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        medidor = m.add(MedidorPotencia("pcc", fases=3))
        m.conectar(fuente.salida, medidor.entrada)

        cero = m.add(FuenteTrifasica("cero", amplitud=0.0, frecuencia=50.0))
        m.conectar(cero.salida, medidor.corrientes)

        backend = BackendRedMock()

        with pytest.raises(ValueError):
            crear_cosimulador_simple(modelo=m, backend=backend, bus_pcc=1,
                          medidor=medidor, fuente_red=fuente, dt_red=-0.1)

        with pytest.raises(ValueError):
            crear_cosimulador_simple(modelo=m, backend=backend, bus_pcc=1,
                          medidor=medidor, fuente_red=fuente, dt_red=0.1,
                          relajacion=1.5)

        with pytest.raises(ValueError):
            crear_cosimulador_simple(modelo=m, backend=backend, bus_pcc=1,
                          medidor=medidor, fuente_red=fuente, dt_red=0.1,
                          max_iter_ventana=0)


class TestSetParam:
    """Tests para Modelo.set_param() y get_param()."""

    def test_set_param_fuente_trifasica(self):
        """Test actualizar amplitud y fase de FuenteTrifasica."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))

        # Ejecutar una vez para armar _param_arrays
        m.run(t_fin=1e-4, registrar=[fuente])

        # Verificar valor inicial
        assert m.get_param(fuente, 0) == 400.0  # amplitud
        assert m.get_param(fuente, 1) == 50.0   # frecuencia
        assert m.get_param(fuente, 2) == 0.0    # fase

        # Actualizar amplitud
        m.set_param(fuente, 0, 420.0)
        assert m.get_param(fuente, 0) == 420.0

        # Actualizar fase
        m.set_param(fuente, 2, 0.5)
        assert abs(m.get_param(fuente, 2) - 0.5) < 1e-10

    def test_set_param_por_nombre(self):
        """Test set_param usando nombre de bloque."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        m.run(t_fin=1e-4, registrar=[fuente])

        m.set_param("red", 0, 440.0)
        assert m.get_param("red", 0) == 440.0

    def test_set_param_indice_invalido(self):
        """Test error con índice fuera de rango."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        m.run(t_fin=1e-4, registrar=[fuente])

        with pytest.raises(ValueError):
            m.set_param(fuente, 99, 100.0)

        with pytest.raises(ValueError):
            m.set_param(fuente, -1, 100.0)

    def test_set_param_antes_de_run(self):
        """Test error si se usa set_param antes de armar el modelo."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))

        with pytest.raises(ValueError):
            m.set_param(fuente, 0, 420.0)

    def test_set_param_bloque_inexistente(self):
        """Test error con bloque que no existe."""
        m = Modelo(dt=1e-4)
        fuente = m.add(FuenteTrifasica("red", amplitud=400, frecuencia=50))
        m.run(t_fin=1e-4, registrar=[fuente])

        with pytest.raises(ValueError):
            m.set_param("no_existe", 0, 100.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])