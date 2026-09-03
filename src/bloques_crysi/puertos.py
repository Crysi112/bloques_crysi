from typing import Optional, Sequence, Union
from ._contexto import _get_modelo_actual
class Puerto:
    def __init__(self, bloque, tipo: str, offset: int, n: int,
                 canales: Optional[Sequence[str]] = None):
        self.bloque = bloque
        self.tipo = tipo
        self.offset = offset
        self.n = n
        self.canales = canales
    def indices(self):
        arr = self.bloque.in_idx if self.tipo == "ent" else self.bloque.out_idx
        return arr[self.offset:self.offset + self.n]
    def __len__(self):
        return self.n
    def __getitem__(self, item):
        if isinstance(item, int):
            if item < 0:
                item += self.n
            if not (0 <= item < self.n):
                raise IndexError("Indice de canal fuera de rango.")
            sub_canales = [self.canales[item]] if self.canales else None
            return Puerto(self.bloque, self.tipo, self.offset + item, 1,
                          canales=sub_canales)
        if isinstance(item, slice):
            start, stop, step = item.indices(self.n)
            if step != 1:
                raise ValueError("Solo se admiten rebanados continuos (step=1).")
            sub_n = max(0, stop - start)
            sub_canales = self.canales[start:stop] if self.canales else None
            return Puerto(self.bloque, self.tipo, self.offset + start, sub_n,
                          canales=sub_canales)
        raise TypeError("Indice debe ser entero o slice.")
    def _obtener_modelo(self):
        from ._contexto import _get_modelo_actual
        m = _get_modelo_actual()
        if m is None:
            raise RuntimeError(
                "Operadores algebraicos requieren contexto 'with modelo:'. "
                "Ejemplo: with modelo: error = v_ref - v_med"
            )
        return m
    def _crear_suma(self, otro, signo_otro: float = 1.0):
        m = self._obtener_modelo()
        if not isinstance(otro, Puerto):
            val = float(otro) * signo_otro
            from .bloques import FuenteConstante, Suma
            c_nombre = f"_c_{abs(val)}"
            c = next((b for b in m.bloques if b.nombre == c_nombre), None)
            if c is None:
                c = m.add(FuenteConstante(c_nombre, abs(val)))
            s = m.add(Suma(f"_sum_c_{len(m.bloques)}", signos=[1.0, 1.0 if val >= 0 else -1.0]))
            m.conectar(self, s.entrada[0])
            m.conectar(c.salida, s.entrada[1])
            return s.salida
        if self.n != otro.n:
            raise ValueError(f"Puertos con distinto tamaño: {self.n} != {otro.n}")
        from .bloques import Suma
        s = m.add(Suma(f"_sum_{len(m.bloques)}", signos=[1.0, signo_otro]))
        m.conectar(self, s.entrada[0])
        m.conectar(otro, s.entrada[1])
        return s.salida
    def __add__(self, otro):
        return self._crear_suma(otro, 1.0)
    def __radd__(self, otro):
        return self._crear_suma(otro, 1.0)
    def __sub__(self, otro):
        return self._crear_suma(otro, -1.0)
    def __rsub__(self, otro):
        if not isinstance(otro, Puerto):
            val = float(otro)
            m = self._obtener_modelo()
            from .bloques import FuenteConstante, Suma
            c_nombre = f"_c_{val}"
            c = next((b for b in m.bloques if b.nombre == c_nombre), None)
            if c is None:
                c = m.add(FuenteConstante(c_nombre, val))
            s = m.add(Suma(f"_rsub_{len(m.bloques)}", signos=[1.0, -1.0]))
            m.conectar(c.salida, s.entrada[0])
            m.conectar(self, s.entrada[1])
            return s.salida
        return otro._crear_suma(self, -1.0)
    def __mul__(self, otro):
        return self._crear_mul(otro)
    def __rmul__(self, otro):
        return self._crear_mul(otro)
    def _crear_mul(self, otro):
        m = self._obtener_modelo()
        if isinstance(otro, (int, float)):
            k = float(otro)
            from .bloques import Ganancia
            g = m.add(Ganancia(f"_k_{k}_{len(m.bloques)}", k))
            m.conectar(self, g.entrada)
            return g.salida
        if not isinstance(otro, Puerto):
            return NotImplemented
        if self.n != otro.n:
            raise ValueError(f"Puertos con distinto tamaño: {self.n} != {otro.n}")
        from .bloques import Multiplicador
        mult = m.add(Multiplicador(f"_mult_{len(m.bloques)}"))
        m.conectar(self, mult.entrada[0])
        m.conectar(otro, mult.entrada[1])
        return mult.salida
    def __truediv__(self, otro):
        m = self._obtener_modelo()
        if isinstance(otro, (int, float)):
            k = float(otro)
            if abs(k) < 1e-12:
                raise ZeroDivisionError("División por cero")
            from .bloques import Ganancia
            g = m.add(Ganancia(f"_div_{k}_{len(m.bloques)}", 1.0 / k))
            m.conectar(self, g.entrada)
            return g.salida
        return NotImplemented
    def __rtruediv__(self, otro):
        return NotImplemented
    def __neg__(self):
        m = self._obtener_modelo()
        from .bloques import Ganancia
        g = m.add(Ganancia(f"_neg_{len(m.bloques)}", -1.0))
        m.conectar(self, g.entrada)
        return g.salida
    def __repr__(self):
        return f"<Puerto {self.tipo}@{self.bloque.nombre}[{self.offset}:{self.offset + self.n}]>"
class Sensor:
    def __init__(self, nombre: str, bloque, tipo: str, offset: int, n: int,
                 canales: Optional[Sequence[str]] = None):
        self.nombre = nombre
        self.bloque = bloque
        self.tipo = tipo
        self.offset = offset
        self.n = n
        self.canales = canales or [f"{nombre}[{i}]" for i in range(n)]
    @property
    def puerto(self) -> Puerto:
        return Puerto(self.bloque, self.tipo, self.offset, self.n,
                      canales=self.canales)
    def indices(self):
        return self.puerto.indices()
    _indices = indices
    def __len__(self):
        return self.n
    def __getitem__(self, item):
        return self.puerto[item]
    def _obtener_modelo(self):
        return self.puerto._obtener_modelo()
    def _crear_suma(self, otro, signo_otro: float = 1.0):
        return self.puerto._crear_suma(otro, signo_otro)
    def __add__(self, otro):
        return self.puerto.__add__(otro)
    def __radd__(self, otro):
        return self.puerto.__radd__(otro)
    def __sub__(self, otro):
        return self.puerto.__sub__(otro)
    def __rsub__(self, otro):
        return self.puerto.__rsub__(otro)
    def __mul__(self, otro):
        return self.puerto.__mul__(otro)
    def __rmul__(self, otro):
        return self.puerto.__rmul__(otro)
    def __truediv__(self, otro):
        return self.puerto.__truediv__(otro)
    def __rtruediv__(self, otro):
        return self.puerto.__rtruediv__(otro)
    def __neg__(self):
        return self.puerto.__neg__()
    def __repr__(self):
        return f"<Sensor {self.nombre!r} {self.n} canales>"
