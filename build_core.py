"""Compila la DLL del núcleo C de bloques_crysi.

Uso:
    python build_core.py            # compila si la DLL no existe o está vieja
    python build_core.py --force    # recompila siempre
"""

import sys

from bloques_crysi import _clib


def main():
    fuerza = "--force" in sys.argv
    dll = _clib.compilar(fuerza=fuerza)
    print(f"DLL lista: {dll}")


if __name__ == "__main__":
    main()