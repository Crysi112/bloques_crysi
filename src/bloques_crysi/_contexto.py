import threading

_modelo_actual = threading.local()

def _get_modelo_actual():
    stack = getattr(_modelo_actual, "stack", None)
    if stack:
        return stack[-1]
    return None

def _set_modelo_actual(m):
    if not hasattr(_modelo_actual, "stack"):
        _modelo_actual.stack = []
    _modelo_actual.stack.append(m)

def _pop_modelo():
    stack = getattr(_modelo_actual, "stack", None)
    if stack:
        stack.pop()
