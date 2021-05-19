def u2f(u):
    e = ((u >> 23) & 0xFF)
    if e == 0: raise NotImplemented
    e -= 127
    s = -1 if (u & 0x80000000) else 1
    m = (u & 0x007FFFFF) | 0x00800000
    f = float(m) / 2 ** (23 - e)
    return f * s

