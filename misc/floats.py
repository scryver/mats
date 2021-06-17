def u2f32(u):
    e = ((u >> 23) & 0xFF)
    if e == 0: raise NotImplemented
    e -= 127
    s = -1 if (u & 0x80000000) else 1
    m = (u & 0x007FFFFF) | 0x00800000
    f = float(m) / 2 ** (23 - e)
    return f * s

def u2f64(u):
    e = ((u >> 52) & 0x3FF)
    if e == 0: raise NotImplemented
    e -= 1023
    s = -1 if (u & 0x8000000000000000) else 1
    m = (u & 0x000FFFFFFFFFFFFF) | 0x0010000000000000
    f = float(m) / 2 ** (52 - e)
    return f * s

