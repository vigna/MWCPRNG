#!/usr/bin/env python3
"""Independent model of the (G)MWC generators with the scrambler X: checks that
the moduli are safe primes, and that the outputs and the results of jump(),
long_jump(), jump_ce(), and jump_n() printed by test.c are correct (usage:
check.py DIR, where DIR contains the outputs of check.sh)."""
import sys, random


def is_prime(n, rounds=32):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    rnd = random.Random(42)
    for _ in range(rounds):
        x = pow(rnd.randrange(2, n - 1), d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

GENS = {  # name: (digit bits, lag, multiplier, e or None, jump and long-jump exponents)
    "MWC128X": (64, 1, 0xff6da236274e252e, None, 64, 96),
    "MWC192X": (64, 2, 0xff4e24dd418f5433, None, 96, 144),
    "MWC256X": (64, 3, 0xff32e7e2bdcb0a33, None, 128, 192),
    "GMWC128X": (64, 1, 0xfff0cfb4ed6366fa, 0x68935f57, 64, 96),
    "GMWC256X": (64, 3, 0xff76e600f25ae046, 0x34dabaf938e45b, 128, 192),
    "MWC64X_32": (32, 1, 0xfff145f8, None, 32, 48),
    "MWC128X_32": (32, 3, 0xffcccc66, None, 64, 96),
}
allok = True
for name, (w, k, a, e, je, le) in GENS.items():
    b = 1 << w; M = b - 1
    m = a * b**k - 1 if e is None else a * b**k + e
    A = pow(b, -1, m)                     # per-step multiplier of the associated MCG
    prime = is_prime(m) and is_prime((m - 1) // 2)
    digits = [v & M for v in [0x0123456789abcdef, 0xfedcba9876543210, 0x0f1e2d3c4b5a6978][:k]]
    c0 = 1 if e is None else 2 * e       # seeds used by test.c
    A0INV = None if e is None else pow(-e % b, -1, b)
    def step(st):
        d, c = st
        t = a * d[0] + c
        if e is None: xn, cn = t & M, t >> w
        else: xn = (A0INV * t) & M; cn = (t + e * xn) >> w
        return (d[1:] + [xn], cn)
    def out(before, after):
        d, c = before
        if e is None: return d[-1] ^ c   # X on the current state
        t = a * d[0] + c                 # X: next digit ^ upper half of a x_{n-k} + c_{n-1}
        return after[0][-1] ^ (t >> w)
    def h(st):                           # state integer
        d, c = st
        X = sum(v << (w * i) for i, v in enumerate(d))
        return (c * b**k + X) % m if e is None else (c * b**k - e * X) % m
    def decode(s):
        if e is None:
            return ([(s >> (w * i)) & M for i in range(k)], s >> (w * k))
        X = s * pow(-e % b**k, -1, b**k) % b**k
        return ([(X >> (w * i)) & M for i in range(k)], (s + e * X) >> (w * k))
    def jump(st, n): return decode(h(st) * pow(A, n, m) % m)
    st = (digits, c0)
    s2 = st
    for _ in range(1000): s2 = step(s2)
    assert jump(st, 1000) == s2, name       # the model agrees with stepping
    outs, s2 = [], st
    for _ in range(8):
        s3 = step(s2); outs.append(out(s2, s3)); s2 = s3
    fmt = lambda tag, s: tag + " " + " ".join("%0*x" % (w // 4, v) for v in s[0] + [s[1]])
    n = 0x9e3779b97f4a7c15 + ((0x0123456789abcdef << 64) if name != "MWC64X_32" else 0)  # jump_n() argument in test.c
    expected = ["out %0*x" % (w // 4, v) for v in outs] + [
        fmt("jump", jump(st, 2**je)), fmt("long_jump", jump(st, 2**le)),
        fmt("jump_ce_3_70", jump(st, 3 << 70)), fmt("jump_ce_12345_0", jump(st, 12345)),
        fmt("jump_n", jump(st, n))]
    got = open(f"{sys.argv[1]}/{name}.out").read().split("\n")[:len(expected)]
    ok = got == expected
    allok &= ok and prime
    print(f"{name}: m and (m-1)/2 prime: {prime}; outputs and jumps match the model: {ok}")
print("ALL OK" if allok else "FAILURE")
sys.exit(0 if allok else 1)
