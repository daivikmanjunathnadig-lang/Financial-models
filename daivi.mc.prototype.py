import argparse
import math
import numpy as np


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2))


def black_scholes(S, K, T, r, sigma, q=0.0, kind="call"):
    d1 = (math.log(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if kind == "call":
        return S * math.exp(-q * T) * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2)
    return K * math.exp(-r * T) * norm_cdf(-d2) - S * math.exp(-q * T) * norm_cdf(-d1)


def make_paths(S, T, r, sigma, q, z):
    steps = z.shape[1]
    dt = T / steps
    incr = (r - q - 0.5 * sigma ** 2) * dt + sigma * math.sqrt(dt) * z
    return S * np.exp(np.cumsum(incr, axis=1))


def vanilla(x, K, kind):
    if kind == "call":
        return np.maximum(x - K, 0)
    return np.maximum(K - x, 0)


def payoff(paths, K, kind, style, barrier, btype):
    final = paths[:, -1]
    if style == "european":
        return vanilla(final, K, kind)
    if style == "asian":
        return vanilla(paths.mean(axis=1), K, kind)
    if style == "barrier":
        if btype.startswith("up"):
            hit = paths.max(axis=1) >= barrier
        else:
            hit = paths.min(axis=1) <= barrier
        alive = ~hit if btype.endswith("out") else hit
        return vanilla(final, K, kind) * alive
    raise ValueError("unknown style: " + style)


def draw_normals(rng, n, steps, antithetic):
    if antithetic:
        z = rng.standard_normal(((n + 1) // 2, steps))
        return np.vstack([z, -z])[:n]
    return rng.standard_normal((n, steps))


def price(S, K, T, r, sigma, q=0.0, kind="call", style="european",
          barrier=None, btype="up-and-out", n=100000, steps=None,
          antithetic=True, control=True, seed=42):
    if steps is None:
        steps = 1 if style == "european" else 252
    rng = np.random.default_rng(seed)
    z = draw_normals(rng, n, steps, antithetic)
    paths = make_paths(S, T, r, sigma, q, z)
    disc = math.exp(-r * T)
    y = disc * payoff(paths, K, kind, style, barrier, btype)

    if control:
        if style == "european":
            x = disc * paths[:, -1]
            x_mean = S * math.exp(-q * T)
        else:
            x = disc * vanilla(paths[:, -1], K, kind)
            x_mean = black_scholes(S, K, T, r, sigma, q, kind)
        if x.var() > 0:
            beta = np.cov(y, x)[0, 1] / x.var(ddof=1)
            y = y - beta * (x - x_mean)

    if antithetic:
        half = n // 2
        y = 0.5 * (y[:half] + y[half:2 * half])
    return y.mean(), y.std(ddof=1) / math.sqrt(len(y))


def greeks(S, K, T, r, sigma, q=0.0, kind="call", style="european",
           barrier=None, btype="up-and-out", n=200000, steps=None, seed=42):
    if steps is None:
        steps = 1 if style == "european" else 252
    rng = np.random.default_rng(seed)
    z = draw_normals(rng, n, steps, True)

    def p(S=S, T=T, r=r, sigma=sigma):
        paths = make_paths(S, T, r, sigma, q, z)
        return math.exp(-r * T) * payoff(paths, K, kind, style, barrier, btype).mean()

    dS = 0.01 * S
    dv = 0.001
    dr = 0.0001
    dT = min(1 / 365, T / 10)

    mid = p()
    up = p(S=S + dS)
    down = p(S=S - dS)
    return {
        "price": mid,
        "delta": (up - down) / (2 * dS),
        "gamma": (up - 2 * mid + down) / dS ** 2,
        "vega": (p(sigma=sigma + dv) - p(sigma=sigma - dv)) / (2 * dv) * 0.01,
        "rho": (p(r=r + dr) - p(r=r - dr)) / (2 * dr) * 0.01,
        "theta": (p(T=T - dT) - mid) / (dT * 365),
    }


def show(label, res):
    mean, se = res
    print(f"{label:<28}{mean:9.4f}  +/- {1.96 * se:.4f}")


def demo():
    S, K, T, r, sigma = 100, 100, 1.0, 0.05, 0.2

    for kind in ("call", "put"):
        show("European " + kind, price(S, K, T, r, sigma, kind=kind, n=200000))
        print(f"{'  black-scholes':<28}{black_scholes(S, K, T, r, sigma, kind=kind):9.4f}")

    print()
    print(f"{'paths':>10}  {'plain':>20}  {'antithetic + control':>22}")
    for n in (1000, 10000, 100000, 1000000):
        a = price(S, K, T, r, sigma, n=n, antithetic=False, control=False, seed=1)
        b = price(S, K, T, r, sigma, n=n, seed=1)
        print(f"{n:>10,}  {a[0]:9.4f} (se {a[1]:.4f})  {b[0]:9.4f} (se {b[1]:.4f})")

    print()
    show("Asian call", price(S, K, T, r, sigma, style="asian"))
    show("Up-and-out call (B=130)", price(S, K, T, r, sigma, style="barrier", barrier=130))

    print()
    for name, val in greeks(S, K, T, r, sigma).items():
        print(f"{name:>8}: {val:.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--S0", type=float)
    ap.add_argument("--K", type=float)
    ap.add_argument("--T", type=float)
    ap.add_argument("--r", type=float, default=0.05)
    ap.add_argument("--sigma", type=float, default=0.2)
    ap.add_argument("--q", type=float, default=0.0)
    ap.add_argument("--type", choices=["call", "put"], default="call")
    ap.add_argument("--style", choices=["european", "asian", "barrier"], default="european")
    ap.add_argument("--barrier", type=float)
    ap.add_argument("--barrier-type", default="up-and-out",
                    choices=["up-and-out", "down-and-out", "up-and-in", "down-and-in"])
    ap.add_argument("--paths", type=int, default=100000)
    ap.add_argument("--steps", type=int)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--greeks", action="store_true")
    a = ap.parse_args()

    if a.S0 is None or a.K is None or a.T is None:
        demo()
        return

    if a.style == "barrier" and a.barrier is None:
        ap.error("--barrier is required for barrier options")

    args = dict(S=a.S0, K=a.K, T=a.T, r=a.r, sigma=a.sigma, q=a.q, kind=a.type,
                style=a.style, barrier=a.barrier, btype=a.barrier_type)
    show("MC price", price(**args, n=a.paths, steps=a.steps, seed=a.seed))
    if a.style == "european":
        print(f"{'Black-Scholes':<28}{black_scholes(a.S0, a.K, a.T, a.r, a.sigma, a.q, a.type):9.4f}")
    if a.greeks:
        for name, val in greeks(**args, n=a.paths, steps=a.steps, seed=a.seed).items():
            print(f"{name:>8}: {val:.4f}")


if __name__ == "__main__":
    main()
