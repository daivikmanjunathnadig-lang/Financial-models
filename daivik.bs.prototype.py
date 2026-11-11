import argparse
import math


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2))


def norm_pdf(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)


def d_values(S, K, T, r, sigma, q=0.0):
    d1 = (math.log(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    return d1, d1 - sigma * math.sqrt(T)


def price(S, K, T, r, sigma, q=0.0, kind="call"):
    if T <= 0:
        return max(S - K, 0) if kind == "call" else max(K - S, 0)
    d1, d2 = d_values(S, K, T, r, sigma, q)
    if kind == "call":
        return S * math.exp(-q * T) * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2)
    return K * math.exp(-r * T) * norm_cdf(-d2) - S * math.exp(-q * T) * norm_cdf(-d1)


def greeks(S, K, T, r, sigma, q=0.0, kind="call"):
    d1, d2 = d_values(S, K, T, r, sigma, q)
    eq = math.exp(-q * T)
    er = math.exp(-r * T)
    sqrt_t = math.sqrt(T)
    pdf = norm_pdf(d1)

    if kind == "call":
        delta = eq * norm_cdf(d1)
        theta = (-S * pdf * sigma * eq / (2 * sqrt_t)
                 - r * K * er * norm_cdf(d2) + q * S * eq * norm_cdf(d1))
        rho = K * T * er * norm_cdf(d2)
    else:
        delta = -eq * norm_cdf(-d1)
        theta = (-S * pdf * sigma * eq / (2 * sqrt_t)
                 + r * K * er * norm_cdf(-d2) - q * S * eq * norm_cdf(-d1))
        rho = -K * T * er * norm_cdf(-d2)

    return {
        "delta": delta,
        "gamma": eq * pdf / (S * sigma * sqrt_t),
        "vega": S * eq * pdf * sqrt_t * 0.01,
        "theta": theta / 365,
        "rho": rho * 0.01,
    }


def implied_vol(market_price, S, K, T, r, q=0.0, kind="call"):
    low, high = 1e-6, 5.0
    if market_price < price(S, K, T, r, low, q, kind) or market_price > price(S, K, T, r, high, q, kind):
        return None
    for _ in range(200):
        mid = 0.5 * (low + high)
        if price(S, K, T, r, mid, q, kind) > market_price:
            high = mid
        else:
            low = mid
        if high - low < 1e-10:
            break
    return 0.5 * (low + high)


def sweep(S, K, T, r, sigma, q, kind, what, lo, hi, points=9):
    rows = []
    for i in range(points):
        x = lo + (hi - lo) * i / (points - 1)
        args = dict(S=S, K=K, T=T, r=r, sigma=sigma)
        args[what] = x
        p = price(args["S"], args["K"], args["T"], args["r"], args["sigma"], q, kind)
        g = greeks(args["S"], args["K"], args["T"], args["r"], args["sigma"], q, kind)
        rows.append((x, p, g["delta"], g["gamma"], g["vega"]))
    return rows


def print_sweep(rows, label):
    print(f"{label:>10} {'price':>10} {'delta':>9} {'gamma':>9} {'vega':>9}")
    for x, p, d, g, v in rows:
        print(f"{x:>10.4f} {p:>10.4f} {d:>9.4f} {g:>9.4f} {v:>9.4f}")


def payoff_table(S, K, T, r, sigma, q, kind):
    premium = price(S, K, T, r, sigma, q, kind)
    print(f"{'spot at expiry':>16} {'payoff':>10} {'profit':>10}")
    for i in range(9):
        st = K * (0.7 + 0.075 * i)
        pay = max(st - K, 0) if kind == "call" else max(K - st, 0)
        print(f"{st:>16.2f} {pay:>10.2f} {pay - premium:>10.2f}")


def parity_check(S, K, T, r, sigma, q):
    c = price(S, K, T, r, sigma, q, "call")
    p = price(S, K, T, r, sigma, q, "put")
    lhs = c - p
    rhs = S * math.exp(-q * T) - K * math.exp(-r * T)
    print(f"call - put           {lhs:.6f}")
    print(f"S*e^-qT - K*e^-rT    {rhs:.6f}")


def demo():
    S, K, T, r, sigma = 100, 100, 1.0, 0.05, 0.2
    for kind in ("call", "put"):
        print(f"{kind}: {price(S, K, T, r, sigma, kind=kind):.4f}")
        for name, val in greeks(S, K, T, r, sigma, kind=kind).items():
            print(f"  {name:>6}: {val:.4f}")
    print()
    parity_check(S, K, T, r, sigma, 0.0)
    print()
    print("price vs spot (call)")
    print_sweep(sweep(S, K, T, r, sigma, 0.0, "call", "S", 80, 120), "spot")
    print()
    print("price vs volatility (call)")
    print_sweep(sweep(S, K, T, r, sigma, 0.0, "call", "sigma", 0.1, 0.5), "vol")
    print()
    iv = implied_vol(10.4506, S, K, T, r)
    print(f"implied vol for a call priced at 10.4506: {iv:.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--S0", type=float)
    ap.add_argument("--K", type=float)
    ap.add_argument("--T", type=float)
    ap.add_argument("--r", type=float, default=0.05)
    ap.add_argument("--sigma", type=float, default=0.2)
    ap.add_argument("--q", type=float, default=0.0)
    ap.add_argument("--type", choices=["call", "put"], default="call")
    ap.add_argument("--iv", type=float)
    ap.add_argument("--sweep", choices=["S", "sigma", "T", "r"])
    ap.add_argument("--lo", type=float)
    ap.add_argument("--hi", type=float)
    ap.add_argument("--payoff", action="store_true")
    a = ap.parse_args()

    if a.S0 is None or a.K is None or a.T is None:
        demo()
        return

    S, K, T, r, sigma, q, kind = a.S0, a.K, a.T, a.r, a.sigma, a.q, a.type

    if a.iv is not None:
        vol = implied_vol(a.iv, S, K, T, r, q, kind)
        if vol is None:
            print("no implied vol exists for that price")
        else:
            print(f"implied vol: {vol:.4f}")
        return

    print(f"{kind} price: {price(S, K, T, r, sigma, q, kind):.4f}")
    for name, val in greeks(S, K, T, r, sigma, q, kind).items():
        print(f"{name:>8}: {val:.4f}")

    if a.sweep:
        if a.lo is None or a.hi is None:
            ap.error("--sweep needs --lo and --hi")
        print()
        print_sweep(sweep(S, K, T, r, sigma, q, kind, a.sweep, a.lo, a.hi), a.sweep)

    if a.payoff:
        print()
        payoff_table(S, K, T, r, sigma, q, kind)


if __name__ == "__main__":
    main()
