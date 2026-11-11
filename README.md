# Financial-models
Monte carlo and black scholes simulator
1. How the simulator works
Simulate the stock. Under risk-neutral pricing the stock drifts at the risk-free rate r, not at its real expected return. One year ahead:
S_T = S₀ · exp[(r − σ²/2)T + σ√T · Z], with Z a standard normal draw. This is exact, so a European option needs only one step.
Compute the payoff on each path. A call pays max(S_T − K, 0). An Asian call pays max(average price − K, 0). A barrier call pays the vanilla payoff only if the barrier condition holds.
Discount and average. The price is e^(−rT) times the mean payoff across all paths.
Report the error. Standard error = (std of payoffs) / √N. The ± figures in the output are 1.96 standard errors, a 95% confidence interval.

A sanity check on the simulated paths: the mean terminal price came out at 105.116 against the theoretical 100·e^0.05 = 105.127, so the paths have the right drift.
2. European call and put
	MC (200k paths)	Black-Scholes
Call	10.4661 ± 0.0345	10.4506
Put	5.5891 ± 0.0345	5.5735

The call differs from the exact value by 0.0155, about half a standard error, so it is inside the confidence interval. The put misses by almost exactly the same amount, 0.0156, because both use the same random draws, so they share one sampling error and put-call parity still holds in the simulation.

About 56.0% of simulated paths finished above the strike, against the exact N(d2) = 55.96%. That confirms the model’s in-the-money probability directly.

3. Convergence and variance reduction

Standard error of the call price at different path counts:

Paths	Plain MC	Antithetic + control variate
1,000	9.61 (se 0.446)	9.68 (se 0.234)
10,000	10.30 (se 0.146)	10.43 (se 0.079)
100,000	10.36 (se 0.046)	10.41 (se 0.025)
1,000,000	10.43 (se 0.015)	10.43 (se 0.008)
Plain Monte Carlo error shrinks like 1/√N, so 10× the paths cuts the error by about 3.2×. The payoff standard deviation is 14.70, which gives 14.70/√1,000,000 = 0.0147.
Antithetic variates pair each draw Z with −Z. A path that ends high is paired with one that ends low, so their errors partly cancel.
Control variate uses a quantity whose exact mean is known: the discounted terminal stock price, whose true mean is S₀. The estimate is adjusted by β·(X − known mean), with β = Cov(Y,X)/Var(X) = 0.674. Payoff and stock price have a correlation of 0.924, so this alone removes about 85% of the variance.
found an inefficiency while checking this. Tested separately at 1M paths:

Setting	Std error
Plain	0.0147
Antithetic only	0.0104
Control variate only	0.0056
Both (the script’s default)	0.0078

Combining them is worse than the control variate alone for this option. The control already removes the linear part of the payoff, which is the part antithetic pairing helps with. What remains is the convex kink at the strike, where antithetic pairs are positively correlated, so pairing makes things slightly worse. The default is still much better than plain Monte Carlo, but for a European option the best setting is antithetic=False.

Two other things in the table. The estimates in a column all sit a little below 10.4506 because each row reuses the same seed, so larger runs contain the smaller runs’ draws and their errors are correlated. And the 1,000-path estimate of 9.68 is a 3-standard-error miss, which is just a small-sample effect: with so few paths, β and the standard error themselves are noisy.

4. Asian call: 5.7873 ± 0.0232

It is far cheaper than the European call (10.45) because averaging over 252 daily prices smooths out the randomness. The average is much less volatile than the final price, roughly σ/√3 ≈ 11.5%, and less volatility means less option value. Asian options have no closed form, which is why Monte Carlo is the natural tool. The control variate here is the vanilla call, which is strongly correlated with the Asian payoff, and it cut the standard error from 0.0254 to 0.0119.

5. Up-and-out barrier call (B = 130): 3.5728 ± 0.0352

The option dies if the stock touches 130 on any trading day. About 21.6% of paths hit that level, and those are the paths that would have paid the most. Removing them takes the value from 10.45 down to 3.57. The control variate helps little here (std error 0.0202 to 0.0179) because the knock-out payoff correlates weakly with the vanilla one. Monitoring only at daily closes misses intraday breaches, so this price is slightly higher than a continuously monitored barrier would give.

1. How the model prices the option

Black-Scholes assumes the stock follows geometric Brownian motion with constant volatility. Under that assumption the option’s fair price has a closed form. It needs two intermediate numbers:

d1 = [ln(S/K) + (r + σ²/2)T] / (σ√T) = (0 + 0.07) / 0.2 = 0.35
d2 = d1 − σ√T = 0.15

N(d1) = 0.6368 and N(d2) = 0.5596, where N is the standard normal CDF. Then:

Call = S·N(d1) − K·e^(−rT)·N(d2) = 63.68 − 95.12 × 0.5596 = 10.4506
Put = K·e^(−rT)·N(−d2) − S·N(−d1) = 5.5735

The call formula reads as the expected value of receiving the stock minus the discounted expected cost of paying the strike. N(d2) is the risk-neutral probability of finishing in the money, about 56% here. e^(−rT) discounts the strike back to today. The Monte Carlo script gave 10.466 ± 0.035, so the simulation agrees with the formula.

Put-call parity check: call − put = 4.8771, and S − K·e^(−rT) = 100 − 95.1229 = 4.8771. They match exactly, as arbitrage-free pricing requires.

2. The Greeks (sensitivities)

Each Greek is a derivative of the price formula with respect to one input.

Greek	Call	Put	Formula (call)	Meaning
Delta	0.6368	−0.3632	N(d1)	Price change per $1 rise in the stock
Gamma	0.0188	0.0188	φ(d1)/(Sσ√T)	How fast delta changes per $1 move
Vega	0.3752	0.3752	S·φ(d1)·√T	Price change per 1 point of volatility
Theta	−0.0176	−0.0045	−Sφ(d1)σ/(2√T) − rKe^(−rT)N(d2)	Price change per day, all else equal
Rho	0.5323	−0.4189	K·T·e^(−rT)·N(d2)	Price change per 1% rise in rates

φ is the normal density, and φ(0.35) = 0.3752.

Delta: the call gains about 64 cents per $1 the stock rises. Put delta equals call delta minus 1 (0.6368 − 1 = −0.3632), which also gives the probability-like reading that the put is in the money about 36% of the time.
Gamma and vega: these are the same for calls and puts, because parity ties the two prices together by a term that doesn’t depend on σ or on curvature in S.
Theta: the call loses about 1.76 cents a day, mostly the −3.75/365 decay in time value, partly offset by interest. The put decays more slowly because higher rates push put value down but time decay is offset by the strike’s discounting.
Rho: calls gain from higher rates, since the discounted strike gets cheaper. Puts lose.
3. Price vs spot (call)
Spot	80	90	100	110	120
Price	1.86	5.09	10.45	17.66	26.17
Delta	0.22	0.43	0.64	0.80	0.90

The price curve is convex, so delta rises as spot rises, climbing toward 1 as the option goes deep in the money. Gamma is that rate of change. It peaks near the money (about 0.022 around spot 90–95 here, because the forward price is above spot) and fades on both sides.

4. Price vs volatility (call)

The price rises almost linearly from 6.81 at 10% vol to 21.79 at 50% vol, about 0.38 per volatility point, which is vega. Vega itself peaks near 30% volatility (0.3794) and then eases. Gamma falls steadily as volatility rises, since higher volatility spreads the price distribution out and flattens delta.

5. Implied volatility

Implied vol runs the model backwards: given a market price, find the σ that makes the formula reproduce it. The formula can’t be solved for σ directly, so the code uses bisection between 0.0001% and 500%. Feeding in 10.4506 recovered 0.2000, the volatility I started with, which confirms the pricing and the solver agree.

A second run, a put with K = 105, T = 0.5 and a price of 8.7008, returned 0.2608. That price came from an earlier example with a 4% rate and 25% vol, while the command used the default 5% rate. The mismatch in rate is what pushed the recovered volatility away from 25%, so implied vol only means something if r is right.

6. Payoff table (put, K = 105, T = 0.5, σ = 25%)

The model premium is 8.396.

At expiry the put pays max(105 − S_T, 0), and profit is that payoff minus the 8.396 paid.
Breakeven is 105 − 8.396 ≈ 96.60. At 97.12 the profit is −0.52, just short of breakeven.
Above 105 the put expires worthless and the loss is capped at the premium, −8.40.
Below breakeven, profit grows about $1 per $1 drop in the stock, up to 31.50 − 8.40 = 23.10 at a stock price of 73.50.
