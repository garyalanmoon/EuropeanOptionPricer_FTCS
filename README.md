This file contains Python code this is designed to price a European option under the standard Black-Scholes framework. These assumptions include:
**(1)** The underlying stock follows a geometric Brownian motion.

**(2)** Regarding the movement of the underlying asset, it is assumed that:

**(2a)** Volatility is constant over the life of the option (volatility is typically measured as standard deviation of returns on the underlying asset or something similar).

**(2b)** The risk-free rate is constant over the life of the option (in the risk-free world, this rate corresponds to the drift of the asset price).

**(3)** The underlying stock pays no dividends (this can be easily relaxed to include continuous dividend payments).

**(4)** All assets and investments are infinitely divisible. Thus, you can borrow or lend any fraction of a dollar you like. Likewise, you can buy or sell any fraction of the underlying stock you like.

**(5)** Markets are frictionless (no transaction costs).

Under these assumptions, the value of the option is governed by a backward parabolic PDE, the Black-Scholes equation. This program prices such an option by numerically solving the Black-Scholes equation.
The PDE is augmented with a terminal condition (it is backwards parabolic and so must be solved backwards in time) given by the payoff of the option at expiry.
Moreover, we have a boundary condition at S = 0: V(0,t) = 0 for a call and V(0,t) = K for a put (K being the strike price).
Using any numerical method requires working on a finite computational domain (as opposed to the half-line, the more natural setting for this problem). This introduces an artificial boundary.
We must also impose a BC at this artificial boundary, denoted by S_max (the maximum considered stock price). There are a number of natural choices for this BC.
We choose to impose a zero Gamma BC to approximate the far-field behavior of the option.
To solve the PDE, we use a forward time centered space (FTCS) finite difference scheme. These schemes are only conditionally stable and so one must choose appropriate time steps in relation to the
resolution of the grid in terms of the underlying asset price.
Further, this code prices the option using the closed-form Black-Scholes formula and returns the error incurred in pricing via the given numerical scheme.
