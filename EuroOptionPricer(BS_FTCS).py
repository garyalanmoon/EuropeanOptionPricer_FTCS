######## This program prices a European option by an explicit (FTCS) finite difference scheme to numerically solve the Black-Scholes equation.
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import os

####Option and Market Parameters
r = 0.05 #risk-free rate of interest (annualized)
sig = 0.2 #market volatility (on annualized basis)
K = 100 #strike price for option
T = 1 #time to expiry (in years)
q = -1 #option type parameter (use q=1 for call and q=-1 for a put)


####Grid Parameters
S_min = 0 #minimum stock price
S_max = 2*K #maximum stock price (approximates infinity on truncated domain)
NS = 200 #number of grid points
dS = (S_max - S_min)/(NS-1) #"spatial" step size
Nt = 8000 #number of time steps (would improve to automatically choose Nt or dt based on dS)
dt = T/Nt #time step size

####Check for Stability
lam = 0.5 * sig**2 * S_max**2 * dt/(dS**2)
print(f'The stability parameter lambda is {lam:.4f}.')
if lam > 0.5:
    response = input('Warning: The FTCS scheme may be unstable (lambda > 0.5). Continue? (y/n): ')
    if response.lower() != 'y':
        print('Terminating program due to stability concerns.')
        os._exit(0)

####Generate the Grid
S = np.linspace(S_min, S_max, NS)
t = np.linspace(0,T,Nt+1)

####Initialize the Value Function V = V(S,t), Impose the Terminal Condition and Impose BCs at time t = T
####The truncation of the domain introduces an artificial boundary. We want this artificial boundary to represent S = +infinity.
####Thus, this will impose a far-field-type BC at the artificial boundary. Some common choices include
####Dirichlet: V(S,t) = S or V(S,t) = S - K*exp(-r(T-t)) as S -> +infinity
####Neumann: Delta = 1 as S -> +infinity
####Our choice: Gamma = 0 as S -> +infinity
V = np.zeros((Nt+1,NS))
V[0,:] = np.maximum( q * ( S - K ), 0 ) #Terminal condition: payoff of the option

##BCs
if q == 1: #Call option
    V[0,0] = 0 #BC at S = 0 (imposed at terminal time): V(0,T) = 0
    V[0,-1] = 2* V[0,-2] - V[0,-3] #Far-field BC (imposed at terminal time): Gamma = 0
else: #put option
    V[0,0] = K #BC at S = 0 (imposed at terminal time): V(0,T) = K
    V[0,-1] = 0 #Far-field BC (imposed at terminal time): V(S,T) = 0 for S >> K.


####March forward (backward, technically) through time via FTCS scheme
####An important note is that Black-Scholes is a backward parabolic equation. If we try to solve forward in time, the output will appear unstable.
####This is not a result of numerical instability. Rather, compare with trying to solve the standard heat equation backward in time.
####The fix is that we essentially change t to T - t. This will turn it to a forward parabolic equation that we can march forward through time.
for n in range(Nt):
    for j in range(1,NS-1):
        delta = (V[n,j+1] - V[n,j-1])/(2*dS) #Approximate first derivative of V wrt S via central difference
        gamma = (V[n,j+1] - 2*V[n,j] + V[n,j-1])/(dS**2) #Approximate the second derivative of V wrt S via central difference

        ##Update interior V using Black-Scholes FTCS
        a = 0.5 * sig**2 * S[j]**2
        b = r * S[j]
        c = -r*V[n,j]
        V[n+1,j] = V[n,j] + dt * ( a * gamma + b * delta + c )

    ##Enforce BCs
    if q == 1: #Call
        V[n+1,0] = 0 #V(0,t) = 0
        V[n+1, -1] = 2 * V[n+1,-2] - V[n+1,-3] #Zero Gamma far-field BC
    else: #Put
        V[n+1,0] = K * np.exp(-r*(T-t[n+1])) #V(0,t) = K * exp(-r*(T-t))
        V[n+1,-1] = 0 #Put is worthless in the far field S >> K


# Interpolate the option price at S = S_0 or any other given price.
S_0 = K #Stock price today (stock price at which we wish to price the option)
V_0 = np.interp(S_0, S, V[-1, :])

print(f'Option Value at S = {S_0}, t = 0: {V_0:.4f}.')

####Price the Option via BS Formula
d_1 = ( np.log(S_0/K) + ( r + 0.5 * sig**2 ) * T )/( sig * np.sqrt(T) )
d_2 = d_1 - sig*np.sqrt(T)
V_BS = q * ( S_0 * norm.cdf(q*d_1) - K * np.exp(-r * T) * norm.cdf(q*d_2) )

####Check Error
error = np.abs( V_0 - V_BS )
print(f'The error is: {error:.5f}.')


####Visualize the Value Function
plt.plot(S, V[-1,:], label="Option Price at t=0")
plt.axvline(K, linestyle="--", color="r", label="Strike Price")
plt.xlabel("Stock Price S")
plt.ylabel("Option Value V(S,0)")
plt.title("European Option Price via FTCS")
plt.legend()
plt.grid()
plt.show()
