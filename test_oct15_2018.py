import scipy.integrate as integrate
import scipy.special as special
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, dblquad
import numpy as np


result = integrate.quad(lambda x: special.jv(2.5, x), 0, 4.5)
print(result)
# (1.1178179380783249, 7.8663172481899801e-09)

I = sqrt(2/pi)*(18.0/27*sqrt(2)*cos(4.5) - 4.0/27*sqrt(2)*sin(4.5) + sqrt(2*pi) * special.fresnel(3/sqrt(pi))[0])
print(I)
# 1.117817938088701
print(abs(result[0]-I))

def integrand(x, a, b):
    return a*x**2 + b
a,b = 2,1
I = quad(integrand, 0, 1, args=(a,b))
print(I)

def integrand(t, n, x):
    return np.exp(-x*t) / t**n
def expint(n, x):
    return quad(integrand, 1, np.inf, args=(n, x))[0]
vec_expint = np.vectorize(expint)
vec_expint(3, np.arange(1.0, 4.0, 0.5))
special.expn(3, np.arange(1.0,4.0,0.5))

result = quad(lambda x: expint(3, x), 0, np.inf)
print(result)
I3 = 1.0/3.0
print(I3)
print(I3 - result[0])

def super_archimedean(ang):
    return (abs(np.cos(m * ang / 4) / a) ** n2 + abs(np.sin(m * ang / 4) / b) ** n3)**(-1 / n1)\
        * (Top_Curve_Ri + Top_Curve_Ratio * ang * 180 / np.pi)
def super_archimedean_deriv(ang):
    pass
def I(n):
    return dblquad(lambda t, x: np.exp(-x*t)/t**n, 0, np.inf, lambda x: 1, lambda x: np.inf)

print(I(4))

area = dblquad(lambda x, y: x*y, 0, 0.5, lambda x: 0, lambda x: 1-2*x)
print(area)

N = 5
def f(t, x):
   return np.exp(-x*t) / t**N

integrate.nquad(f, [[1, np.inf],[0, np.inf]])

def f(x, y):
    return x*y

def bounds_y():
    return [0, 0.5]

def bounds_x(y):
    return [0, 1-2*y]

integrate.nquad(f, [bounds_x, bounds_y])