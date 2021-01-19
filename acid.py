import sympy
from sympy import symbols,Eq,nsolve,solve,re
import numpy
from numpy import log10

Ka = 2E-05
Kw = 1e-14

cNaOH0 = 0.1
cHA0 = 0.1
cA0 = 0.05

VHA0 = 25/1000

VNaOH0_range = numpy.linspace(0,50/1000,100)

pH_list = []
VNaOH0_list = []

for VNaOH0 in VNaOH0_range:
    z = cNaOH0*VNaOH0/(VNaOH0+VHA0)
    x,y = symbols('x,y')

    cHA = cHA0*VHA0/(VNaOH0+VHA0)
    HA = cHA - x
    A = x + cA0
    H3O = x + y
    OH = z + y
    Na = z

    acid = Eq(Ka,H3O*A/HA)
    water = Eq(Kw,H3O*OH)

    yf = solve(acid,y)[0]
    xsol = solve(water.subs(y,yf),x)
    ysol = [yf.subs(x,xi) for xi in xsol]

    t = 0

    for xi in xsol:
        for yi in ysol:
            HAi = re(cHA - xi)
            Ai = re(xi)
            H3Oi = re(xi + yi)
            OHi = re(z + yi)
            if HAi > 0 and Ai > 0 and H3Oi > 0 and OHi > 0:
                HA = HAi
                A = Ai
                H3O = H3Oi
                OH = OHi
                t = 1
    print(t)
    #solution = nsolve((water,acid),(x,y),((Ka*cHA)**.5,-1e-7))

    #HA = HA.subs(x,solution[0]).subs(y,solution[1])
    #A = A.subs(x,solution[0]).subs(y,solution[1])
    #H3O = H3O.subs(x,solution[0]).subs(y,solution[1])
    #OH = OH.subs(x,solution[0]).subs(y,solution[1])

    pH = -log10(float(H3O))

    VNaOH0_list.append(VNaOH0)
    pH_list.append(pH)

import matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(VNaOH0_list,pH_list)
fig.savefig("test.png")
