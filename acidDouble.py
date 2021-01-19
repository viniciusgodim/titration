import sympy
from sympy import symbols,Eq,nsolve,solve,re
import numpy
from numpy import log10

Ka1 = 1.75e-05
Ka2 = 0
Kw = 1e-14

cNaOH0 = 0.1
cHA10 = 0.1
cHA20 = 0
VHA0 = 25/1000
cA10 = 0
cA20 = 0

VNaOH0_range = numpy.linspace(0,50/1000,100)

pH_list = []
VNaOH0_list = []

for VNaOH0 in VNaOH0_range:
    z = cNaOH0*VNaOH0/(VNaOH0+VHA0)
    x1,x2,y = symbols('x1,x2,y')

    cHA1 = cHA10*VHA0/(VNaOH0+VHA0)
    cHA2 = cHA20*VHA0/(VNaOH0+VHA0)
    HA1 = cHA1 - x1
    HA2 = cHA2 - x2
    A1 = x1
    A2 = x2
    H3O = x1 + x2 + y
    OH = z + y
    Na = z

    acid1 = Eq(Ka1,H3O*A1/HA1)
    acid2 = Eq(Ka2,H3O*A2/HA2)
    water = Eq(Kw,H3O*OH)

    yf1 = solve(acid1,y)[0]
    yf2 = solve(acid2,y)[0]
    x1f = solve(Eq(yf1,yf2),x1)[0]

    x2sol = solve(water.subs(y,yf1).subs(x1,x1f),x2)

    t = 0

    for x2i in x2sol:
        x1i = re(x1f.subs(x2,x2i))
        yi = re(yf1.subs(x1,x1i).subs(x2,x2i))
        HA1i = re(cHA1 - x1i)
        HA2i = re(cHA2 - x2i)
        A1i = re(x1i)
        A2i = re(x2i)
        H3Oi = re(x1i + x2i + yi)
        OHi = re(z + yi)
        if HA1i > 0 and HA2i > 0 and A1i > 0 and A2i > 0 and H3Oi > 0 and OHi > 0:
            HA1 = HA1i
            HA2 = HA2i
            A1 = A1i
            A2 = A2i
            H3O = H3Oi
            OH = OHi
            t = t + 1
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
