import uncertainties as uc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.odr import ODR, Model, RealData
from uncertainties import unumpy
from uncertainties.umath import *
def linear_model(params, x):
    k, b = params
    return k * x + b
n = 125
df1 = pd.read_csv("eng.csv")
df2 = pd.read_csv("vel.csv")
df3 = pd.read_csv("msd.csv")
Vx = df2['vx'].to_numpy()
Vy = df2['vy'].to_numpy()
Vz = df2['vz'].to_numpy()
px = np.array([0.0]*(len(Vx)//n))
for i in range(0, len(px)):
    px[i] = sum(Vx[n*i:n*(i+1)])
V = np.sqrt(Vx**2+Vy**2+Vz**2)
msd = df3['msd'].to_numpy()
msdN  = df3['N'].to_numpy()
F = df1["F"].to_numpy()
K = df1["K"].to_numpy()
P = df1["P"].to_numpy()
T = K/n*2/3
data = RealData(msdN*0.001, msd)
model = Model(linear_model)
odr = ODR(data, model, beta0=[1.0, 0.0])
output = odr.run()
k = output.beta[0]
b = output.beta[1]
plt.plot(msdN*0.001, msd)
plt.plot(msdN*0.001, msdN*0.001*k+b)
print("В = ", k/6)
plt.show()
def func1(a, m):
    d = round(np.sqrt(m))
    b = np.array([0.0]*d)
    m1 = min(a)
    m2 = max(a)
    e = (m2 - m1) / d
    S = 0
    for i in range(0, d):
        for j in range(0, len(a)):
            if(m1 + e * i <= a[j] <= m1 + e * (i+1)):
                b[i] += 1
    for i in range(0, d):
        S += b[i] * e
    return b / S
def func2(a, m):
    d = round(np.sqrt(m))
    b = np.array([0.0]*d)
    m1 = min(a)
    m2 = max(a)
    e = (m2 - m1) / d
    for i in range(0, d):
        b[i] = m1 + e * (i+1)
    return b
def func3(a, m):
    d = round(np.sqrt(m))
    m1 = min(a)
    m2 = max(a)
    e = (m2 - m1) / d
    return e
def func4(v, t):
    return np.exp(-v*v/t/2)/np.sqrt(2*np.pi*t)
def func5(v, t):
    return 4*np.pi*v*v*np.exp(-v*v/t/2)/np.pow(2*np.pi*t, 1.5)
def func6(f, t):
    return -2*t*np.log(f*np.sqrt(2*np.pi*t))
def func7(f, t, v):
    return -2*t*np.log(f*np.pow(2*np.pi*t, 1.5)/4/np.pi/(v*v))
plt.plot(F)
plt.plot(K)
plt.plot(P)
plt.show()

plt.plot(T)
plt.plot([np.average(T)]*len(T))
plt.show()

b = func1(V[len(V)-n:len(V)], n)
nb = func2(V[len(V)-n:len(V)], n)
wb = func3(V[len(V)-n:len(V)], n)
b1 = func1(Vx[len(Vx)-n:len(Vx)], n)
nb1 = func2(Vx[len(Vx)-n:len(Vx)], n)
wb1 = func3(Vx[len(Vx)-n:len(Vx)], n)
b2 = func1(Vy[len(Vy)-n:len(Vy)], n)
nb2 = func2(Vy[len(Vy)-n:len(Vy)], n)
wb2 = func3(Vy[len(Vy)-n:len(Vy)], n)
b3 = func1(Vz[len(Vz)-n:len(Vz)], n)
nb3 = func2(Vz[len(Vz)-n:len(Vz)], n)
wb3 = func3(Vz[len(Vz)-n:len(Vz)], n)


plt.bar(nb, b, align = "edge", width = -wb)
plt.plot(nb, func5(nb, T[-1]))
plt.show()
plt.bar(nb**2, func7(b, T[-1], nb), align = "edge")
plt.plot(nb**2, nb**2)
plt.show()
plt.bar(nb1, b1, align = "edge", width = -wb1)
plt.plot(nb1, func4(nb1, T[-1]))
plt.show()
plt.bar(nb1**2, func6(b1, T[-1]), align = "edge")
plt.plot(nb1**2, nb1**2)
plt.show()
plt.bar(nb2, b2, align = "edge", width = -wb2)
plt.plot(nb2, func4(nb2, T[-1]))
plt.show()
plt.bar(nb2**2, func6(b2, T[-1]), align = "edge")
plt.plot(nb2**2, nb2**2)
plt.show()
plt.bar(nb3, b3, align = "edge", width = -wb3)
plt.plot(nb3, func4(nb3, T[-1]))
plt.show()
plt.bar(nb3**2, func6(b3, T[-1]), align = "edge")
plt.plot(nb3**2, nb3**2)
plt.show()

plt.plot(px)
plt.show()


