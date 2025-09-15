import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("data.csv")

N = 10000
n = 27
dt = 0.001

time = np.array([i*dt for i in range(N)])
X = df['X'].to_numpy()
Y = df['Y'].to_numpy()
Z = df['Z'].to_numpy()
msd = np.array([0.0]*N)
for i in range(N):
    msd_temp = 0.0
    for j in range(n):
        msd_temp += (X[j+i*n]-X[j])**2+(Y[j+i*n]-Y[j])**2+(Z[j+i*n]-Z[j])**2
    msd[i] = msd_temp/n
plt.plot(time, msd)
print(msd[-1]/time[-1]/6)
plt.show()