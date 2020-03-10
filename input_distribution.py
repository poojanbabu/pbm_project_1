import random
import matplotlib.pyplot as plt
import collections
import math
import numpy as np
from scipy.interpolate import make_interp_spline

alpha = 0.6
beta = 0.3

iterations = 1000
timesteps = 1000
freq1 = {}
freq2 = {}
dt = 0.1
for j in range(iterations):
    X = 10
    for i in range(timesteps):
        if X == 0:
            beta_dt = 0
            alpha_dt = alpha
        if X > 0:
            stop = 1/X
            beta_dt = beta * X * dt
            alpha_dt = alpha * dt
        p = round(random.uniform(0, 1), 2)
        if p <= alpha_dt:
            print("Birth")
            X = X + 1
        elif alpha_dt < p <= (alpha_dt + beta_dt):
            print("Death")
            X = X - 1
        else:
            print("Do nothing: Stay")
        if timesteps-200 < i < timesteps-100:
            if X in freq1:
                freq1[X] += 1
            else:
                freq1[X] = 1
        elif timesteps-100 < i < timesteps:
            if X in freq2:
                freq2[X] += 1
            else:
                freq2[X] = 1
freq_ord1 = collections.OrderedDict(sorted(freq1.items()))
freq_ord2 = collections.OrderedDict(sorted(freq2.items()))
max_frequency = max(freq2.values())
for key in freq_ord2:
    freq_ord2[key] /= max_frequency
print(freq_ord2)

n = 0
e = 2.71828

mu = alpha/beta

y = []
x = []

while n <= 10:

    C = 3.5 * (((mu**n)*(e**-mu))/(math.factorial(n)))

    y.append(C)
    x.append(n)

    n += 1

ax = np.array(x)
ay = np.array(y)

# Spline Just for Visual Guide
xnew = np.linspace(ax.min(), ax.max(), 1000)
spl = make_interp_spline(ax, ay, k=3)
y_smooth = spl(xnew)

# Plot on Graph
plt.style.use('seaborn-darkgrid')
plt.title('Distribution of input molecules')
plt.xlabel('X')
plt.ylabel('Probability')
plt.bar(freq_ord2.keys(), freq_ord2.values())
plt.scatter(x, y)
plt.plot(xnew, y_smooth, 'red')
plt.show()
