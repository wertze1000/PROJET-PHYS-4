import numpy as np
import matplotlib.pyplot as plt

def test(x):
    result = []
    for i in x:
        if i == 0:
            result.append(10)
        elif i == max(x):
            result.append(10)
            break
        else:
            result.append(0)
    return result

print(test(np.arange(0,10,0.1)))

def splitInterval(l, n):
    w = (l[1] - l[0]) // n
    return [[l[0]+i*w, l[0]+(i+1)*w] for i in range(n)]

x_min = 0.01
x_max = 2

xRange = np.linspace(x_min, x_max, 1000)
intervals = np.array_split(xRange, 2)
print(np.linspace(x_min, x_max, 1000), intervals[0])

plt.plot(np.arange(0,10,0.1), test(np.arange(0,10,0.1)))
plt.show()