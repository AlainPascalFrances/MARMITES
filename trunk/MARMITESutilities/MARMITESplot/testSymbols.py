import matplotlib.pyplot as plt
import matplotlib.dates as mdates

x = list(range(1,100,1))
y = list(range(1,200,2))[0:99]

#plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
plt.plot(x, y, 'bo')
plt.plot(y, 'r+')
plt.plot(x, y, 'r-')
plt.plot(y, x, 'b--')
plt.plot(y, x, 'bs')
plt.plot(y, x, 'g^')

plt.show()