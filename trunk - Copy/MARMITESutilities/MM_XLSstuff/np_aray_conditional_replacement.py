import timeit
import numpy as np

t = timeit.Timer('''
import numpy as np
data = np.clip(np.random.randn(2.5E6), -1, 1)
np.where(data <= 0.0, -10, data)
''')
print t.timeit(2)

t = timeit.Timer('''
import numpy as np
data = np.clip(np.random.randn(2.5E6), -1, 1)
for i, d in enumerate(data):
    if d <= 0.0:
       data[i] = -10
''')
print t.timeit(2)

