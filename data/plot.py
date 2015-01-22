import matplotlib.pyplot as pl
import numpy as np

x, y, z = np.loadtxt("solution.txt").transpose()
temp = z.reshape((101,101))

pl.imshow(temp, vmin = 0, vmax = 1)
pl.colorbar()
pl.show()