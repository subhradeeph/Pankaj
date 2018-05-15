'''
    Testing fftw wrapper for python (2.7). Needs the following packages:

        pyfftw (can be installed with pip in linux/mac/windows-bash, but
                not with anaconda in windows) (probably requires
                libfftw-dev as dependency)
        numpy (for other arrays)
        matplotlib (for plotting)
'''
import pyfftw as pfw
import numpy as np
import matplotlib.pylab as plt

nx = 128 # size of the array

# creating fftw-type array
a = pfw.empty_aligned(nx, dtype = 'complex128', n = 16)

# defining x-domain grid points
x = np.linspace(0, 2 * np.pi, nx)

# populating array with a cosine function
a[:] = np.cos(x)

# forward transform
b = pfw.interfaces.numpy_fft.fft(a)

# inverse transform and normalization
c = pfw.interfaces.numpy_fft.fft(b) / nx

# plotting

plt.figure(0)
plt.plot(x, a, x, c)
plt.figure(1)
plt.plot(x, b.real)
plt.show()
