"""
Plot 1-level Haar transforms of the following functions-sampled uniformly over [0,1) using 1024 points.
(a) f(x) = (x^2)(1-x)
(b) f(x) = (x^4)[(1-x)^6][cos (64𝜋𝑥)]
(c) f(x) =
        −3 for 0.4 < x < 0.5
        2 for 0.5 ≤ x < 0.8
        0 for others
(d) f(x) = sgn(sin12*pi*x)
"""


import numpy as np
import matplotlib.pyplot as plt


def haar_transform(y):
    n = len(y)
    scaling = [(y[2*i] + y[2*i+1]) / np.sqrt(2) for i in range(n // 2)]
    wavelet = [(y[2*i] - y[2*i+1]) / np.sqrt(2) for i in range(n // 2)]
    return scaling, wavelet


plt.figure(figsize=(14, 11))


# (a)
def f_a(x):
    return x**2 * (1 - x)


x_a_values = np.linspace(0, 1, 1024)
y_a_values = f_a(x_a_values)
scaling_a, wavelet_a = haar_transform(y_a_values)

plt.subplot2grid((2, 2), (0, 0))
plt.plot(scaling_a, color='green', label='Scaling Coefficients')
plt.plot(wavelet_a, color='red', label='Wavelet Coefficients')
plt.title('1-level Haar Transform of (a)')
plt.xlabel('Index')
plt.ylabel('Coefficient Value')
plt.legend()
plt.grid(True)


# (b)
def f_b(x):
    return (x**4)*((1 - x)**6)*(np.cos(64*np.pi*x))


x_b_values = np.linspace(0, 1, 1024)
y_b_values = f_b(x_b_values)
scaling_b, wavelet_b = haar_transform(y_b_values)

plt.subplot2grid((2, 2), (0, 1))
plt.plot(scaling_b, color='green', label='Scaling Coefficients')
plt.plot(wavelet_b, color='red', label='Wavelet Coefficients')
plt.title('1-level Haar Transform of (b)')
plt.xlabel('Index')
plt.ylabel('Coefficient Value')
plt.legend()
plt.grid(True)


# (c)
def f_c(x):
    if 0.4 < x < 0.5:
        return -3
    elif 0.5 <= x < 0.8:
        return 2
    else:
        return 0


f_c_vectorized = np.vectorize(f_c)

x_c_values = np.linspace(0, 1, 1024)
y_c_values = f_c_vectorized(x_c_values)
scaling_c, wavelet_c = haar_transform(y_c_values)

plt.subplot2grid((2, 2), (1, 0))
plt.plot(scaling_c, color='green', label='Scaling Coefficients')
plt.plot(wavelet_c, color='red', label='Wavelet Coefficients')
plt.title('1-level Haar Transform of (c)')
plt.xlabel('Index')
plt.ylabel('Coefficient Value')
plt.legend()
plt.grid(True)


# (d)
def f_d(x):
    return np.sign(np.sin(12*np.pi*x))


x_d_values = np.linspace(0, 1, 1024)
y_d_values = f_d(x_d_values)
scaling_d, wavelet_d = haar_transform(y_d_values)

plt.subplot2grid((2, 2), (1, 1))
plt.plot(scaling_d, color='green', label='Scaling Coefficients')
plt.plot(wavelet_d, color='red', label='Wavelet Coefficients')
plt.title('1-level Haar Transform of (d)')
plt.xlabel('Index')
plt.ylabel('Coefficient Value')
plt.legend()
plt.grid(True)


plt.show()


