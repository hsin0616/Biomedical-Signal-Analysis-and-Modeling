"""
f(t) = 3*cos(4*pi*t) + 5*sin(12*pi*t) + 6*cos(44*pi*t) + 2*sin(52*pi*t)
g(t) = f(2*t - 0.1)

Find Fourier transform (FT) of f(t) and g(t)
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq


def f(t):
    return 3 * np.cos(4 * np.pi * t) + 5 * np.sin(12 * np.pi * t) + 6 * np.cos(44 * np.pi * t) + 2 * np.sin(52 * np.pi * t)


def g(t):
    return f(2 * t - 0.1)


N = 512
t = np.linspace(0, 1, N, endpoint=False)
dt = t[1] - t[0]

plt.figure(figsize=(12, 8))

# FT of f(t)
f_transform = fft(f(t))
f_freq = fftfreq(N, dt)

plt.subplot(2, 1, 1)
plt.plot(f_freq, np.abs(f_transform), label='FT of f(t)', color='blue')
plt.title('Fourier Transform of f(t)')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')

# FT of g(t)
g_transform = fft(g(t))
g_freq = fftfreq(N, dt)

plt.subplot(2, 1, 2)
plt.plot(g_freq, np.abs(g_transform), label='FT of g(t)', color='red')
plt.title('Fourier Transform of g(t)')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')

plt.tight_layout()
plt.show()
