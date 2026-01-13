"""
f(t) =
    cos(10*pi*t) for 0s <= t < 0.3s
    cos(50*pi*t) for 0.3s <= t < 0.6s
    cos(100*pi*t) for 0.6s <= t <= 1s

Find FT and STFT of f(t)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal import spectrogram


def f(t):
    condition1 = (t >= 0) & (t < 0.3)
    condition2 = (t >= 0.3) & (t < 0.6)
    condition3 = (t >= 0.6) & (t <= 1)

    result = np.zeros_like(t)
    result[condition1] = np.cos(10 * np.pi * t[condition1])
    result[condition2] = np.cos(50 * np.pi * t[condition2])
    result[condition3] = np.cos(100 * np.pi * t[condition3])

    return result


t = np.linspace(0, 1, 1000)  # 1 second duration with 1000 samples

plt.figure(figsize=(12, 6))

# FT
ft = fft(f(t))

plt.subplot(2, 1, 1)
plt.plot(np.abs(ft))
plt.title('Fourier Transform of f(t)')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')

# STFT
frequencies, times, Sxx = spectrogram(f(t), fs=1000)  # fs is the sampling frequency, here assumed as 1000 Hz

plt.subplot(2, 1, 2)
plt.pcolormesh(times, frequencies, np.log(Sxx))
plt.title('Short-Time Fourier Transform of f(t)')
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.colorbar(label='Log Amplitude')
plt.tight_layout()

plt.show()
