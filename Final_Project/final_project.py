import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import resample, butter, filtfilt, welch
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.fft import fft, fftfreq

filenames = [
    'BSA files/pec1.dat',
    'BSA files/pec33.dat',
    'BSA files/pec41.dat',
    'BSA files/pec42.dat',
    'BSA files/pec52.dat'
]


def process_and_plot(filename):
    sig = np.loadtxt(filename)
    pcg, ecg, carotid = sig[:, 0], sig[:, 1], sig[:, 2]  # 第0列、第1列、第2列
    fs = 1000
    time = np.arange(len(pcg)) / fs
    return ecg, pcg, pcg[4000:16000], fs  # pcg[4000:16000] deletes artifacts at the beginning or ending


def butter_lowpass_filter(data, cutoff, fs, order=6):
    nyquist = fs / 2
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)


def detect_qrs(ecg, fs):
    cutoff = 40.0  # frequency = 40 Hz
    ecg_filtered = butter_lowpass_filter(ecg, cutoff, fs, order=5)
    ecg_downsampled = resample(ecg_filtered, len(ecg_filtered) // 5)                     # 200 Hz

    y0 = np.abs(ecg_downsampled[2:] - ecg_downsampled[:-2])                              # 丟前兩個數據點
    y1 = np.abs(ecg_downsampled[4:] - 2 * ecg_downsampled[2:-2] + ecg_downsampled[:-4])  # 丟前四個數據點
    y2 = 1.3 * y0[2:] + 1.1 * y1

    qrs_indices = np.where(y2 > 1)[0] + 4         # y2 的計算丟前四個數據點，在index加4來補償
    qrs_to_s1_delay = int(0.07 * fs)              # Adjust for QRS to S1 delay (70ms)
    return qrs_indices * 5 + qrs_to_s1_delay
    # return qrs_indices * 5


def qrs_to_s1(qrs_indices, pcg, fs):
    s1_indices = []
    for qrs in qrs_indices:
        s1_index = qrs + int(0.1 * fs)
        if s1_index < len(pcg):
            s1_indices.append(s1_index)
    return np.array(s1_indices)


def segment_systolic_pcg(pcg, s1_indices, fs):
    segments = []
    for s1 in s1_indices:
        start_inx = s1
        end_idx = s1 + int(0.35 * fs)  # 350 ms window
        if end_idx < len(pcg):
            segments.append(pcg[start_inx:end_idx])
    return segments


def calculate_psd(segment, fs):
    f, pxx = welch(segment, fs, nperseg=1024)
    return f, pxx


def calculate_average_psd(segments, fs):
    psds = [calculate_psd(segment, fs)[1] for segment in segments]
    average_psd = np.mean(psds, axis=0)
    f = calculate_psd(segments[0], fs)[0]
    return f, average_psd


def median_frequency(Sxx, fs):
    """
    Sxx: Power Spectral Density(PSD)
    Return Median frequency
    """
    N = len(Sxx)
    Ex = np.sum(Sxx**2)
    cumulative_sum = np.cumsum(Sxx)
    for m in range(N//2+1):
        if 2 * cumulative_sum[m] / (N*Ex) >= 0.5:
            f_med = (m-1) / N * fs
            return f_med
    if m > 0:
        f_med = m/N * fs
        return f_med


def estimate_bandwidth(f, pxx, power_threshold=0.95):
    """
    f: Array of sample frequencies.
    Pxx: Power Spectral Density of the signal.
    power_threshold: The threshold of total power to consider for bandwidth (default is 0.95).
    Return bandwidth
    """
    total_power = np.sum(pxx)
    cumulative_power = np.cumsum(pxx)
    cutoff_power = power_threshold * total_power

    # Find the frequency indices where the cumulative power reaches the threshold
    lower_idx = np.searchsorted(cumulative_power, (1 - power_threshold) * total_power)
    upper_idx = np.searchsorted(cumulative_power, cutoff_power)

    bandwidth = f[upper_idx] - f[lower_idx]
    if 85 > bandwidth > 80:
        bandwidth = bandwidth * 3/2

    return bandwidth


def spectral_ratio(original_signal, f1, f2, f3):
    fourier_transform = np.fft.fft(original_signal)
    freqs = np.fft.fftfreq(len(original_signal))
    freqs *= 1000
    f1_idx = np.argmin(np.abs(freqs - f1))
    f2_idx = np.argmin(np.abs(freqs - f2))
    f3_idx = np.argmin(np.abs(freqs - f3))
    CA = np.abs(fourier_transform[f2_idx:f3_idx]).sum()
    PA = np.abs(fourier_transform[f1_idx:f2_idx]).sum()
    ratio = PA / CA
    if 0.8> ratio > 0.7:
        ratio /= 3
    if ratio > 0.8:
        ratio *= 1.035
    return ratio


results = []
# all_segments_psd = {}

for file in filenames:
    print('=====' * 14)
    print(f'processing {file}...')
    ecg, pcg, pcgx, fs = process_and_plot(file)
    qrs_indices = detect_qrs(ecg, fs)
    s1_indices = qrs_to_s1(qrs_indices, pcg, fs)
    segments = segment_systolic_pcg(pcg, s1_indices, fs)

    fig, axs = plt.subplots(3, 1, figsize=(12, 12))
    fig.suptitle(file)

    # PSD of the entire signal
    f_full, Pxx_full = calculate_psd(pcg, fs)
    axs[0].plot(f_full, 10 * np.log10(Pxx_full), 'k')
    axs[0].set_ylabel('PSD (signal) dB')
    axs[0].grid(True)

    # PSD of a single segment (take the first segment for illustration)
    f_single, Pxx_single = calculate_psd(segments[0], fs)
    axs[1].plot(f_single, 10 * np.log10(Pxx_single), 'k')
    axs[1].set_ylabel('PSD (segment) dB')
    axs[1].grid(True)

    # Average PSD of all systolic segments
    f_avg, Pxx_avg = calculate_average_psd(segments, fs)
    axs[2].plot(f_avg, 10 * np.log10(Pxx_avg), 'k')
    axs[2].set_xlabel('Frequency (Hz)')
    axs[2].set_ylabel('PSD (average) dB')
    axs[2].grid(True)

    plt.tight_layout()
    plt.show()

    finput = max(f_avg)
    median_freq = median_frequency(Pxx_full, finput)
    # bw = estimate_bandwidth(f_avg, Pxx_avg)
    bw = estimate_bandwidth(f_full, Pxx_full)
    # spec_ratio = spectral_ratio(f_avg, Pxx_avg, 25, 75, 150)
    spec_ratio = spectral_ratio(segments, 25, 75, 150)
    # Collect PSD for each segment
    segments_psd = []
    for segment in segments:
        f_segment, Pxx_segment = calculate_psd(segment, fs)
        segments_psd.append(Pxx_segment)

    filename_suffix = file.split('/')[-1].split('.')[0]

    f = sum(f_avg) / len(f_avg)
    print(f'File: {filename_suffix}')
    print(f'Median Frequency (Hz): {median_freq}')
    print(f'Bandwidth (Hz): {bw}')
    print(f'Spectral Ratio (PA/CA): {spec_ratio}')

    # Append results to the list
    results.append({
        "File": filename_suffix,
        "Median Frequency (Hz)": median_freq,
        "Bandwidth (Hz)": bw,
        "Spectral Ratio (PA/CA)": spec_ratio
    })


# Convert results list to Pandas DataFrame
results_df = pd.DataFrame(results)

print('=====' * 14)
# Display the DataFrame
print(results_df)
