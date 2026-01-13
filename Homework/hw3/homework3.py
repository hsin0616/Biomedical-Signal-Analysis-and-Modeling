import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter


# Pan-Tompkins algorithm
def pan_tompkins(ecg_signal, fs):
    # Step 1: Bandpass filter
    b_low = [1, -2, 1, 0, 0, 0, -2, 1, 0, 0, 0, 1, -2, 1]
    a_low = [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    b_high = [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]
    a_high = [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

    filtered_signal_low = lfilter(b_low, a_low, ecg_signal)
    filtered_signal = lfilter(b_high, a_high, filtered_signal_low)

    # Step 2: Differentiator
    diff_signal = np.zeros_like(filtered_signal)
    diff_signal[0] = (2 * filtered_signal[0] + filtered_signal[0] - filtered_signal[2] - 2 * filtered_signal[3]) / 8
    for i in range(1, len(diff_signal)):
        diff_signal[i] = (2 * filtered_signal[i] + filtered_signal[i - 1] - filtered_signal[i - 3] - 2 *
                          filtered_signal[i - 4]) / 8

    # Step 3: Squaring operation
    squared_signal = diff_signal ** 2

    # Step 4: Moving-window integrator
    window_width = int(0.150 * fs)  # width of moving window, about 150ms
    mov_avg_signal = np.convolve(squared_signal, np.ones(window_width) / window_width, mode='same')

    # Find QRS peak positions
    peak_indices = np.where((mov_avg_signal[1:-1] > mov_avg_signal[:-2]) & (mov_avg_signal[1:-1] > mov_avg_signal[2:]))[
                       0] + 1

    detected_beats = len(peak_indices)
    expected_beats = len(ecg_signal) // fs  # Expected number of beats based on signal duration and sampling rate
    missed_beats = detected_beats - expected_beats
    return peak_indices, missed_beats


# Function to calculate activity, mobility, and form factor (FF)
def calculate_activity_mobility_ff(signal):
    # Calculate activity (variance of the signal segment)
    activity = np.var(signal)

    # Calculate first derivative and its activity
    first_derivative = np.gradient(signal)
    first_derivative_activity = np.var(first_derivative)

    # Calculate second derivative and its activity
    second_derivative = np.gradient(first_derivative)
    second_derivative_activity = np.var(second_derivative)

    # Calculate mobility
    mobility = np.sqrt(first_derivative_activity / activity)

    # Calculate form factor (FF)
    FF = np.sqrt(second_derivative_activity / first_derivative_activity)

    return activity, mobility, FF


with open('ecgpvc.dat', 'r') as file:
    ecg_signal = np.array([float(line.strip()) for line in file.readlines()])

fs = 2048  # sampling rate

# Detected QRS peaks using Pan-Tompkins
qrs_indices, missed_beats = pan_tompkins(ecg_signal, fs)

# Calculate RR intervals
RR_intervals = np.diff(qrs_indices)

# Normal beats and PVCs
normal_indices = [i for i, rr in enumerate(RR_intervals) if rr < 0.92 * np.mean(RR_intervals)]
PVC_indices = [i for i, rr in enumerate(RR_intervals) if rr > 1.08 * np.mean(RR_intervals)]

# Mean and standard deviation of RR intervals for normal beats and PVCs
mean_RR_normal = np.mean(RR_intervals[normal_indices])
mean_RR_PVC = np.mean(RR_intervals[PVC_indices])
std_RR_normal = np.std(RR_intervals[normal_indices])
std_RR_PVC = np.std(RR_intervals[PVC_indices])

# Calculate FF for each QRS complex
FF_values = []
for qrs_index in qrs_indices:
    qrs_t_part = ecg_signal[qrs_index:qrs_index + 80]  # QRS-T, 80 samples
    activity, mobility, FF = calculate_activity_mobility_ff(qrs_t_part)
    FF_values.append(FF)

# Mean and standard deviation of FF values for normal beats and PVCs
mean_FF_normal = np.mean([FF_values[i] for i in normal_indices])
mean_FF_PVC = np.mean([FF_values[i] for i in PVC_indices])
std_FF_normal = np.std([FF_values[i] for i in normal_indices])
std_FF_PVC = np.std([FF_values[i] for i in PVC_indices])

# Variation of RR intervals and FF values between normal and PVC beats
RR_mean_difference = abs(mean_RR_normal - mean_RR_PVC)
RR_std_difference = abs(std_RR_PVC - std_RR_normal)
FF_mean_difference = abs(mean_FF_normal - mean_FF_PVC)
FF_std_difference = abs(std_FF_normal - std_FF_PVC)

print('===================================================')
print("Number of beats missed by detection procedure:", missed_beats)
print('===================================================')
print("Mean RR interval for normal heartbeats:", mean_RR_normal)
print("Standard deviation of RR interval for normal heartbeats:", std_RR_normal)
print("Mean RR interval for PVCs:", mean_RR_PVC)
print("Standard deviation of RR interval for PVCs:", std_RR_PVC)
print('===================================================')
print("Mean FF for normal heartbeats:", mean_FF_normal)
print("Standard deviation of FF for normal heartbeats:", std_FF_normal)
print("Mean FF for PVCs:", mean_FF_PVC)
print("Standard deviation of FF for PVCs:", std_FF_PVC)
print('===================================================')
print("Difference in mean RR interval between normal and PVCs:", RR_mean_difference)
print("Difference in standard deviation of RR interval between normal and PVCs:", RR_std_difference)
print("Difference in mean FF between normal and PVCs:", FF_mean_difference)
print("Difference in standard deviation of FF between normal and PVCs:", FF_std_difference)


# Split QRS indices into groups of 10 heartbeats
num_periods = len(qrs_indices) // 10
grouped_qrs_indices = [qrs_indices[i * 10:(i + 1) * 10] for i in range(num_periods)]
# Plot each group of 10 heartbeats in a new figure
for i, group in enumerate(grouped_qrs_indices):
    plt.figure(figsize=(12, 6))
    start_index = group[0] - 400  # Start 400 samples before the first QRS peak
    end_index = group[-1] + 400  # End 400 samples after the last QRS peak
    if start_index < 0:
        start_index = 0
    if end_index > len(ecg_signal):
        end_index = len(ecg_signal)
    time = np.arange(start_index, end_index) / fs
    signal_segment = ecg_signal[start_index:end_index]
    plt.plot(time, signal_segment, label='ECG signal')
    r_points_indices = [peak for peak in group if end_index >= peak >= start_index]
    # The highest amplitude of the first cycle
    max_amp_index_first = np.argmax(signal_segment)
    max_amp_time_first = time[max_amp_index_first]
    # The first R peak
    plt.scatter(max_amp_time_first, signal_segment[max_amp_index_first], color='red', label='R peak')

    # Calculate the highest point in every 0.05 second interval
    time_points = np.arange(0.05, time[-1], 0.05)
    for t in time_points:
        # Find the index of the closest time point in the segment
        closest_index = np.argmin(np.abs(time - t))
        # Find the maximum amplitude within +/- 0.025 seconds from the closest time point
        window_indices = np.where((time >= (t - 0.025)) & (time <= (t + 0.025)))[0]
        if len(window_indices) > 0:
            max_amp_index = window_indices[np.argmax(signal_segment[window_indices])]
            max_amp_time = time[max_amp_index]

            # Plot the red point only if amplitude > 2500
            max_amplitude = signal_segment[max_amp_index]
            if max_amplitude > 2500:
                plt.scatter(max_amp_time, max_amplitude, color='red')

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'ECG Signal - Group {i+1}')
    plt.grid(True)
    plt.legend()
    plt.show()

    if i == 28:
        break
