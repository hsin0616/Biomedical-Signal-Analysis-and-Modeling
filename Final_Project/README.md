# Biomedical Signal Analysis and Modeling
## Final Project
### Objectives
1) Segmentation of phonocardiographic (PCG) signals using the ECG as a reference. 
2) Computation of an averaged power spectral density (PSD) of the systolic segments of 
a PCG signal. 
3) Parametric representation of PCG signals. 
4) Auditory analysis and classification of heart sounds. 

### Specific Tasks
#### 1. Data collection + basic preprocessing (PCG/ECG/carotid)
- Use pec1.dat, pec33.dat, pec41.dat, pec42.dat, pec52.dat and pec_sound.m (from final.zip → “BSA files”).
- Each **.dat** file contains three channels: ECG, PCG, carotid pulse.
- Sampling rate per channel: fs = 1000 Hz.
- Use pec_sound.m to (1) read the file, (2) separate channels, (3) plot signals, (4) listen to the PCG.
- Some records contain artifacts at the beginning/end → remove those portions before analysis.
- Subject types (1) Normal: pec1, pec52 (2) Ventricular septal defect (VSD) → systolic murmur: pec33, pec42 (3) Aortic stenosis → systolic murmur: pec41.

#### 2. QRS detection in ECG (Balda et al. method) + mapping to PCG S1 start
- Detect QRS complexes using Balda et al. algorithm (Textbook Section 4.3.1).
- The algorithm expects 200 Hz sampling, but ECG here is 1000 Hz.
- Apply a Butterworth low-pass filter to the ECG channel only (anti-aliasing).
- Downsample by factor of 5 (1000 Hz → 200 Hz).
- Run the QRS detection on the 200 Hz ECG.
- After detecting QRS times, convert detected QRS indices back to the original timing (scaling factor), and transfer each detected QRS time to the corresponding start of S1 in the PCG channel by applying the needed correction/offset (alignment between ECG electrical event and PCG mechanical sound).

#### 3. Segment systolic PCG + PSD + synchronized averaging
- For each cardiac cycle, segment the systolic portion of PCG using a window: 300–350 ms starting at the beginning of each QRS complex.
- Choose window length so that it including the beginning of S1, systolic murmur (if present), but excludes S2.
- Compute PSD (Power Spectral Density) for: one individual systolic segment, and the averaged systolic signal/PSD across many beats.
- Perform synchronized averaging: use as many clean cycles as possible, align all segments using QRS timing, compute the average PSD (more stable spectral estimate).
- Plot comparisons:
PSD of a single systolic segment vs PSD of the averaged systolic segments.
#### 4. Frequency-domain features + spectral ratio + results table
- From the averaged systolic PSD for each subject, compute median frequency and estimate bandwidth of the PCG from the PSD shape (approx frequency range containing most energy).
- Design a spectral ratio feature (similar to Eq. 6.49) to pick suitable frequency limits f1, f2, f3 and define ratio of energy in different bands to separate normal PCG vs systolic murmur cases.
- Prepare a table for all 5 signals including: median frequency (Hz), estimated bandwidth (Hz), spectral ratio (unitless), any other PSD-based parameters you compute.
- Include correct units: frequency-related → Hz、PSD units depend on scaling; if not calibrated, label as AU (arbitrary units).
- Analyze and discuss how parameters differ between normal and murmur cases, which features best separate the classes.
#### 5. Listening test + relate sound to PSD
- Listen to each PCG signal and write auditory observations: S1/S2 clarity, presence of systolic murmur (e.g., harsh/noisy/“whooshing” during systole), loudness differences, timing consistency, etc.
- Relate auditory impressions to frequency results: murmurs typically show more high-frequency energy and broader bandwidth in systole.
- For listening comfort, resample PCG from 1 kHz → 8 kHz using interp, and remove artifact segments first (loud strange sounds can bias perception).
- Export audio using the approach in pec_sound.m: (1) convert to .au audio file, (2)optionally export to .wav using wavwrite,listen using audio tools.
#### Plot/Reporting requirements (applies to all parts)
- Use figures to show: signals before/after filtering, QRS detection marks, systolic segmentation windows, PSD comparisons, feature distributions (normal vs murmur).
- Always label axes with correct units.
- If units are not calibrated/unknown → label as “AU”.