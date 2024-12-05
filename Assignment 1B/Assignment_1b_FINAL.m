%% MSD 2024 Assignment Part 1: Open-loop System Identification
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% Define parameters for chirp signal
ts = 30e-6; % Sampling Time (in seconds)
fmin = 1; % Start Frequency (in Hz)
fmax = 100000; % Max Frequency in Chirp (in Hz)

totaltime = 4; % in seconds
t = 0:ts:totaltime; % Time vector

% Define the chirp signal
tmax = t(end); % Time of chirp signal
% If tmax < t(end), append the signals to generate multiple chirp signals
% Total time for appended signal should be the specified totaltime (variable)

% Define chirp method
method = 'linear'; % or 'logarithmic'
u = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal

%% Define the white noise (realisitic with normal distribution)
% vr = 10; %Define variance
% u = sqrt(vr)*randn(length(t),1); % Generate white noise with variance = vr

%% Call the obfuscated Plant.p function
y = Plant(u, t); % Process the chirp signal through the system

% Plot the input and output signals
figure(1); clf(1);
subplot(2, 1, 1);
plot(t, u);
title('Input Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, y);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

%% System Identification
% This section to be run after generating the input-output signals

ts = 30e-6; % Sampling Time (seconds)
fs = 1/ts; % Sampling Frequency (Hz)

input = u';
output = y;

ft = logspace(-1,4,1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()
L = length(u);
wind = rectwin(L/2); % hann(), or rectwin() etc.
[T,f] = tfestimate(input, output, wind, [], ft, fs); 
[C,f] = mscohere(input, output, wind, [], ft, fs);

figure(2); clf(2);
subplot(3,1,1);
semilogx(f, mag2db(abs(T)), 'b'); grid on; hold on;
title('Magnitude of Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

subplot(3,1,2);
semilogx(f, rad2deg(unwrap(angle(T))), 'b'); grid on; hold on;
title('Phase of Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

subplot(3,1,3);
semilogx(f, C, 'b'); grid on; hold on;
title('Coherence');
xlabel('Frequency (Hz)');
ylabel('Coherence');

%% Obtain Continuous Time Transfer Function

% Generate data file to be used in tfest() function
% data = iddata(output, input, ts); % Make Input-Output Data
% OR
data = frd(T, f, ts); % Make FRD Data

% Use tfest() function to obtain transfer function from data
% General Configuration: sys = tfest(data, np, nz, iodelay);
% Important to model delay for controller design

np = 10; % Tune number of poles
nz = 7; % Tune number of zeros
iodelay = 1150e-6; % Tune delay 
sys = tfest(data, np, nz, iodelay);
Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump, Pdenp)

[mag, phase, f] = bode(Ptf, f); % Plant Transfer Function from Identification

% Add Bode plot of the identified transfer function to figure(2)
subplot(3,1,1);
semilogx(f, mag2db(squeeze(mag)), 'r--'); % Add magnitude in red dotted line
legend('Estimated T', 'System Identification', Location='southwest'); % Add legend

% 
subplot(3,1,2);
semilogx(f, squeeze(phase), 'r--'); % Add phase in red dotted line
legend('Estimated T', 'System Identification', Location='southwest'); % Add legend

subplot(3,1,3);
legend('Estimated T', Location='southwest'); % Only one curve here

figure(3);clf(3);
bode(Ptf)