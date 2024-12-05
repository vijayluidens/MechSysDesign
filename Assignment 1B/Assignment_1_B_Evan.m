%% MSD 2024 Assignment Part 1: Open-loop System Identification
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% Define parameters for chirp signal
clc;clear
ts = 30e-6; % Sampling Time (in seconds)
fmin = 1; % Start Frequency (in Hz)
fmax = 10000; % Max Frequency in Chirp (in Hz)

totaltime = 10; % in seconds
t = 0:ts:totaltime; % Time vector

% Define the chirp signal
tmax = t(end); % Time of chirp signal
% If tmax < t(end), append the signals to generate multiple chirp signals
% Total time for appended signal should be the specified totaltime (variable)

% Define chirp method
method = 'linear'; % 'linear' or 'logarithmic'
u = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal

%% Define the white noise (realisitic with normal distribution)
% vr = 10 % Define variance
% u = sqrt(vr)*randn(length(t),1) % Generate white noise with variance = vr

%% Call the obfuscated Plant.p function
y = Plant(u, t); % Process the chirp signal through the system

% Plot the input and output signals
figure(1);clf(1);
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

ts = 30e-6; %Sampling Time (seconds)
fs = 1/ts; %Sampling Frequency (Hz)

input = u';
output = y;

ft = logspace(-1,4,1000); % Define frequency vector data set

% Estimate Frequency Response using tfestimate()
L = length(u);
window =[]; % hann(), or rectwin(), hamming() etc.
[H,f] = tfestimate(input,output,window,[],ft,fs); 
[C,f] = mscohere(input,output,window,[],ft,fs);

[H2,f2] = tfestimate(input,output,hann(L),[],ft,fs); 
[C2,f2] = mscohere(input,output,hann(L),[],ft,fs);

[H3,f3] = tfestimate(input,output,rectwin(L),[],ft,fs); 
[C3,f3] = mscohere(input,output,rectwin(L),[],ft,fs);

figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(H))); grid on; hold on;
ylabel('|G| dB');
subplot(3,1,2);semilogx(f,rad2deg(angle(H))); grid on; hold on;
ylabel('Phase G(dB)');
subplot(3,1,3);semilogx(f,C); grid on; hold on;
ylabel('Coherence');

figure(7);clf(7);
subplot(3,1,1);semilogx(f2,mag2db(abs(H2))); grid on; hold on;
ylabel('|G| dB'); xlim([1e2 1e5])
subplot(3,1,2);semilogx(f2,rad2deg(angle(H2))); grid on; hold on;
ylabel('Phase G(dB)');
subplot(3,1,3);semilogx(f2,C2); grid on; hold on;
ylabel('Coherence');

figure(8);clf(8);
subplot(3,1,1);semilogx(f3,mag2db(abs(H3))); grid on; hold on;
ylabel('|G| dB');
subplot(3,1,2);semilogx(f3,rad2deg(angle(H3))); grid on; hold on;
ylabel('Phase G(dB)');
subplot(3,1,3);semilogx(f3,C3); grid on; hold on;
ylabel('Coherence');

% figure(4);clf(4);
% subplot(3,1,1);semilogx(f,mag2db(abs(H))); grid on; hold on;
% ylabel('|G| dB'); xlim([10 1e4]);
% subplot(3,1,2);semilogx(f,rad2deg(angle(H))); grid on; hold on;
% ylabel('Phase G(dB)'); xlim([10 1e4]);
% subplot(3,1,3);semilogx(f,C); grid on; hold on;
% ylabel('Coherence'); xlim([10 1e4]);

%% Obtain Continuous Time Transfer Function

% Generate data file to be used in tfest() function
%data = iddata(output,input,ts); %Make Input-Output Data
%OR
data = frd(H,f,ts); %Make FRD Data

% Use tfest() function to obtain transfer function from data
% General Configuration: sys = tfest(data,np,nz,iodelay);
% Important to model delay for controller design

np = 10; % Tune number of poles
nz = 8; % Tune number of zeros
iodelay = 1150e-6; % Tune delay 
sys = tfest(data,np,nz,iodelay);
figure('Name', 'p12z8d0', 'NumberTitle', 'off');
compare(data, sys);
Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump,Pdenp);

