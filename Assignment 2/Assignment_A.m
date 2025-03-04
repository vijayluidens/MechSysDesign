clear; clc; close all;

%% Load Plant
load MSD2024_P2_Plant.mat  % Load plant transfer function
s = tf('s'); 

%% Define Bode Options
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.XLim= [0.1 2000];
opts.PhaseVisible = 'on';
opts.PhaseWrapping = 'off';
opts.Grid = 'on';
opts.XLabel.String = 'Frequency';
opts.XLabel.FontSize = 16;
opts.XLabel.FontWeight = 'bold';
opts.YLabel.String = {'Magnitude','Phase'};
opts.YLabel.FontSize = 16;
opts.YLabel.FontWeight = 'bold';
opts.TickLabel.FontSize = 16;

%% Plot Plant Response
figure(1);clf(1);
bode(G, opts)
title('Transfer Function G')

%% Design PID Controller
wc = 500 * 2 * pi; % Desired crossover frequency
a = 3; % Differentiator band
kp = (1 / a) * (1 / abs(freqresp(G, wc)));
wi = wc / 10; 
wd = wc / a; 
wt = wc * a; 

C = (kp / a) * (1 + wi / s) * (s / wd + 1) / (s / wt + 1);

%% Plot Controller Response
figure(2);clf(2);
bode(C, opts);
title('PID Controller')
grid on; grid minor;

%% Open-Loop Response
L = G * C; 
figure(3);clf(3);
bode(L, opts);
title('Open-Loop')
grid on; grid minor;

figure(4); 
margin(L); % Evaluate Gain and Phase Margins

%% Sensitivity Function
S = 1 / (1 + L);
figure(5);
bode(S, opts);
title('Sensitivity Function S(s)');
grid on; grid minor;

%% Closed-Loop Step Response
T = feedback(G * C, 1); % Closed-loop transfer function
figure(6);
step(T);
title('Step Response of Closed-Loop System');
grid on; grid minor;

% %% Design a Notch Filter
% wnf = 536*2*pi; % Notch at second resonance
% Q1 = 20; Q2 = 10; % Change both ratio (Q1/Q2) and absolute value (|Q1|,|Q2|)
% N = ((s/wnf)^2 + (s/(Q1*wnf)) + 1)/((s/wnf)^2 + (s/(Q2*wnf)) + 1);
% 
% figure(5);clf(5);
% bode(N,opts) % Observe notch depth and width
% 
% wc = 120*2*pi; % Increase from 80 to 120 Hz
% a = 3; % Differentiator band
% kp = (1/a)*(1/abs(freqresp(G,wc)));
% wi = wc/10; 
% wd = wc/a; 
% wt = wc*a; 
% 
% C2 = (kp/a)*(1 + wi/s)*((s/wd + 1)/(s/wt + 1))*N; 
% 
% figure(2);clf(2);
% bode(C,C2,opts)
% legend('No Notch','With Notch',Location='southwest')
% grid on; grid minor;
% 
% % Open-Loop
% L2 = P*C2; 
% figure(2);clf(2);
% bode(L,L2,opts);
% title('Open-Loop')
% legend('No Notch','With Notch',Location='southwest')
% grid on; grid minor;
% 
% figure(3); margin(L2); % Evaluate Gain and Phase Margins 