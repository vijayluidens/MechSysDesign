%% Assignment 1A: Question 6


% Define bode options
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.XLim= [1 1000];
opts.PhaseVisible = 'on';
opts.PhaseWrapping = 'on';
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 0;
opts.Grid = 'on';
opts.XLabel.String = 'Frequency';
opts.XLabel.FontSize = 16;
opts.XLabel.FontWeight = 'bold';
opts.YLabel.String = {'Magnitude','Phase'};
opts.YLabel.FontSize = 16;
opts.YLabel.FontWeight = 'bold';
opts.TickLabel.FontSize = 16;

% System Parameters
m1 = 2;    % [kg] mass of the base
m2 = 0.2;  % [kg] mass of the manipulator
m3 = 0.05; % [kg] mass of the parastic part

k1 = 1e4;  % [N/m] stiffness coefficient between the ground and the base
k2 = 3e4;  % [N/m] stiffness coefficient between the base and the manipulator
k3 = 4e4;  % [N/m] stiffness coefficient between the manipulator and parasitic part

c1 = 0.1e0; % [Ns/m] damping coefficient between the ground and the base
c2 = 0.1e0; % [Ns/m] damping coefficient between the base and the manipulator
c3 = 0.1e0; % [Ns/m] damping coefficient between the manipulator and parasitic part

% Mass matrix
M = [m1 0 0;
    0 m2 0;
    0 0 m3];

% Damping matrix
C = [c1+c2 -c2 0;
    -c2 c2+c3 -c3;
    0 -c3 c3];

%Stiffness matrix
K = [k1+k2 -k2 0;
    -k2 k2+k3 -k3
    0 -k3 k3];

%% Modal Analysis

% To compute the eigenvalues and eigenmodes with the modal analysis we
% will ignore damping. 

% Eigenvector and Eigenvalues computation
[phi, E] = eig(M^-1*K);

phi1 = phi(:,1);   % eigenvector of mode 1
phi2 = phi(:,2);   % eigenvector of mode 2
phi3 = phi(:,3);   % eigenvector of mode 3

MM1 = phi1'*M*phi1; % Modal Mass 1
MM2 = phi2'*M*phi2; % Modal Mass 2
MM3 = phi3'*M*phi3; % Modal Mass 3 

MK1=phi1'*K*phi1; % Modal Stiffness 1
MK2=phi2'*K*phi2; % Modal Stiffness 2
MK3=phi3'*K*phi3; % Modal Stiffness 3

MC1=phi1'*C*phi1; % Modal Damping 1
MC2=phi2'*C*phi2; % Modal Damping 2
MC3=phi3'*C*phi3; % Modal Damping 3

% Eigenfrequencies
f01 = sqrt(MM1\MK1)/(2*pi); % Eigenfrequency 1   
f02 = sqrt(MM2\MK2)/(2*pi); % Eigenfrequency 2
f03 = sqrt(MM3\MK3)/(2*pi); % Eigenfrequency 3

s=tf('s'); w=logspace(-1,4,20000)*2*pi;

%% Transfer Functions of Modal Analysis

% Transfer Function of x1/F1 (l=1, k=1)
G11=phi1(1)*phi1(1)/(MM1*s^2+MK1)+phi2(1)*phi2(1)/(MM2*s^2+MK2)+phi3(1)*phi3(1)/(MM3*s^2+MK3); 
%close all
figure(1); clf(1); bode(G11,w,opts); grid on; hold on
legend

%% Other DOF wrt F1
% Transfer Function of x2/F1 (l=2, k=1)
G21 = phi1(1)*phi1(2)/(MM1*s^2+MK1)+phi2(1)*phi2(2)/(MM2*s^2+MK2)+phi3(1)*phi3(2)/(MM3*s^2+MK3);
bode(G21,w,opts); grid on; hold on
legend 
% Transfer Function of x3/F1 (l=3, k=1)
G31=phi1(1)*phi1(3)/(MM1*s^2+MK1)+phi2(1)*phi2(3)/(MM2*s^2+MK2)+phi3(1)*phi3(3)/(MM3*s^2+MK3);
bode(G31,w,opts); grid on
legend

%% Transfer Functions wrt F2
% Transfer Function of x1/F2 (l=1, k=2)
l=1;k=2;
G12=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

G12f = G12 - G11; 
figure(2); clf(2); 
bode(G12f,w,opts); grid on; hold on;
legend

% Transfer Function of x2/F2 (l=2, k=2)
l=2;k=2;
G22=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
G22f = G22 - G21;
bode(G22f,w,opts); grid on; hold on

% Transfer Function of x3/F2 (l=3, k=2)
l=3;k=2;
G32=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
G32f = G32 - G31;
bode(G32f,w,opts); grid on

%% Transfer Functions wrt F3
% Transfer Function of x1/F3 (l=1, k=3)
l=1;k=3;
G13=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

G13f = G13 - G12; 
figure(3); clf(3); 
bode(G13f,w,opts); grid on; hold on;
legend

% Transfer Function of x2/F3 (l=2, k=3)
l=2;k=3;
G23=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

G23f = G23 - G22; 
bode(G23f,w,opts); grid on; hold on;
legend

% Transfer Function of x3/F3 (l=3, k=3)
l=3; k=3;
G33=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

G33f = G33 - G32;
bode(G33f, w,opts); grid on


%% Create System Transfer Function Matrix

%[x] = [MA1]*[Fnet] :- Relatting MA1 to input (Actuation Forces on Body) and output vectors
MA1 = [G11, G12, G13;
       G21, G22, G23;
       G31, G32, G33];

%[x] = [MA2]*[F] :- Relatting  MA2 to input (Actuator Forces) and output vectors
MA2 = [G11, G12f, G13f;
       G21, G22f, G23f;
       G31, G32f, G33f];

opts.PhaseVisible = 'off';

figure(4);clf(4);
bode(MA1,'k',MA2,'r--',opts)
legend

%% Compare Matrices from EOM and Modal Analysis

figure(5);clf(5);
load("data_H.mat", 'H'); %Loading H variable from Assignment A5 file
bode(H,'k',MA2,'r:',opts)
legend Location southwest
% Customize figure appearance
set(gcf, 'Color', 'w');
fig = gcf; 
fig.Color = 'w';
% Maximize the figure window
set(fig, 'WindowState', 'maximized');