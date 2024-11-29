%% 
clear
clc

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

%% Symbolically compute EoM

%Create s

syms s x1 x2 x3 dx1 dx2 dx3 ddx1 ddx2 ddx3 m1 m2 m3 c1 c2 c3 k1 k2 k3 F1 F2 F3

% System matrices

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

%Force matrix
F = [F1 - F2;
    F2; 
    0];

X = [x1; x2; x3];  % Displacement
Xd = [dx1; dx2; dx3];  % Velocity
Xdd = [ddx1; ddx2; ddx3];  % Acceleration

% General Form of Equations of Motion (EOM)
EOM = M*Xdd + C*Xd + K*X - F;
EOM = subs(EOM, [dx1 dx2 dx3 ddx1 ddx2 ddx3], [s*x1 s*x2 s*x3 s^2*x1 s^2*x2 s^2*x3]);

%% Defining Parameters for the transfer function x/F1 (wrt 1st Actuator Force)
c = 0.1;
EOMtf = subs(EOM, [m1 m2 m3 c1 c2 c3 k1 k2 k3], [2 0.2 0.05 c c c 1e4 3e4 4e4]);
EOMtf = subs(EOMtf, F2, 0);  % F2 is set to 0

% The system is a set of equations for x1, x2, and x3.
% Create the system of equations
eqs = EOMtf;

% Solve for x1, x2, and x3 in terms of F1
sol = solve(eqs, [x1, x2, x3]);

% Extract x in terms of F1
x1_sol = sol.x1;
x2_sol = sol.x2;
x3_sol = sol.x3;

% Now, compute the transfer function H(s) = x1/F1
H_x1F1 = simplify(x1_sol / F1);
H_x2F1 = simplify(x2_sol / F1);
H_x3F1 = simplify(x3_sol / F1);

% Display the transfer function
disp('Transfer Function x1/F1:');
disp(H_x1F1);

% Extract the coefficients of the numerator and denominator:
[num, den] = numden(simplify(x1_sol / F1));
[num2, den2] = numden(simplify(x2_sol / F1));
[num3, den3] = numden(simplify(x3_sol / F1));

% Now, you can create the transfer function using the coefficients:
H_x1F1 = tf([sym2poly(num)], [sym2poly(den)]);
H_x2F1 = tf([sym2poly(num2)], [sym2poly(den2)]);
H_x3F1 = tf([sym2poly(num3)], [sym2poly(den3)]);

% Display the transfer function
disp('Transfer Function x1/F1:');
disp(H_x1F1);

figure(1);clf(1);
bode(H_x1F1,H_x2F1,H_x3F1,opts)
legend

%% Defining Parameters for the transfer function x/F2 (wrt 2nd Actuator Force)
c = 0.1e0;
EOMtf = subs(EOM, [m1 m2 m3 c1 c2 c3 k1 k2 k3], [2 0.2 0.05 c c c 1e4 3e4 4e4]);
EOMtf = subs(EOMtf, F1, 0);  % F1 is set to 0 

% Create the system of equations
eqs = EOMtf;

% Solve for x1, x2, and x3 in terms of F2
sol = solve(eqs, [x1, x2, x3]);

% Extract x1 in terms of F1 (assuming we want x1 only)
x1_sol = sol.x1;
x2_sol = sol.x2;
x3_sol = sol.x3;

% Now, compute the transfer function H(s) = x1/F1
H_x1F2 = simplify(x1_sol / F2);
H_x2F2 = simplify(x2_sol / F2);
H_x3F2 = simplify(x3_sol / F2);

% Extract the coefficients of the numerator and denominator:
[num, den] = numden(simplify(x1_sol / F2));
[num2, den2] = numden(simplify(x2_sol / F2));
[num3, den3] = numden(simplify(x3_sol / F2));

% Now, you can create the transfer function using the coefficients:
H_x1F2 = tf([sym2poly(num)], [sym2poly(den)]);
H_x2F2 = tf([sym2poly(num2)], [sym2poly(den2)]);
H_x3F2 = tf([sym2poly(num3)], [sym2poly(den3)]);

figure(2);clf(2);
bode(H_x1F2,H_x2F2,H_x3F2,opts)
legend

%% Create System Transfer Function Matrix

%[x] = [H]*[F] :- Relatting H to input and output vectors
H = [H_x1F1, H_x1F2;
     H_x2F1, H_x2F2;
     H_x3F1, H_x3F2];

% Input-Output Names
H.u{1} = 'F_1'; H.u{2} = 'F_2';
H.y{1} = 'x_1'; H.y{2} = 'x_2'; H.y{3} = 'x_3';

% Plot the System Matrix
figure(3);clf(3);
bode(H,opts)

% Customize figure appearance
set(gcf, 'color', 'w');
fig = gcf; 
fig.Color = 'w';
% Maximize the figure window
set(fig, 'WindowState', 'maximized');

%% Question A.3

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

%% Transfer Functions from Modal Analysis
% General Form
% xl/Fk (for i^th mode) = phi_i(l)*phi_i(k)/(MM_i*s^2 + MK_i)

% Transfer Function of x1/F1 (l=1, k=1)
xF11=phi1(1)*phi1(1)/(MM1*s^2+MK1)+phi2(1)*phi2(1)/(MM2*s^2+MK2)+phi3(1)*phi3(1)/(MM3*s^2+MK3); % Damping Excluded
%close all
figure(1); clf(1); bode(xF11,w,opts); grid on
legend

%% Other DOF wrt F1
% Transfer Function of x2/F1 (l=2, k=1)
xF21 = phi1(1)*phi1(2)/(MM1*s^2+MK1)+phi2(1)*phi2(2)/(MM2*s^2+MK2)+phi3(1)*phi3(2)/(MM3*s^2+MK3);
bode(xF21,w,opts); grid on
legend 
% Transfer Function of x3/F1 (l=3, k=1)
xF31=phi1(1)*phi1(3)/(MM1*s^2+MK1)+phi2(1)*phi2(3)/(MM2*s^2+MK2)+phi3(1)*phi3(3)/(MM3*s^2+MK3);
bode(xF31,w,opts); grid on
legend

%% Transfer Functions wrt F2
% Transfer Function of x1/F2 (l=1, k=2)
l=1;k=2;
xF12=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

xF12f = xF12 - xF11; % Net Force on m1 is difference of two actuator forces (F1-F2)
figure(2); clf(2); 
bode(xF12f,w,opts); grid on; hold on;
legend

% Transfer Function of x2/F2 (l=2, k=2)
l=2;k=2;
xF22=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
xF22f = xF22 - xF21;
bode(xF22f,w,opts); grid on

% Transfer Function of x3/F2 (l=3, k=2)
l=3;k=2;
xF32=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
xF32f = xF32 - xF31;
bode(xF32f,w,opts); grid on


%% %% Create System Transfer Function Matrix

%[x] = [MA1]*[Fnet] :- Relatting MA1 to input (Actuation Forces on Body) and output vectors
MA1 = [xF11, xF12;
       xF21, xF22;
       xF31, xF32];

%[x] = [MA2]*[F] :- Relatting  MA2to input (Actuator Forces) and output vectors
MA2 = [xF11, xF12f;
       xF21, xF22f;
       xF31, xF32f];

opts.PhaseVisible = 'off';

figure(3);clf(3);
bode(MA1,'k',MA2,'r--',opts)
legend


%% Modal Analysis with Damping 

% Transfer Function of x1/F1 (l=1, k=1)
xF11d=phi1(1)*phi1(1)/(MM1*s^2+MC1*s+MK1)+phi2(1)*phi2(1)/(MM2*s^2+MC2*s+MK2)+phi3(1)*phi3(1)/(MM3*s^2+MC3*s+MK3); % Damping Included

% Transfer function of x2/F1 (l=2, k=1)
xF21d = phi1(1)*phi1(2)/(MM1*s^2 + MC1*s + MK1)+phi2(1)*phi2(2)/(MM2*s^2+MC2*s+MK2)+phi3(1)*phi3(2)/(MM3*s^2+MC3*s+MK3); % Damping included

% Transfer Function of x3/F1 (l=3, k=1)
xF31d = phi1(1)*phi1(3)/(MM1*s^2 + MC1*s + MK1)+phi2(1)*phi2(3)/(MM2*s^2+MC2*s+MK2)+phi3(1)*phi3(3)/(MM3*s^2+MC3*s+MK3);