% System Parameters
m1 = 2; %kg
m2 = 0.2; %kg
m3 = 0.05; %kg
k1 = 10e3; %N/m
k2 = 30e3; %N/m
k3 = 40e3; %N/m
c1 = 0.1; %Ns/m
c2 = 0.1; %Ns/m
c3 = 0.1; %Ns/m


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

%Eigenmatrix

% This is an updated version



