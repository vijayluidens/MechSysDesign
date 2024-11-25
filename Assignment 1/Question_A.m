% Given Parameters
m1 = 2; m2 = 0.2; m3 = 0.05; % Mass [kg]
k1 = 10e3; k2 = 30e3; k3 = 40e3; % Spring Constant [N/m]
c1 = 0.1; c2 = 0.1; c3 = 0.1; % Damping constant [Ns/m]

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

% Frequency range for plotting 
f = logspace(-1, 3, 1001); % Linear Frequency
w = f * 2 * pi; % Angular Frequency

% Initialize transfer function H(s) = X(s) / F(s)
H_s = zeros(3, 2, length(w));

for i = 1:length(w)
    s = 1j * w(i); % Complex frequency
    A = M * s^2 + C * s + K; % System dynamics in Laplace
    A_inv = inv(A); % Inverse of system dynamics
    H_s(:, :, i) = A_inv(:, 1:2);  % Extract the first two columns (F1 and F2)
end

% Magnitudes of transfer functions for plotting
H_magnitudes = abs(H_s);

% Plot transfer functions (magnitude vs frequency in Hz)
figure;
for i = 1:3  % Rows for x1, x2, x3
    for j = 1:2  % Columns for F1, F2
        subplot(3, 2, (i - 1) * 2 + j);
        loglog(f, squeeze(H_magnitudes(i, j, :)));
        xlabel('Frequency (Hz)');
        ylabel('|H(j\omega)|');
        title(['Transfer Function H(', num2str(i), ',', num2str(j), ')']);
        grid on;
    end
end

sgtitle('Transfer Function Matrix H(s) (Frequency in Hz)');

%% Question A.3

% In order to include damping in our system, 
% we must convert our system to state space

%State space matrix A
zero = zeros(size(M)); % Zero matrix
I = eye(size(M));      % Identity matrix
A = [zero, I;
    -inv(M)*K, -inv(M)*C];

% Eigenvalue computation
[eigenvectors, eigenvalues] = eig(A);

% Retrieve eigenvalues
eigenvalues = diag(eigenvalues); % Convert diagonal matrix to column vector

