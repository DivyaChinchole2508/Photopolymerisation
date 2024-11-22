% Parameters
f = 0.8;            % Initiation efficiency (fraction)
k_d = 1.0e-3;       % Decomposition rate constant (s^-1)
k_p = 1.0;          % Propagation rate constant (L/mol/s)
k_t = 1.0e5;        % Termination rate constant (L/mol/s)
I0 = 1.0;           % Initial photoinitiator concentration (mol/L)
M0 = 1.0;           % Initial monomer concentration (mol/L)
R0 = 0;             % Initial radical concentration (mol/L)
t_final = 5000;    % Total time (s)
dt = 1;             % Time step (s)

% Time vector
t = 0:dt:t_final;
N = length(t);

% Initialize concentration, rate, and conversion vectors
I = zeros(1, N);    % Photoinitiator concentration
M = zeros(1, N);    % Monomer concentration
R = zeros(1, N);    % Radical concentration
R_i = zeros(1, N);  % Initiation rate
R_p = zeros(1, N);  % Propagation rate
R_t = zeros(1, N);  % Termination rate
X = zeros(1, N);    % Conversion

% Set initial concentrations
I(1) = I0;
M(1) = M0;
R(1) = R0;

% Finite difference solution
for i = 2:N
    % Initiation rate
    R_i(i-1) = 2 * f * k_d * I(i-1);
    
    % Update photoinitiator concentration
    dI_dt = -k_d * I(i-1); 
    I(i) = I(i-1) + dI_dt * dt;
    if I(i) < 0, I(i) = 0; end  % Ensure non-negative concentration
    
    % Update radical concentration
    R_t(i-1) = k_t * R(i-1)^2;  % Termination rate
    dR_dt = R_i(i-1) - R_t(i-1); 
    R(i) = R(i-1) + dR_dt * dt;
    if R(i) < 0, R(i) = 0; end  % Ensure non-negative concentration
    
    % Propagation rate
    R_p(i-1) = k_p * M(i-1) * R(i-1);
    
    % Update monomer concentration
    dM_dt = -R_p(i-1);  
    M(i) = M(i-1) + dM_dt * dt;
    if M(i) < 0, M(i) = 0; end  % Ensure non-negative concentration
    
    % Calculate conversion
    X(i) = (M0 - M(i)) / M0;
end

% Calculate final rates and conversion
R_i(N) = 2 * f * k_d * I(N);
R_t(N) = k_t * R(N)^2;
R_p(N) = k_p * M(N) * R(N);
X(N) = (M0 - M(N)) / M0;

% Plot the results
figure;

% Plot initiation rate
subplot(3,2,1);
plot(t, R_i, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Initiation Rate, R_i (mol/L/s)');
title('Initiation Rate vs. Time');
grid on;

% Plot propagation rate
subplot(3,2,2);
plot(t, R_p, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Propagation Rate, R_p (mol/L/s)');
title('Propagation Rate vs. Time');
grid on;

% Plot termination rate
subplot(3,2,3);
plot(t, R_t, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Termination Rate, R_t (mol/L/s)');
title('Termination Rate vs. Time');
grid on;

% Plot radical concentration
subplot(3,2,4);
plot(t, R, 'g', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Radical Concentration, [R] (mol/L)');
title('Radical Concentration vs. Time');
grid on;

% Plot monomer concentration
subplot(3,2,5);
plot(t, M, 'm', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Monomer Concentration, [M] (mol/L)');
title('Monomer Concentration vs. Time');
grid on;

% Plot conversion
subplot(3,2,6);
plot(t, X, 'c', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Conversion, X');
title('Conversion vs. Time');
grid on;
