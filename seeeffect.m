clear; clc; close all;

% Parameters
L = 100;          % Domain size (arbitrary units)
N = 128;          % Number of grid points (NxN grid)
dx = L / N;       % Spatial resolution
dt = 0.01;        % Time step (arbitrary units)
T_end = 50;       % Total simulation time
phi_init = 1.5;   % Initial monomer concentration
k_p = 0.05;       % Photopolymerization rate constant

% Grid initialization
x = linspace(0, L, N);
y = x;
[xx, yy] = meshgrid(x, y);

% Define a spatially varying light intensity pattern (Gaussian)
[I_x, I_y] = meshgrid(x, y);
I = 0.5 + 0.5 * exp(-((I_x - L/2).^2 + (I_y - L/2).^2) / (L/10)^2); % Gaussian intensity

% Parameter ranges to test
M_values = [0.1, 0.25, 0.5];      % Mobility coefficient
kappa_values = [0.1, 0.5, 1.0];   % Gradient energy coefficient
chi_values = [0.5, 1.0, 1.5];     % Interaction parameter

% Fourier domain setup
kx = (-N/2:N/2-1) * (2*pi/L);    % Wavenumber in x-direction
ky = kx;                        % Wavenumber in y-direction

% Create separate figures for phase separation and Fourier transforms
figure('Name', 'Phase Separation Patterns', 'Position', [100, 100, 1200, 800]);
subplot_idx = 1;

% Run simulations for each combination of parameters
for M_idx = 1:length(M_values)
    for kappa_idx = 1:length(kappa_values)
        for chi_idx = 1:length(chi_values)
            % Extract current parameter values
            M = M_values(M_idx);
            kappa = kappa_values(kappa_idx);
            chi = chi_values(chi_idx);
            
            % Initialize concentration field
            phi = phi_init + 0.1 * randn(N, N); 
            phi(phi > 1) = 1; phi(phi < 0) = 0; % Clamp concentration
            
            % Simulation loop
            time_steps = ceil(T_end / dt);
            for t = 1:time_steps
                % Compute chemical potential
                f_phi = log(phi + 1e-10) - log(1 - phi + 1e-10) + chi * (1 - 2 * phi);
                mu = f_phi - kappa * del2(phi, dx, dx); % Chemical potential
                
                % Update phase field (Cahn-Hilliard dynamics)
                lap_mu = del2(mu, dx, dx);          % Laplacian of chemical potential
                phi = phi + dt * M * lap_mu;        % Diffusion term
                
                % Monomer consumption due to photopolymerization
                phi = phi - dt * k_p * I .* phi;    % Reaction term modified by light intensity
                
                % Enforce physical bounds
                phi(phi > 1) = 1; phi(phi < 0) = 0; % Clamp concentration
            end
            
            % Plot phase separation pattern
            subplot(3, 3, subplot_idx);
            imagesc(x, y, phi); axis equal; axis tight; colormap(jet);
            title(['M=', num2str(M), ', \kappa=', num2str(kappa), ', \chi=', num2str(chi)]);
            xlabel('x'); ylabel('y');
            
            % Increment subplot index for phase separation figure
            subplot_idx = subplot_idx + 1;
        end
    end
end

% Create separate figure for Fourier Transforms
figure('Name', 'Fourier Transforms', 'Position', [100, 100, 1200, 800]);
subplot_idx = 1;

% Run simulations for each combination of parameters to calculate Fourier Transforms
for M_idx = 1:length(M_values)
    for kappa_idx = 1:length(kappa_values)
        for chi_idx = 1:length(chi_values)
            % Extract current parameter values
            M = M_values(M_idx);
            kappa = kappa_values(kappa_idx);
            chi = chi_values(chi_idx);
            
            % Initialize concentration field
            phi = phi_init + 0.1 * randn(N, N); 
            phi(phi > 1) = 1; phi(phi < 0) = 0; % Clamp concentration
            
            % Simulation loop
            time_steps = ceil(T_end / dt);
            for t = 1:time_steps
                % Compute chemical potential
                f_phi = log(phi + 1e-10) - log(1 - phi + 1e-10) + chi * (1 - 2 * phi);
                mu = f_phi - kappa * del2(phi, dx, dx); % Chemical potential
                
                % Update phase field (Cahn-Hilliard dynamics)
                lap_mu = del2(mu, dx, dx);          % Laplacian of chemical potential
                phi = phi + dt * M * lap_mu;        % Diffusion term
                
                % Monomer consumption due to photopolymerization
                phi = phi - dt * k_p * I .* phi;    % Reaction term modified by light intensity
                
                % Enforce physical bounds
                phi(phi > 1) = 1; phi(phi < 0) = 0; % Clamp concentration
            end
            
            % Compute Fourier Transform of the final phase pattern
            FT_phi = abs(fftshift(fft2(phi))).^2; % Power spectrum
            
            % Plot Fourier Transform
            subplot(3, 3, subplot_idx);
            imagesc(kx, ky, log(FT_phi + 1e-10)); axis equal; axis tight; colormap(jet);
            title(['FT: M=', num2str(M), ', \kappa=', num2str(kappa), ', \chi=', num2str(chi)]);
            xlabel('k_x'); ylabel('k_y');
            
            % Increment subplot index for Fourier transform figure
            subplot_idx = subplot_idx + 1;
        end
    end
end
