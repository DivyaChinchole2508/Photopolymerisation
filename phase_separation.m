% Modeling Phase Separation Dynamics with Varying Light Intensity
clear; clc; close all;

% Parameters
L = 100;          % Domain size (arbitrary units)
N = 128;          % Number of grid points (NxN grid)
dx = L / N;       % Spatial resolution
dt = 0.01;        % Time step (arbitrary units)
T_end = 50;       % Total simulation time
save_freq = 50;   % Visualization save frequency

% Material parameters
M = 0.25;            % Mobility coefficient
kappa = 0.5;        % Gradient energy coefficient
chi = 1;          % Interaction parameter (Flory-Huggins)
k_p = 0.05;       % Photopolymerization rate constant
phi_init = 1.5;   % Initial monomer concentration

% Grid initialization
x = linspace(0, L, N); y = x;
[xx, yy] = meshgrid(x, y);
phi = phi_init + 0.1 * randn(N, N); % Initial concentration with random noise
phi(phi > 1) = 1; phi(phi < 0) = 0; % Bound concentration between 0 and 1

% Define a spatially varying light intensity pattern (e.g., Gaussian)
[I_x, I_y] = meshgrid(x, y);
I = 0.5 + 0.5 * exp(-((I_x - L/2).^2 + (I_y - L/2).^2) / (L/10)^2); % Gaussian intensity centered in the domain

% Visualization setup
figure; colormap(jet);

% Simulation loop
time_steps = ceil(T_end / dt);
for t = 1:time_steps
    % Compute chemical potential
    f_phi = log(phi + 1e-10) - log(1 - phi + 1e-10) + chi * (1 - 2 * phi);
    mu = f_phi - kappa * del2(phi, dx, dx); % Chemical potential with gradient energy
    
    % Update phase field (Cahn-Hilliard dynamics)
    lap_mu = del2(mu, dx, dx);           % Laplacian of chemical potential
    phi = phi + dt * M * lap_mu;         % Diffusion term
    
    % Monomer consumption due to photopolymerization, using spatial light intensity
    phi = phi - dt * k_p * I .* phi;     % Reaction term modified by light intensity
    
    % Enforce physical bounds
    phi(phi > 1) = 1; phi(phi < 0) = 0;  % Clamp concentration
    
    % Visualization
    if mod(t, save_freq) == 0 || t == 1
        imagesc(x, y, phi); axis equal; axis tight;
        colorbar;
        
        % Title
        title(['Phase Separation, Time = ', num2str(t * dt)]);
        xlabel('x'); ylabel('y');
        drawnow;
    end
end

% Final pattern analysis
% Fourier transform for pattern quantification
FT_phi = abs(fftshift(fft2(phi))).^2; % Power spectrum
kx = (-N/2:N/2-1) * (2*pi/L);
ky = kx;

% Find the dominant frequency
[max_val, max_idx] = max(FT_phi(:));  % Find the peak in the Fourier Transform
[kx_idx, ky_idx] = ind2sub(size(FT_phi), max_idx); % Convert index to 2D subscripts
kx_dominant = kx(kx_idx);             % Dominant wavenumber in x
ky_dominant = ky(ky_idx);             % Dominant wavenumber in y
wavelength_dominant = 2 * pi / sqrt(kx_dominant^2 + ky_dominant^2); % Characteristic wavelength

% Display results
disp(['Dominant Wavenumber: kx = ', num2str(kx_dominant), ', ky = ', num2str(ky_dominant)]);
disp(['Characteristic Wavelength: ', num2str(wavelength_dominant)]);

figure;
imagesc(kx, ky, log(FT_phi + 1e-10)); % Log-scale for better visualization
axis equal; axis tight; colormap(jet); colorbar;
title('Fourier Transform of Final Pattern');
xlabel('k_x'); ylabel('k_y');
