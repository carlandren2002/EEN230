% --- Parameters based on personal identification numbers ---
p4 = 4;
p6 = 4;
p8 = 6;

% --- Calculate geometry parameters (in meters) ---
w = (24 + 2*p8)/1000;             % Horizontal distance between wires [m]
h1 = (12 + p6)/1000;              % Height of wire #1 [m]
h2 = (12 + p4)/1000;              % Height of wire #2 [m]

% --- Current parameters ---
i1 = 3;                           % Current in wire #1 [A]
i2 = -3;                          % Current in wire #2 [A]
mu0 = 4*pi*1e-7;                  % Permeability of free space [H/m]

% --- Grid for visualization ---
xmin = -w;
xmax = 2*w;
ymin = 0;                         % Start from y=0 (ground plane)
ymax = 2*w;

xVtr = linspace(xmin, xmax, 25);  % Reduced resolution for quiver plot
yVtr = linspace(ymin, ymax, 25);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% --- Initialize magnetic field components ---
Bx = zeros(size(xMtx));
By = zeros(size(xMtx));

% --- Small offset to avoid division by zero ---
small_offset = 1e-12;

% --- Calculate B-field contributions from each wire and its image ---
% Wire 1 at (0, h1) and its image at (0, -h1)
r1_sq = xMtx.^2 + (yMtx - h1).^2 + small_offset;
r1prime_sq = xMtx.^2 + (yMtx + h1).^2 + small_offset;

Bx = Bx + (mu0/(2*pi)) * i1 * (-(yMtx - h1) ./ r1_sq + (yMtx + h1) ./ r1prime_sq);
By = By + (mu0/(2*pi)) * i1 * (xMtx ./ r1_sq - xMtx ./ r1prime_sq);

% Wire 2 at (w, h2) and its image at (w, -h2)
r2_sq = (xMtx - w).^2 + (yMtx - h2).^2 + small_offset;
r2prime_sq = (xMtx - w).^2 + (yMtx + h2).^2 + small_offset;

Bx = Bx + (mu0/(2*pi)) * i2 * (-(yMtx - h2) ./ r2_sq + (yMtx + h2) ./ r2prime_sq);
By = By + (mu0/(2*pi)) * i2 * ((xMtx - w) ./ r2_sq - (xMtx - w) ./ r2prime_sq);

% --- Calculate magnitude of B-field for color map ---
B_magnitude = sqrt(Bx.^2 + By.^2);

% --- Plotting ---
figure;
pcolor(xMtx, yMtx, B_magnitude);
shading interp;
hold on;

% Quiver plot for vector field
quiver(xMtx, yMtx, Bx, By, 1, 'k', 'linewidth', 2)

% Plot wire positions
plot(0, h1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Wire 1
plot(w, h2, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Wire 2

% Plot ground plane
plot([xmin xmax], [0 0], 'k-', 'LineWidth', 2);

% Add colorbar and labels
colorbar;
title('Magnetic Flux Density');
xlabel('Position x [m]');
ylabel('Position y [m]');
axis equal;
axis([xmin xmax ymin ymax]);

hold off;