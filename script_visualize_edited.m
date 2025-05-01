close all
clear all

% --- Parameters based on personal identification numbers ---
p4 = 4;
p6 = 4;
p8 = 6;

% --- Calculate geometry parameters (in meters) ---
a = (0.4 + 0.2*p4)/1000;          % Wire radius [m]
w = (24 + 2*p8)/1000;             % Horizontal distance between wires [m]
h1 = (12 + p6)/1000;              % Height of wire #1 [m]
h2 = (12 + p4)/1000;              % Height of wire #2 [m]
l = (450 + 20*p4)/1000;           % Length [m]

% --- Charge densities ---
rho_l1 = 7e-9;  % Line charge density of wire #1 [C/m]
rho_l2 = -7e-9; % Line charge density of wire #2 [C/m]

% --- Common parameters ---
xmin = -w;
xmax = 2*w;
ymin = 0;
ymax = 2*w;

% Constants
epsilon0 = 8.85418782e-12; % Vacuum permittivity [F/m]
k = 1/(2*pi*epsilon0);     % Coulomb constant for line charges

% --- Visualize scalar field ---
% Define grid points
xVtr = linspace(xmin, xmax, 90);
yVtr = linspace(ymin, ymax, 60);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% Initialize electric field components
Ex = zeros(size(xMtx));
Ey = zeros(size(xMtx));

% Calculate electric field due to wire #1 at (0,h1) and its image at (0,-h1)
% Real wire #1
r1_squared = (xMtx).^2 + (yMtx - h1).^2;
Ex = Ex + k * rho_l1 * (xMtx) ./ r1_squared;
Ey = Ey + k * rho_l1 * (yMtx - h1) ./ r1_squared;

% Image of wire #1
r1_image_squared = (xMtx).^2 + (yMtx + h1).^2;
Ex = Ex + k * (-rho_l1) * (xMtx) ./ r1_image_squared;
Ey = Ey + k * (-rho_l1) * (yMtx + h1) ./ r1_image_squared;

% Calculate electric field due to wire #2 at (w,h2) and its image at (w,-h2)
% Real wire #2
r2_squared = (xMtx - w).^2 + (yMtx - h2).^2;
Ex = Ex + k * rho_l2 * (xMtx - w) ./ r2_squared;
Ey = Ey + k * rho_l2 * (yMtx - h2) ./ r2_squared;

% Image of wire #2
r2_image_squared = (xMtx - w).^2 + (yMtx + h2).^2;
Ex = Ex + k * (-rho_l2) * (xMtx - w) ./ r2_image_squared;
Ey = Ey + k * (-rho_l2) * (yMtx + h2) ./ r2_image_squared;

% Calculate electric field magnitude
E_magnitude = sqrt(Ex.^2 + Ey.^2);

% Concatenate scalar field (remove values that are too large)
sIdx = find(E_magnitude > 3e6); E_magnitude(sIdx) = NaN*sIdx;

% Visualize scalar field
figure(1), clf
pcolor(xMtx, yMtx, E_magnitude), hold on
shading interp
colorbar

% --- visualize vector field ---
% Define grid points for vector field
xVtr = linspace(xmin, xmax, 18);
yVtr = linspace(ymin, ymax, 12);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% Initialize electric field components
vxFld = zeros(size(xMtx));
vyFld = zeros(size(xMtx));

% Add a small offset to avoid division by zero
small_offset = 1e-6;

% Calculate field due to wire #1 and its image
r1_squared = (xMtx).^2 + (yMtx - h1).^2 + small_offset;
vxFld = vxFld + k * rho_l1 * (xMtx) ./ r1_squared;
vyFld = vyFld + k * rho_l1 * (yMtx - h1) ./ r1_squared;

r1_image_squared = (xMtx).^2 + (yMtx + h1).^2 + small_offset;
vxFld = vxFld + k * (-rho_l1) * (xMtx) ./ r1_image_squared;
vyFld = vyFld + k * (-rho_l1) * (yMtx + h1) ./ r1_image_squared;

% Calculate field due to wire #2 and its image
r2_squared = (xMtx - w).^2 + (yMtx - h2).^2 + small_offset;
vxFld = vxFld + k * rho_l2 * (xMtx - w) ./ r2_squared;
vyFld = vyFld + k * rho_l2 * (yMtx - h2) ./ r2_squared;

r2_image_squared = (xMtx - w).^2 + (yMtx + h2).^2 + small_offset;
vxFld = vxFld + k * (-rho_l2) * (xMtx - w) ./ r2_image_squared;
vyFld = vyFld + k * (-rho_l2) * (yMtx + h2) ./ r2_image_squared;

% Calculate vector magnitudes
v_magnitude = sqrt(vxFld.^2 + vyFld.^2);

% Normalize vectors for better visualization
vxFld_norm = vxFld ./ v_magnitude;
vyFld_norm = vyFld ./ v_magnitude;

% Handle potential infinity/NaN values
vxFld_norm(isnan(vxFld_norm) | isinf(vxFld_norm)) = 0;
vyFld_norm(isnan(vyFld_norm) | isinf(vyFld_norm)) = 0;

% Visualize vector field with normalized vectors
figure(1)
quiver(xMtx, yMtx, vxFld_norm, vyFld_norm, 0.5, 'k', 'linewidth', 2)

% Add wire positions
plot(0, h1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % Wire #1
plot(w, h2, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');  % Wire #2
plot([xmin xmax], [0 0], 'k-', 'LineWidth', 2);  % Ground plane

axis equal
% --- font sizes ---
set(gca, 'fontsize', 16)
axis([xmin xmax ymin ymax])