close all
clear all

% --- Common parameters ---
xmin = -3;
xmax = 6;
ymin = 0;
ymax = 6;

% --- Visualize scalar field ---

% Define grid points
xVtr = linspace(xmin, xmax, 90);
yVtr = linspace(ymin, ymax, 60);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% Compute scalar field
sFld = 1./sqrt(xMtx.^2 + (yMtx-2).^2);

% Concatenate scalar field (remove values that are too large)
sIdx = find(abs(sFld) > 1); sFld(sIdx) = NaN*sIdx;

% Visualize scalar field
figure(1), clf
pcolor(xMtx, yMtx, sFld), hold on
shading interp
colorbar


% --- visualize vector field ---

% Define grid points
xVtr = linspace(xmin, xmax, 18);
yVtr = linspace(ymin, ymax, 12);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% Compute vector field
vxFld = xMtx./(xMtx.^2 + (yMtx-2).^2);
vyFld = (yMtx-2)./(xMtx.^2 + (yMtx-2).^2);

% Concatenate vector field (remove values that are too large)
vIdx = find(sqrt(vxFld.^2 + vyFld.^2) > 1); 
vxFld(vIdx) = NaN*vIdx;
vyFld(vIdx) = NaN*vIdx;

% Visualize vector field
figure(1)
quiver(xMtx, yMtx, vxFld, vyFld, 0, 'k', 'linewidth', 2)
axis equal


% --- font sizes ---
set(gca, 'fontsize', 16)
axis([xmin xmax ymin ymax])