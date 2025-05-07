% Clear workspace and close figures
close all;
clear all;
clc;

% --- Define Constants ---
mu0 = 4*pi*1e-7; % Permeability of free space (H/m)

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

% --- Define Currents (from 3.3.4a) ---
i1 = 3;  % Current in wire #1 (A)
i2 = -3; % Current in wire #2 (A)

% --- Define Plotting Region (from 3.3.4a) ---
xmin = -w;
xmax = 2*w;
ymin = 0;     % Ground plane is at y=0
ymax = 2*w;   % Upper limit for y

% --- Visualize Magnetic Vector Potential Az ---
% Define grid points
numX = 1000; % Number of points in x
numY = 1000; % Number of points in y
xVtr = linspace(xmin, xmax, numX);
yVtr = linspace(ymin, ymax, numY);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% Wire 1 and its image
r1_sq = xMtx.^2 + (yMtx - h1).^2;
r1_img_sq = xMtx.^2 + (yMtx + h1).^2;
% Wire 2 and its image
r2_sq = (xMtx - w).^2 + (yMtx - h2).^2;
r2_img_sq = (xMtx - w).^2 + (yMtx + h2).^2;

small_offset = 1e-12; % Mycket litet värde
r1_sq = r1_sq + small_offset;
r2_sq = r2_sq + small_offset;


% Compute magnetic vector potential Az using superposition
Az = (-mu0/(2*pi)) * i1 * log(sqrt(r1_img_sq)./sqrt(r1_sq)) + ...
     (-mu0/(2*pi)) * i2 * log(sqrt(r2_img_sq)./sqrt(r2_sq));

% --- Visualize Az ---
figure(1), clf; % Create figure 1 and clear it
pcolor(xMtx, yMtx, Az); % Use pcolor for filled cells
hold on; % Hold the plot for potential overlays

% Add wire locations for context (optional)
plot(0, h1, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); % Wire 1 location
plot(w, h2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8); % Wire 2 location

shading interp; % Interpolate colors for a smooth look
colorbar; % Show the color scale
axis equal; % Use equal scaling for x and y axes
axis([xmin xmax ymin ymax]); % Set axis limits to the defined region

% --- Add Labels and Title ---
xlabel('x (m)');
ylabel('y (m)');
title(['Magnetic Vector Potential A_z']);
set(gca, 'fontsize', 14); % Set font size for the axes

hold off; % Release the plot hold
