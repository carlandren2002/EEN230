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

rho_l1 = 7e-9;  
rho_l2 = -7e-9; 

epsilon0 = 8.85418782e-12; % Vakuumpermittivitet [F/m]
k = 1/(2*pi*epsilon0);     % Konstant för potential från linjeladdning

% --- Grid för visualisering ---
xmin = -w;
xmax = 2*w;
ymin = 0; % Starta från y=0 (jordplan)
ymax = 2*w;

xVtr = linspace(xmin, xmax, 180); 
yVtr = linspace(ymin, ymax, 120);
[xMtx, yMtx] = meshgrid(xVtr, yVtr);

% V(x,y) = (k/2)*rho_l1*log( (x^2+(y+h1)^2) / (x^2+(y-h1)^2) ) + ...
%          (k/2)*rho_l2*log( ((x-w)^2+(y+h2)^2) / ((x-w)^2+(y-h2)^2) )

r1_sq = xMtx.^2 + (yMtx - h1).^2;
r1prime_sq = xMtx.^2 + (yMtx + h1).^2;
r2_sq = (xMtx - w).^2 + (yMtx - h2).^2;
r2prime_sq = (xMtx - w).^2 + (yMtx + h2).^2;

% Lägg till liten offset för att undvika log(0) eller division med noll
% exakt vid trådpositionerna (y=h1, x=0) och (y=h2, x=w)
small_offset = 1e-12; % Mycket litet värde
r1_sq = r1_sq + small_offset;
r2_sq = r2_sq + small_offset;

% Beräkna termerna i potentialformeln
term1_V = (k/2) * rho_l1 * log(r1prime_sq ./ r1_sq);
term2_V = (k/2) * rho_l2 * log(r2prime_sq ./ r2_sq);

% Total potential
V = term1_V + term2_V;

pcolor(xMtx, yMtx, V);
shading interp;
hold on;

contour(xMtx, yMtx, V, 20, 'k'); 

colorbar;

plot(0, h1, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
plot(w, h2, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 

plot([xmin xmax], [0 0], 'k-', 'LineWidth', 2);

title('Elektrisk Potential V(x,y)');
xlabel('Position x [m]');
ylabel('Position y [m]');
axis equal; 

axis([xmin xmax ymin ymax]);

hold off;
