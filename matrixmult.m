% --- Parameters based on personal identification numbers ---
p4 = 4;
p6 = 4;
p8 = 6;
% --- Calculate geometry parameters (in meters) ---
a = (0.4 + 0.2*p4)/1000; % Wire radius [m]
w = (24 + 2*p8)/1000; % Horizontal distance between wires [m]
h1 = (12 + p6)/1000; % Height of wire #1 [m]
h2 = (12 + p4)/1000; % Height of wire #2 [m]
l = (450 + 20*p4)/1000; % Length [m]
% --- Given line charge densities ---
rho_l1 = 7e-9; % Line charge density on wire #1 [C/m]
rho_l2 = -7e-9; % Line charge density on wire #2 [C/m]
epsilon0 = 8.85418782e-12;
k = 1/(2*pi*epsilon0);
d = sqrt((h1-h2)^2 + w^2);

rho_l = [rho_l1; rho_l2]; 


P = [
    k*log(2*h1/a), k*log(d/w);
    k*log(d/w), k*log(2*h2/a)
];

C = inv(P);

fprintf('P\n');
disp(P);
fprintf('\nC\n');
disp(C);