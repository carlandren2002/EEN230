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

% Definiera x-vektor för plotten
x_plot = linspace(-w, 2*w, 500); % 500 punkter i intervallet [-w, 2w]

% Beräkna rho_s(x) enligt formeln 
% rho_s(x) = -(1/pi) * [ (rho_l1 * h1) / (x^2 + h1^2) + (rho_l2 * h2) / ((x-w)^2 + h2^2) ]
term1 = (rho_l1 * h1) ./ (x_plot.^2 + h1^2);
term2 = (rho_l2 * h2) ./ ((x_plot - w).^2 + h2^2);
rho_s = -(1/pi) * (term1 + term2);

figure(); 

plot(x_plot, rho_s, 'b-', 'LineWidth', 2);
hold on; 

plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
plot(w, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 

title('Inducerad ytladdningstäthet på jordplanet');
xlabel('Position x [m]');
ylabel('\rho_s [C/m^2]');
grid on; 

hold off;
