% --- Parameters based on personal identification numbers ---
p4 = 4;
p6 = 4;
p8 = 6;

% --- Calculate geometry parameters (in meters) ---
a = (0.4 + 0.2*p4)/1000;          % Wire radius [m]
w = (24 + 2*p8)/1000;             % Horizontal distance between wires [m]
h1 = (12 + p6)/1000;              % Height of wire #1 [m]
h2 = (12 + p4)/1000;              % Height of wire #2 [m]
h = h1;
l = (450 + 20*p4)/1000;           % Length [m]  
I1 = 3;
I2 = -3;

rho_l1 = 7e-9;  
rho_l2 = -7e-9; 

epsilon0 = 8.85418782e-12; % Vakuumpermittivitet [F/m]

% Definiera x-vektor för plotten
x_plot = linspace(-w, 2*w, 500); 

J_s = (h / pi) * (I1 ./ (x_plot.^2 + h^2) + (I2 ./ ((x_plot + w).^2 + h^2)));

figure(); 

plot(x_plot, J_s, 'b-', 'LineWidth', 2);
hold on; 

plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
plot(w, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 

title('induced surface current density');
xlabel('Position x [m]');
ylabel('J_s');
grid on; 

hold off;
