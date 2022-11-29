% This script sizes and outputs a nozzle contour for a rocket thrust
% chamber given expansion ratio and throat area/radius, based on a Rao
% nozzle contour. The sizing guide used is available at
clc; 
clear all; 

% initial inputs from CEA or other source 
p_chamber = 300; % [lb/ft^2]
p_exit = 14.7; % [lb/ft^2]
exp_ratio = 8; % [N/A]
r_throat = 1; % [N/A]
theta_n = 25; % [deg]
theta_e = 8; % [deg]

% throat calculations 
theta = -135:1:-90;
x_t = 1.5 * r_throat * cosd(theta); % Eq. 4
y_t = (1.5 * r_throat * sind(theta)) + (1.5 * r_throat) + r_throat; % Eq. 4

% exit section calculations up to bell nozzle
theta = -90:1:(theta_n-90);
x_e = 0.382 * r_throat * cosd(theta); % Eq. 5
y_e = (0.382 * r_throat * sind(theta)) + (.382 * r_throat) + r_throat; % Eq. 5

% Bezier calculations for bell nozzle
t = 0:.05:1;
Nx = x_e(end);
Ny = y_e(end);
nozzle_length = .8 * (((sqrt(exp_ratio)-1) * r_throat)/ tand(15));
r_exit = sqrt(exp_ratio) * r_throat;
Ex = nozzle_length;
Ey = r_exit;
m1 = tand(theta_n); % Eq. 8
m2 = tand(theta_e); % Eq. 8
c1 = Ny - (m1 * Nx); % Eq. 9
c2 = Ey - (m2 * Ex); % Eq. 9
Qx = (c2 - c1) / (m1 - m2); % Eq. 10
Qy = ((m1 * c2) - (m2 * c1)) / (m1 - m2); % Eq. 10
x = ((1-t).^2 .* Nx) + (2 .* (1 - t) .* t .* Qx) + (t.^2 .* Ex); % Eq. 6
y = ((1-t).^2 .* Ny) + (2 .* (1 - t) .* t .* Qy) + (t.^2 .* Ey); % Eq. 6

% plotting nozzle curve 
plot(x_t,y_t,LineWidth=2,color="b")
hold on
plot(x_t,(y_t.*-1),LineWidth=2,color="b")
axis equal 
plot(x_e,y_e,LineWidth=2,color="r");
plot(x_e,(y_e.*-1),LineWidth=2,color="r");
plot(x,y,LineWidth=2,color="r")
plot(x,(y.*-1),LineWidth=2,color="r")
legend("converging section", "","diverging section")
grid on

% formatting for output to text file 
output_array = [x_t,x_e,x];
output_array(2,:) = [y_t,y_e,y];
output_array(3,:) = zeros(1,length(output_array));
output_array = transpose(output_array);
writematrix(output_array,'nozzle_contour.txt','Delimiter','space')