% This script aims to model the orbit of a spacecraft as it passes through
% the atmosphere. The calculation is fairly simply; drag is calculated from
% the spacecraft's current velocity and density, which is found from the
% "getVenusDensity()" function. Acceleration and incremental position change
% are then calculated. Trajectory is assumed to follow the intitial
% hyperbolic trajectory of the spacecraft; this assumption was made for
% ease of calculations, and should not differ from the ideal performance
% significantly. 
% 
% Author: Jacob Bell (bell339@purdue.edu)
% Designed for use in AAE 251 Introduction to Aerospace Design
clc; 
clear all; 

mu = 3.2486 * 10^14; % [m^3/s^2] gravitational parameter, Venus
semi_major = -44348000; % [m] semi-major axis of the initial trajectory
venus_rad = 6051.8 * 10^3; % [m]
rp = venus_rad + 85000; % [m] closest approach, perigee
c = semi_major - rp; % [N/A] 
e0 = c / semi_major; % [N/A] eccentricity
h0 = 300000; % [m] % initial height above Venus
r0 = h0 + venus_rad; % initial radius [m] from center of Venus
v0 = sqrt(((2 * mu) / r0) - mu / semi_major); % [m/s] initial velocity 
rho0 = getVenusDensity(h0); % [kg/m^3] initial density 
cd = .1; % [N/A] Engineering toolbox estimate for cd
mass = 7500; % [kg] spacecraft mass
S = 5; % [m^2] reference area

entrance_angle = abs(acosd((rp + (rp * e0) - r0) / (r0 * e0))); % [deg.]
timestep = 1; % [sec]

% creating matrices 
time = 0:timestep:1000;

% pre-allocating
force = zeros(1,length(time));
acc = zeros(1,length(time));
rho = zeros(1,length(time));
velo = zeros(1,length(time));
radii = zeros(1,length(time));
theta = zeros(1,length(time));
work = zeros(1,length(time));
velo2 = zeros(1,length(time));

% initial calculations and setting initial conditions
velo(1) = v0;
radii(1) = r0;
theta(1) = entrance_angle;
velo2(1) = sqrt(((2 * mu) / r0) - mu / semi_major);

% this loop iterates with a timestep to calculate instaneous acceleration
% and velocity at each time interval
lcv = 2;
while lcv < length(time)-1
        % New theta and radii values from previous instantaneous velocity, 
        % radius, and theta values
        theta(lcv) = theta(lcv-1) - (velo(lcv-1) ./ radii(lcv-1)) * (180 * timestep / pi);
        radii(lcv) = (rp + rp*e0) ./ (1 + (e0 .* cosd(theta(lcv))));

        % force calcs
        rho(lcv) = getVenusDensity(radii(lcv)-venus_rad);
        force(lcv) = (.5 * rho(lcv) * velo(lcv-1)^2 * S * cd);
    
        % acceleration
        acc(lcv) = force(lcv) / mass;

        % updating next value for velocity
        velo(lcv) = velo(lcv-1) - (acc(lcv) * timestep);
        velo2(lcv) = sqrt(((2 * mu) / radii(lcv)) - mu / semi_major);
        
        if radii(lcv) > venus_rad + 300000
            break
        end

    %fprintf("%d",lcv)
    lcv = lcv + 1;
end 

% removing empty space in arrays
logical = radii==0;
force(logical) = [];
acc(logical) = [];
rho(logical) = [];
velo(logical) = [];
radii(logical) = [];
theta(logical) = [];
work(logical) = [];
velo2(logical) = [];

% plotting the modified values with drag 
figure 
plot(radii(theta~=0).*cosd(theta(theta~=0)),radii(theta~=0).*sind(theta(theta~=0)),'.')
hold on 

% Creating uniform circles for Venus and the atmosphere plotting
circle = zeros(5,360);
circles(1,:) = 1:1:360;
circles(2,:) = venus_rad; % venus itself
circles(5,:) = venus_rad + 300000; % upper layer in the atmosphere
c = [.25,.25,.1];
fill(circles(5,:).*cosd(circles(1,:)),circles(5,:).*sind(circles(1,:)),c)
alpha(.5)
c = [1,.647,0];
fill(circles(2,:).*cosd(circles(1,:)),circles(2,:).*sind(circles(1,:)),c) % venus plot

% actual hyperbola without drag
hyperbola(1,:) = -80:1:80;
hyperbola(2,:) = (rp + rp*e0) ./ (1 + (e0 .* cosd(hyperbola(1,:)))); % equation of an hyperbola

% finding the final elliptical orbit after calculating drag
dv_loss = velo2(end) - velo(end);
semi_major_final = mu / ((2 * mu) / radii(end) - velo(end).^2);
c_final = semi_major_final - rp;
e_final = c_final / semi_major_final;
final_ellipse(1,:) = -180:1:0;
final_ellipse(2,:) = (rp + rp * e_final) ./ (1 + (e_final .* cosd(final_ellipse(1,:))));

% finding the angular difference between where both trajectories exit the
% atmosphere and rotating the elliptical final orbit to match the exit
% point
[minval,indexmin] = min(abs(repmat(300000+venus_rad,1,length(final_ellipse)) - final_ellipse(2,:)));
final_ellipse(3,:) = final_ellipse(1,:) - (final_ellipse(1,indexmin)-theta(end));
logical = final_ellipse(2,:) >= venus_rad + 300000;
plot(hyperbola(2,:).*cosd(hyperbola(1,:)),hyperbola(2,:).*sind(hyperbola(1,:)),'r',"LineWidth",1) % actual hyperbola
plot((venus_rad+60000).*cosd(circles(1,:)),(venus_rad+60000).*sind(circles(1,:)),'g') % 60 km
plot(final_ellipse(2,logical).*cosd(final_ellipse(3,logical)),final_ellipse(2,logical).*sind(final_ellipse(3,logical)),'b',"LineWidth",1) 
grid on
axis equal
legend("Actual Path with Aerodynamic Drag","Atmosphere at 300 km","Venus Surface", ...
    "Actual Trajectory Without Drag","60 KM Above Surface","Resultant Elliptical Orbit ")
title("Ellptical Trajectory with Aerodynamic Drag")
xlabel("Location (m)")
ylabel("Location (m)")
xlim([-venus_rad (venus_rad+400000)])

% plotting velocity values
figure
plot(theta,velo,'r')
hold on
plot(theta,velo2,'b')
title("Orbital Velocity vs. True Anomaly")
legend("Velocity After Drag", "Velocity Without Drag")
xlabel("True Anomaly (deg.)")
ylabel("Velocity (m/s)")