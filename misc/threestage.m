% This script takes an input for the total delta v required for a
% 3 stage maneuver of some kind, as well as the three specific impuse values of the
% stages, and finds the optimal delta v fraction for the manuever.
% AAE 251 Fall 2022
% PM7 Team R301 
% Rocket staging optimization
% Authors: Jacob Bell
% Collaborators: Caroline Gee, Sydney Brown, Thomas Brannon, Justin
% Armstrong, Lu Larest, Carter Doby, Abby Woodburry 
% 
% This code creates a 3D mesh visualizing the GLOW for various combinations of
% 2nd and 3rd stage delta v split (1st stage delta v fraction is found
% from 2nd and 3rd). A minimum value is found for GLOW corresponding to 
% optimal staging for a given maneuver
clc;
clear all; 

% inputs/known variables 
dv_total = 11.4 * 10^3; % [m/s] total delta-v required 
payload_mass = 10000; % [kg] payload mass
empty_frac1 = .06; % [N/A] empty mass fraction stage 1
empty_frac2 = .09; % [N/A] empty mass fraction stage 2
empty_frac3 = .115; % [N/A] empty mass fraction stage 3
Isp1  = 289.1; % [s] Stage 1 specific impulse 
Isp2 = 315.4; % [s] Stage 2 specific impulse
Isp3 = 325.1; % [s] Stage 3 specific impulse
g0 = 9.81; % [m/s^2] gravitational acceleration

% characteristic exhaust velocity calculation 
c1 = Isp1 * g0;
c2 = Isp2 * g0;
c3 = Isp3 * g0;

% creating an array of values for delta v fractions
alpha = .2:.01:.6;
[X,Y,Z] = meshgrid(alpha,alpha,alpha);
[X1,Y1] = meshgrid(alpha,alpha);
howdy = (X + Y + Z) == 1;
X = X .* howdy;
Y = Y .* howdy;
Z = Z .* howdy;

% 3rd stage calcs
prop_mass3 = (payload_mass .* (exp(dv_total .* X./c3) - 1) .* (1 - empty_frac3) ./ ...
(1 - empty_frac3 .* exp(dv_total .* X./c3)));
inert_mass3 = -1 .* prop_mass3 .* empty_frac3 ./ (empty_frac3 - 1);
tot_mass3 = (prop_mass3 + inert_mass3 + payload_mass) .* howdy;

 % 2nd stage calcs
prop_mass2 = (tot_mass3 .* (exp(dv_total .* Y/c2) - 1) .* (1 - empty_frac2) ./ ...
(1 - empty_frac2 .* exp(dv_total .* Y./c2)));
inert_mass2 = -1 .* prop_mass2 .* empty_frac2 ./ (empty_frac2 - 1);
tot_mass2 = (prop_mass2 + inert_mass2 + tot_mass3) .* howdy;

% 1st stage calcs
prop_mass1 = (tot_mass2 .* (exp(dv_total .* Z/c1) - 1) .* (1 - empty_frac3) ./ ...
(1 - empty_frac3 .* exp(dv_total .* Z/c1)));
inert_mass1 = -1 .* prop_mass1 .* empty_frac3 ./ (empty_frac3 - 1);
tot_mass1 = (prop_mass1 + inert_mass1 + tot_mass2) .* howdy;

% Reformating total mass values to plot better, removing low and high vals
total_mass = sum(tot_mass1(:,:,:),3);
[GLOW_min,indexmin] = min(total_mass(total_mass>0));
howdy = (total_mass <= 0);
total_mass(howdy) = GLOW_min;
GLOW_mean = mean(mean(total_mass));
howdy = (total_mass >= GLOW_mean);
total_mass(howdy) = GLOW_mean;

% creating mesh 
mesh(X1,Y1,total_mass);

ylabel("2nd Stage \DeltaV fraction")
xlabel("3rd Stage \DeltaV fraction")
zlabel("GLOW (kg)")
title("GLOW/3 Stage Total Mass vs. 3rd Stage \DeltaV fraction")
c = colorbar;
c.Label.String = 'GLOW (kg)';

[GLOW_min1,indexmin] = min(tot_mass1(tot_mass1>0));
frac3 = X1(indexmin);
frac2 = Y1(indexmin);
frac1 = 1 - (frac3 + frac2);


