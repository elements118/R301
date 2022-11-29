% AAE 251 Fall 2022
% PM7 Team R301 
% Rocket staging optimization
% Authors: Jacob Bell
% Collaborators: Caroline Gee, Sydney Brown, Thomas Brannon, Justin
% Armstrong, Lu Larest, Carter Doby, Abby Woodburry 
% This code plots various Delta-V fractions and their respective GLOW
% values for a theoretical 3 stage launch vehicle. The code applies for any
% 3 stage dv split on a vehicle, however, and "GLOW" is simply the mass of
% all three stages
clc;
clear all; 

% inputs/known variables 
dv_total = 15 * 10^3; % [m/s] total delta-v required 
payload_mass = 100; % [kg] payload mass
empty_frac = .115; % [N/A] empty mass fraction 
Isp1  = 362; % [s] Stage 1 specific impulse 
Isp2 = 465; % [s] Stage 2 specific impulse
Isp3 = 350; % [s] Stage 3 specific impulse
g0 = 9.81; % [m/s^2] gravitational acceleration

% characteristic exhaust velocity calculation 
c1 = Isp1 * g0;
c2 = Isp2 * g0;
c3 = Isp3 * g0;

% creating an array of values for delta v fractions
alpha = .2:.01:.6;

GLOW = [];
lcv = 1;
for i = alpha
    for j = alpha
        for k = alpha
            if (i + j + k) == 1
                % 3rd stage calcs
                prop_mass3 = (payload_mass * (exp(dv_total * i/c3) - 1) * (1 - empty_frac) / ...
                    (1 - empty_frac * exp(dv_total * i/c3)));
                inert_mass3 = -1 * prop_mass3 * empty_frac / (empty_frac - 1);
                tot_mass3 = prop_mass3 + inert_mass3 + payload_mass;

                % 2nd stage calcs
                prop_mass2 = (tot_mass3 * (exp(dv_total * j/c2) - 1) * (1 - empty_frac) / ...
                    (1 - empty_frac * exp(dv_total * j/c2)));
                inert_mass2 = -1 * prop_mass2 * empty_frac / (empty_frac - 1);
                tot_mass2 = prop_mass2 + inert_mass2 + tot_mass3;

                % 1st stage calcs
                prop_mass1 = (tot_mass2 * (exp(dv_total * k/c1) - 1) * (1 - empty_frac) / ...
                    (1 - empty_frac * exp(dv_total * k/c1)));
                inert_mass1 = -1 * prop_mass1 * empty_frac / (empty_frac - 1);
                tot_mass1 = prop_mass1 + inert_mass1 + tot_mass2;
                GLOW(lcv) = tot_mass1;
                fracs(1,lcv) = i;
                fracs(2,lcv) = j;
                fracs(3,lcv) = k;
                lcv = lcv + 1;
            end
        end
    end
end

[GLOW_min,i] = min(GLOW(GLOW>0));
stage1_dv = fracs(1,i);
stage2_dv = fracs(2,i);
stage3_dv = fracs(3,i);

figure
plot3(fracs(1,:),fracs(2,:),GLOW,'.');
hold on 
plot3(stage3_dv,stage2_dv,GLOW_min,"Marker","pentagram","markersize",12,"MarkerFaceColor","r")
ylabel("2nd Stage \DeltaV fraction")
xlabel("3rd Stage \DeltaV fraction")
zlabel("GLOW (kg)")
title("GLOW/3 Stage Total Mass vs. 3rd Stage \DeltaV Fraction and 2nd Stage \DeltaV Fraction")
grid
legend("GLOW (kg)", "Minimum GLOW")

