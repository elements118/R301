function [mach_num, density, temperature, pressure] = atmosphere(altitude,units,increment)
% This function takes inputs for altitude and desired unit system, as well as altitude increment, and
% returns atmospheric conditions (pressure, density, temperature, pressure,
% mach number. Each returned value will be a 1xm array for each altitude increment from sea level to the maximum input altitude.
% Units = 1 is metric, 2 is english

% Constants and surface values
earth_rad = 6371000; % [m] radius of the earth
g = 9.8; % [m/s^2] acceleration due to gravity on surface
p_surf = 1.01*10^5; % [N/m^2] pressure on surface
rho_surf = 1.2250; % [kg/m^3] density on surface
t_surf = 288.16; % [K] temperature on surface
R = 287; % [] gas constant
gamma = 1.4;% [N/A] ratio of specific heats

% conversion factors
m2f = 3.28084; % [N/A] meters to feet conversion factor
km2f = m2f * 1000; % [N/A] kilometers to feet conversion factor
kpm2rpf = .54864; % [N/A] K/m to R/ft conversion factor
k2r = 1.8; % [N/A] degrees kelvin to degrees rankine
npm2lbpf = .020885; % [N/A] N/m^2 to lbm/ft^2
kgpm2slgpf = 0.00194032; % [N/A] kg/m^3 to slug/ft^3
mps2fps = 3.28084; % [N/A] m/s to f/s 
ft2m = .3048; % [N/A] feet to meters

if units == 2
    altitude = altitude * ft2m; 
end
if (altitude > 100000) 
    altitude = 100000; 
end


% Creating matrices
geom_alt = 0:increment:100000;
outputs = zeros(7,length(geom_alt));
outputs(1,:) = geom_alt;
outputs(2,:) = (earth_rad ./ (earth_rad + ...
    geom_alt)) .* geom_alt;

% predetermined a values
a1 = -6.5 * 10^-3; % [K/m] 
a2 = 3 * 10^-3; % [K/m]
a3 = -4.5 * 10^-3; % [K/m]
a4 = 4 * 10^-3; % [K/m]

% predetermined altitude values 
h1 = 11 * 1000; % [m]
h2 = 25 * 1000; % [m]
h3 = 47 * 1000; % [m]
h4 = 53 * 1000; % [m]
h5 = 79 * 1000; % [m]
h6 = 90 * 1000; % [m]
h7 = 105 * 1000; % [m]

% predetermined temperature values
t1 = 216.66; % [K]
t2 = 282.66; % [K]
t3 = 165.66; % [K]
t4 = 225.66; % [K]

% Time to loop!
    outputs(2,:) = (earth_rad ./ (earth_rad + ...
    geom_alt(:))) .* geom_alt(:); % [m] geopotential altitude
    outputs(3,:) = t_surf + a1 .* (outputs(1,:));
    outputs(4,:) = p_surf .* (outputs(3,:) ./ t_surf) .^ ((-1 * g) / (a1 * R));
    outputs(5,:) = rho_surf .* (outputs(3,:) ./ t_surf) .^ (-1 * ((g / (a1 * R)) + 1));

% base constants first isothermal
p1 = p_surf * (t1 / t_surf) ^ ((-1 * g) / (a1 * R));
rho1 = rho_surf * (t1 / t_surf) ^ (-1 * ((g / (a1 * R)) + 1));

% First isothermal region
tempmat = outputs(1,:) >= h1;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t1;
    outputs(4,tempmat) = p1 .* exp(-1 .* (g / (R * t1)) .* (outputs(2,tempmat) - h1));
    outputs(5,tempmat) = rho1 .* exp(-1 .* (g / (R * t1)) .* (outputs(2,tempmat) - h1));

% base constants 2nd gradient
p2 = p1 * exp(-1 * ((g / (R * t1)) * (h2 - h1)));
rho2 = rho1 * exp(-1 * (g / (R * t1)) * (h2 - h1));

% 2nd gradient region
tempmat = outputs(1,:) >= h2;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t1 + a2 .* (outputs(1,tempmat) - h2);
    outputs(4,tempmat) = p2 .* (outputs(3,tempmat) ./ t1) .^ ((-1 * g) / (a2 * R));
    outputs(5,tempmat) = rho2 .* (outputs(3,tempmat) ./ t1) .^ (-1 * ((g / (a2 * R)) + 1));

% 2nd isothermal region base
p3 = p2 * (t2 / t1) ^ ((-1 * g) / (a2 * R));
rho3 = rho2 * (t2 / t1) ^ (-1 * ((g / (a2 * R)) + 1));

% 2nd isothermal region loop
tempmat = outputs(1,:) >= h3;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t2;
    outputs(4,tempmat) = p3 .* exp(-1 * (g / (R * t2)) * (h4-h3));
    outputs(5,tempmat) = rho3 .* exp(-1 * (g / (R * t2)) * (h4-h3));

% base constants 3rd gradient
p4 = p3 * exp(-1 * ((g / (R * t2)) * (h4-h3)));
rho4 = rho3 * exp(-1 * (g / (R * t2)) * (h4-h3));

% 3rd gradient region
tempmat = outputs(1,:) >= h4;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t2 + a3 .* (outputs(1,tempmat) - h4);
    outputs(4,tempmat) = p4 .* (outputs(3,tempmat) ./ t2) .^ ((-1 * g) / (a3 * R));
    outputs(5,tempmat) = rho4 .* (outputs(3,tempmat) ./ t2) .^ (-1 * ((g / (a3 * R)) + 1));

% base constants third isothermal
p5 = p4 * (t3 / t2) ^ ((-1 * g) / (a3 * R));
rho5 = rho4 * (t3 / t2) ^ (-1 * ((g / (a3 * R)) + 1));

% Third isothermal region
tempmat = outputs(1,:) >= h5;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t3;
    outputs(4,tempmat) = p5 .* exp(-1 * (g / (R * t3)) .* (outputs(2,tempmat) - h6));
    outputs(5,tempmat) = rho5 .* exp(-1 * (g / (R * t3)) .* (outputs(2,tempmat) - h6));

% base constants 4th gradient
p6 = p5 * exp(-1 * ((g / (R * t3)) * (h7-h6)));
rho6 = rho5 * exp(-1 * (g / (R * t3)) * (h7-h6));

% 4th gradient region
tempmat = outputs(1,:) >= h6;
    outputs(2,tempmat) = (earth_rad ./ (earth_rad + ...
    geom_alt(tempmat))) .* geom_alt(tempmat); % [ft] geopotential altitude
    outputs(3,tempmat) = t4 + a4 .* (outputs(1,tempmat) - h7);
    outputs(4,tempmat) = p6 .* (outputs(3,tempmat) / t4) .^ ((-1 * g) / (a4 * R));
    outputs(5,tempmat) = rho6 .* (outputs(3,tempmat) / t4) .^ (-1 * ((g / (a4 * R)) + 1));

% Speed of sound 
outputs(7,:) = sqrt(gamma * R * outputs(3,:));

% Create matrix of english values
eng_outputs = repmat(outputs,1);
eng_outputs(1:2,:) = eng_outputs(1:2,:) * m2f;
eng_outputs(3,:) = eng_outputs(3,:) * k2r;
eng_outputs(4,:) = eng_outputs(4,:) * npm2lbpf;
eng_outputs(5,:) = eng_outputs(5,:) * kgpm2slgpf;
eng_outputs(7,:) = eng_outputs(7,:) * mps2fps;

[~,indexMax] = min(abs(geom_alt - altitude));

if units == 2
     mach_num = eng_outputs(7,1:indexMax);
     density = eng_outputs(5,1:indexMax);
     temperature = eng_outputs(3,1:indexMax);
     pressure = eng_outputs(4,1:indexMax);
elseif units == 1
     mach_num = outputs(7,1:indexMax);
     density = outputs(5,1:indexMax);
     temperature = outputs(3,1:indexMax);
     pressure = outputs(4,1:indexMax);
end

