% -------------------------------------
%  CEA Wrapper for Engine Sizing Code
% Author: Kamon Blong (kblong@purdue.edu)
% First Created: 7/17/2022
% Last Updated: 
% -------------------------------------

function [cstar, isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son] = RunCEA(P_c, P_e, fuel, fuel_temp, oxidizer, oxidizer_temp, OF, sub, sup, file_name)

%{ 
Description: Based off of cea_rocket_example.m and previous code by Tim
    Kayser. Sends engine performance and fuel properties to NASA CEA where
    a bunch of voodoo stoichiometry magic is performed to return relevant
    values for combustion properties of propellants.

Inputs:
- 

Outputs: 
- 
%}

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;

inp('type') = 'eq';                   % Sets the type of CEA calculation
inp('p') = P_c;                       % Chamber pressure
inp('p_unit') = 'psi';                % Chamber pressure units
inp('o/f') = OF;                      % Mixture ratio
if sub(1) || sup(1) ~= 0
    inp('sub') = sub;                     % subsonic area ratios
    inp('sup') = sup;                   % supersonic area ratios
else
    inp('pip') = P_c / P_e;                 % Pressure ratios
end
inp('fuel') = fuel;             % Fuel name from thermo.inp
inp('fuel_t') = fuel_temp;                % Fuel inlet temperature
inp('ox') = oxidizer;               % Ox name from thermo.inp
inp('ox_t') = oxidizer_temp;                  % Ox inlet temperature
inp('file_name') = file_name;   % Input/output file name
if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').
data_eq = data('eq');

% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.
isp = squeeze(data_eq('isp'));
cstar = squeeze(data_eq('cstar'));
exp_ratio = data_eq('ae/at');
M = squeeze(data_eq('mach'));
gamma = squeeze(data_eq('gammas'));
P = squeeze(data_eq('p'));
T = squeeze(data_eq('t'));
rho = squeeze(data_eq('rho'));
mu = squeeze(data_eq('visc'));
Pr = squeeze(data_eq('prandtl'));
Mw = squeeze(data_eq('m'));
k = squeeze(data_eq('k'));
son = squeeze(data_eq('son'));

% convert and correct units
%exp_ratio = 1 / M(end) * ((2 + (gamma(1) - 1 ) * M(end) ^ 2) / (gamma(1) + 1)) ^ ((gamma(1) + 1) / (2 * (gamma(1) - 1))); % calculate expansion ratio manually
P = P.* 1.45038E-4; % convert [Pa] to [psi]
isp = isp(end) / 9.8067; % normalize isp to seconds
exp_ratio = exp_ratio(end); % select expansion ratio
cstar = cstar(1) * 3.281; % convert [m/s] to [ft/s]
T = T .* 1.8; % convert [K] to [R]
rho = rho * 3.613E-5; % convert [kg/m^3] to [lbm/in^3]
mu = mu * 1.450E-4; % convert [Pa-s] to [psi-s]
Mw = Mw * 2.205; % convert [kg/kmol] to [lbm/kmol]
k = k * 0.5782; % convert [W/m-K] to [Btu/hr-ft-R]
end