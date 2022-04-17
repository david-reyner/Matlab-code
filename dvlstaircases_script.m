%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Devil's staircases %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we calculate the devil's staircases emerging when the PING/ING
% oscillator is externally driven by a sinusoidal (type_prt = 1) or a Von
% Mises (type_prt = 2) stimulus. We compute the rotation number of the
% phase equation as a function of the relation between the periods of the
% perturbation and the oscillator. The T/T*-discretization consists of 1500
% points per single period T* (i.e up to T/T* <= 1), each one up to 750
% iterates of the stroboscopic map.

% The next set of parameters must be initialized before running:
type
coord
type_prt
k
A
Tpini
Tpend
Ie_ext

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi or 3 for both)
    type_prt --> Type of perturbation (1 for Sinusoidal, 2 for Von Mises)
    k --> Von Mises parameter (input concentration factor)
    A --> Amplitude of the perturbation
    Tpini --> Beginning of the period discretization (Tpini <= T/T*)
    Tpend --> End of the period discretization (T/T* <= Tpend)
    Ie_ext --> Constant current to exc. neurons
%}

format long;

% Integration tolerances (stroboscopic map)
options_ode = odeset('AbsTol', 1e-12, 'RelTol', 1e-10, 'InitialStep', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% Auxiliar strings within the file names to load
auxstr = ['_Ie_ext', num2str(Ie_ext)];
auxstr1 = ['_Ie_ext', num2str(Ie_ext), '_V'];

% --> Oscillator period
if type == 1
    name_file = ['period_ping', auxstr, '.txt'];
else
    name_file = ['period_ing', auxstr, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'r');
T = fscanf(file, '%f');
fclose(file);

% --> Infinitesimal phase response curve
if type == 1
    name_file = ['iPRC_ping', auxstr1, '.txt'];
else
    name_file = ['iPRC_ing', auxstr1, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'r');
formatSpec = '%f %f %f'; sizeZ = [3 Inf];
res = fscanf(file, formatSpec, sizeZ); res = res';
t = res(:,1); Z = res(:,2:3);
fclose(file);

% --> Infinitesimal phase response curve (at equidistant points)
if type == 1
    name_file = ['equi_iPRC_ping', auxstr1, '.txt'];
else
    name_file = ['equi_iPRC_ing', auxstr1, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'r');
formatSpec = '%f %f %f'; sizeZ = [3 Inf];
res = fscanf(file, formatSpec, sizeZ); res = res';
t_equi = res(:,1); Z_equi = res(:,2:3);
fclose(file);

% Sum of the iPRC-V's (only when the perturbation is applied to both variables)
Z(:,3) = Z(:,1) + Z(:,2);
Z_equi(:,3) = Z_equi(:,1) + Z_equi(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case the original iPRC Z (at non-equidistant points) is used, then
% some kind of interpolation must be considered. Two possible options are
% available: using the Matlab functions interp1q or griddedInterpolant.

% Linear interpolation of the component coord of the iPRC
% lin_intpol = @(l) interp1q(t, Z(:,coord), l); % Quick 1D linear interpolation (not recommended)
% lin_intpol = griddedInterpolant(t, Z(:,coord), 'pchip'); % GriddedInterpolant function

% Periodic perturbation with period Tp
if type_prt == 1
    % Sinusoidal input
    p = @(t,Tp) 1 + cos(2*pi*t/Tp);
else
    % Von Mises input
    p = @(t,Tp) vonmises_dist(t, 0, k, Tp);
end

% Devil's staircases discretization
nPeriod = 1500;
rot_number_sort = zeros(nPeriod,1);
if Tpend > 1
    Tps = linspace(Tpini*T,Tpend*T,nPeriod);
else
    Tps = linspace(0.01,Tpend*T,nPeriod);
end

% Create file where data will be saved
if type == 1
    name_file = ['dvlstaircases_ping_coord_', num2str(coord)];
else
    name_file = ['dvlstaircases_ing_coord_', num2str(coord)];
end
name_file = strcat(name_file, 'eqlinint');
if type_prt == 2
    name_file = strcat(name_file, ['_k', num2str(k)]);
end
name_file = strcat(name_file, ['_A', num2str(A), '_numTp', num2str(nPeriod), ...
    '_upto', num2str(Tpend), 'T', auxstr, '.txt']); % file name
fclose(fopen(fullfile(folder, name_file), 'w+'));

try

for i = 1:length(Tps)
    i
    T_prime = Tps(i); % Period of the perturbation
    
    %%%%%%%%%%%%%%%%%%%%%% Computing rotation number %%%%%%%%%%%%%%%%%%%%%%
    % Initial condition: phase zero
    theta_0 = 0;
    
    % Stroboscopic map: flow of the phase equation after a period of the input

    % First iterate (uncomment next line if original Z is to be used)
%     [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%         lin_intpol, p(l,T_prime), A, T), [0 T_prime], theta_0, options_ode);
    [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
        t_equi, Z_equi(:,coord), p(l,T_prime), A, T), [0 T_prime], theta_0, options_ode);

    % Up to maxiter iterates
    maxiter = 750; % Max. number of iterations
    laps = zeros(maxiter,1); % Number of completed laps
    phases = zeros(maxiter,1); % Phase after each iteration
    laps(1) = floor(theta(end)/T); phases(1) = mod(theta(end),T);
    for iter = 2:maxiter
        % Autonomous phase equation (uncomment next line if original Z is to be used)
%         [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%             lin_intpol, p(l,T_prime), A, T), [0 T_prime], phases(iter-1), options_ode);
        [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
            t_equi, Z_equi(:,coord), p(l,T_prime), A, T), [0 T_prime], phases(iter-1), options_ode);
        laps(iter) = laps(iter-1) + floor(theta(end)/T); % Number of completed laps
        phases(iter) = mod(theta(end),T); % Phase after each iteration
    end
    
    % Alternative way of computing the rotation number
    [~, ind] = sort(phases);
    aux_rot_number1 = zeros(length(phases),1); j1 = 1;
    aux_rot_number2 = zeros(length(phases),1); j2 = 1;
    for j = 1:length(phases)-1
        ind1 = ind(j); ind2 = ind(j+1);
        if ind1 < ind2 % rho min
            aux_rot_number1(j1) = (laps(ind2) - laps(ind1))/(ind2-ind1); j1 = j1 + 1;
        else % rho max
            aux_rot_number2(j2) = (laps(ind1) - laps(ind2))/(ind1-ind2); j2 = j2 + 1;
        end
    end
    aux_rot_number1 = aux_rot_number1(1:j1-1);
    aux_rot_number2 = aux_rot_number2(1:j2-1);
    if isempty(aux_rot_number1) == 0
        rot_number_sort(i) = max(aux_rot_number1);
    else
        rot_number_sort(i) = min(aux_rot_number2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end

% To save data
file = fopen(fullfile(folder, name_file), 'a');
for i = 1:length(Tps)
    fprintf(file, '%16.15f %16.15f\r\n', [(1/T)*Tps(i) rot_number_sort(i)]);
end
fclose(file);

catch ME
    file = fopen('errorFileEixam_DevilStaircases.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%