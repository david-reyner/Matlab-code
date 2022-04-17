%%
%%%%%%%%%%%%%%%%% Entrainment - Two identical competitors %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Rotation number and Sync. index %%%%%%%%%%%%%%%%%%%%%

% Here we determine whether the oscillator is still entrained by the
% initial Von Mises stimulus or, on the contrary, the 1:1 phase-locking is
% broken by the presence of a 2nd identical stimulus. We compute the
% rotation number of the phase equation and the synchronization index as a
% function of the concentration factor k2 and the periods relation T/T*.
% The k2-discretization consists of 90 points from k2 = 0 to k2 = 20 while
% the T/T*-discretization is made up of 8 points along the 1:1 relation.
% For each point (k2, T/T*), 1500 iterates of the stroboscopic map have
% been used.

% The next set of parameters must be initialized before running:
type
coord
k1
A
Ie_ext

% Arguments from command window (only octave)
% args = argv()
% type = args{1};
% coord = args{2};
% k1 = args{3};
% A = args{4};
% Ie_ext = args{5};

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi or 3 for both)
    k1 --> Von Mises parameter (input concentration factor)
    A --> Amplitude of the perturbation
    Ie_ext --> Constant current to exc. neurons
%}

format long;

% Integration tolerances (Stroboscopic map)
options_ode = odeset('AbsTol', 1e-14, 'RelTol', 1e-11, 'InitialStep', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% Auxiliar string for non-zero values of Ie_ext
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
% lin_intpol = griddedInterpolant(t, Z(:,coord)); % GriddedInterpolant function

% Concentration factor of the second input (finer non-equidistant discretization)
k = fliplr(linspace(0,20^(1/3),82));
k = k.^3; k(1) = 20;
k = [k, 10 5 2 1 0.8 0.5 0.25 0.1];
k = sort(k);

% Relation between periods fulfilling a 1:1 resonant condition
numTpT = 8;
if type == 1
    if coord == 1
%         TpT = linspace(0.855, 1, numTpT); % PING-Ve (normalized von mises)
%         TpT = linspace(0.71, 1, numTpT); % PING-Ve
%         TpT = linspace(0.8, 1.01, numTpT); % for k1 = 20
        TpT = [0.85, 0.865, 0.88, 0.915, 0.93, 0.95, 0.97, 0.99]; % for k1 = 2
    elseif coord == 2
%         TpT = linspace(0.965, 1, numTpT); % PING-Vi (normalized von mises)
        TpT = linspace(0.91, 1.05, numTpT); % PING-Vi
    else
        TpT = [0.845 0.86 0.89 0.91 0.93 0.95 0.975 1]; % PING-VeVi (for k1 = 2)
    end
else
%     TpT = linspace(0.932, 1, numTpT); % ING-Vi (normalized von mises)
    TpT = linspace(0.92, 1, numTpT); % ING-Vi
end

% Create file where data will be saved
if type == 1
    name_file = ['rotnumber_sync_2inputs_ping_coord_', num2str(coord)];
else
    name_file = ['rotnumber_sync_2inputs_ing_coord_', num2str(coord)];
end
name_file = strcat(name_file, ['eqlinint_k1', num2str(k1)]);
if A ~= 1
    name_file = strcat(name_file, ['_A', num2str(A)]);
end
name_file = strcat(name_file, ['_k', num2str(length(k)), ...
    '_numTpT', num2str(numTpT), auxstr, '.txt']); % filename
fclose(fopen(fullfile(folder, name_file), 'w+'));

try

for ii = 1:length(TpT)
    
    p = TpT(ii) % Relation between periods of the 1st input and the oscillator

    % Periodic perturbation with period Tp (Von Mises distribution)
    Tp = p*T; % Fulfilling a 1:1 Phase-locking condition
    pvm1 = @(t) vonmises_dist(t, 0, k1, Tp);
    pvm2 = @(t,k2) vonmises_dist(t, Tp/2, k2, Tp);
    
    % Rotation number (Entrainment)
    rot_number_sort = zeros(length(k),1);
    
    % Synchronization index
    r = zeros(length(k),1);
    mean_phase = zeros(length(k),1);

    for i = 1:length(k)
        k2 = k(i);
        
        %%%%%%%%%%%%%%%%%%%% Computing rotation number %%%%%%%%%%%%%%%%%%%%
        % Initial condition: phase zero
        theta_0 = 0;

        % Stroboscopic map: flow of the phase equation up to period of the forcer

        % First iterate (uncomment next line if original Z is to be used)
%         [l, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%             lin_intpol, pvm1(l)+pvm2(l,k2), A, T), [0 Tp], theta_0, options_ode);
        [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
            t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l,k2), A, T), [0 Tp], theta_0, options_ode);

        % Up to maxiter iterates
        maxiter = 1500; % Max. number of iterations
        laps = zeros(maxiter,1); % Number of completed laps
        phases = zeros(maxiter,1); % Phase after each iteration
        laps(1) = floor(theta(end)/T); phases(1) = mod(theta(end),T);
        for iter = 2:maxiter
            % Autonomous phase equation (uncomment next line if original Z is to be used)
%             [l, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%                 lin_intpol, pvm1(l)+pvm2(l,k2), A, T), [0 Tp], phases(iter-1), options_ode);
            [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
                t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l,k2), A, T), [0 Tp], phases(iter-1), options_ode);
            laps(iter) = laps(iter-1) + floor(theta(end)/T); % Number of completed laps
            phases(iter) = mod(theta(end),T); % Phase after each iteration     
        end
        
        % Measuring the degree of phase-locking (Synchronization index)
        mean_vector = [sum(cos(2*pi*phases./T)); sum(sin(2*pi*phases./T))];
        mean_phase(i) = atan2(mean_vector(2), mean_vector(1)); % Mean phase
        if mean_phase(i) < 0
            mean_phase(i) = mean_phase(i) + 2*pi;
        end
        r(i) = 1/maxiter*sqrt(mean_vector(1)^2 + mean_vector(2)^2); % Sync. index
        
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % To save data
    file = fopen(fullfile(folder, name_file), 'a');
    for i = 1:length(k)
        fprintf(file, '%16.15f %16.15f %16.15f %16.15f %16.15f\r\n', ...
            [p k(i) rot_number_sort(i) r(i) mean_phase(i)]);
    end
    fclose(file);
end

catch ME
    file = fopen('errorFileEixam_rotnumber_sync_2inputs.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%