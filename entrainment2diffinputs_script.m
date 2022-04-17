%%
%%%%%%%%%%%%%%% Entrainment - Two non-identical competitors %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Synchronization index %%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we determine whether the oscillator is still synchronized by the
% initial Von Mises stimulus or, on the contrary, the 1:1 phase-locking is
% broken by the presence of a 2nd non-identical stimulus. We compute only 
% the synchronization index as a function of the concentration factor k2 
% and the periods relation between both inputs, i.e. T2/T1. We fix the
% periods relation T/T* at three different points: left, middle and right
% of the 1:1 plateau. The k2-discretization consists of 9 points from 
% k2 = 0.1 to k2 = 20. For each point (k2, T2/T1), 1500 iterates of the 
% stroboscopic map have been used.

% The next set of parameters must be initialized before running:
type
coord
k1
A
p
strbmap
Ie_ext

% Arguments from command window (only octave)
% args = argv()
% type = args{1};
% coord = args{2};
% k1 = args{3};
% A = args{4};
% p = args{5};
% strbmap = args{6};
% Ie_ext = args{7};

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi or 3 for both)
    k1 --> Von Mises parameter (input concentration factor)
    A --> Amplitude of the perturbation
    p --> Ratio between periods of the 1st input and oscillator
    strbmap --> Time of stroboscopic map (1 for up to T1, 2 for up to T2)
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
lin_intpol = griddedInterpolant(t, Z(:,coord), 'pchip'); % GriddedInterpolant function

% Concentration factor of the second input
k = [20 10 5 2 1 0.8 0.5 0.25 0.1];

% A finer non-equidistant k2-discretization (more points in the [0,1] interval)
% k = fliplr(linspace(0,20^(1/3),82));
% k = k.^3; k(1) = 20;
% k = [k, 10 5 2 1 0.8 0.5 0.25 0.1];
% k = sort(k);

% Relation between input periods
T1 = p*T; % First input oscillates at p times the frequency of the oscillator
numT2 = 10; % Number of points in the T2-discretization
T2vec_dec = linspace(T1/2, T1, numT2/2); % T1/2 <= T2
T2vec_inc = linspace(T1, 3*T1/2, numT2/2 + 1); % T2 <= 3*T1/2
T2vec = [fliplr(T2vec_dec) T2vec_inc(2:end)];

% Create file where data will be saved
if type == 1
    name_file = ['syncindex_2diffinputs_ping_coord_', num2str(coord)];
else
    name_file = ['syncindex_2diffinputs_ing_coord_', num2str(coord)];
end
name_file = strcat(name_file, ['eqlinint_k1', num2str(k1)]);
if A ~= 1
    name_file = strcat(name_file, ['_A', num2str(A)]);
end
name_file = strcat(name_file, ['_p', num2str(p), '_k', num2str(length(k)), ...
    '_numT2', num2str(numT2)]);
if strbmap == 2
    name_file = strcat(name_file, '_strbT2');
end
name_file = strcat(name_file, [auxstr, '.txt']); % filename
fclose(fopen(fullfile(folder, name_file), 'w+'));

try
    
for ii = 1:length(T2vec)
    
    T2 = T2vec(ii) % Second input period
    
    % Phase-shift
    mu1 = 0;
    if T1 > T2
        mu2 = T2/2;
    else
        mu2 = T1/2;
    end

    % Periodic perturbation (Von Mises distribution)
    pvm1 = @(t) vonmises_dist(t, mu1, k1, T1); % 1st input (k1 and T1 on the 1:1 plateau)
    pvm2 = @(t,k2) vonmises_dist(t, mu2, k2, T2); % 2nd input
    
    % Synchronization index
    r = zeros(length(k),1);
    mean_phase = zeros(length(k),1);
    
    % Defining stroboscopic map up to T1 or to T2
    if strbmap == 1
        Ts = T1;
    else
        Ts = T2;
    end

    for i = 1:length(k)
        
        k2 = k(i);
        
        %%%%%%%%%%%%%%%%%%%% Computing rotation number %%%%%%%%%%%%%%%%%%%%
        % Initial condition: phase zero
        theta_0 = 0;
        
        % Stroboscopic map: flow of the phase equation up to period of the forcer
        
        % First iterate (uncomment next line if Z at equidistant points is to be used)
        [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
            lin_intpol, pvm1(l)+pvm2(l,k2), A, T), [0 Ts], theta_0, options_ode);
%         [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
%             t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l,k2), A, T), [0 Ts], theta_0, options_ode);
        
        % Up to maxiter iterates
        maxiter = 1500; % Max. number of iterations
        laps = zeros(maxiter,1); % Number of completed laps
        phases = zeros(maxiter,1); % Phase after each iteration
        laps(1) = floor(theta(end)/T); phases(1) = mod(theta(end),T);
        for iter = 2:maxiter
            % Non-autonomous phase equation (uncomment next line if Z at equidistant points is to be used)
            [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
                lin_intpol, pvm1(l)+pvm2(l,k2), A, T), [(iter-1)*Ts iter*Ts], phases(iter-1), options_ode);
%             [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
%                 t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l,k2), A, T), [(iter-1)*Ts iter*Ts], phases(iter-1), options_ode);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % To save data
    file = fopen(fullfile(folder, name_file), 'a');
    for i = 1:length(k)
        fprintf(file, '%16.15f %16.15f %16.15f %16.15f\r\n', ...
            [T2/T1 k(i) r(i) mean_phase(i)]);
    end
    fclose(file);
end

catch ME
    file = fopen('errorFileEixam_syncindex_2diffinputs.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%