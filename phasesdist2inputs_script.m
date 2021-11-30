%%
%%%%%%%%%% Phases distribution - Two (non)-identical competitors %%%%%%%%%%

% Here we calculate the phases of the (T1 or T2)-stroboscopic map resulting
% from the interference of a second (non)-identical stimulus. Given the
% periods ratios p and q (T1/T* and T2/T1, resp.) and the concentration
% factors of each input (k1 and k2), we iterate the stroboscopic map up to
% 10000 times and save the resulting phase distribution. Parameters p and 
% k1 are such that there is a 1:1 phase-locking between the first input
% and the oscillator.

type
coord
A
k1
p
k2
q
strbmap
Ie_ext

% Arguments from command window (only octave)
% args = argv()
% type = args{1};
% coord = args{2};
% A = args{3};
% k1 = args{4};
% p = args{5};
% k2 = args{6};
% q = args{7};
% strbmap = args{8};
% Ie_ext = args{9};

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi)
    A --> Amplitude of the perturbation
    k1 --> Von Mises parameter (input concentration factor)
    p --> Ratio between periods of the 1st input and oscillator
    k2 --> Von Mises parameter (second input concentration factor)
    q --> Ratio between periods of the 2nd and 1st inputs
    strbmap --> Time of stroboscopic map (1 for up to T1, 2 for up to T2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case the original iPRC Z (at non-equidistant points) is used, then
% some kind of interpolation must be considered. Two possible options are
% available: using the Matlab functions interp1q or griddedInterpolant.

% Linear interpolation of the component coord of the iPRC
% lin_intpol = @(l) interp1q(t, Z(:,coord), l); % Quick 1D linear interpolation (not recommended)
% lin_intpol = griddedInterpolant(t, Z(:,coord), 'pchip'); % GriddedInterpolant function

% Relation between the periods of the inputs
T1 = p*T; % First input oscillates at p times the frequency of the oscillator
T2 = q*T1; % Second input oscillates at q times the frequency of the 1st one

% Determine phase-shift of the second perturbation
mu1 = 0;
if T1 > T2
    mu2 = T2/2;
else
    mu2 = T1/2;
end

% Periodic perturbation (Von Mises distribution)
pvm1 = @(t) vonmises_dist(t, mu1, k1, T1); % 1st input (k1 and T1 on the 1:1 plateau)
pvm2 = @(t) vonmises_dist(t, mu2, k2, T2); % 2nd input

% Defining stroboscopic map up to T1 or to T2
if strbmap == 1
    Ts = T1;
else
    Ts = T2;
end

% Create file where data will be saved
if type == 1
    name_file = ['phasesdist_2inputs_ping_coord_', num2str(coord)];
else
    name_file = ['phasesdist_2inputs_ing_coord_', num2str(coord)];
end
name_file = strcat(name_file, ['eqlinint_k1', num2str(k1), ...
    '_A', num2str(A), '_p', num2str(p), '_k2', num2str(k2), ...
    '_q', num2str(q)]);
if strbmap == 1
    name_file = strcat(name_file, '_strbT1');
else
    name_file = strcat(name_file, '_strbT2');
end
name_file = strcat(name_file, auxstr, '.txt'); % file name
fclose(fopen(fullfile(folder, name_file), 'w+'));

try

%%%%%%%%%%%%%%%%% Computing phases of the stroboscopic map %%%%%%%%%%%%%%%%
% Initial condition: phase zero
theta_0 = 0;

% Stroboscopic map: flow of the phase equation after a period of the input

% First iterate (uncomment next line if original Z is to be used)
% [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%     lin_intpol, pvm1(l)+pvm2(l), A, T), [0 Ts], theta_0, options_ode);
[~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
    t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l), A, T), [0 Ts], theta_0, options_ode);

% Up to maxiter iterates
maxiter = 10000; % Max. number of iterations
laps = zeros(maxiter,1); % Number of completed laps
phases = zeros(maxiter,1); % Phase after each iteration
laps(1) = floor(theta(end)/T); phases(1) = mod(theta(end),T);
for iter = 2:maxiter
    % Non-autonomous phase equation (uncomment next line if original Z is to be used)
%     [~, theta] = ode45(@(l, theta) phase_equation(l, theta, ...
%         lin_intpol, pvm1(l)+pvm2(l), A, T), [(iter-1)*Ts iter*Ts], phases(iter-1), options_ode);
    [~, theta] = ode45(@(l, theta) phase_equation_eqlinint(l, theta, ...
        t_equi, Z_equi(:,coord), pvm1(l)+pvm2(l), A, T), [(iter-1)*Ts iter*Ts], phases(iter-1), options_ode);
    laps(iter) = laps(iter-1) + floor(theta(end)/T); % Number of completed laps
    phases(iter) = mod(theta(end),T); % Phase after each iteration
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To save data
file = fopen(fullfile(folder, name_file), 'a');
for i = 1:length(phases)
    fprintf(file, '%d %16.15f\r\n', [i phases(i)]);
end
fclose(file);

catch ME
    file = fopen('errorFileEixam_phases_dist_2diffinputs.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%