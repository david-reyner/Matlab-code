%%
%%%%%%%%%%%%%%%%%%%%%%%%% Adjoint method and iPRC %%%%%%%%%%%%%%%%%%%%%%%%%

% Here we compute the infinitesimal Phase Response Curves using backward
% integration of the Adjoint method. Ordinary backward integration of the 
% adjoint equation together with the original system is not possible since 
% the limit cycle becomes unstable. Therefore, to integrate the adjoint
% equation, some numerical interpolation of the Jacobian matrix (and thus 
% also of the variables along the periodic orbit) is required. Given the 
% external currents to the excitatory and inhibitory populations, this
% script first computes the emerging periodic orbit and its period, and
% second, the infinitesimal Phase Response Curve, iPRC, (i.e. a periodic 
% solution of the adjoint equation) and its first and second derivatives 
% (using a fast Fourier transform algorithm).

% Arguments from command window (only octave)
% args = argv()
% Ie_ext = args{1};
% Ii_ext = args{2};

Ie_ext
Ii_ext

% Script arguments:
%{
    Ie_ext --> Constant current to exc. neurons
    Ii_ext --> Constant current to inh. neurons
%}

format long;

options_ode = odeset('AbsTol', 1e-18, 'RelTol', 5e-14);

%%%%%%%%%%%%%%%% 1ST STEP: Computation of a periodic orbit %%%%%%%%%%%%%%%%

% Initial/Final integration time
t0 = 0; tmax = 500;

% Initial condition
x00 = [0.1; -1.5; 0.5; 0.5;
       0.1; -1.5; 0.5; 0.5];
   
% PING VS ING
type = 1; % (1 PING, 2 ING)

if type == 1
    % Parameters (PING)
    tau_e = 8; tau_i = 8; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 5;
    Jee = 0; Jei = 13; Jie = 13; Jii = 0;
else
    % Parameters (ING)
    tau_e = 10; tau_i = 10; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 1;
    Jee = 0; Jei = 0; Jie = 0; Jii = 15; % <-----------------------------------------------------------------
    
    % External excitatory and inhibitory inputs
%     Ie_ext = 25; Ii_ext = 25;
end
parm = [tau_e; tau_i; delta_e; delta_i; eta_e; eta_i;
        tau_se; tau_si; Jee; Jei; Jie; Jii; Ie_ext; Ii_ext];

x0 = x00; % initial condition
err = 1; % initial error
while err >= 1e-12 % until difference between first and last integration points is small
    % Oscillatory behaviour arises from the parameters choice
    [~, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
        [t0 tmax], x0, options_ode);

    % Initial condition lying on a periodic orbit (PO)
    x0 = x(end,:);

    %-----------> Computing the period of the periodic orbit <------------%
    % Poincare Section (PS): V = 0.5*(min(Ve) + max(Ve))
    Vps = round((min(x(:,2)) + max(x(:,2)))/2, 1);
    if type == 1
        coord = 2;
    else
        coord = 6;
    end

    % Remark: In general initial condition x0 does not lie on the Poincare 
    % Section (PS). Thus, we perform a preliminary integration to obtain an
    % initial condition on the PS and also on the PO

    % EVENTS handle function: solves and finds where functions of (t,x) are 0
    myfunc = @(t, x) periodPO(t, x, Vps, coord);
    options_ode2 = odeset('AbsTol', 1e-18, 'RelTol', 5e-14, 'Events', myfunc);

    [~, ~, ~, xe, ~] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, ...
        x, parm), [t0 tmax], x0, options_ode2);

    % Initial Condition on the Poincare Section PS
    x0 = xe'; y0 = x0;

    % Computing the period of the periodic orbit
    n_crossing = 2; ind = 1; T = 0;
    while ind <= n_crossing
        [~, y, te, ~, ~] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, ...
            x, parm), [t0 tmax], y0, options_ode2);
        y0 = y(end,:);
        T = T + te;
        if abs(te) >= 1e-12
            ind = ind + 1;
        end
    end

    % Checkig period
    [~, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
        [t0 T], x0, options_ode);
    err = norm(x0(coord-1:end) - x(end,coord-1:end)') % update error
    norm(x0(coord-1:end) - x(end,coord-1:end)') < 5e-14 % 13 significant digits
    norm(x0(coord-1:end) - x(end,coord-1:end)') < 5e-13 % 12 significant digits
    %---------------------------------------------------------------------%
end

% Setting phase zero at the voltage peak
[~, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
    [t0 T], x0, options_ode);
[~, ind] = max(x(:,coord));
x0 = x(ind,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try

%%%%%%%%%%%%%%%%%%% 2ND STEP: Adjoint Equation and iPRC %%%%%%%%%%%%%%%%%%%

% Interpolating variables along a cycle
[t, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
    [t0 T], x0, options_ode);

sample = {t, 1:size(x,2)};
Fint = griddedInterpolant(sample, x, 'spline'); % cubic spline interpolation

% Initial condition (normalized so that Z·F = 1)
F0 = full_synaptic_firing_rate_equations(0, x0, parm);
zaux1 = (1 - F0(2) - F0(5)/2 - F0(6))/F0(1);
z0 = [zaux1 1 0 0 0.5 1 0 0]';
dot(z0,F0)

%-----> Convergence to the infinitesimal Phase Response Curve (iPRC) <----%
conv = false; % initializing convergence indicator
count = 0; % counting number of iterations to converge
inf_loop = false; % detecting infinite loop (periodic orbit does no longer exist)
phase0 = z0; % initial condition to find the periodic solution
t0 = 0; % initial time
first = true;
abserr = 1;
while ~conv && ~inf_loop
    % Transition to the T-periodic solution of the Adjoint equation
    [t, x] = ode45(@(t, x) adjoint_method_interp(t, x, parm, T, Fint), [t0 t0-T], ...
       phase0, options_ode);
   
    % Absolute difference of errors between two consecutive iterations
    if first
        first = false;
        resn = norm(x(end,:)' - x(1,:)');
    else
        resn1 = norm(x(end,:)' - x(1,:)');
        abserr = abs(resn1 - resn)
        resn = resn1;
    end

    res = norm(x(end,:)' - x(1,:)')
    count = count + 1

    t0 = t(end); % initial time for the next iteration
    phase0 = x(end,:)'; % initial condition for the next iteration
    inf_loop = count >= 100 && res > 1e-3; % detection of infinite loop
    conv = abserr < 5e-11 || res < 5e-10; % checking convergence
end
%-------------------------------------------------------------------------%

% Infinitesimal phase response curve (iPRC)
Z = flipud(x(:,:));

% Extracting the iPRC at equidistant points
t_aux = flipud(t + count*T);
t_equi = linspace(0, T, length(t_aux));
sample = {t_aux, 1:size(Z,2)};
Zequi = griddedInterpolant(sample, Z);

% Z1_equi = interp1(t_aux, Z(:, 1), t_equi);
% Z2_equi = interp1(t_aux, Z(:, 2), t_equi);
% Z3_equi = interp1(t_aux, Z(:, 3), t_equi);
% Z4_equi = interp1(t_aux, Z(:, 4), t_equi);
% Z5_equi = interp1(t_aux, Z(:, 5), t_equi);
% Z6_equi = interp1(t_aux, Z(:, 6), t_equi);
% Z7_equi = interp1(t_aux, Z(:, 7), t_equi);
% Z8_equi = interp1(t_aux, Z(:, 8), t_equi);

% Derivatives of the infinitesimal PRC (using fast Fourier transform)

% Linear interpolation of the components Ve and Vi of the iPRC
PRCVe = griddedInterpolant(t_aux, Z(:,2)); % griddedInterpolant function
PRCVi = griddedInterpolant(t_aux, Z(:,6)); % griddedInterpolant function

N = 600; % number of points
dd = linspace(0,T,N+1); aa = dd(end); dd(end) = [];
iiVe = PRCVe(dd')'; iiVi = PRCVi(dd')';
cVe = fft(iiVe); cVi = fft(iiVi); % coeffs. discrete fourier transform (fft algo.)
kk = (1/T)*[0:N/2-1, 0, -N/2+1:-1];
dcVe = 2*pi*1i*kk.*cVe; dcVi = 2*pi*1i*kk.*cVi; % coeffs. 1st derivative of DFT
ddcVe = 2*pi*1i*kk.*dcVe; ddcVi = 2*pi*1i*kk.*dcVi; % coeffs. 2nd derivative of DFT
dd(end+1) = aa;

% Getting the derivatives of the iPRCs and imposing the periodicity condition
dZVe = ifft(dcVe); dZVe(end+1) = dZVe(1);
dZVi = ifft(dcVi); dZVi(end+1) = dZVi(1);
ddZVe = ifft(ddcVe); ddZVe(end+1) = ddZVe(1);
ddZVi = ifft(ddcVi); ddZVi(end+1) = ddZVi(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Saving data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% Auxiliar name for non-default values of I_ext
auxstr = ['_Ie_ext', num2str(Ie_ext)];
if Ii_ext ~= 0
    auxstr = strcat(auxstr, ['_Ii_ext', num2str(Ii_ext)]);
end

% Initial condition
if type == 1
    name_file = ['initial_condition_ping', auxstr, '.txt'];
else
    name_file = ['initial_condition_ing', auxstr, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'w+');
for i = 1:length(x0)
    fprintf(file, '%16.15f\r\n', x0(i));
end
fclose(file);

% Period
if type == 1
    name_file = ['period_ping', auxstr, '.txt'];
else
    name_file = ['period_ing', auxstr, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'w+');
fprintf(file, '%16.15f\n', T);
fclose(file);

% Infinitesimal PRC
str = {'r', 'V', 'Se', 'Si'};
for ii = 1:length(str)
    if type == 1
        name_file = ['iPRC_ping', auxstr];
    else
        name_file = ['iPRC_ing', auxstr];
    end
    name_file = strcat(name_file, ['_', str{ii}, '.txt']);
    fclose(fopen(fullfile(folder, name_file), 'w+'));
    file = fopen(fullfile(folder, name_file), 'a');
    for i = 1:length(Z)
        fprintf(file, '%16.15f %16.15f %16.15f\r\n', [t_aux(i) Z(i,ii) Z(i,ii+4)]);
    end
    fclose(file);
end

% Infinitesimal PRC (at equidistant points)
for ii = 1:length(str)
    if type == 1
        name_file = ['equi_iPRC_ping', auxstr];
    else
        name_file = ['equi_iPRC_ing', auxstr];
    end
    name_file = strcat(name_file, ['_', str{ii}, '.txt']);
    fclose(fopen(fullfile(folder, name_file), 'w+'));
    file = fopen(fullfile(folder, name_file), 'a');
    Z_equi = Zequi({t_equi,[ii,ii+4]});
    for i = 1:length(t_equi)
        fprintf(file, '%16.15f %16.15f %16.15f\r\n', [t_equi(i) Z_equi(i,1) Z_equi(i,2)]);
    end
    fclose(file);
end

% Derivatives of the iPRCs
if type == 1
    name_file_d1 = ['dPRC_ping', auxstr, '.txt'];
    name_file_d2 = ['ddPRC_ping', auxstr, '.txt'];
else
    name_file_d1 = ['dPRC_ing', auxstr, '.txt'];
    name_file_d2 = ['ddPRC_ing', auxstr, '.txt'];
end
fclose(fopen(fullfile(folder, name_file_d1), 'w+'));
fclose(fopen(fullfile(folder, name_file_d2), 'w+'));

% 1st Derivatives
file = fopen(fullfile(folder, name_file_d1), 'a');
for i = 1:length(dZVe)
    fprintf(file, '%16.15f %16.15f %16.15f\r\n', [dd(i) dZVe(i) dZVi(i)]);
end
fclose(file);

% 2nd Derivatives
file = fopen(fullfile(folder, name_file_d2), 'a');
for i = 1:length(ddZVe)
    fprintf(file, '%16.15f %16.15f %16.15f\r\n', [dd(i) ddZVe(i) ddZVi(i)]);
end
fclose(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catch ME
    file = fopen('errorFileEixam_adjoint_method_iPRC.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return