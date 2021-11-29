%%
%%%%%%%%%%%%%%%%%%%%% Communication through coherence %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% (Von Mises or sinusoidal stimulus) %%%%%%%%%%%%%%%%%%%

% Here we measure effective communication (in the CTC framework) between
% the PING/ING oscillator and the sinusoidal (type_prt = 1) or the Von
% Mises (type_prt = 2) external drive. These factors consists in the timing
% of the inhibition with respect to the perturbation and the increase in
% the excitatory response due to that perturbation. These factors are
% computed along amplitude sections A = ctt of the corresponding Arnold 
% tongues. The T/T*-discretization along the sections increases according
% to the amplitude.

% The next set of parameters must be initialized before running:
type
coord
type_prt
p
q
k
Asec
Ie_ext

% Arguments from command window (only octave)
% args = argv()
% type = args{1};
% coord = args{2};
% type_prt = args{3};
% p = args{4};
% q = args{5};
% k = args{6};
% Asec = args{7};
% Ie_ext = args{8};

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi)
    type_prt --> Type of perturbation (1 for Sinusoidal, 2 for Von Mises)
    p --> Revolutions completed by the oscillator
    q --> Revolutions completed by the input
    k --> Von Mises parameter (input concentration factor)
    Asec --> Set of amplitude sections
    Ie_ext --> Constant current to exc. neurons
%}

format long;

% Parameters
if type == 1
    % Parameters (PING)
    tau_e = 8; tau_i = 8; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 5;
    Jee = 0; Jei = 13; Jie = 13; Jii = 0;
         
    % External inhibitory input
    Ii_ext = 0; % <-------------------------------------------
else
    % Parameters (ING)
    tau_e = 10; tau_i = 10; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 1;
    Jee = 0; Jei = 10; Jie = 0; Jii = 15; % <------------------------------

    % External excitatory and inhibitory inputs
    Ie_ext = 25; Ii_ext = 25;
end
parm = [tau_e; tau_i; delta_e; delta_i; eta_e; eta_i; 
        tau_se; tau_si; Jee; Jei; Jie; Jii; Ie_ext; Ii_ext];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% Auxiliar string for non-zero values of Ie_ext
auxstr = ['_Ie_ext', num2str(Ie_ext)];

% --> Oscillator period
if type == 1
    name_file = ['period_ping', auxstr, '.txt'];
else
    name_file = ['period_ing', auxstr, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'r');
T = fscanf(file, '%f');
fclose(file);

% --> Initial condition (on the oscillator)
if type == 1
    name_file = ['initial_condition_ping', auxstr, '.txt'];
else
    name_file = ['initial_condition_ing', auxstr, '.txt'];
end
file = fopen(fullfile(folder, name_file), 'r');
formatSpec = '%f %f %f %f %f %f %f %f'; sizeA = [8 Inf];
res = fscanf(file, formatSpec, sizeA); x0 = res;
fclose(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Periodic perturbation with period Tp
if type_prt == 1
    % Sinusoidal input
    pt = @(t,Tp) 1 + cos(2*pi*t/Tp);
else
    % Von Mises input
    pt = @(t,Tp) vonmises_dist(t, 0, k, Tp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Devil staircases %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Ie_ext ~= 10
    Tpend = 2; % Period discretization up to T/T^* <= Tpend
    nPeriod = 1500; % Number of points in a single period discretization

    % File/s to be loaded
    t_period = []; rot_number = t_period;
    while Tpend >= 1
        % Filename
        if type == 1
            name_file = ['dvlstaircases_ping_coord_', num2str(coord)];
        else
            name_file = ['dvlstaircases_ing_coord_', num2str(coord)];
        end
        name_file = strcat(name_file, 'eqlinint');
        if type_prt == 2
            name_file = strcat(name_file, ['_k', num2str(k)]);
        end
        name_file = strcat(name_file, ['_A', num2str(Asec), '_numTp', num2str(nPeriod), ...
            '_upto', num2str(Tpend), 'T', auxstr, '.txt']);

        % Loading data
        file = fopen(fullfile(folder, name_file), 'r');
        formatSpec = '%f %f'; size_rot = [2 Inf];
        res = fscanf(file, formatSpec, size_rot); res = res';
        t_period = [t_period; flipud(res(:,1))];
        rot_number = [rot_number; flipud(res(:,2))];
        fclose(file);

        Tpend = Tpend - 1;
    end
    
    t_period = flipud(t_period);
    rot_number = flipud(rot_number);

    % Starting and ending points of the p:q plateau
    rt = p/q;
    ind1 = find(rot_number(1:end-1) ~= rt & rot_number(2:end) == rt);
    ind2 = find(rot_number(1:end-1) == rt & rot_number(2:end) ~= rt);
    ind1 = ind1 + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinate where the perturbation is applied
if coord == 1
    coord = 2;
else
    coord = 6;
end

%%%%%%%%%%%%%%%%%%%%%%%% Arnold tongues boundaries %%%%%%%%%%%%%%%%%%%%%%%%
% Files to be loaded
if type == 1
    aux_name = ['_', num2str(p), num2str(q), '_ping_coord_'];
    if coord == 2
        aux_name = strcat(aux_name, 'Ve');
    else
        aux_name = strcat(aux_name, 'Vi');
    end
else
    aux_name = ['_', num2str(p), num2str(q), '_ing_coord_Vi'];
end
if type_prt == 2
    aux_name = strcat(aux_name, ['_k', num2str(k)]);
end
name_left = strcat('left_branch', aux_name, auxstr, '.txt');
name_right = strcat('right_branch', aux_name, auxstr, '.txt');

% Loading left branch
file = fopen(fullfile(folder, name_left), 'r');
formatSpec = '%f %f %f %d %d %f'; sizeA = [6 Inf];
res = fscanf(file, formatSpec, sizeA); res = res';
T_prime = res(:,2); T_prime(end+1) = p*T/q;
Al = res(:,3); Al(end+1) = 0; [Al_ord, ind] = sort(Al);
Tl_prime = T_prime(ind)./T;
fclose(file);

% Loading right branch
file = fopen(fullfile(folder, name_right), 'r');
formatSpec = '%f %f %f %d %d %f'; sizeA = [6 Inf];
res = fscanf(file, formatSpec, sizeA); res = res';
T_prime = res(:,2); T_prime(end+1) = p*T/q;
Ar = res(:,3); Ar(end+1) = 0; [Ar_ord, ind] = sort(Ar);
Tr_prime = T_prime(ind)./T;
fclose(file);

% Sample points in interpolation must be unique and in ascending order
Al_ord = unique(Al_ord); Ar_ord = unique(Ar_ord);
Tl_prime = unique(Tl_prime); Tr_prime = unique(Tr_prime);

% Linear interpolation of left and right branches
% left_branch = @(l) interp1q(Al_ord, Tl_prime, l);
% right_branch = @(l) interp1q(Ar_ord, Tr_prime, l);
left_branch = griddedInterpolant(Al_ord, Tl_prime);
if Tl_prime(1) ~= p/q
    left_branch = griddedInterpolant(Al_ord, flipud(Tl_prime));
end

right_branch = griddedInterpolant(Ar_ord, Tr_prime);
if Tr_prime(1) ~= p/q
    right_branch = griddedInterpolant(Ar_ord, flipud(Tr_prime));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phase zero
phase0 = x0; t0 = 0;

% Integration tolerances (unperturbed and perturbed firing rate model)
options_ode = odeset('AbsTol', 1e-14, 'RelTol', 1e-13);

%%%%%%%%%%%%%%%%% Periodic orbit of the unperturbed system %%%%%%%%%%%%%%%%
[t_unprt, x_unprt] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
    [t0 T], phase0, options_ode);

% Maximum of the excitatory firing rate
Re0 = max(x_unprt(:,1));

% Time-average of the excitatory firing rate (integral with trapezoidal rule)
Re0_mean = (1/T)*trapz(t_unprt, x_unprt(:,1));

if type == 2
    % Maximum of the inhibitory firing rate
    Ri0 = max(x_unprt(:,5));

    % Time-average of the inhibitory firing rate (integral with trapezoidal rule)
    Ri0_mean = (1/T)*trapz(t_unprt, x_unprt(:,5));
end

norm(phase0 - x_unprt(end,:)')
T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating file/s where data will be saved
name_file = 'ctc_factors_eixam';
name_file_theta2 = strcat(name_file, aux_name);
name_file_tau2 = strcat(name_file, aux_name);
name_file = strcat(name_file, aux_name);

A_name = num2str(Asec); A_name(A_name == '.') = [];
name_file = strcat(name_file, '_A', A_name);
name_file_tau2 = strcat(name_file_tau2, '_A', A_name, '_sndpeak');
name_file_theta2 = strcat(name_file_theta2, '_A', A_name, '_sndEpeak');
if type == 2
    name_file_ri = strcat(name_file, '_Ri.txt');
    name_file = strcat(name_file, '_Re');    
    fclose(fopen(fullfile(folder, name_file_ri), 'w'));
end
name_file = strcat(name_file, auxstr, '.txt'); % to save data from the (single) 1st peak
fclose(fopen(fullfile(folder, name_file), 'w'));

name_file_theta2 = strcat(name_file_theta2, auxstr, '.txt'); % to save data from the 2nd exc. peak
fclose(fopen(fullfile(folder, name_file_theta2), 'w'));

name_file_tau2 = strcat(name_file_tau2, auxstr, '.txt'); % to save data from the 2nd peak
fclose(fopen(fullfile(folder, name_file_tau2), 'w'));

try

for ii = 1:length(Asec)

    A = Asec(ii); % Amplitude section
    
    % Period discretization
    if Ie_ext == 10
        lt = left_branch(A);
        rt = right_branch(A);
    else
        lt = t_period(ind1);
        rt = t_period(ind2);
    end
    
    if q == 1
        Tdis = linspace(lt, rt, ceil(A*2400));
    else
        Tdis = linspace(lt, rt, ceil(A*3000));
    end
    
    if A == 0
        Tdis = Tdis(1);
    end
    
    % Initial condition for the first discretization point is the phase zero
    if ii == 1
        phaseA = phase0;
    end
    first_phaseA = true;
        
    for i = 1:length(Tdis)
        i
        T_prime = Tdis(i)*T;

        % Initial condition for the next case is the periodic solution from
        % the previous one (if convergence was achieved)
        if i == 1
            phaseT = phaseA;
        end

        %%%%%% Convergence to periodic orbit of the perturbed system %%%%%%
        conv = false; % initializing convergence indicator
        count = 0; % counting number of iterations to converge
        inf_loop = false; % detecting infinite loop (periodic orbit does no longer exist)
        while ~conv && ~inf_loop
            % Transition to (if exists) periodic orbit of the perturbed system
            [~, x_prt] = ode45(@(t, x) perturbed_full_synaptic_firing_rate_equations(t, ...
                x, parm, A, pt(t,T_prime), coord), [0 q*T_prime], phaseT, options_ode);

            res = norm(x_prt(end,:)' - x_prt(1,:)');
            count = count + 1;

            phaseT = x_prt(end,:)'; % initial condition for the next iteration
            inf_loop = count >= 750 && res > 1e-3; % detection of infinite loop
            conv = res < 5e-11; % updating convergence criterion (10 significant digits)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if res < 5e-9
            % Initial condition for the first discretization point of the
            % next amplitude section (if any)
            if first_phaseA
                first_phaseA = false;
                phaseA = phaseT;
            end
            
            % Periodic orbit of the perturbed system
            [t, x] = ode45(@(t, x) perturbed_full_synaptic_firing_rate_equations(t, ...
                x, parm, A, pt(t,T_prime), coord), [0 q*T_prime], phaseT, options_ode);
            
            % Times where the maxima of inhibitory/excitatory firing rates take place
            [~, ti_ind] = max(x(:,5)); ti = t(ti_ind);
            [~, te_ind] = max(x(:,1)); te = t(te_ind);
            
            % Determining number of inhibitory/excitatory peaks in a q*T_prime cycle
            xre = x(:,1); xri = x(:,5);
            peaks_indE = find(xre(2:end-1) >= xre(1:end-2) & xre(2:end-1) >= xre(3:end)) + 1;
            peaks_indI = find(xri(2:end-1) >= xri(1:end-2) & xri(2:end-1) >= xri(3:end)) + 1;
            te_aux = t(peaks_indE); ti_aux = t(peaks_indI);
            
            % Check if any local maxima is reached at the endpoints
            if xre(end) >= xre(end-1) && xre(end) >= xre(2)
                peaks_indE = [peaks_indE; length(xre)];
                te_aux = [te_aux; t(end)];
            end

            % Perturbation peaks
            tp = 0; % First peak
            if q > 1
                tp2 = T_prime; % Second peak (middle one)
            end
            
            % Time differences between maximum of inhibition and maxima of perturbation/excitation
            delta_tau = (1/T_prime)*(ti - tp);
            delta_theta = (1/T_prime)*(ti - te);
            
            ti_aux
            te_aux
            if length(ti_aux) > 1 % two excitatory/inhibitory peaks
                % Check if variable ti_last exists or not
                if ~exist('ti_last', 'var')
                    ti_last = ti_aux(1);
                end
                
                % Index of the original inhibition peak
                snd_ti = ismembertol(ti_aux, ti_last, 0.05);
                if ~any(snd_ti) % the original inhibition peak occurs near the final time
                    snd_ti = ismembertol(q*T_prime - ti_aux, ti_last, 0.05);
                end
                ti = ti_aux(snd_ti == 1); % original inhibitory peak
                ti_last = ti; % update time of the original inhibition for the next iteration
                delta_tau = (1/T_prime)*(ti - tp); % delta tau for the original peak
                
                % Find the time of the excitatory peak before the original inhibition
                te = te_aux(te_aux - ti <= 0);
                if isempty(te) || length(te) > 1
                    te = te_aux(2); % choose the furthest excitatory peak
                end
                delta_theta = (1/T_prime)*(ti - te);
            elseif length(te_aux) > 1 % two excitatory peaks per one of inhibition
                % Check if variable te_last exists or not
                if ~exist('te_last', 'var')
                    te_last = te_aux(1);
                end
                
                % Indices of the excitatory peaks
                snd_te = ismembertol(te_aux, te_last, 0.05);
                snd_te1 = find(snd_te == 1);
                snd_te2 = find(snd_te == 0);
                
                te = te_aux(snd_te1); % original excitatory peak
                te_last = te; % update time of the original excitation for the next iteration
                ti_last = ti; % update time of the max. inhibition for the next iteration
                delta_theta = (1/T_prime)*(ti - te); % delta theta for the original peak
                
                % Inhibition-excitation timing for the second excitation peak
                delta_theta2 = (1/T_prime)*(ti - te_aux(snd_te2));
                if te_aux(snd_te2) > ti
                    % Take the complementary distance to avoid a jump in delta_theta2
                    delta_theta2 = (1/T_prime)*(q*T_prime - abs(ti - te_aux(snd_te2)));
                end
            else % one-to-one relation between excitation and inhibition
                te_last = te; % time in which excitation reaches a maximum
                ti_last = ti; % time in which inhibition reaches a maximum
            end
            
            % -- If excitation occurs after inhibition due to the periodicity condition
            if te > ti
                % Take the complementary distance to avoid a jump in delta_theta
                delta_theta = (1/T_prime)*(T_prime - abs(ti - te));
            end
            
            % -- If inhibition is found during the second cycle of the perturbation (q = 2)
            if q > 1 && ti/T_prime > q/2
                delta_tau = (1/T_prime)*(ti - tp2); % Definition from 0 to 1
            end
            
            % Time differences between maximum of inhibition and maxima of 
            % perturbation and excitation (second peak if any)
            if length(ti_aux) > 1
                
                % Index of the inhibition's second peak
                snd_ti = ismembertol(ti_aux, ti, 1e-5);
                snd_indI = find(snd_ti == 0);
                
                % Second pair inhibition-perturbation
                delta_tau2 = (1/T_prime)*(ti_aux(snd_indI) - tp);
                
                % -- If second inhibition is found during the second cycle of the perturbation
                if q > 1 && ti_aux(snd_indI)/T_prime > q/2 
                    delta_tau2 = (1/T_prime)*(ti_aux(snd_indI) - tp2); % Definition from 0 to 1
                end
                
                if length(te_aux) > 1
                    
                    % Index of the excitation's second peak
                    snd_te = ismembertol(te_aux, te, 1e-5);
                    snd_indE = find(snd_te == 0);
                    
                    % Second pair inhibition-excitation
                    delta_theta2 = (1/T_prime)*(ti_aux(snd_indI) - te_aux(snd_indE));
                    if te_aux(snd_indE) > ti_aux(snd_indI)
                        % Take the complementary distance to avoid a jump in delta_theta2
                        delta_theta2 = (1/T_prime)*(T_prime - abs(ti_aux(snd_indI) - te_aux(snd_indE)));
                    end
                    
                    % Second response of the excitatory receiving population to the external perturbation (p > 1)
                    alpha2 = xre(peaks_indE(snd_indE))/Re0;
                    
                    % -- If second exc. peak is the maximum one, re-assigning the minimum one
                    if all(xre(peaks_indE(snd_indE)) >= xre(peaks_indE))
                        all(xre(peaks_indE(snd_indE)) >= xre(peaks_indE))
                        snd_indE = find(snd_te == 1);
                        alpha2 = xre(peaks_indE(snd_indE))/Re0;
                    end
                end
            end
            
            if type == 2
                % ING
                RiA = max(x(:,5));
                RiA_mean = (1/(q*T_prime))*trapz(t, x(:,5));

                % Response of the inhibitory receiving population to the external perturbation (maximum form)
                alpha = RiA/Ri0;

                % Response of the inhibitory receiving population to the external perturbation (integral form)
                alpha_mean = RiA_mean/Ri0_mean;
                
                % Saving data (Inhibitory population)
                file = fopen(fullfile(folder, name_file_ri), 'a');
                fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, delta_tau, alpha, alpha_mean, delta_theta);
                fclose(file);
            end
            
            % PING (or excitatory population in the ING)
            ReA = max(x(:,1));
            ReA_mean = (1/(q*T_prime))*trapz(t, x(:,1));

            % Response of the excitatory receiving population to the external perturbation (maximum form)
            alpha = ReA/Re0;

            % Response of the excitatory receiving population to the external perturbation (integral form)
            alpha_mean = ReA_mean/Re0_mean;

            % Saving data (excitatory population PING/ING)
            file = fopen(fullfile(folder, name_file), 'a');
            fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, delta_tau, alpha, alpha_mean, delta_theta);
            fclose(file);
            
            % Saving second factor tau/theta (if defined)
            if length(ti_aux) > 1 % two pairs of excitatory-inhibitory peaks
                file = fopen(fullfile(folder, name_file_tau2), 'a');
                if length(te_aux) > 1
                    fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, delta_tau2, delta_theta2, alpha2);
                else
                    fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, delta_tau2, -1, -1);
                end
                fclose(file);
            elseif length(te_aux) > 1 % two excitatory peaks
                file = fopen(fullfile(folder, name_file_theta2), 'a');
                fprintf(file, '%d %d %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, delta_theta2);
                fclose(file);
            end
        end
    end
end

catch ME
    file = fopen('errorFileEixam.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = ME.message;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%