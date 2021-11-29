%%
%%%%%%%%%%%%%%%%%%%%% Communication through coherence %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% (Von Mises or sinusoidal stimulus) %%%%%%%%%%%%%%%%%%%

% Here we measure the increase of the half-width in the excitatory and
% inhibitory firing rates (PING/ING oscillators) due to the sinusoidal
% (type_prt = 1) or Von Mises (tyep_prt = 2) perturbation. The half-width
% is the time difference of those points that are half of the maximum
% amplitude in a waveform (existence guaranteed by the intermediate value
% theorem). This factor is computed along amplitude sections A = ctt of the 
% corresponding Arnold tongues. The T/T*-discretization along the sections 
% increases according to the amplitude.

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
    Asec --> Amplitude section
    Ie_ext --> Constant current to exc. neurons
%}

format long;

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
    Tpend = 2; %input('\nPeriod discretization up to (T/T^* < q): q = ');
    nPeriod = 1500; %input('\nNumber of points in a single period discretization: ');

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
            '_upto', num2str(Tpend), 'T']);
        if type == 1 && Ie_ext ~= 10
            name_file = strcat(name_file, ['_Ie_ext', num2str(Ie_ext)]);
        end
        name_file = strcat(name_file, '.txt');

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

% Parameters
if type == 1
    % Parameters (PING)
    tau_e = 8; tau_i = 8; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 5;
    Jee = 0; Jei = 13; Jie = 13; Jii = 0;
         
    % External inhibitory input
    Ii_ext = 0;
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

% Phase zero
phase0 = x0; t0 = 0;

% Integration tolerances (unperturbed and perturbed firing rate model)
options_ode = odeset('AbsTol', 1e-14, 'RelTol', 1e-13);

%%%%%%%%%%%%%%%%% Periodic orbit of the unperturbed system %%%%%%%%%%%%%%%%
[t_unprt, x_unprt] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
    [t0 T], phase0, options_ode);

%-----------------> Half-width of excitatory firing rate <----------------%
[Re0, indmax] = max(x_unprt(:,1)); Te_max = t_unprt(indmax)/T;
[Re1, indmin] = min(x_unprt(:,1)); Te_min = t_unprt(indmin)/T;

% Half of the maximum of the excitatory firing rate
Re_mid = (Re0 + Re1)/2;

Re_top = x_unprt(:,1) > Re_mid;
Re2mn = find(Re_top(1:end-1) == 1 & Re_top(2:end) == 0); % road to minimum
Re2mx = find(Re_top(1:end-1) == 0 & Re_top(2:end) == 1); % road to maximum

% Compute half-width of the excitatory peak of the unperturbed oscillator
Hwe0 = hw(t_unprt, x_unprt, Re2mn, Re2mx, Re_mid, 1, 1, T);
%-------------------------------------------------------------------------%

%-----------------> Half-width of inhibitory firing rate <----------------%  
[Ri0, indmax] = max(x_unprt(:,5)); Ti_max = t_unprt(indmax)/T;
[Ri1, indmin] = min(x_unprt(:,5)); Ti_min = t_unprt(indmin)/T;

% Half of the maximum of the inhibitory firing rate
Ri_mid = (Ri0 + Ri1)/2;

Ri_top = x_unprt(:,5) > Ri_mid;
Ri2mn = find(Ri_top(1:end-1) == 1 & Ri_top(2:end) == 0); % road to minimum
Ri2mx = find(Ri_top(1:end-1) == 0 & Ri_top(2:end) == 1); % road to maximum

% Compute half-width of the inhibitory peak of the unperturbed oscillator
Hwi0 = hw(t_unprt, x_unprt, Ri2mn, Ri2mx, Ri_mid, 5, 1, T);
%-------------------------------------------------------------------------%

norm(phase0 - x_unprt(end,:)')
T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating file/s where data will be saved
name_file = 'ctc_factors_halfwidth_eixam';
name_file = strcat(name_file, aux_name);

A_name = num2str(Asec); A_name(A_name == '.') = [];
name_file = strcat(name_file, '_A', A_name);
name_file_snd = strcat(name_file, '_sndpeak');
name_file_sndE = strcat(name_file, '_sndEpeak');

name_file = strcat(name_file, auxstr, '.txt'); % to save data from the (single) 1st peak
fclose(fopen(fullfile(folder, name_file), 'w'));

name_file_sndE = strcat(name_file_sndE, auxstr, '.txt'); % to save data from the 2nd exc. peak
fclose(fopen(fullfile(folder, name_file_sndE), 'w'));

name_file_snd = strcat(name_file_snd, auxstr, '.txt'); % to save data from the 2nd peak
fclose(fopen(fullfile(folder, name_file_snd), 'w'));


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
            [ri_max, ti_indmax] = max(x(:,5)); ti_max = t(ti_indmax);
            [re_max, te_indmax] = max(x(:,1)); te_max = t(te_indmax);
            
            % Times where the minima of inhibitory/excitatory firing rates take place
            [ri_min, ti_indmin] = min(x(:,5)); ti_min = t(ti_indmin);
            [re_min, te_indmin] = min(x(:,1)); te_min = t(te_indmin);
            
            % Determining number of inhibitory/excitatory peaks in a q*T_prime cycle
            xre = x(:,1); xri = x(:,5);
            peaks_indE = find(xre(2:end-1) >= xre(1:end-2) & xre(2:end-1) >= xre(3:end)) + 1;
            peaks_indI = find(xri(2:end-1) >= xri(1:end-2) & xri(2:end-1) >= xri(3:end)) + 1;
            te_aux = t(peaks_indE); ti_aux = t(peaks_indI);
            
            % Determining number of inhibitory/excitatory valleys in a q*T_prime cycle
            peaks_indE_min = find(xre(2:end-1) <= xre(1:end-2) & xre(2:end-1) <= xre(3:end)) + 1;
            peaks_indI_min = find(xri(2:end-1) <= xri(1:end-2) & xri(2:end-1) <= xri(3:end)) + 1;
            te_aux_min = t(peaks_indE_min); ti_aux_min = t(peaks_indI_min);
            if p > 1 && length(peaks_indE_min) == 1
                if xre(1) <= xre(2) && xre(1) <= xre(end)
                    peaks_indE_min = [1 peaks_indE_min];
                    te_aux_min = [t(1) te_aux_min];
                elseif xre(end) <= xre(end-1) && xre(end) <= xre(1)
                    peaks_indE_min = [peaks_indE_min length(xre)];
                    te_aux_min = [te_aux_min t(end)];
                end
            end
            
            % Half of the maximum value for the exc. and inh. firing rates
            re_mid = (re_max + re_min)/2;
            ri_mid = (ri_max + ri_min)/2;
            
            % Half-width for the second peak of the exc. and/or the inh. firing rates
            if length(peaks_indE) > 1
                % Second maximum peak for the exc. firing rate
                re_max_snd = xre(peaks_indE(2));
                if abs(te_max - te_aux(2)) < 0.1
                    re_max_snd = xre(peaks_indE(1));
                end
                
                % Second minimum peak for the exc. firing rate
                re_min_snd = xre(peaks_indE_min(2));
                if abs(te_min - te_aux_min(2)) < 0.1
                    re_min_snd = xre(peaks_indE_min(1));
                end
                
                % Half of the second maximum value for the exc. and inh. firing rates
                re_mid_snd = (re_max_snd + re_min_snd)/2;
                
                % Excitatory firing rate start above or below the half of the 2nd max
                re_ini = x(1,1) > re_mid_snd;
                
                %----Computing half-width of the second excitatory peak---%
                re_top_snd = x(:,1) > re_mid_snd;
                re2mn_snd = find(re_top_snd(1:end-1) == 1 & re_top_snd(2:end) == 0); % road to minimum
                re2mx_snd = find(re_top_snd(1:end-1) == 0 & re_top_snd(2:end) == 1); % road to maximum
                
                if re_ini % Excitatory activity starts above the half-maximum of the 2nd peak of the firing rate
                    if te_indmax < re2mn_snd(1) || te_indmax > re2mx_snd(2)
                        re2mn_snd = re2mn_snd(2); re2mx_snd = re2mx_snd(1);
                    else
                        re2mn_snd = re2mn_snd(1); re2mx_snd = re2mx_snd(2);
                    end
                else % Excitatory activity starts below the half-maximum of the 2nd peak of the firing rate
                    if te_indmax > re2mx_snd(1) && te_indmax < re2mn_snd(1)
                        re2mn_snd = re2mn_snd(2); re2mx_snd = re2mx_snd(2);
                    else
                        re2mn_snd = re2mn_snd(1); re2mx_snd = re2mx_snd(1);
                    end
                end
                
                % Compute half-width of the second excitatory peak
                HweA = hw(t, x, re2mn_snd, re2mx_snd, re_mid_snd, 1, q, T_prime);
                %---------------------------------------------------------%
                
                if length(peaks_indI) > 1
                    % Second maximum peak for the inh. firing rate
                    ri_max_snd = xri(peaks_indI(2));
                    if abs(ti_max - ti_aux(2)) < 0.1
                        ri_max_snd = xri(peaks_indI(1));
                    end
                
                    % Second minimum peak for the inh. firing rate
                    ri_min_snd = xri(peaks_indI_min(2));
                    if abs(ti_min - ti_aux_min(2)) < 0.1
                        ri_min_snd = xri(peaks_indI_min(1));
                    end
                
                    % Half of the second maximum value for the inh. firing rate
                    ri_mid_snd = (ri_max_snd + ri_min_snd)/2;
                
                    % Inhibitory firing rate start above or below the half of the 2nd max
                    ri_ini = x(1,5) > ri_mid_snd;
                    
                    %--Computing half-width of the second inhibitory peak-%
                    ri_top_snd = x(:,5) > ri_mid_snd;
                    ri2mn_snd = find(ri_top_snd(1:end-1) == 1 & ri_top_snd(2:end) == 0); % road to minimum
                    ri2mx_snd = find(ri_top_snd(1:end-1) == 0 & ri_top_snd(2:end) == 1); % road to maximum

                    if ri_ini % Inhibitory activity starts above the half-maximum of the firing rate
                        if ti_indmax < ri2mn_snd(1) || ti_indmax > ri2mx_snd(2)
                            ri2mn_snd = ri2mn_snd(2); ri2mx_snd = ri2mx_snd(1);
                        else
                            ri2mn_snd = ri2mn_snd(1); ri2mx_snd = ri2mx_snd(2);
                        end
                    else % Inhibitory activity starts below the half-maximum of the firing rate
                        if ti_indmax > ri2mx_snd(1) && ti_indmax < ri2mn_snd(1)
                            ri2mn_snd = ri2mn_snd(2); ri2mx_snd = ri2mx_snd(2);
                        else
                            ri2mn_snd = ri2mn_snd(1); ri2mx_snd = ri2mx_snd(1);
                        end
                    end

                    % Compute half-width of the second inhibitory peak
                    HwiA = hw(t, x, ri2mn_snd, ri2mx_snd, ri_mid_snd, 5, q, T_prime);
                    %-----------------------------------------------------%
                    
                    % Increase of the second half-width of the excitatory/inhibitory firing rates to perturbations
                    gammaE = HweA/Hwe0;
                    gammaI = HwiA/Hwi0;

                    % Saving data for the second excitatory and inhibitory peaks
                    file = fopen(fullfile(folder, name_file_snd), 'a');
                    fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, gammaE, gammaI);
                    fclose(file);
                else
                    % Increase of the second half-width of the excitatory firing rate to perturbations
                    gammaE = HweA/Hwe0;

                    % Saving data only for the second excitatory peak
                    file = fopen(fullfile(folder, name_file_sndE), 'a');
                    fprintf(file, '%d %d %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, gammaE);
                    fclose(file);
                end
            end
            
            % Excitatory/inhibitory firing rates start above or below the half maximum
            re_ini = x(1,1) > re_mid;
            ri_ini = x(1,5) > ri_mid;
            
            %------------------> Excitatory population <------------------%
            re_top = [x(:,1) > re_mid; x(1,1) > re_mid];
            re2mn = find(re_top(1:end-1) == 1 & re_top(2:end) == 0); % road to minimum
            re2mx = find(re_top(1:end-1) == 0 & re_top(2:end) == 1); % road to maximum
            xre(peaks_indE)
            length(re2mn)
            length(re2mx)
            if length(re2mn) > 1
                % If there is more than one peak higher than the half of
                % the maximum value re_mid, then take the halfwidth of the
                % maximum peak
                if re_ini % Excitatory activity starts above the half-maximum of the firing rate
                    if te_indmax < re2mn(1) || te_indmax > re2mx(2)
                        re2mn = re2mn(1); re2mx = re2mx(2);
                    else
                        re2mn = re2mn(2); re2mx = re2mx(1);
                    end
                else % Excitatory activity starts below the half-maximum of the firing rate
                    if te_indmax > re2mx(1) && te_indmax < re2mn(1)
                        re2mn = re2mn(1); re2mx = re2mx(1);
                    else
                        re2mn = re2mn(2); re2mx = re2mx(2);
                    end
                end
            end
            
            % Compute half-width of the excitatory peak
            HweA = hw(t, x, re2mn, re2mx, re_mid, 1, q, T_prime);
            %-------------------------------------------------------------%
            
            %------------------> Inhibitory population <------------------%
            ri_top = [x(:,5) > ri_mid; x(1,5) > ri_mid];
            ri2mn = find(ri_top(1:end-1) == 1 & ri_top(2:end) == 0); % road to minimum
            ri2mx = find(ri_top(1:end-1) == 0 & ri_top(2:end) == 1); % road to maximum
            xri(peaks_indI)
            length(ri2mn)
            length(ri2mx)
            if length(ri2mn) > 1
                % If there is more than one peak higher than the half of
                % the maximum value ri_mid, then take the halfwidth of the 
                % maximum peak
                if ri_ini % Inhibitory activity starts above the half-maximum of the firing rate
                    if ti_indmax < ri2mn(1) || ti_indmax > ri2mx(2)
                        ri2mn = ri2mn(1); ri2mx = ri2mx(2);
                    else
                        ri2mn = ri2mn(2); ri2mx = ri2mx(1);
                    end
                else % Inhibitory activity starts below the half-maximum of the firing rate
                    if ti_indmax > ri2mx(1) && ti_indmax < ri2mn(1)
                        ri2mn = ri2mn(1); ri2mx = ri2mx(1);
                    else
                        ri2mn = ri2mn(2); ri2mx = ri2mx(2);
                    end
                end
            end
            
            % Compute half-width of the excitatory peak
            HwiA = hw(t, x, ri2mn, ri2mx, ri_mid, 5, q, T_prime);
            %-------------------------------------------------------------%

            % Increase of the half-widths of the excitatory/inhibitory firing rates to perturbations
            gammaE = HweA/Hwe0;
            gammaI = HwiA/Hwi0;

            % Saving data
            file = fopen(fullfile(folder, name_file), 'a');
            fprintf(file, '%d %d %16.15f %16.15f %16.15f %16.15f\r\n', p, q, T_prime/T, A, gammaE, gammaI);
            fclose(file);
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

%_________________________________________________________________________%
function HwA = hw(t, x, r2mn, r2mx, r_mid, coord, q, T_prime)
    % Local function to compute the width at half the maximum amplitude of 
    % a waveform. Given the points before and after intersecting the
    % half-maximum line r_mid, we find the intersection points (linear
    % interpolation) and compute their (normalized) time difference.
    
    % Linear interpolation of the time at which the half of its
    % maximum value is attained (road to minimum)
    t0 = t(r2mn); t1 = t(r2mn+1);
    x0 = x(r2mn,coord); x1 = x(r2mn+1,coord);
    t2mn = t0*(1 - (r_mid - x0)/(x1 - x0)) + t1*((r_mid - x0)/(x1 - x0));
    t2mn = t2mn/(q*T_prime);

    % Linear interpolation of the time at which the half of its
    % maximum value is attained (road to maximum)
    t0 = t(r2mx); t1 = t(r2mx+1);
    x0 = x(r2mx,coord); x1 = x(r2mx+1,coord);
    t2mx = t0*(1 - (r_mid - x0)/(x1 - x0)) + t1*((r_mid - x0)/(x1 - x0));
    t2mx = t2mx/(q*T_prime);

    % Half-width of the excitatory/inhibitory firing rate
    if t2mx < t2mn
        HwA = t2mn - t2mx;
    else % take the complementary distance
        HwA = 1 - abs(t2mn - t2mx);
    end
end
%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%