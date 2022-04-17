%%
%%%%%%%%%%% Firing rate and synaptic variables - Time evolution %%%%%%%%%%%

% Here we display how the excitatory firing rate and the inhibitory synapse
% of the single-perturbed system are modified by the interference of a 
% second identical/non-identical Von Mises input. The unperturbed system is
% first integrated with a single Von Mises perturbation for a sufficiently 
% long time. We then add a second perturbation and integrate the resulting
% system for another long period.

close all; clear all; clc;

type = input('\nGamma mechanism (1 for PING or 2 for ING): '); % Gamma mechanism
coord = input('\nCoordinate (1 for Ve, 2 for Vi or 3 for both): '); % Component where perturbation is applied
if coord == 1
    coord = 2;
elseif coord == 2
    coord = 6;
else
    coord = [2, 6];
end

% External constant current (oscillations closer/further to the Hopf bifurcation)
Ie_ext = input('\nConstant current to exc. neurons: ');
auxstr = ['_Ie_ext', num2str(Ie_ext)];

% Amplitude, intrinsic parameters and frequency relationships
A = input('\nAmplitude of the perturbation: ');
k1 = input('\nConcentration factor of the first input: ');
p = input('\nRelation between periods of the 1st input and the oscillator (i.e. p = T1/T*): ');
k2 = input('\nConcentration factor of the second input: ');
qs = input('\nRelations between periods of the 2nd and the 1st inputs (i.e. q = T2/T1): ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

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

% Parameters
if type == 1
    % Parameters (PING)
    tau_e = 8; tau_i = 8; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 5;
    Jee = 0; Jei = 13; Jie = 13; Jii = 0;

    % External excitatory and inhibitory inputs
    Ii_ext = 0;
else
    % Parameters (ING)
    tau_e = 10; tau_i = 10; delta_e = 1; delta_i = 1;
    eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 1;
    Jee = 0; Jei = 10; Jie = 0; Jii = 15;
    
    % External excitatory and inhibitory inputs
    Ie_ext = 25; Ii_ext = 25;
end
parm = [tau_e; tau_i; delta_e; delta_i; eta_e; eta_i; 
        tau_se; tau_si; Jee; Jei; Jie; Jii; Ie_ext; Ii_ext];

% Initial condition lying on the oscillator of the unperturbed case
phase0 = x0;

% Integration tolerances (Perturbed firing rate model)
options_ode = odeset('AbsTol', 1e-14, 'RelTol', 1e-13);

%%% Digression: Impact of perturbation on firing rate (Ve VS Vi VS VeVi) %%
% Frequency-identical periodic von Mises distributions
pt1 = @(t) vonmises_dist(t, 0, k1, p*T);
pt2 = @(t) vonmises_dist(t, p*T/2, k2, p*T);

% Integrating system perturbed by the two von Mises inputs
[t_prt, x_prt] = ode45(@(t, x) perturbed_full_synaptic_firing_rate_equations(t, ...
    x, parm, A, pt1(t) + pt2(t), coord), [0 500], phase0, options_ode);

figure; hold on; box on; set(gca, 'Fontsize', 15); axfr = gca;
xlabel('Time (ms)', 'Fontsize', 18, 'FontWeight', 'bold');
ylabel('Firing rate (kHz)', 'Fontsize', 18, 'FontWeight', 'bold');
plot(t_prt, x_prt(:,1), 'r-', 'Linewidth', 1.5); % exc. firing rate
plot(t_prt, x_prt(:,5), 'b-', 'Linewidth', 1.5); % inh. firing rate
plot(t_prt, A*pt1(t_prt), '-k', 'Linewidth', 1.5); % 1st input
plot(t_prt, A*pt2(t_prt), '--k', 'Linewidth', 1.5); % 2nd input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap = hsv(6);

% Loop on the given set qs of relations T2/T1
for i = 1:length(qs)
    q = qs(i); % Period relation between the 2nd and 1st inputs
    
    % Periods of the periodic perturbations
    T1 = p*T;
    T2 = q*T1;
    T3 = max(T1,T2);

    % Phase-shift
    if T1 > T2
        mu2 = T2/2;
    else
        mu2 = T1/2;
    end

    % First periodic perturbation with period T1 (Von Mises distribution)
    pt1 = @(t) vonmises_dist(t, 0, k1, T1);

    % Onset of the second perturbation
    t_snd = 25*T1;

    % Second periodic perturbation with period T2 (Von Mises distribution)
    pt2 = @(t) vonmises_dist(t-t_snd, mu2, k2, T2);

    % Integrating perturbed system
    [t_prt, x_prt] = ode45(@(t, x) perturbed_full_synaptic_firing_rate_equations(t, ...
        x, parm, A, pt1(t), coord), [0 t_snd], phase0, options_ode);

    [t, x] = ode45(@(t, x) perturbed_full_synaptic_firing_rate_equations(t, ...
        x, parm, A, pt1(t)+pt2(t), coord), [t_prt(end) t_prt(end)+14*T3], ...
        x_prt(end,:), options_ode);

    % Axes and labels fontsize
    afs = 15; lfs = 17;

    % Creating figure
    figure; subplot(2,1,1); hold on; box on;
    set(gca, 'Fontsize', afs, 'XLabel', [], 'XTickLabel', []); ax = gca;
    xlim([-T1 t_prt(end)+5*T3-t_snd]);
    
    yyaxis left;
    ylabel('r_e (kHz)', 'Fontsize', lfs, 'FontWeight', 'bold');
    plot(t_prt-t_snd, x_prt(:,1), '-', 'Color', cmap(1,:), ...
        'Linewidth', 1.5); % exc. firing rate bef. 2nd input
    
    yyaxis right;
    ylabel('S_{ei}', 'Fontsize', lfs, 'FontWeight', 'bold');
    plot(t_prt-t_snd, x_prt(:,4), '-', 'Color', cmap(3,:), ...
        'Linewidth', 1.5); % inh. synapse bef. 2nd input

    yyaxis left;
    plot(t-t_snd, x(:,1), '-', 'Color', cmap(1,:), ...
        'Linewidth', 1.5); % exc. firing rate aft. 2nd input
    
    yyaxis right;
    plot(t-t_snd, x(:,4), '-', 'Color', cmap(3,:), ...
        'Linewidth', 1.5); % inh. synapse aft. 2nd input
    
    % YAxis color
    ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
    
    % Plotting 1st and 2nd inputs
    subplot(2,1,2); hold on; box on;
    set(gca, 'Fontsize', afs); ax1 = gca;
    xlim([-T1 t_prt(end)+5*T3-t_snd]);
    xlabel('Time (ms)', 'Fontsize', lfs, 'FontWeight', 'bold');
    ylabel('A p(t)', 'Fontsize', lfs, 'FontWeight', 'bold');
    plot(t_prt-t_snd, A*pt1(t_prt), '-k', 'Linewidth', 1.5);
    plot(t-t_snd, A*pt1(t), '-k', 'Linewidth', 1.5);
    plot(t-t_snd, A*pt2(t), '--k', 'Linewidth', 1.5);
    
    % Axes positions
    ax.Position = [ax.Position(1) ax.Position(2)-0.225 ...
                   ax.Position(3) ax.Position(4)+0.275];
    ax1.Position = [ax1.Position(1)+0.017 ax1.Position(2)+0.035 ...
                    ax1.Position(3)-0.06 ax1.Position(4)-0.18];
    ax.Position(1) = ax1.Position(1);
    ax.Position(3) = ax1.Position(3);
    
    % Saving figure
    folder_fig = strcat(pwd, '\Firing Rate - 2 competitors'); % folder
    name_fig = 'firingrate_2inputs';
    if type == 1
        name_fig = strcat(name_fig, ['_ping_coord_', strrep(num2str(coord), ' ', '')]);
    else
        name_fig = strcat(name_fig, ['_ing_coord_', strrep(num2str(coord), ' ', '')]);
    end
    name_fig = strcat(name_fig, ['_A', num2str(A), '_k1', num2str(k1), ...
        '_p', num2str(p), '_k2', num2str(k2), '_q', num2str(q)]); % filename
    
    % PDF format
    set(ax, 'Fontsize', 36); set(ax1, 'Fontsize', 36);
    set(get(ax, 'XLabel'), 'FontSize', 38); set(get(ax1, 'XLabel'), 'FontSize', 38);
    set(get(ax, 'YLabel'), 'FontSize', 38); set(get(ax1, 'YLabel'), 'FontSize', 38);
    set(gcf, 'PaperPositionMode', 'manual');   
    set(gcf, 'PaperOrientation', 'landscape');
    name_pdf = strcat(name_fig, '.pdf');
    fprintf('\nSaving figure: ''%s''\n', name_pdf);
    print(gcf, '-dpdf', '-fillpage', fullfile(folder_fig, name_pdf));

    % EPS format
    set(ax, 'Fontsize', afs); set(ax1, 'Fontsize', afs);
    set(get(ax, 'XLabel'), 'FontSize', lfs); set(get(ax1, 'XLabel'), 'FontSize', lfs);
    set(get(ax, 'YLabel'), 'FontSize', lfs); set(get(ax1, 'YLabel'), 'FontSize', lfs);
    name_eps = strcat(name_fig, '.eps');
    fprintf('\nSaving figure: ''%s''\n', name_eps);
    print(gcf, '-depsc', '-tiff', fullfile(folder_fig, name_eps));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%