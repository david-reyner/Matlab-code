%%
%%%%%%%%%%%%%%%%%% Characterization of gamma oscillations %%%%%%%%%%%%%%%%%

% Here we aim to study the evolution of the oscillation period, the maximum
% and integral mean values (time-average) of the excitatory and  inhibitory
% firing rates and their time and phase differences of the emerging 
% oscillations when the excitatory or inhibitory driving currents are 
% modulated. The discretization of the external current to modulate
% consists of 500 points. To find the periodic orbit we integrate the
% firing rate model for a sufficiently large time to ensure that the final
% integration point lies on the periodic orbit (if any, for that value of
% the external current). For weakly attracting periodic orbits we integrate
% the model for a larger time until the norm difference between the initial
% and final points in a single integration of time T is small enough.

% The next parameter must be initialized before running:
driv_curr

% Arguments from command window (only octave)
% args = argv()
% driv_curr = args{1};

% Script arguments
%{
    driv_curr --> Applying a variable driving current (1 to exc. pop. or 2 to inh. pop.) 
%}

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% Integration tolerances (firing rate model)
options_ode = odeset('AbsTol', 1e-18, 'RelTol', 1e-13);

% Initial/Final integration time
t0 = 0; tmax = 500;

% Initial condition
x00 = [0.1; -1.5; 0.5; 0.5;
       0.1; -1.5; 0.5; 0.5];

% Parameters (PING)
tau_e = 8; tau_i = 8; delta_e = 1; delta_i = 1;
eta_e = -5; eta_i = -5; tau_se = 1; tau_si = 5;
Jee = 0; Jei = 13; Jie = 13; Jii = 0;
parm_ping = [tau_e; tau_i; delta_e; delta_i; eta_e; eta_i; 
             tau_se; tau_si; Jee; Jei; Jie; Jii];

% External excitatory and inhibitory inputs
if driv_curr == 1
%     I_ext = linspace (5, 15, 500); Ii_ext = 0;
    I_ext = linspace (6, 15, 500); Ii_ext = 0;
else
%     Ie_ext = 15; I_ext = linspace(10, -5, 500);
    Ie_ext = 12; I_ext = linspace(6, -1, 500);
end

% Period of oscillations
Tosc = zeros(length(I_ext), 3);

% Time and phase difference
diff_phases = zeros(length(I_ext), 1);
diff_phases_norm = zeros(length(I_ext), 1);

% Maximum firing rates (as a function of the driving current)
re = zeros(length(I_ext), 1);
ri = zeros(length(I_ext), 1);

% Mean firing rates (as a function of the driving current)
re_mean = zeros(length(I_ext), 1);
ri_mean = zeros(length(I_ext), 1);

% Poincare Section (PS): Ve = 0
Vps = 0; coord = 2;    

% EVENTS handle function: solves and finds where functions of (t,x) are 0
myfunc = @(t, x) periodPO(t, x, Vps, coord);
options_ode2 = odeset('AbsTol', 1e-18, 'RelTol', 1e-13, 'Events', myfunc);

% Files where data will be saved
if driv_curr == 1
    name_file = 'oscillation_frequency_time_difference_exc_population';
else
    name_file = 'oscillation_frequency_time_difference_inh_population';
end
if tau_si ~= 1
    name_file = strcat(name_file, ['_tau_si', num2str(tau_si)]);
end
fclose(fopen(fullfile(folder, strcat(name_file, '.txt')), 'w'));

% Loop varying the excitatory/inhibitory driving current
for i = 1:length(I_ext)
    % Temporal variable that integrates the exact firing rate model for a
    % larger time in case of weakly attracting periodic orbits near the
    % Hopf bifurcation
    add_time = 0;
    
    period = false;
    while ~period % until fully convergence to the periodic orbit (if any)
        i
        if driv_curr == 1
            parm = [parm_ping; I_ext(i); Ii_ext];
        else
            parm = [parm_ping; Ie_ext; I_ext(i)];
        end
        [t, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
            [t0 tmax + add_time], x00, options_ode);

        % Last integration point (initial condition lying on a periodic orbit, if any)
        x0 = x(end,:);

        % Initial seed for the next iteration (closer to the periodic solution, if any)
        x00 = x0; 

        % In general the initial condition x0 on the PO (if any) does not 
        % lie on the Poincare Section (PS). Thus, we perform a preliminary 
        % integration to obtain an initial condition on the PS and also on 
        % the PO. If no orbit exists, the preliminary integration will not
        % cut the Poincare section and we proceed with the next driving
        % current value
        [~, ~, ~, xe, ~] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
            [t0 tmax], x0, options_ode2);

        % Initial condition also on the Poincare Section PS (if xe not empty)
        x0 = xe'; y0 = x0;

        if isempty(y0) == 0
            % Computing the period of the periodic orbit
            n_crossing = 2; ind = 1; T = 0;
            while ind <= n_crossing
                [~, y, te, ~, ~] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
                    [t0 tmax], y0, options_ode2);
                y0 = y(end,:);
                T = T + te(end);
                ind = ind + 1;
            end

            % Period of the periodic orbit
            Tosc(i,1:2) = [i T];

            % Checkig period
            [t, x] = ode45(@(t, x) full_synaptic_firing_rate_equations(t, x, parm), ...
                [t0 T], x0, options_ode);

            % Firing rates maximum values
            [re(i), indre] = max(x(:,1));
            [ri(i), indri] = max(x(:,5));
            
            % Firing rate integral mean values
            re_mean(i) = (1/T)*trapz(t, x(:,1));
            ri_mean(i) = (1/T)*trapz(t, x(:,5));
           
            % Measuring time and phase differences
            diff_phases(i) = (t(indri) - t(indre));
            diff_phases_norm(i) = (1/T)*(t(indri) - t(indre));
            
            % Error between final and initial point after integration of T
            Tosc(i,3) = norm(x0(coord-1:end) - x(end,coord-1:end)');
            norm(x0(coord-1:end) - x(end,coord-1:end)')
            norm(x0(coord-1:end) - x(end,coord-1:end)') < 5e-14 % 13 significant digits
            norm(x0(coord-1:end) - x(end,coord-1:end)') < 5e-13 % 12 significant digits
            
            period = Tosc(i,3) < 1e-12;
            add_time = add_time + 100; % increase integration time
        else % If there is no intersection, we assume there is no periodic orbit
            period = true; % Next value of the driving current I_ext
        end
    end    
end

% Saving data
file = fopen(fullfile(folder, strcat(name_file, '.txt')), 'a');
formatdata = '%16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f %16.15f\r\n';
for i = 1:length(I_ext)
    fprintf(file, formatdata, I_ext(i), Tosc(i,2), re(i), ri(i), ...
        re_mean(i), ri_mean(i), diff_phases(i), diff_phases_norm(i));
end
fclose(file);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%% Reading and interpretation of data %%%%%%%%%%%%%%%%%%%

% In this second part of the script we proceed to read the saved files,
% plot the corresponding variables computed in the previous part of the
% script and save the resulting figures.

format long;
close all; clc;

% Folder from where to read the corresponding files
folder = fullfile(pwd, 'Oscillation frequency');

% Excitatory or inhibitory driving current
driv_curr = input('\nApplied variable driving current ([1] to Exc. pop. or [2] to Inh. pop.): ');
if driv_curr == 1
    name_file = 'oscillation_frequency_time_difference_exc_population';
else
    name_file = 'oscillation_frequency_time_difference_inh_population';
end

% Inhibitory-to-inhibitory synapse
Jii = input('\nSynaptic strength of the inhibitory synapse Sii: ');
if Jii ~= 0
    name_file = strcat(name_file, ['_Jii', num2str(Jii)]);
end

% Time constant of inhibitory synapse
tau_si = input('\nTime constant of the inhibitory synapses tau_si: ');
if tau_si ~= 1
    name_file = strcat(name_file, ['_tau_si', num2str(tau_si)]);
end
name_file = strcat(name_file, '_tau8(1).txt')

file = fopen(fullfile(folder, name_file), 'r');
formatdata = '%f %f %f %f %f %f'; size_file = [8 Inf];
res = fscanf(file, formatdata, size_file); res = res';
fclose(file);

% Unpacking variables
I_ext = res(:,1); Tosc = res(:,2);
re = res(:,3); ri = res(:,4); % Firing rates maximum values
re_mean = res(:,5); ri_mean = res(:,6); % Firing rate integral mean values
diff_phases = res(:,7); diff_phases_norm = res(:,8); % Time and phase differences

pos_osc = find(Tosc ~= 0); % Discarding I_ext values with no periodic orbit
mn = I_ext(pos_osc(1)); mx = I_ext(pos_osc(end)); % First & last values with PO data

%%%%%%%%%%%%%%% Frequency, maximum and integral mean values %%%%%%%%%%%%%%%
figure; box on; hold on; set(gca, 'Fontsize', 16); ax = gca;
yyaxis left;
ax.YAxis(1).Color = 'k';
p1 = plot(I_ext(pos_osc), 1000./Tosc(pos_osc), 'Color', [0 0.7 0], 'Linewidth', 1.5); % Frequency
p4 = plot(I_ext(pos_osc), 1000*re_mean(pos_osc), 'r--', 'Linewidth', 1.5); % Exc. firing rate integral mean value
p5 = plot(I_ext(pos_osc), 1000*ri_mean(pos_osc), 'b--', 'Linewidth', 1.5); % Inh. firing rate integral mean value
ylabel('Rate/Frequency (Hz)', 'Fontsize', 18, 'Fontweight', 'bold');

yyaxis right;
ax.YAxis(2).Color = 'k';
p2 = plot(I_ext(pos_osc), re(pos_osc), 'r-', 'Linewidth', 1.5); % Exc. firing rate max value
p3 = plot(I_ext(pos_osc), ri(pos_osc), 'b-', 'Linewidth', 1.5); % Inh. firing rate max value
ylabel('Max. firing rates (kHz)', 'Fontsize', 18, 'Fontweight', 'bold');
if driv_curr == 1
    xlabel('Excitatory driving current, I_e^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
else
%     set(gca, 'XDir', 'Reverse');
    xlabel('Inhibitory driving current, I_i^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
end
xlim([min(mn, mx)-0.5 max(mn, mx)+0.5]);

% hleg = legend([p1 p2 p3 p4 p5], {'Frequency', 'Max. r_e', 'Max. r_i', ...
%     'Mean value r_e', 'Mean value r_i'}, 'Location', 'Northwest', ...
%     'Fontsize', 16, 'Fontweight', 'bold');
% 
% % -- Linewidth of the displayed graphic line in each LegendEntry's Icon
% % (http://undocumentedmatlab.com/articles/plot-legend-customization)
% for i = 1:length(hleg.PlotChildren)
%     hLegendEntry = hleg.EntryContainer.NodeChildren(i);
%     hLegendIconLine = hLegendEntry.Icon.Transform.Children.Children;
%     hLegendIconLine.LineWidth = 2;
% end

% Saving figure
name_fig = 'max_rate_frequency';
if driv_curr == 1
    name_fig = strcat(name_fig, '_excitatory');
else
    name_fig = strcat(name_fig, '_inhibitory');
end
if Jii ~= 0
    name_fig = [name_fig '_Jii', num2str(Jii)];
end
if tau_si ~= 1
    name_fig = [name_fig '_tau_si', num2str(tau_si)];
end
print(gcf, '-depsc', '-tiff', fullfile(folder, [name_fig, '_tau8.eps']));

figure; box on; hold on; set(gca, 'Fontsize', 16);
p1 = plot(I_ext(pos_osc), 1000./Tosc(pos_osc), 'Color', [0 0.7 0], 'Linewidth', 1.5); % Frequency
p4 = plot(I_ext(pos_osc), 1000*re_mean(pos_osc), 'r--', 'Linewidth', 1.5); % Exc. firing rate integral mean value
p5 = plot(I_ext(pos_osc), 1000*ri_mean(pos_osc), 'b--', 'Linewidth', 1.5); % Inh. firing rate integral mean value
ylabel('Rate/Frequency (Hz)', 'Fontsize', 18, 'Fontweight', 'bold');
if driv_curr == 1
    xlabel('Excitatory driving current, I_e^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
    posLeg = 'Southeast';
else
%     set(gca, 'XDir', 'Reverse');
    xlabel('Inhibitory driving current, I_i^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
    posLeg = 'Southwest';
end
xlim([min(mn, mx)-0.5 max(mn, mx)+0.5]);

hleg = legend([p1 p4 p5], {'Frequency', 'Mean value r_e', ...
    'Mean value r_i'}, 'Location', posLeg, 'Fontsize', 16, ...
    'Fontweight', 'bold');

% -- Linewidth of the displayed graphic line in each LegendEntry's Icon
% (http://undocumentedmatlab.com/articles/plot-legend-customization)
for i = 1:length(hleg.PlotChildren)
    hLegendEntry = hleg.EntryContainer.NodeChildren(i);
    hLegendIconLine = hLegendEntry.Icon.Transform.Children.Children;
    hLegendIconLine.LineWidth = 2;
end

% Saving figure
name_fig = 'rate_frequency';
if driv_curr == 1
    name_fig = strcat(name_fig, '_excitatory');
else
    name_fig = strcat(name_fig, '_inhibitory');
end
if Jii ~= 0
    name_fig = [name_fig '_Jii', num2str(Jii)];
end
if tau_si ~= 1
    name_fig = [name_fig '_tau_si', num2str(tau_si)];
end
print(gcf, '-depsc', '-tiff', fullfile(folder, [name_fig, '_tau8.eps']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Time and phase differences %%%%%%%%%%%%%%%%%%%%%%%
figure; box on; hold on; set(gca, 'Fontsize', 16);
yyaxis left
plot(I_ext(pos_osc), diff_phases(pos_osc), 'Linewidth', 1.5); % Time difference
ylabel('Time difference (ms)', 'Fontsize', 18, 'Fontweight', 'bold');

yyaxis right
plot(I_ext(pos_osc), diff_phases_norm(pos_osc), 'Linewidth', 1.5); % Phase difference
ylabel('Phase difference', 'Fontsize', 18, 'Fontweight', 'bold');
if driv_curr == 1
%     \bf Excitatory driving current, $\bf \bar I_e^{ext}$
    xlabel('Excitatory driving current, I_e^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
else
%     set(gca, 'XDir', 'Reverse');
    xlabel('Inhibitory driving current, I_i^{ext}', 'Fontsize', 18, 'Fontweight', 'bold');
end
xlim([min(mn, mx)-0.5 max(mn, mx)+0.5]);

% Saving figure
name_fig = 'phase_difference';
if driv_curr == 1
    name_fig = strcat(name_fig, '_excitatory');
else
    name_fig = strcat(name_fig, '_inhibitory');
end
if Jii ~= 0
    name_fig = [name_fig '_Jii', num2str(Jii)];
end
if tau_si ~= 1
    name_fig = [name_fig '_tau_si', num2str(tau_si)];
end
print(gcf, '-depsc', '-tiff', fullfile(folder, [name_fig, '_tau8.eps']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%