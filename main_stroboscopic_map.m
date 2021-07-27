%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Stroboscopic map and fixed points %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% for Von Mises distributions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script computes the composition of the stroboscopic map P q times in
% order to spot its q-periodic points (fixed points if q = 1). After that,
% such fixed/periodic points are determined using a bisection method.

format long;

% program_invocation_name (only in octave)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Poincare Phase Map %%%%%%%%%%%%%%%%%%%%%%%%%%%
%_________________________________________________________________________%

% Integration tolerances (phase equation)
options_ode_str = odeset('AbsTol', 1e-12, 'RelTol', 1e-10, 'InitialStep', 1e-3);

% Period relation between the forcing and the oscillator
T_prime = p*T;

% Periodic perturbation with period Tp (Von Mises distribution)
pt = @(t) vonmises_dist(t,0,k,T_prime);

% Linear interpolation of the component coord of the iPRC (interp1q or griddedInterpolant)
% lin_intpol = @(l) interp1q(t, Z(:,coord), l); % Quick 1D linear interpolation (not recommended)
lin_intpol = griddedInterpolant(t, Z(:,coord)); % GriddedInterpolant

% Poincare Phase Map
num = 7500;
theta_n = linspace(0, T, num)';
theta_n1 = zeros(num,1);

n = q; % Fixed or q-periodic points

try

for i = 1:length(theta_n)
    % Stroboscopic map: Flow of the phase equation after n periods of the forcer
    [~, theta] = ode45(@(l, theta) phase_equation(l, theta, lin_intpol, ...
        pt(l), A, T), [0 n*T_prime], theta_n(i), options_ode_str);
    theta_n1(i) = mod(theta(end),T);
end

catch
    file = fopen('errorFileEixam_StroboscopicMap.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = lasterr;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

% Saving data
if type == 1
    name_file = ['stroboscopic_map_ping_coord_', num2str(coord)];
else
    name_file = ['stroboscopic_map_ing_coord_', num2str(coord)];
end
name_file = strcat(name_file, ['_A', num2str(A), '_k', num2str(k), ...
    '_p', num2str(round(1000*p)/1000), 'approx_', num2str(n), ...
    'periodicpoints.txt']);

file = fopen(fullfile(folder, name_file), 'w');
for i = 1:length(theta_n)
    fprintf(file, '%16.15f %16.15f\r\n', [theta_n(i) theta_n1(i)]);
end
fclose(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bisection method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_________________________________________________________________________%

% Integration tolerances (stroboscopic map)
options_ode_bis = odeset('AbsTol', 1e-16, 'RelTol', 5e-14, 'InitialStep', 1e-3);

% Restricting to interval containing fixed/periodic points of the Poincare Phase Map
dif = theta_n1 - theta_n;
sta_fxd_pts = find(dif(1:end-1) >= 0 & dif(2:end) <= 0); % Change of sign: from positive to negative --> stability
sta_fxd_pts = sta_fxd_pts(abs(dif(sta_fxd_pts) - dif(sta_fxd_pts+1)) < 1); % Removing the sign change due to the modulus
usta_fxd_pts = find(dif(1:end-1) <= 0 & dif(2:end) >= 0); % Change of sign: from negative to positive --> instability
usta_fxd_pts = usta_fxd_pts(abs(dif(usta_fxd_pts) - dif(usta_fxd_pts+1)) < 1); % Removing the sign change due to the modulus

try

% Bisection method to find fixed/periodic points of the Poincare Phase Map
% options_ode_bis.MaxStep = 0.1*abs(n*T_prime);
F = @(x) phase_equation_flow(lin_intpol, pt, A, T, T_prime, n, x, options_ode_bis) - x; % Function to find its zeros
tol = 5e-12; % Bisection method tolerance

% Files where stable/unstable fixed/periodic points will be stored
if type == 1
    name_file1 = ['bisect_mth_stb_ping_coord_', num2str(coord), ...
        '_A', num2str(A), '_k', num2str(k), '_p', num2str(round(1000*p)/1000), ...
        'approx_', num2str(n), 'periodicpoints.txt'];
    name_file2 = ['bisect_mth_ustb_ping_coord_', num2str(coord), ...
        '_A', num2str(A), '_k', num2str(k), '_p', num2str(round(1000*p)/1000), ...
        'approx_', num2str(n), 'periodicpoints.txt'];
else
    name_file1 = ['bisect_mth_stb_ing_coord_', num2str(coord), ...
        '_A', num2str(A), '_k', num2str(k), '_p', num2str(round(1000*p)/1000), ...
        'approx_', num2str(n), 'periodicpoints.txt'];
    name_file2 = ['bisect_mth_ustb_ing_coord_', num2str(coord), ...
        '_A', num2str(A), '_k', num2str(k), '_p', num2str(round(1000*p)/1000), ...
        'approx_', num2str(n), 'periodicpoints.txt'];
end

file1 = fopen(fullfile(folder, name_file1), 'w'); % To save stable fixed/periodic points
file2 = fopen(fullfile(folder, name_file2), 'w'); % To save unstable fixed/periodic points

i = 1;
xs = zeros(length(sta_fxd_pts),1);
xu = zeros(length(usta_fxd_pts),1);
while i <= length(sta_fxd_pts)
    % Stable fixed/periodic points
    x0 = theta_n(sta_fxd_pts(i)); % Left endpoint
    x1 = theta_n(sta_fxd_pts(i)+1); % Right endpoint
    xs(i) = bisection_method(F, x0, x1, tol)
    F(xs(i))
    
    % Unstable fixed/periodics points
    x0 = theta_n(usta_fxd_pts(i)); % Left endpoint
    x1 = theta_n(usta_fxd_pts(i)+1); % Right endpoint
    xu(i) = bisection_method(F, x0, x1, tol)
    F(xu(i))
    
    % Saving stable/unstable points
    fprintf(file1, '%16.15f %16.15f\r\n', [xs(i), F(xs(i))+xs(i)]);
    fprintf(file2, '%16.15f %16.15f\r\n', [xu(i), F(xu(i))+xu(i)]);
    i = i + 1;
end

fclose(file1);
fclose(file2);

catch
    file = fopen('errorFileEixam_BisectionMethod.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = lasterr;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
