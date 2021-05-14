%%
%%%%%%%%%%%%%%%%%%% Arnold tongues (Von Mises stimulus) %%%%%%%%%%%%%%%%%%%

% Here we compute the boundaries of the Arnold tongues arising from a
% PING/ING oscillator externally forced by a Von Mises stimulus. To do so
% we determine the left/right saddle-node bifurcations of the stroboscopic 
% map and continuate them using a corrector-predictor method, based on a 
% Newton method and a tangent approximation. The initial condition consists
% of one of the stable fixed/periodic points of the stroboscopic map (found
% using a bisection method), together with the periods relation T/T* and 
% the amplitude A.

% The next set of parameters must be initialized before running:
type
coord
k
p_ord
q
branch

% Arguments from command window (in case Octave is used)
% args = argv()
% type = args{1};
% coord = args{2};
% k = args{3};
% p_ord = args{4};
% q = args{5};
% branch = args{6};

% Script arguments:
%{
    type --> Gamma mechanism (1 for PING, 2 for ING)
    coord --> Component where perturbation is applied (1 for Ve, 2 for Vi)
    k --> Von Mises parameter (input concentration factor)
    p_ord --> Oscillator's revolutions (order of the Arnold tongue)
    q --> Period of q-periodic points
    branch --> Branch to compute (1 for the left, 2 for the right)
%}

format long;

% Changing tolerances
options_ode = odeset('AbsTol', 1e-14, 'RelTol', 5e-14, 'InitialStep', 1e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd; % Current folder

% --> Oscillator period
if type == 1
    name_file = 'period_ping.txt';
else
    name_file = 'period_ing.txt';
end
file = fopen(fullfile(folder, name_file), 'r');
formatSpec = '%f';
T = fscanf(file, formatSpec);
fclose(file);

% --> Infinitesimal phase response curve
if type == 1
    name_file = 'iPRC_ping.txt';
else
    name_file = 'iPRC_ing.txt';
end
file = fopen(fullfile(folder, name_file), 'r');
formatSpec = '%f %f %f'; sizeZ = [3 Inf];
res = fscanf(file, formatSpec, sizeZ); res = res';
t = res(:,1); Z = res(:,2:3);
fclose(file);

% --> First and second derivatives of the iPRC
if type == 1
    name_file_d1 = 'dPRC_ping.txt';
    name_file_d2 = 'ddPRC_ping.txt';
else
    name_file_d1 = 'dPRC_ing.txt';
    name_file_d2 = 'ddPRC_ing.txt';
end

file = fopen(fullfile(folder, name_file_d1), 'r');
res = fscanf(file, formatSpec, sizeZ); res = res';
dd1 = res(:,1); dZ = res(:,2:3); % First derivative
fclose(file);

file = fopen(fullfile(folder, name_file_d2), 'r');
res = fscanf(file, formatSpec, sizeZ); res = res';
dd2 = res(:,1); ddZ = res(:,2:3); % Second derivative
fclose(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear interpolation of the component coord of the iPRC (interp1q or griddedInterpolant)
% PRC = @(l) interp1q(t, Z(:,coord), l); % Quick 1D linear interpolation (not recommended)
% dPRC = @(l) interp1q(dd1, dZ(:,coord), l); % Quick 1D linear interpolation (not recommended)
% ddPRC = @(l) interp1q(dd2, ddZ(:,coord), l); % Quick 1D linear interpolation (not recommended)
PRC = griddedInterpolant(t, Z(:,coord)); % GriddedInterpolant function
dPRC = griddedInterpolant(dd1, dZ(:,coord)); % GriddedInterpolant function
ddPRC = griddedInterpolant(dd2, ddZ(:,coord)); % GriddedInterpolant function

% Initial conditions (PING - Ve, PING - Vi & ING - Vi cases)
A = 1;
if branch == 1
    
    % Choose a point p = Tp/T lying on left edge of the p_ord : q plateau
    p = iniconp(type, coord, k, p_ord, q, branch);
    if type == 1
        name_file = ['bisect_mth_stb_ping_coord_', num2str(coord)];
    else
        name_file = ['bisect_mth_stb_ing_coord_', num2str(coord)];
    end
    name_file = strcat(name_file, ['_A', num2str(A), '_k', num2str(k), ...
        '_p', num2str(round(1000*p)/1000), 'approx_', num2str(q), ...
        'periodicpoints.txt']);
    
    % Check if the file with fixed/periodic points exists or not
    if exist(name_file, 'file') ~= 0
        % Loading fixed/periodic points of the corresponding stroboscopic map
        file = fopen(fullfile(folder, name_file), 'r');
        formatSpec = '%f %f'; size_stb = [2 Inf];
        res1 = fscanf(file, formatSpec, size_stb); res1 = res1';
        xs = res1(:,1);
        fclose(file);
    else
        % Computing fixed/periodic points of the corresponding stroboscopic map
        main_stroboscopic_map
    end
    
    [~, ind] = min(abs(xs-T/2));
    xl = xs(ind); % Choose stable fixed/periodic point closer to T/2
    Tl = p*T; % Choose the leftmost point on the "p_ord : q plateau"
    
    % Initial condition for the left branch
    inc_left = [xl; Tl; A];
else
    
    % Choose a point p = Tp/T lying on right edge of the p_ord : q plateau
    p = iniconp(type, coord, k, p_ord, q, branch);
    if type == 1
        name_file = ['bisect_mth_stb_ping_coord_', num2str(coord)];
    else
        name_file = ['bisect_mth_stb_ing_coord_', num2str(coord)];
    end
    name_file = strcat(name_file, ['_A', num2str(A), '_k', num2str(k), ...
        '_p', num2str(round(1000*p)/1000), 'approx_', num2str(q), ...
        'periodicpoints.txt']);
    
    % Check if the file with fixed/periodic points exists or not
    if exist(name_file, 'file') ~= 0
        % Loading fixed/periodic points of the corresponding stroboscopic map
        file = fopen(fullfile(folder, name_file), 'r');
        formatSpec = '%f %f'; size_stb = [2 Inf];
        res1 = fscanf(file, formatSpec, size_stb); res1 = res1';
        xs = res1(:,1);
        fclose(file);
    else
        % Computing fixed/periodic points of the corresponding stroboscopic map
        main_stroboscopic_map
    end
    
    [~, ind] = min(abs(xs-T/2));
    xr = xs(ind); % Choose stable fixed/periodic point closer to T/2
    Tr = p*T; % Choose the rightmost point on the "p_ord : q plateau"
    
    % Initial condition for the right branch
    inc_right = [xr; Tr; A];
end

% Arc step
delta_s = 0.025;

% Building filename where data will be stored
if branch == 1
    % File to save the left branch
    name_file = ['left_branch_', num2str(p_ord), num2str(q)];
else
    % File to save the right branch
    name_file = ['right_branch_', num2str(p_ord), num2str(q)];
end

if type == 1
    name_file = strcat(name_file, '_ping');
else
    name_file = strcat(name_file, '_ing');
end

name_file = strcat(name_file, '_coord_');
if coord == 1
    name_file = strcat(name_file, 'Ve');
else
    name_file = strcat(name_file, 'Vi');
end

% Opening/creating file (discarding existing contents)
name_file = strcat(name_file, ['_k', num2str(k), '.txt']);
fclose(fopen(fullfile(folder, name_file), 'w'));

try

j = 1;
while j <= 2
    % Initial point near the real left/right Arnold tongue's boundary
    if branch == 1
        w0 = inc_left;
    else
        w0 = inc_right;
    end

    % Number of continuation steps
    arnold_iter = 1; max_arnold_iter = 3000;
    
    % Number of reductions of the step along the tangent line
    reduction = 0;

    while arnold_iter <= max_arnold_iter && ((j==1)*(w0(3) > 1e-6) || (j==2)*(w0(3) < 2.11))
        fprintf('\nContinuation iteration: %d\n', arnold_iter);
        
        % CONTINUATION: Correction step (Modified Newton Method)
        tol = 5e-11; % Tolerance
        stopp = false; % Stopping criterion
        control_norm = true; % Control norm
        iter = 1; maxiter = 40; % Maximum number of iterations
        [G, DG] = arnold_boundaryq_vonmises(w0, PRC, dPRC, ddPRC, 0, k, q, T, options_ode);

        while iter <= maxiter && ~stopp && control_norm
            inc_w = -DG'*((DG*DG')\G); % Correction
            w1 = w0 + inc_w; % Newton step
            [G, DG] = arnold_boundaryq_vonmises(w1, PRC, dPRC, ddPRC, 0, k, q, T, options_ode);
            fprintf('  \nIteration: %d \tError norm: %e\n', iter, norm(G));
            stopp = norm(G) < tol; % Update stopping criterion
            control_norm = norm(G) < 1; % Control norm: Prevent uncontrolled growth
            w0 = w1; % Update initial point
            iter = iter + 1;
        end

        if control_norm % No uncontrolled growth
            if ~stopp
                fprintf('\nNo convergence (Exceeded max. number of iterations)\n\n');
            else
                fprintf('\nConvergence at iteration: %d\n', iter-1);
                fprintf('Norm of the error: %e\n', norm(G));
            end

            % Saving data
            file = fopen(fullfile(folder, name_file), 'a');
            fprintf(file, '%16.15f %16.15f %16.15f %d %d %16.15f\r\n', w0, stopp, iter-1, norm(G));
            fclose(file);
            
            % CONTINUATION: Prediction step (Tangent approximation)
            aux = w0; % In case Newton fails to control the norm of G
            reduction = 0; % Number of reductions on the delta_s step
            tg = tangent_vector(DG); % Tangent vector along Arnold tongue's boundary
            if j == 1 
                % Lower branch
                aux_w = w0 - delta_s*tg;
                if aux_w(3) > w0(3)
                    aux_w = w0 + delta_s*tg;
                end
            else
                % Upper branch
                aux_w = w0 + delta_s*tg;
                if aux_w(3) < w0(3)
                    aux_w = w0 - delta_s*tg;
                end
            end
            w0 = aux_w; % Update initial seed

            arnold_iter = arnold_iter + 1; % Increment continuation steps counter
        else % Uncontrolled growth: Reduce step delta_s
            reduction = reduction + 1; % Counting the number of delta_s reductions
            if reduction > 10
                fprintf('\nStep delta_s too small: %f\n', delta_s/(2^reduction));
                return
            end
            w0 = aux; % Before doing a prediction (so that we can control better the step along the tangent line)
            new_delta = delta_s/(2^reduction);
            fprintf('\nStep reduction. From %f to %f\n', delta_s/(2^(reduction-1)), new_delta);
            if j == 1 
                % Lower branch
                aux_w = w0 - new_delta*tg;
                if aux_w(3) > w0(3)
                    aux_w = w0 + new_delta*tg;
                end
            else
                % Upper branch
                aux_w = w0 + new_delta*tg;
                if aux_w(3) < w0(3)
                    aux_w = w0 - new_delta*tg;
                end
            end
            w0 = aux_w; % Update initial seed
        end

    end
    
    j = j + 1;

end

catch
    file = fopen('errorFileEixam_ArnoldTongues.txt', 'a');
    datenow = datestr(now);
    datenow = strcat(datenow, ': ');
    lastmsg = lasterr;
    fprintf(file, '%s %s\r\n', datenow, lastmsg);
    fclose(file);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
