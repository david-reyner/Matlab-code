function [G, DG] = arnold_boundaryq_vonmises(w, Z, dZ, ddZ, mu, k, q, T, options_ode)
    
    % Variables
    theta = w(1); T_prime = w(2); A = w(3);
    
    % Integration of reduced set of variationals equations + phase equation
    inc = [theta; 1; 0; 0; 0; 0; 0];    
    [t, x] = ode45(@(t, x) phase_equation_variationals_vonmises(t, x, ...
        Z, dZ, ddZ, mu, k, A, T, T_prime), [0 q*T_prime], inc, options_ode);
    
    y = x(end,:); % Solution of variational equations of phase equation
    y(1) = mod(y(1),T);

    % Periodic point condition
    dx1 = y(1) - theta; % Phi_{qT}(\theta) mod T* - \theta = 0
    if abs(dx1) > T/2
        if dx1 > 0
            dx1 = dx1 - T;
        else
            dx1 = dx1 + T;
        end
    end

    % Saddle-node condition (on maps)
    dx2 = y(2) - 1; % \frac{\partial \Phi_{qT}}{\partial \theta} - 1 = 0
    
    % Arnold tongue's boundary system
    G = [dx1; dx2];
    
    % Von Mises distribution
    vonmises = vonmises_dist(t(end),mu,k,T_prime);
    
    % Differential Arnold tongue's boundary system
    F1 = 1 + A*Z(y(1))*vonmises; % F_1 evaluated at t = qT;
    DF_theta = A*dZ(y(1))*vonmises; % DF_{\theta} evaluated at t = qT;
    DG = [y(2)-1      y(3)+q*F1          y(4);
          y(5)    y(6)+q*DF_theta*y(2)   y(7)];

end
