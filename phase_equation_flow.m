function theta_n1 = phase_equation_flow(lin_intpol, p, A, T, T_prime, q, inc, options_ode)

    % Stroboscopic map: Flow of the phase equation after q periods of the forcer
    [~, theta] = ode45(@(l, theta) phase_equation(l, theta, lin_intpol, ...
        p(l), A, T), [0 q*T_prime], inc, options_ode);
    theta_n1 = mod(theta(end),T);

end

