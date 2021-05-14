function dtheta = phase_equation(~, x, func, p, A, T)

    % Phase variable
    theta = x(1);
    
    % Linear interpolation of iPRC
    theta = mod(theta,T);
    lin_intpol = func(theta);
    
    % Phase equation
    dtheta = 1 + A*lin_intpol*p;

end

