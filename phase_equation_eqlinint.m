function dtheta = phase_equation_eqlinint(~, x, t, Z, p, A, T)

    % Phase variable
    theta = x(1);
    
    % Linear interpolation of iPRC
    theta = mod(theta,T);
    len = length(Z)-1;
    wdth = T/len;
    ip1 = floor(theta/wdth);
    ip2 = ceil(theta/wdth);
    if ip2 == ip1
        wght = 0;
    else
        t0 = t(ip1+1);
        t1 = t(ip2+1);
        wght = (theta - t0)/(t1 - t0);
    end
    Z0 = Z(ip1+1);
    Z1 = Z(ip2+1);
    Ztheta = Z0*(1 - wght) + Z1*wght;

    % Phase equation
    dtheta = 1 + A*Ztheta*p;

end

