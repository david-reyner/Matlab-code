function dx = phase_equation_variationals_vonmises(t, x, Z, dZ, ddZ, mu, k, A, T, T_prime)

    % Variables
    theta = x(1); % \theta
    Phi_theta = x(2); % \frac{\partial \varphi_1}{\partial \theta}
    Phi_T = x(3); % \frac{\partial \varphi_1}{\partial T}
    Phi_A = x(4); % \frac{\partial \varphi_1}{\partial A}
    Phi_theta2 = x(5); % \frac{\partial^2 \varphi_1}{\partial \theta^2}
    Phi_Ttheta = x(6); % \frac{\partial^2 \varphi_1}{\partial T \partial \theta}
    Phi_Atheta = x(7); % \frac{\partial^2 \varphi_1}{\partial A \partial \theta}
    
    % mod T
    theta = mod(theta,T);

    % Modified Bessel function of first kind (order 0)
    I0k = besseli(0, k);
    
    % Von Mises distribution
    vonmises = vonmises_dist(t,mu,k,T_prime);
    
    % Von Mises' partial derivative with respect to T_prime
    % --> If normalizing by T
%     dvonmises = (1/(T_prime^2*I0k))*exp(k*cos(2*pi*(t-mu)/T_prime))*((2*pi*k*(t-mu)/T_prime)*sin(2*pi*(t-mu)/T_prime)-1);
    
    % --> If normalizing by 2\pi
%     dvonmises = (1/(T_prime^2*I0k))*exp(k*cos(2*pi*(t-mu)/T_prime))*(k*(t-mu)*sin(2*pi*(t-mu)/T_prime));
    
    % --> Without normalizing
    dvonmises = (1/(T_prime^2*I0k))*exp(k*cos(2*pi*(t-mu)/T_prime))*(2*pi*k*(t-mu)*sin(2*pi*(t-mu)/T_prime));
    
    % Partial derivatives of the phase equation's vector field
    DF_theta = A*dZ(theta)*vonmises;
    DF_T = A*Z(theta)*dvonmises;
    DF_A = Z(theta)*vonmises;
    DF_theta2 = A*ddZ(theta)*vonmises;
    DF_Ttheta = A*dZ(theta)*dvonmises;
    DF_Atheta = dZ(theta)*vonmises;
    
    % System of equations
    dx1 = 1 + A*Z(theta)*vonmises;
    dx2 = DF_theta*Phi_theta;
    dx3 = DF_theta*Phi_T + DF_T;
    dx4 = DF_theta*Phi_A + DF_A;
    dx5 = DF_theta2*(Phi_theta)^2 + DF_theta*Phi_theta2;
    dx6 = DF_theta2*Phi_T*Phi_theta + DF_Ttheta*Phi_theta + DF_theta*Phi_Ttheta;
    dx7 = DF_theta2*Phi_A*Phi_theta + DF_Atheta*Phi_theta + DF_theta*Phi_Atheta;
    
    dx = [dx1; dx2; dx3; dx4; dx5; dx6; dx7;];

end