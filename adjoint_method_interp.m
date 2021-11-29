function dx = adjoint_method_interp(t, x, parm, T, F)

    % Parameters
    tau_e = parm(1); tau_i = parm(2); delta_e = parm(3); delta_i = parm(4);
    eta_e = parm(5); eta_i = parm(6); tau_se = parm(7); tau_si = parm(8);
    Jee = parm(9); Jei = parm(10); Jie = parm(11); Jii = parm(12);
    Ie_ext = parm(13); Ii_ext = parm(14);
    
    % Time
    tau = mod(t,T);
    
    % Interpolating variables along the limit cycle
    xinterp = F({tau,1:8});
    re = xinterp(1); Ve = xinterp(2); See = xinterp(3); Sei = xinterp(4);
    ri = xinterp(5); Vi = xinterp(6); Sie = xinterp(7); Sii = xinterp(8);
    
    % Infinitesimal-PRC variables
    z1 = x(1); z2 = x(2); z3 = x(3); z4 = x(4);
    z5 = x(5); z6 = x(6); z7 = x(7); z8 = x(8);
    
    Z = [z1 z2 z3 z4 z5 z6 z7 z8]';
    
    % Linearization of vector field around limit cycle
    DF = [2*Ve/tau_e 2*re/tau_e 0 0 0 0 0 0;
          -2*(tau_e*pi)^2*re/tau_e 2*Ve/tau_e 1 -1 0 0 0 0;
          1/tau_se*Jee 0 -1/tau_se 0 0 0 0 0;
          0 0 0 -1/tau_si 1/tau_si*Jei 0 0 0;
          0 0 0 0 2*Vi/tau_i 2*ri/tau_i 0 0;
          0 0 0 0 -2*(tau_i*pi)^2*ri/tau_i 2*Vi/tau_i 1 -1;
          1/tau_se*Jie 0 0 0 0 0 -1/tau_se 0;
          0 0 0 0 1/tau_si*Jii 0 0 -1/tau_si];
    
    % Adjoint equations
    dx = -DF'*Z;

end
