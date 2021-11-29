function dx = full_synaptic_firing_rate_equations(~, x, parm)
    % G. Dumont and B. Gutkin research article (PLOS Computational Biology)
    
    % Variables
    re = x(1); Ve = x(2); See = x(3); Sei = x(4); 
    ri = x(5); Vi = x(6); Sie = x(7); Sii = x(8);
    
    % Parameters
    tau_e = parm(1); tau_i = parm(2); delta_e = parm(3); delta_i = parm(4);
    eta_e = parm(5); eta_i = parm(6); tau_se = parm(7); tau_si = parm(8);
    Jee = parm(9); Jei = parm(10); Jie = parm(11); Jii = parm(12);
    Ie_ext = parm(13); Ii_ext = parm(14);
    
    %%%%%%%%%%%%%%%%%%%%%%% Exact firing rate model %%%%%%%%%%%%%%%%%%%%%%%
    
    % Synaptic interactions
    Ie = Ie_ext + tau_e*See - tau_e*Sei;
    Ii = Ii_ext + tau_i*Sie - tau_i*Sii;
    
    % Excitatory Cells
    dre = 1/tau_e*(delta_e/(pi*tau_e) + 2*re*Ve);
    dVe = 1/tau_e*(Ve^2 + eta_e - (tau_e*pi*re)^2 + Ie);
    dSee = 1/tau_se*(-See + Jee*re);
    dSei = 1/tau_si*(-Sei + Jei*ri);
    
    % Inhibitory Cells
    dri = 1/tau_i*(delta_i/(pi*tau_i) + 2*ri*Vi);
    dVi = 1/tau_i*(Vi^2 + eta_i - (tau_i*pi*ri)^2 + Ii);
    dSie = 1/tau_se*(-Sie + Jie*re);
    dSii = 1/tau_si*(-Sii + Jii*ri);
    
    dx = [dre; dVe; dSee; dSei;
          dri; dVi; dSie; dSii];

end