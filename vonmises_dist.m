function f = vonmises_dist(x, mu, k, T)
    
    % T-periodic Von Mises distribution normalized by 2\pi
    I0k = besseli(0, k); % Modified Bessel function of first kind (order 0)
    f = exp(k*cos(2*pi*(x-mu)/T))./(2*pi*I0k);

end

