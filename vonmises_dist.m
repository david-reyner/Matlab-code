function f = vonmises_dist(x, mu, k, T)
    % T-periodic Von Mises distribution (without normalizing)
    I0k = besseli(0, k); % Modified Bessel function of first kind (order 0)
    f = exp(k*cos(2*pi*(x-mu)/T))./I0k;
end

