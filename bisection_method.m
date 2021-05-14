function x_mid = bisection_method(F, x0, x1, tol)
    
    iter = 1; % Current iteration
    maxiter = 80; % Max. number of iterations
    stopp = false; % Stopping criterion: no solution found
    x_mid = (x0+x1)/2; % Mid point
    while iter <= maxiter && abs(F(x_mid)) >= tol && ~stopp
        % Bisecting interval
        if F(x_mid)*F(x0) < 0
            x1 = x_mid;
        else
            x0 = x_mid;
        end
        
        % Updating mid point
        x_mid = (x0+x1)/2;
        
        % Increment iteration counter
        iter = iter + 1;
        
        % Stop if interval's length is less than tol
        if abs((x1-x0)/2) < tol
            stopp = true;
        end
    end
end

