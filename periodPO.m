function [value, isterminal, direction] = periodPO(t, y, val, coord)
    
    value = y(coord) - val; % function to find a zero
    isterminal = 1; % 1 if to terminate at a zero of event function and
                    % 0 otherwise
    direction = 0; % 1 if function is increasing
                   % -1 if function is decreasing

end