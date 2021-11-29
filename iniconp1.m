function p = iniconp1(~, ~, k, p_ord, q_ord, branch)
    % This function returns the point p = T/T* lying at one of the ends of
    % the p_ord : q_ord plateau of the Devil's staircases corresponding to
    % factor k. Branch determines which end we are looking upon: left or 
    % right. The cases correspond to a PING setting when both excitatory
    % and inhibitory populations receive the same external input
    
    if branch == 1 % Left branch
        if q_ord == 1
            if p_ord == 1 % 1:1 plateau
                if k == 20
                    p = 0.785278763082547;
                elseif k == 2
                    p = 0.829956784428725;
                else % k = 0.5
                    p = 0.883970511727836;
                end
            else % 2:1 plateau
                if k == 20
                    p = 1.677118078719146;
                elseif k == 2
                    p = 1.797198132088059;
                else % k = 0.5
                    p = 1.825216811207472;
                end
            end
        else % 1:2 plateau
            if k == 20
                p = 0.439857642824036;
            elseif k == 2
                p = 0.445859168079493;
            else % k = 0.5
                p = 0.4536;
            end
        end
    else % Right branch
        if q_ord == 1
            if p_ord == 1 % 1:1 plateau
                if k == 20
                    p = 1.035356904603069;
                elseif k == 2
                    p = 1.0071;
                else % k = 0.5
                    p = 0.944652600421899;
                end
            else % 2:1 plateau
                if k == 20
                    p = 2.064709806537691;
                elseif k == 2
                    p = 1.909272848565710;
                else % k = 0.5
                    p = 1.834556370913943;
                end
            end
        else % 1:2 plateau
            if k == 20
                p = 0.465864252264349;
            elseif k == 2
                p = 0.464530579985359;
            else % k = 0.5
                p = 0.459862727008892;
            end
        end
    end
end