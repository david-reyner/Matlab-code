function p = iniconp(type, coord, k, p_ord, q_ord, branch)
    % This function returns the point p = T/T* lying at one of the ends of
    % the p_ord : q_ord plateau of the Devil's staircases corresponding to
    % factor k. Branch determines which end we are looking upon: left or 
    % right. Type and coord determine which case we are exploring: PING-Ve,
    % PING-Vi or ING-Vi.
    
    if branch == 1 % Left branch
        if type == 1
            if coord == 1 % PING - Ve
                if q_ord == 1
                    if p_ord == 1 % 1:1 plateau
                        if k == 20
                            p = 0.701042922662681;
                        elseif k == 2
                            p = 0.750371394416885;
                        else % k = 0.5
                            p = 0.811450734080839;
                        end
                    else % 2:1 plateau
                        if k == 20
                            p = 1.567423230974633;
                        elseif k == 2
                            p = 1.683578104138852;
                        else % k = 0.5
                            p = 1.692923898531375;
                        end
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.402085845325362;
                    elseif k == 2
                        p = 0.409787685203752;
                    else % k = 0.5
                        p = 0.419082367326527;
                    end
                end
            else % PING - Vi
                if q_ord == 1 % 1:1 plateau
                    if k == 20
                        p = 0.901676339009060;
                    elseif k == 2
                        p = 0.9786;
                    else % k = 0.5
                        p = 0.991369223743137;
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.483136430736813;
                    elseif k == 2
                        p = 0.492112012576907;
                    else % k = 0.5
                        p = 0.4958;
                    end
                end
            end
        else % ING - Vi
            if q_ord == 1 % 1:1 plateau
                if k == 20
                    p = 0.915555404454237;
                elseif k == 2
                    p = 0.943951572399698;
                else
                    p = 0.961260645629203;
                end
            else % 1:2 plateau
                if k == 20
                    p = 0.478264709164509;
                elseif k == 2
                    p = 0.467958381989373;
                else % k = 0.5
                    p = 0.4696;
                end
            end
        end
    else % Right branch
        if type == 1
            if coord == 1 % PING - Ve
                if q_ord == 1
                    if p_ord == 1 % 1:1 plateau
                        if k == 20
                            p = 1.001575048066344;
                        elseif k == 2
                            p = 0.964812989106635;
                        else % k = 0.5
                            p = 0.883816473465305;
                        end
                    else % 2:1 plateau
                        if k == 20
                            p = 1.995994659546062;
                        elseif k == 2
                            p = 1.781041388518024;
                        else % k = 0.5
                            p = 1.699599465954606;
                        end
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.433974600241343;
                    elseif k == 2
                        p = 0.432360484644778;
                    else % k = 0.5
                        p = 0.426385331851565;
                    end
                end
            else % PING - Vi
                if q_ord == 1 % 1:1 plateau
                    if k == 20
                        p = 1.052069425901202;
                    elseif k == 2
                        p = 1.0115;
                    else % k = 0.5
                        p = 1.0025;
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.507052996923798;
                    elseif k == 2
                        p = 0.502070600565595;
                    else % k = 0.5
                        p = 0.4991;
                    end
                end
            end
        else % ING - Vi
            if q_ord == 1 % 1:1 plateau
                if k == 20
                    p = 1.002670226969292;
                elseif k == 2
                    p = 0.994230308923498;
                else
                    p = 0.977745477276351;
                end
            else % 1:2 plateau
                if k == 20
                    p = 0.481596352145578;
                elseif k == 2
                    p = 0.471914081751533;
                else
                    p = 0.4709;
                end
            end
        end
    end
end
