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
                            p = 0.795948141314470;
                            % p = 0.775261585047620; % tau_si = 5 % p = 0.701042922662681; % tau_se = tau_si
                        elseif k == 2
                            p = 0.843293507218629;
                            % p = 0.825277552766993; % tau_si = 5 % p = 0.750371394416885; % tau_se = tau_si
                        else % k = 0.5
                            p = 0.893306217680769;
                            % p = 0.877294159195140; % tau_si = 5 % p = 0.811450734080839; % tau_se = tau_si
                        end
                    else % 2:1 plateau
                        if k == 20
                            p = 1.703802535023349;
                            % p = 1.677785190126751; % tau_si = 5 % p = 1.567423230974633; % tau_se = tau_si
                        elseif k == 2
                            p = 1.811874583055370;
                            % p = 1.787191460973983; % tau_si = 5 % p = 1.683578104138852; % tau_se = tau_si
                        else % k = 0.5
                            p = 1.835890593729153;
                            % p = 1.809206137424950; % tau_si = 5 % p = 1.692923898531375; % tau_se = tau_si
                        end
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.437190298266056;
                        % p = 0.429151088429564; % tau_si = 5 % p = 0.402085845325362; % tau_se = tau_si
                    elseif k == 2
                        p = 0.444525495800503;
                        % p = 0.437153643264664; % tau_si = 5 % p = 0.409787685203752; % tau_se = tau_si
                    else % k = 0.5
                        p = 0.454528037892931;
                        % p = 0.447823716378130; % tau_si = 5 % p = 0.419082367326527; % tau_se = tau_si
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
                            p = 1.018012008005337;
                            % p = 1.010673782521681; % tau_si = 5 % p = 1.001575048066344; % tau_se = tau_si
                        elseif k == 2
                            p = 0.9944;
                            % p = 0.9876; % tau_si = 5 % p = 0.964812989106635; % tau_se = tau_si
                        else % k = 0.5
                            p = 0.943985764282404;
                            % p = 0.933312043040837; % tau_si = 5 % p = 0.883816473465305; % tau_se = tau_si
                        end
                    else % 2:1 plateau
                        if k == 20
                            p = 2.025350233488993;
                            % p = 2.011340893929286; % tau_si = 5 % p = 1.995994659546062; % tau_se = tau_si
                        elseif k == 2
                            p = 1.905270180120080;
                            % p = 1.881921280853903; % tau_si = 5 % p = 1.781041388518024; % tau_se = tau_si
                        else % k = 0.5
                            p = 1.843895930620413;
                            % p = 1.816544362908606; % tau_si = 5 % p = 1.699599465954606; % tau_se = tau_si
                        end
                    end
                else % 1:2 plateau
                    if k == 20
                        p = 0.479867811193748;
                        % p = 0.473832019592203; % tau_si = 5 % p = 0.433974600241343; % tau_se = tau_si
                    elseif k == 2
                        p = 0.473866285938291;
                        % p = 0.467830103465879; % tau_si = 5 % p = 0.432360484644778; % tau_se = tau_si
                    else % k = 0.5
                        p = 0.464530579985359;
                        % p = 0.457826909922004; % tau_si = 5 % p = 0.426385331851565; % tau_se = tau_si
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
