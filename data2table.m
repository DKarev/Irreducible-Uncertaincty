function tbl = data2table(data,standardize,no_timeouts)
    
    if nargin < 2; standardize = 0; end
    if nargin < 3; no_timeouts = 1; end

    % take care of timeouts
    %
    if no_timeouts
        % get rid of them
        F = fieldnames(data);
        for s = 1:length(data)
            which_trials = ~data(s).timeout;
            for i = 1:length(F)
                if size(data(s).(F{i}), 1) == length(which_trials)
                    data(s).(F{i}) = data(s).(F{i})(which_trials,:);
                end
            end
        end
    end
    
    RS = []; SS = []; C = []; V = []; S = []; RU = []; TU = []; rt = []; risky = []; cond = []; TU_irr = []; RU_irr = [];
    for s = 1:length(data)
        latents = kalman_filter(data(s));
        V = [V; latents.m(:,1) - latents.m(:,2)];
        risky = [risky; double((data(s).cond==1&data(s).choice==1)|(data(s).cond==2&data(s).choice==2)) - double((data(s).cond==1&data(s).choice==2)|(data(s).cond==2&data(s).choice==1))];
        RS = [RS; double(data(s).cond==1) - double(data(s).cond==2)];
        SS = [SS; double(data(s).cond==4) - double(data(s).cond==3)];
        cond = [cond; data(s).cond];
        C = [C; double(data(s).choice==1)];
        N = length(data(s).choice);
        S = [S; zeros(N,1)+s];
        TU = [TU; sqrt(latents.s(:,1) + latents.s(:,2))];
        RU = [RU; sqrt(latents.s(:,1)) - sqrt(latents.s(:,2))];
        TU_irr = [TU_irr; sqrt(latents.tau(:,1) + latents.tau(:,2)) ];
        RU_irr = [RU_irr; sqrt(latents.tau(:,1)) - sqrt(latents.tau(:,2))];
        rt = [rt; data(s).RT];
    end
    
    SSV = SS.*V;
    VTU = V./TU;
    VTU_irr = V./TU_irr;
    rt = log(rt);
    cond = categorical(cond);
    
    % orthogonalize RU_irr w.r.t RU
    tmp = GramSchmidt([RU RU_irr]);
    RU_irr = tmp(:,2);
    tmp = GramSchmidt([VTU VTU_irr]);
    VTU_irr = tmp(:,2);
    
    if standardize==1
        VTU = zscore(VTU);
        V = zscore(V);
        RU = zscore(RU);
        VTU_irr = zscore(VTU_irr);
        RU_irr = zscore(RU_irr);
    elseif standardize == 2
        VTU = VTU / norm(VTU);
        V = V / norm(V);
        RU = RU / norm(RU);
        VTU_irr = VTU_irr / norm(VTU_irr);
        RU_irr = RU_irr / norm(RU_irr);
    end
    
    tbl = table(RS,SS,SSV,C,S,RU,VTU,V,TU,RU_irr,VTU_irr,rt,risky,cond);
