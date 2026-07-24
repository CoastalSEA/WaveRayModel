function stats = wrm_interval_statistics(binvar)
%
%-------function help------------------------------------------------------
% NAME
%   wrm_interval_statistics.m
% PURPOSE
%   Compute Peclet ratio, directional bias, skewness, and classification
%   for binvar
% USAGE
%   stats = wrm_interval_statistics(binvar);
% INPUT
%   binvar - numeric data to used to compute statistics
% OUTPUT
%   stats - struct with the following fields:
%           mean - Mean
%           stdev - Standard deviation
%           ptle95 - 95th Percentile
%           ptle05 - 5th Percentile
%           Pe - Peclet ratio
%           dirBias - Directional bias index
%           netQ - Net transport index
%           skewness - Skewness
%           classification - Classification
%           index - Classification index (+/- 0-6)
%           nPoints - No. of points
% NOTES
%   Classification index definitions: 
%           0 - Mixed-transitional
%           1 - Intermittent bursts
%           2 - Noisy bidirectional
%           3 - Competing processes
%           4 - Intermittent-diffusive
%           5 - Ongoing advection
%           6 - Persistent advection
% SEE ALSO
%   called in wrm_transport_plots.m
%
% Author: Ian Townend & Copilot
% CoastalSEA (c)July 2026
%----------------------------------------------------------------------
%
    stats = struct('mean',[],'stdev',[],'ptle95',[],'ptle05',[],...
                   'Pe',[],'dirBias',[],'netQ',[],'skewness',[],...
                   'classification','','index',[],'nPoints',[]);
    if isempty(binvar)
        stats = struct('mean',NaN,'stdev',NaN,'ptle95',NaN,'ptle05',NaN,...
                       'Pe',NaN,'dirBias',NaN,'netQ',NaN,'skewness',NaN,...
                       'classification','None','index',NaN,'nPoints',0);
        return
    end

    % --- Peclet-like ratio ---
    mu  = mean(binvar,'omitnan');
    sig = std(binvar,'omitnan');
    Pe  = mu / sig;
    
    % --- Directional bias ---
    fpos = mean(binvar > 0);
    fneg = mean(binvar < 0);
    dirBias = (fpos - fneg) / (fpos + fneg);
    
    % --- Skewness ---
    sk = skewness(binvar); %needs statistics and machine learning toolbox
    
    % --- Net transport index netQ ---
    Mpos = sum(binvar(binvar > 0));
    Mneg = sum(abs(binvar(binvar < 0)));
    netQ = (Mpos - Mneg) / (Mpos + Mneg);

    % ---Percentile stats ------

    % --- Classification ---
        if Pe>1 && dirBias>0.8
            classification = 'Persistent advection';
            index = 6*sign(Pe);
        elseif Pe<1 && dirBias>0.8
            classification = 'Ongoing advection';
            index = 5*sign(Pe);
        elseif Pe<1 && abs(dirBias)<0.5 && abs(sk)<0.5
            classification = 'Intermittent-diffusive'; 
            index = 4*sign(Pe);
        elseif abs(Pe - 1)<0.2 && dirBias>=0.5 && dirBias<=0.8 && abs(sk)>= 0.5 && abs(sk)<= 1.5
            classification = 'Competing processes';
            index = 3*sign(Pe);
        elseif Pe<1 &&dirBias>=0.3 && dirBias<=0.6 && abs(sk)>=0.5 && abs(sk)<=1.5
            classification = 'Noisy bidirectional';
            index = 2*sign(Pe);
        elseif Pe>1 && dirBias<0.5 && abs(sk)>1.5
            classification = 'Intermittent bursts';
            index = 1*sign(Pe);
        else
            classification = 'Mixed-transitional';
            index = 0*sign(Pe);
        end

    % --- Store ---
    stats.mean = mu;
    stats.stdev = sig;
    stats.ptle95 = prctile(binvar,95);
    stats.ptle05 = prctile(binvar,05);
    stats.Pe = Pe;
    stats.dirBias = dirBias;
    stats.netQ = netQ;
    stats.skewness = sk;
    stats.classification = classification;
    stats.index = index;
    stats.nPoints = numel(binvar);
end