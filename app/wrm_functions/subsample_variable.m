function [binvar,bintime] = subsample_variable(dst,varsel,binsel)    
%
%-------function help------------------------------------------------------
% NAME
%   subsample_variable.m
% PURPOSE
%   
% USAGE
%   vari = subsample_variable(obj) 
% INPUTS
%   dst - 
%   varsel - struct containing variable selection, isdrift flag, calms
%            threshold and peclet threshold
%   binsel - [1x4] cell array defining user selection for:
%            1 - bin interval; 2 - bin period; 3 - Season; 4 - Direction; 
%            Bins - {'All','Year', 'Quarter', 'Month', 'Week', 'Day', 'Hour'};
%            Season - {'All','Winters','Summers'};
%            Direction - {'All','+ve only','-ve only'};
% OUTPUT
%   vari - 
% SEE ALSO
%   binned_variable.m and get_var_sampling.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2026
%--------------------------------------------------------------------------
%
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    mtime = dst.(pntnames{1}).RowNames;
    varAllPnts = zeros(npnts,numel(mtime));
    for i=1:npnts
        varAllPnts(i,:) = dst.(pntnames{i}).(varsel.name);
        varAllPnts(abs(varAllPnts)<varsel.calms.value) = NaN; %remove near zero values
    end

    if binsel{1}==1         %all years case no subsampling required
        binvar = mat2cell(varAllPnts,ones(1,npnts));
        bintime.intstart = mtime;
        %apply direction subsampling if selected
        if binsel{4}>1
            [binvar,bintime] = getDirectionSample(binvar,bintime,binsel);
        end
        return;
    end

    bins = {'All','year', 'quarter', 'month', 'week', 'day', 'hour'};
    bininterval = bins{binsel{1}};
    binperiod = bins{binsel{2}};
    %loop over each point and subsample for seas and direction
    %binvar = cell(npnts,x,y);
    for i=1:npnts
        vari = varAllPnts(i,:);
        %apply seasonal subsampling if selected
        if binsel{3}>1                  %summer or winter selected
            [vari,bintime] = getAllSummerWinter(vari,mtime,binsel);
        else
            [~,vari,bintime] = binned_variable(vari,mtime,bininterval,binperiod);
        end
        
        %apply direction subsampling if selected
        if binsel{4}>1
            [vari,bintime] = getDirectionSample(vari,bintime,binsel);
        end
        binvar(i,:,:) = vari;
    end
end

%%
function [binvar,bintime] = getAllSummerWinter(vari,mtime,binsel)
    %subsample the data for winters or summers (All years or for each year)
    [~,binvar,bintime] = binned_variable(vari,mtime,'month','year');
    % Find indices of months in the desired range
    idsel = find(bintime.intervals >= 4 & bintime.intervals <= 9);
    if binsel{3}==2                           %'Winters'
        binvar = binvar(:,[1:3,10:12]);
        bintime.intervals(idsel) = [];
        bintime.intstart(idsel) = [];
    else                                      %'Summers'
        binvar = binvar(:,4:9);
        bintime.intervals = bintime.intervals(idsel);
        bintime.intstart = bintime.intstart(idsel);
    end
    %
    if binsel{1}==1
        binvar = vertcat(binvar{:});          %concatenate if all years
    end
end

%%
function [binvar,bintime] = getDirectionSample(binvar,bintime,binsel)
    %subsample the data for all positive or all negative values
    if binsel{4}==2           %'+ve only
        binvar = cellfun(@(x) x .* (x >= 0), binvar, 'UniformOutput', false);
    elseif binsel{4}==3       %'-ve only
        binvar = cellfun(@(x) x .* (x < 0), binvar, 'UniformOutput', false);
    end
end