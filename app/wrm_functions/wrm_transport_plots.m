function wrm_transport_plots(obj,mobj,option)
%                                                                          **************************
%-------function help------------------------------------------------------Currently under development
% NAME
%   wrm_transport_plots.m
% PURPOSE
%   Use sediment transport results for a set of points along the coast to 
%   examine drift rates, the divergence of drift and the Peclet number
%   (indicates balance of advection and diffusion).
% USAGE
%   wrm_transport_plotsl(objm mobj)
% INPUT
%   obj - WRM_SedimentTransport class instance 
%   mobj - CoastalTools class instance (currently not used)
%   option - selected plot options. Currently the options include:
%           1.	Mean drift: at each point the mean and standard deviation 
%               of the selected variable is computed for all years and the 
%               results of the mean and +/- st.dev. are plotted as a function 
%               of point number (non-dimensional equivalent to distance 
%               along the shore). 
%           2.	Summary point drift: for each point plot the monthly and 
%               annual drift values and the positive and negative contributions 
%               (uses littoraldriftstats from the Derive Output function library).
%           3.	Summary shore drift: creates surface plots of a selected 
%               variable downsampled to months or years by applying a suitable
%               function (mean, stdec, sum, min, max, etc) and plotted as a 
%               function of point number (non-dimensional equivalent to distance
%               along the shore) and time.
%           4.	Divergence and Peclet: 
% OUTPUT 
%   plot options for drift at multiple points as detailed above
% NOTES
%    called as part of WaveRayModel from WRM_SedimentTranport
% SEE ALSO
%   Kahl, et al (2024). Characterizing longshore transport potential and 
%   divergence of drift to inform beach loss trends. Coastal Engineering, 
%   189, 104473. https://doi.org/10.1016/j.coastaleng.2024.104473 
%
% Author: Ian Townend
% CoastalSEA (c)May 2025
%--------------------------------------------------------------------------
%  
    switch option
        case 'Annual Mean Drift'
            an_mean_drift(obj,mobj);
        case 'Monthly Mean Drift'
            mn_mean_drift(obj,mobj);
        case 'Summary Point Drift'
            summary_point_drift(obj);
        case 'Summary Shore Drift'
            summary_shore_drift(obj);
        case 'Monthly Peclet Ratio'
            drift_peclet(obj,mobj);
        case 'Cluster Peclet Ratio'
            cluster_peclet(obj,mobj);
        case 'Wave-Drift Tables'
            wavedrift_table(obj,mobj);
    end
end

%%
function an_mean_drift(obj,mobj)
    %plot the annual mean drift and standard deviation for all years
    calmsthreshold = calmsThreshold(mobj);
    if isnan(calmsthreshold), return; end  %user cancelled

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable();
    if isempty(var), return; end
    meanVar = zeros(1,npnts); upper = meanVar; lower = upper; peclet = upper;
    for i=1:npnts
        Var = dst.(pntnames{i}).(var.name);
        Var(abs(Var)<calmsthreshold) = NaN; %remove near zero values
        meanVar(i) = mean(Var,'omitnan');
        stdVar = std(Var,'omitnan');
        upper(i) = meanVar(i)+stdVar;
        lower(i) = meanVar(i)-stdVar;
        peclet(i) = meanVar(i)./stdVar;
    end
    downcoast = meanVar; upcoast = meanVar;
    downcoast(peclet>-1) = NaN;     %downcoast advection fo Pe<-1
    upcoast(peclet<1) = NaN;        %upcoast advection fo Pe>1

    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);
    loc = 1:npnts;
    grey = mcolor('light grey');
    fill(ax,[loc, fliplr(loc)],[upper,fliplr(lower)],grey,...
        'FaceAlpha',0.5,'EdgeColor',grey,'DisplayName',sprintf('Stdev %s',var.name));
    hold on
    plot(ax,loc,meanVar,'-k','DisplayName',sprintf('Mean %s',var.name));
    if any(~isnan(downcoast)) || any(~isnan(upcoast))
        plot(ax,loc,downcoast,'-og','DisplayName','Pe<-1','LineWidth',0.8,'MarkerSize',4);
        plot(ax,loc,upcoast,'-ob','DisplayName','Pe>1','LineWidth',0.8,'MarkerSize',4);
    end
    hold off
    xlabel('Position along shore')
    ylabel(var.desc)
    title(sprintf('Case: %s',dst.(pntnames{1}).Description));
    subtitle('Mean and Standard deviation for all years')
    legend
end

%%
function mn_mean_drift(obj,mobj)
    %plot the monthly mean drift for each year at a selected point
    calmsthreshold = calmsThreshold(mobj);
    if isnan(calmsthreshold), return; end  %user cancelled

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable();
    mtime = dst.(pntnames{1}).RowNames;
    
    for i=1:npnts
        Var = dst.(pntnames{i}).(var.name);        
        Var(abs(Var)<calmsthreshold) = NaN; %remove near zero values
        [~,binvar,bintime] = binned_variable(Var,mtime,'month','year');
        nint = size(binvar,2);
        nper = size(binvar,1);
        for j=1:nper
            for k=1:nint
                meanVar = mean(binvar{j,k},'omitnan');
                stdVar = std(binvar{j,k},'omitnan');
                peclet = meanVar./stdVar;
                monthlyMean(i,j,k) = meanVar;
                % monthlyStdev(i,j,k) = stdVar;
                monthlyPeclet(i,j,k) = peclet;
            end
        end
    end
    downcoast = monthlyMean; upcoast = monthlyMean;
    downcoast(monthlyPeclet>-1) = NaN;     %downcoast advection for Pe<-1
    upcoast(monthlyPeclet<1) = NaN;        %upcoast advection for Pe>1

    ok = 1;
    while ok>0
        [sel,ok] = listdlg('Name','Plot profile', ...
                                 'PromptString','Select variable', ...
                                 'ListSize',[200,300], ...
                                 'SelectionMode','single', ...
                                 'ListString',pntnames);  
        if ok==0, continue; end 
        hf = figure('Name','SedTrans','Tag','PlotFig');
        ax = axes(hf); %#ok<LAXES>
        pointvar = squeeze(monthlyMean(sel,:,:));
        pointdown = squeeze(downcoast(sel,:,:));
        pointup = squeeze(upcoast(sel,:,:));
        nyear = size(pointvar,1);
        mntime = unique(bintime.intervals);
        antime = string(bintime.periods);
        
        hold on
        plot(ax,mntime,pointvar(1,:),'DisplayName','Year','ButtonDownFcn',@godisplay);
        plot(ax,mntime,pointdown(1,:),'-og','DisplayName','Pe<-1',...
                'LineWidth',0.8,'MarkerSize',4,'ButtonDownFcn',@godisplay);
        plot(ax,mntime,pointup(1,:),'-ob','DisplayName','Pe>1',...
                'LineWidth',0.8,'MarkerSize',4,'ButtonDownFcn',@godisplay);
        for i=2:nyear
            p1 = plot(ax,mntime,pointvar(i,:),'DisplayName',antime(i),'ButtonDownFcn',@godisplay);
            p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            p1 = plot(ax,mntime,pointdown(i,:),'-og','LineWidth',0.8,...
                'MarkerSize',4,'DisplayName',antime(i),'ButtonDownFcn',@godisplay);
            p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            p1 = plot(ax,mntime,pointup(i,:),'-ob','LineWidth',0.8,...
                'MarkerSize',4,'DisplayName',antime(i),'ButtonDownFcn',@godisplay);
            p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        end
        hold off
        xlabel('Month')
        ylabel(var.desc)
        title(sprintf('Drift potential for %s (%s)',...
                       dst.(pntnames{sel}).Description,pntnames{sel}));
        subtitle('Monthly Means for each year of data set')
        legend
    end
end
%%
function summary_point_drift(obj)
   %summary plot of monthly and annual drift at a point (littoraldriftstats)
   dst = obj.Data;
   pntnames = fieldnames(dst);
   [sel,ok] = listdlg('Name','Plot profile', ...
                                 'PromptString','Select variable', ...
                                 'ListSize',[200,300], ...
                                 'SelectionMode','single', ...
                                 'ListString',pntnames);  
   if ok==0, return; end 
   Qs = dst.(pntnames{sel}).Qs;
   mtime = dst.(pntnames{sel}).RowNames;
   littoraldriftstats(Qs,mtime,'month',false);
   hf = gcf;
   ht = findobj(hf.Children,'String','Drift potential');
   ht.String = sprintf('Drift potential for %s (%s)',...
                       dst.(pntnames{sel}).Description,pntnames{sel});
end

%%
function summary_shore_drift(obj)
    %summary plot of selected statistical propoerty and period for all 
    %alongshore points
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable();
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;

    sample = setDownsampleSettings();
    if isempty(sample), return; end

    nfunc = ['nan',sample.func];
    for i=1:npnts
        Var = dst.(pntnames{i}).(var.name);
        if strcmp(sample.func,'sum')
            Var = Var*seconds(mtime(2)-mtime(1));
        end
        if strcmp(sample.type,'+ve')
            Var(Var<0) = NaN;     %mask all negative values
        elseif strcmp(sample.type,'-ve')
            Var(Var>0) = NaN;     %mask all positive values
        end
        [mt,mvar(:,i)] = downsample(Var,mtime,sample.period,nfunc); 
    end

    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);
    [X,Y] = meshgrid(1:npnts,datenum(mt)); %#ok<DATNM>
    contourf(ax,X,Y,mvar);
    % Format the axes to display datetime
    if strcmp(sample.period,'month')
        datetick('y', 'mmm-yy', 'keepticks'); %#ok<DATIC>
    else
        datetick('y', 'yyyy', 'keepticks'); %#ok<DATIC>
    end
    axis tight

    hc = colorbar;
    if strcmp(sample.func,'sum')
        hc.Label.String = 'Volume (m^3)';
    else
        hc.Label.String = 'Volume (m^3/s)';
    end
    xlabel('Point number')
    ylabel(sprintf('Time (%s)',sample.period))
    title(sprintf('Case: %s',dst.(pntnames{1}).Description));
    subtitle(sprintf('Downsampled %s %s using %s(%s)',sample.type,...
                                      var.name,sample.func,sample.period));
end

%%
function drift_peclet(obj,mobj)
    %plots to examine Peclet ratio using monthly/annual sampling(see Kahl, et al, 2024)
    calmsthreshold = calmsThreshold(mobj);
    if isnan(calmsthreshold), return; end  %user cancelled

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable(1);
    if isempty(var), return; end

    mtime = dst.(pntnames{1}).RowNames;
    % monthlyMean = zeros(npnts,nper*nint); 
    % monthlyStdev = monthlyMean; monthlyPeclet = monthlyMean;
    % annualMean = zerod(npnts,nper); 
    % annualStdev = annualMean; annualPeclet = annualMean;
    for i=1:npnts
        Var = dst.(pntnames{i}).(var.name);        
        Var(abs(Var)<calmsthreshold) = NaN; %remove near zero values
        [~,binvar,bintime] = binned_variable(Var,mtime,'month','year');
        nint = size(binvar,2);
        nper = size(binvar,1);
        nyr = 0;
        for j=1:nper
            for k=1:nint
                meanVar = mean(binvar{j,k},'omitnan');
                stdVar = std(binvar{j,k},'omitnan');
                peclet = meanVar./stdVar;
                if peclet>-1 && peclet<1 %#ok<BDSCI>
                    peclet = 0;
                end
                monthlyMean(i,nyr+k) = meanVar;
                monthlyStdev(i,nyr+k) = stdVar;
                monthlyPeclet(i,nyr+k) = peclet;
            end
            annualMean(i,j) = mean(monthlyMean(i,(nyr+1:nyr+12)));
            annualStdev(i,j) = mean(monthlyStdev(i,(nyr+1:nyr+12)))*sqrt(12);
            peclet = annualMean(i,j)/annualStdev(i,j);
            if isnan(peclet) || isinf(peclet)
                peclet = 0;
            % elseif peclet>-1 && peclet<1
            %     peclet = 0;
            end
            annualPeclet(i,j) = peclet;
            nyr = nyr+nint;
        end
    end
    %plot monthly mean peclet ratio as a surface plot (position,time)
    desc = struct('case',dst.(pntnames{1}).Description,'var','Peclet ratio (<-1 or >1)');
    startdate = datetime(year(bintime.periods(1)),1,1);  %force full year
    enddate = datetime(year(bintime.periods(end)),12,1); %to match variable
    bins = startdate:calmonths(1):enddate;
    axm = plotPeclet(monthlyPeclet,bins,desc);  
    axm.CLim = [-2,2];
    subtitle(axm,'Monthly: Downdrift advection (Pe>1); Updrift advection (Pe<-1)')
    %annual mean peclet ratio as a surface plot (position,time)
    axa = plotPeclet(annualPeclet,bintime.periods,desc);
    subtitle(axa,'Annual: Downdrift advection (Pe>1); Updrift advection (Pe<-1)') 
end

%%
function cluster_peclet(obj,mobj)
     %plots to Peclet ratio using cluster sampling (see Kahl, et al, 2024)
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable(1);  %selects Qs without prompting user
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;

    %set a minimum threshold to remove zero values
    calmsthreshold = calmsThreshold(mobj);
    if isempty(calmsthreshold), return; end  %user cancelled

    ans0 = questdlg('Use absolute values of drift or +/- values?','Clusters',...
                                            'abs(Qs)','+/-(Qs)','abs(Qs)');

    %NB: options defines the variables used in get clusters and includes 
    %additional variables used in mergeSelection for posnegClusters
    options = setClusterOptions(dst.(pntnames{1}).(var.name));
    if strcmp(ans0,'abs(Qs)')        
        [cluster,options] = absClusters(options,dst,calmsthreshold);
    else
        [cluster,options] = posnegClusters(options,dst,calmsthreshold);        
    end

    clustpoint = []; clustpec = []; clustints = [];
    for i=1:npnts
        nint = length([cluster.Ints{i,:}]);
        clustpoint = [clustpoint;repmat(i,nint,1)];
        clustints = [clustints;[cluster.Ints{i,:}]'];
        clustpec = [clustpec;[cluster.Peclet{i,:}]']; 
    end
  
    %plot cluster mean peclet ratio as a surface plot (position,time)
    desc = struct('case',dst.(pntnames{1}).Description,'var','Peclet ratio (<-0.8 or >0.8)');
    bintime = mtime(1):caldays(1):mtime(end);
    x = 1:npnts;
    % Define a grid for interpolation
    [xq, yq] = meshgrid(x, datenum(bintime));    
    % Interpolate scattered data onto the grid
    zq = griddata(clustpoint,datenum(clustints),clustpec, xq, yq, 'linear'); % 'linear', 'nearest', or 'cubic'
    axc = plotPeclet(zq',bintime,desc);    
    axc.CLim = [-2,2];
    % save('clusterplot',"xq","yq","zq","bintime","desc");
    txt1 = 'Clusters: Downdrift advection (Pe>1); Updrift advection (Pe<-1)';
    if strcmp(ans0,'abs(Qs)') 
        txt2 = sprintf('Absolute - Threshold: %0.4f; Interval: %dd; Minimum duration: %dd',...
                    options.threshold,options.clint,options.mincluster);
    else
        txt2 = sprintf('Pos/Neg - Threshold: %0.4f/%0.4f; Interval: %dd/%dd; Minimum duration: %dd/%dd',...
                    options.pos.threshold,options.neg.threshold,...
                    options.pos.clint, options.neg.clint,...
                    options.pos.mincluster, options.neg.mincluster);
    end
    subtitle(axc,sprintf('%s\n%s',txt1,txt2))   
end
%%
function wavedrift_table(obj,mobj)
    %
    calmsthreshold = calmsThreshold(mobj);
    if isnan(calmsthreshold), return; end  %user cancelled

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable(1);
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;




end

%%
function calmsthreshold = calmsThreshold(mobj)
    %set the calms threshold to apply to the data
    calmsthreshold = 100;  %"calms" are drift rates less than threshold
                           % 100m^3/yr ~= 3e-6 m^3/s; 
    answer = inputdlg('Calms threshold (Qs (m^3/yr):','Drift stats',...
                                        1,{num2str(calmsthreshold)});
    if isempty(answer), calmsthreshold = []; return; end
    calmsthreshold = str2double(answer{1})/mobj.Constants.y2s;
end

%%
function var = getVariable(sel)
    %select a variable to use in the plot
    varname = {'Qs','dQdx','Qx'};
    vardesc = {'Alongshore drift rate potential (m^3/s)',...
               'Alongshore drift gradient (m^3/s/m)'...
               'Cross-shore transport rate (m^3/s)'};
    if nargin<1
        [sel,ok] = listdlg('Name','Plot profile', ...
                                     'PromptString','Select variable', ...
                                     'ListSize',[200,80], ...
                                     'SelectionMode','single', ...
                                     'ListString',vardesc);    
        if ok==0, var = []; return; end 
    end
    var.name = varname{sel};
    var.desc = vardesc{sel};
end

%%
function sample = setDownsampleSettings()
    %check that the current sit parameters settings are correct
    %modifications used to update RunParams stored with Case
    type = 'all';
    period = 'year';
    func = 'sum';
    promptxt = {'Variable sampling (all, +ve, -ve)',...
                'Sub-sampling interval (year or month)',...
                'Function to apply (mean, min, max, etc)'};
    defaults = {type,period,func};
    data = inputdlg(promptxt,'Drift settings',1,defaults);
    if isempty(data), sample = []; return; end  %no change to default settings
    sample.type = data{1};
    sample.period = data{2};
    sample.func = data{3}; 
end

%%
function [cluster,userops] = absClusters(options,dst,calmsthreshold)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable(1);  %selects Qs without prompting user
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;

    ans2 = questdlg('Check settings for selected points?','Clusters','Yes','No','Quit','Yes');
    if strcmp(ans2,'Yes')
        promptxt = {'Accept figures are used to adjust the cluster selection';...
                    'The plot shows absolute drift values (all +ve)';...
                    'Use Threshold and Time between clusters to adjust number of clusters'};
        getdialog(promptxt,[],4);
    
        ok = 1;
        while ok>0
               [sel,ok] = listdlg('Name','Plot profile', ...
                                     'PromptString','Select variable', ...
                                     'ListSize',[200,300], ...
                                     'SelectionMode','single', ...
                                     'ListString',pntnames);  
                if ok==0, continue; end 
    
                Var = dst.(pntnames{sel}).(var.name); 
                vardst = getDSTable(dst.(pntnames{sel}),'VariableNames',var.name);
                vardst.(var.name) = abs(vardst.(var.name));
                [idcls,userops] = getclusters(vardst,options);
                userops.mincluster = options.mincluster;
                userops.isplot = options.isplot;
                mergeAbsClusters(Var,mtime,idcls,userops);
    
                ans1 = questdlg('Use selected options or examine another point?',...
                         'Clusters','Use selected','New point','Use selected');
                if strcmp(ans1,'Use selected'), ok = 0; end
        end
    
        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');        
    elseif  strcmp(ans2,'No')
        userops = options;
    end

    if strcmp(ans2,'Quit'), return; end

    userops.isplot = false; %supress plots in for loop
    for i=1:npnts        
        Var = dst.(pntnames{i}).(var.name);        
        Var(abs(Var)<calmsthreshold) = NaN; %remove near zero values
        vardst = getDSTable(dst.(pntnames{i}),'VariableNames',var.name);
        vardst.(var.name) = abs(vardst.(var.name)); %use absolute values for intervals

        % find clusters based on results from peak selection
        idcls = getVarClusters(vardst,userops);        
        %merge any overlaps to define intervals to be used
        medges = mergeAbsClusters(Var,mtime,idcls,userops);
        if isempty(medges)
            medges = [mtime(1),mtime(end)];
            % warndlg(sprintf('No clusters found for point %d\nTry changing the threshold',i));
            % return;
        else
            medges = [mtime(1),medges,mtime(end)];
        end
        
        %find the indices of the variable within each interval
        [intervals,intstart] = discretize(mtime,medges);
        nint = length(intstart);
        %compute the statistics for the point over each interval 
        for j=1:nint
            idint = intervals==j;
            binvar = Var(idint);
            meanVar = mean(binvar,'omitnan');
            stdVar = std(binvar,'omitnan');
            peclet = meanVar./stdVar;

            %set diffusion values to 0
            if isnan(peclet) || isinf(peclet)
                peclet = 0;
            elseif peclet>-0.8 && peclet<0.8 %#ok<BDSCI>
                peclet = 0;
            end
            % %limit the maximum advection values
            % if peclet<-2
            %     peclet = -2;
            % elseif peclet>2
            %     peclet = 2;
            % end

            cluster.Ints{i,j} = intstart(j); %#ok<*AGROW>
            cluster.Mean{i,j} = meanVar;
            cluster.Stdev{i,j} = stdVar;
            cluster.Peclet{i,j} = peclet;
        end
    end
end

%%
function [cluster,userops]  = posnegClusters(options,dst,calmsthreshold)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    var = getVariable(1);  %selects Qs without prompting user
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;

    %set a minimum trheold to remove zero values
    ans2 = questdlg('Check settings for selected points?','Clusters','Yes','No','Quit','Yes');
    if strcmp(ans2,'Yes')
        promptxt = {'Accept figures are used to adjust the threshold selection';...
                    'The first plot sets the positive threshold';...
                    'The second plot sets the negative threshold';...
                    '(NB: variable is inverted in the 2nd plot)'};
        getdialog(promptxt,[],5);

        ok = 1;
        while ok>0
               [sel,ok] = listdlg('Name','Plot profile', ...
                                     'PromptString','Select variable', ...
                                     'ListSize',[200,300], ...
                                     'SelectionMode','single', ...
                                     'ListString',pntnames);  
                if ok==0, continue; end 

                Var = dst.(pntnames{sel}).(var.name); 
                posdst = getDSTable(dst.(pntnames{sel}),'VariableNames',var.name);
                [idpos,posops] = getclusters(posdst,options);
                %posops.mincluster = options.mincluster;
                %posops.isplot = options.isplot;
                negdst = posdst;
                negdst.(var.name) = negdst.(var.name)*-1;
                [idneg,negops] = getclusters(negdst,options);
                %negops.mincluster = options.mincluster;
                %negops.isplot = options.isplot;
                mergePosNegClusters(Var,mtime,idpos,idneg,options);

                ans1 = questdlg('Use selected options or examine another point?',...
                         'Clusters','Use selected','New point','Use selected');
                if strcmp(ans1,'Use selected'), ok = 0; end
        end

        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');
    elseif strcmp(ans2,'No')
        posops = options; negops = options; %use same options for postive and nrgstive drift
    end

    if strcmp(ans2,'Quit'), return; end

    options.isplot = false;
    for i=1:npnts        
        Var = dst.(pntnames{i}).(var.name);        
        Var(abs(Var)<calmsthreshold) = NaN; %remove near zero values
        posdst = getDSTable(dst.(pntnames{i}),'VariableNames',var.name);

        % find clusters based on results from peak selection
        idposcls = getVarClusters(posdst,posops);
        negdst = posdst;
        negdst.(var.name) = negdst.(var.name)*-1;
        idnegcls = getVarClusters(negdst,negops);

        %merge any overlaps to define intervals to be used
        medges = mergePosNegClusters(Var,mtime,idposcls,idnegcls,options); %only uses mincluster field in options
        if isempty(medges)
            warndlg(sprintf('No clusters found for point %d\nTry changing the threshold',i));
            return;
        end
        medges = [mtime(1),medges,mtime(end)];
        %find the indices of the variable within each interval
        [intervals,intstart] = discretize(mtime,medges);
        nint = length(intstart);
        %compute the statistics for the point over each interval 
        for j=1:nint
            idint = intervals==j;
            binvar = Var(idint);
            meanVar = mean(binvar,'omitnan');
            stdVar = std(binvar,'omitnan');
            peclet = meanVar./stdVar;

            %set diffusion values to 0
            if isnan(peclet) || isinf(peclet)
                peclet = 0;
            elseif peclet>-0.8 && peclet<0.8 %#ok<BDSCI>
                peclet = 0;
            end
            % %limit the maximum advection values
            % if peclet<-2
            %     peclet = -2;
            % elseif peclet>2
            %     peclet = 2;
            % end

            cluster.Ints{i,j} = intstart(j); %#ok<*AGROW>
            cluster.Mean{i,j} = meanVar;
            cluster.Stdev{i,j} = stdVar;
            cluster.Peclet{i,j} = peclet;
        end
    end
    userops.pos = posops;
    userops.neg = negops;
end

%%
function options = setClusterOptions(data,opts)
    %define the options used in a peaks and cluster data selection
    if nargin<2 || isempty(opts)
        default = {num2str(mean(data,'omitnan')+std(data,'omitnan')),...
                   '1','0','15','5'};
    else
        default{1} = num2str(opts.threshold);
        default{2} = num2str(opts.method);
        default{3} = num2str(opts.tint);
        default{4} = num2str(opts.clint);
        default{5} = num2str(opts.mincluster); 
    end
    prompt = {'Threshold for peaks:','Selection method (1-4)',...
        'Time between peaks (hours)','Time between clusters (days)',...
        'Minimum duration of a cluster (days)'};
    title = 'Cluster Statistics';
    numlines = 1;
    
    answer = inputdlg(prompt,title,numlines,default);
    if isempty(answer), options = []; return; end
    threshold = str2double(answer{1});   %variable threshold 
    method = str2double(answer{2});      %peak selection method (see peaks.m)
    tint = str2double(answer{3});        %time interval between independent peaks (h)
    clint = str2double(answer{4});       %time interval for clusters (d) 
    mincluster = str2double(answer{5});  %minimum length of a cluster (d)

    options = struct('threshold',threshold,'method',method,'tint',tint,...
                     'clint',clint,'mincluster',mincluster,'isplot',true);
end

% %%
% function options = setPosNegClusterOptions(data,opts)
%     %define the options used in a peaks and cluster data selection
%     if nargin<2 || isempty(opts)
%         default = {num2str(mean(data,'omitnan')+std(data,'omitnan')),...
%                    '1','0','15','5','10'};
%     else
%         default{1} = num2str(opts.threshold);
%         default{2} = num2str(opts.method);
%         default{3} = num2str(opts.tint);
%         default{4} = num2str(opts.clint);
%         default{5} = num2str(opts.mincluster); 
%         default{6} = num2str(opts.splitdur);
%     end
%     prompt = {'Threshold for peaks:','Selection method (1-4)',...
%         'Time between peaks (hours)','Time between clusters (days)',...
%         'Minimum duration of a cluster (days)','Duration of overlaps to split (days)'};
%     title = 'Cluster Statistics';
%     numlines = 1;
% 
%     answer = inputdlg(prompt,title,numlines,default);
%     if isempty(answer), options = []; return; end
%     threshold = str2double(answer{1});   %variable threshold 
%     method = str2double(answer{2});      %peak selection method (see peaks.m)
%     tint = str2double(answer{3});        %time interval between independent peaks (h)
%     clint = str2double(answer{4});       %time interval for clusters (d) 
%     mincluster = str2double(answer{5});  %minimum length of a cluster (d)
%     splitdur = str2double(answer{6});    %maximum duration to split (d)
% 
%     options = struct('threshold',threshold,'method',method,'tint',tint,...
%                      'clint',clint,'mincluster',mincluster,...
%                      'splitdur',splitdur,'isplot',true);
% end

%%
function medges = mergeAbsClusters(var,mtime,idpos,opts)
    %merge  absolute cluster selections to a single set of edges
    mincluster = opts.mincluster*24;      %min length of a cluster (h)
    dt = mode(diff(mtime));
    mincls = floor(mincluster/hours(dt));  
    func = @(x) length(x)<mincls;
    postimes = {idpos(:).date};
    posshort = cellfun(func,postimes,"UniformOutput",false);
    postimes([posshort{:}]) = [];

    mdates.posstart = cellfun(@(x) x(1),postimes);
    mdates.posend = cellfun(@(x) x(end),postimes);
    medges = sort(unique([mdates.posstart,mdates.posend]));

    if opts.isplot
        plotEdges(var,mtime,medges,'Intervals to used for statistics');
    end
end

%%
function medges = mergePosNegClusters(var,mtime,idpos,idneg,opts)
    %merge positive and negative cluster selections to a single set of edges
    mincluster = opts.mincluster*24;      %min length of a cluster (h)
    dt = mode(diff(mtime));
    mincls = floor(mincluster/hours(dt));  
    func = @(x) length(x)<mincls;
    postimes = {idpos(:).date};
    posshort = cellfun(func,postimes,"UniformOutput",false);
    postimes([posshort{:}]) = [];
    negtimes = {idneg(:).date};
    negshort = cellfun(func,negtimes,"UniformOutput",false);
    negtimes([negshort{:}]) = [];

    mdates.posstart = cellfun(@(x) x(1),postimes);
    mdates.posend = cellfun(@(x) x(end),postimes);
    mdates.negstart = cellfun(@(x) x(1),negtimes);
    mdates.negend = cellfun(@(x) x(end),negtimes);

    %find any overlaps and merge any that are short
    mdates = mergeOverlaps(var,mtime,mdates,opts);    
    if isempty(mdates.posstart) && isempty(mdates.negstart)
        medges = []; return;
    elseif isempty(mdates.posstart) 
        medges = sort(unique([mdates.negstart,mdates.negend]));
    elseif isempty(mdates.negstart)
        medges = sort(unique([mdates.posstart,mdates.posend]));
    else
        medges = sort(unique([mdates.posstart,mdates.posend,mdates.negstart,mdates.negend]));
    end
    medges = sort(unique([mdates.posstart,mdates.posend]));

    if opts.isplot
        plotEdges(var,mtime,medges,'Intervals to used for statistics');
    end
end

%%
function mdates = mergeOverlaps(var,mtime,mdates,opts)
    %find any overlaps and merge any that are short
    posstart = mdates.posstart;
    posend = mdates.posend;
    negstart = mdates.negstart;
    negend = mdates.negend;

    overlaps = findOverlaps(posstart,posend,negstart,negend);
    if any(overlaps,'all')
        
        [row, col] = find(overlaps);
        %plotMergedVar(var,mtime,posstart(row),posend(row),negstart(col),negend(col),'Unmerged cluster overlaps');
        for k = 1:length(row)
            overlap_start = max(posstart(row(k)), negstart(col(k)));
            overlap_end = min(posend(row(k)), negend(col(k)));
            overlap_length = hours(overlap_end-overlap_start);
            if overlap_length<opts.mincluster*24    %split between the two
                if posstart(row(k))<negstart(col(k)) && posend(row(k))<negend(col(k))
                    posend(row(k)) = posend(row(k))-hours(overlap_length/2+0.0);
                    negstart(col(k)) = negstart(col(k))+hours(overlap_length/2+0.0);
                elseif posstart(row(k))>negstart(col(k)) && posend(row(k))>negend(col(k))
                    posstart(row(k)) = posstart(row(k))+hours(overlap_length/2+0.0);
                    negend(col(k)) = negend(col(k))-hours(overlap_length/2+0.0);
                end
            else
                fprintf('Overlap %d: %s to %s\n', k, overlap_start, overlap_end);
            end            
        end
        overlaps = findOverlaps(posstart,posend,negstart,negend);
        [row, col] = find(overlaps);
        if ~isempty(row)
            if opts.isplot
                plotMergedVar(var,mtime,posstart(row),posend(row),negstart(col),negend(col),'Cluster overlaps to be subdivided');
            end
            fprintf('%d overlaps have been subdivided into discrete intervals\n', length(row));
        end
        %update struct with merged intervals
        mdates = struct('posstart',posstart,'posend',posend,'negstart',negstart,'negend',negend);
    end

    %-nested function------------------------------------------------------
    function overlaps = findOverlaps(posstart,posend,negstart,negend)
        % Initialize a logical matrix to store overlaps
        overlaps = false(length(posstart), length(negstart));
        
        % Check for overlaps between intervals
        for i = 1:length(posstart)
            for j = 1:length(negstart)
                overlaps(i, j) = (posstart(i) <= negend(j)) && (negstart(j) <= posend(i));
            end
        end
    end
end

%%
function plotMergedVar(var,mtime,pstart,pend,nstart,nend,titxt)
    %plot the merged selection
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf); 
    plot(ax,mtime,var,'Color',[0.75,0.75,0.75],'LineWidth',0.2)
    yy = ylim;
    posy = [0,yy(2)];
    negy = [0,yy(1)];
    hold on
    if length(pstart)==2 
         %needed if there are only 2 points to avoid plotting diagonal
        for i=1:2
            plot([pstart(i), pstart(i)],posy,'-','Color',"#77AC30",'LineWidth',0.8);%#7E2F8E
            plot([pend(i),pend(i)],posy,'--','Color',"#77AC30",'LineWidth',0.8);
            plot([nstart(i),nstart(i)],negy,'-','Color','#A2142F','LineWidth',0.8);
            plot([nend(i),nend(i)],negy,'--','Color','#A2142F','LineWidth',0.8); 
        end
    else
        plot([pstart', pstart'],posy,'-','Color',"#77AC30",'LineWidth',0.8);%#7E2F8E
        plot([pend',pend'],posy,'--','Color',"#77AC30",'LineWidth',0.8);
        plot([nstart',nstart'],negy,'-','Color','#A2142F','LineWidth',0.8);
        plot([nend',nend'],negy,'--','Color','#A2142F','LineWidth',0.8);  
    end
    hold off
    xlabel('Time')
    ylabel('Selected drift variable')
    title(titxt)
end

%%
function plotEdges(var,mtime,medges,titxt)
    %plot the merged edges to be used to compute the statitics
    mvar = max(abs(var),[],'omitnan')/4;
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf); 
    plot(ax,mtime,var,'Color',[0.75,0.75,0.75],'LineWidth',0.2)
    hold on 
        plot([medges',medges'],[-mvar,mvar],'-','Color',"#0072BD")
        plot(medges,0,'.','Color',"#0072BD",'MarkerSize',4)
    hold off
    xlabel('Time')
    ylabel('Selected drift variable')
    title(titxt)
end

%%
function idcls = getVarClusters(dst,opts)
    %extract the clusters for a variable in the dstable
    var = dst.(dst.VariableNames{1});
    mtime = dst.RowNames;
    % find peaks (method 1:all peaks; 2:independent crossings; 3:timing
    % seperation of tint)
    returnflag = 0; %0:returns indices of peaks; 1:returns values       
    idpks = peaksoverthreshold(var,opts.threshold,opts.method,...
                                        mtime,hours(opts.tint),returnflag);
    % find clusters based on results from peak selection
    pk_date = mtime(idpks);    %datetime of peak
    pk_vals = var(idpks);  %value of peak
    idcls = clusters(pk_date,pk_vals,days(opts.clint));
end

%%
function ax = plotPeclet(var,bintime,desc)
    %monthly or cluster mean peclet ratio as a surface plot (position,time)
    npnts = size(var,1);
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);    
    grid on
    [X,Y] = meshgrid(1:npnts,bintime); 
    %contourf(ax,X,datenum(Y),var');
    surf(ax,X,Y,var');
    shading interp
    view(2)
    axis tight
    ax.Layer = 'top';
    hc = colorbar;
    hc.Label.String = desc.var;
    datetick('y', 'yyyy'); %#ok<DATIC>
    xlabel('Position along shore')
    ylabel('Year')
    title(sprintf('Case: %s',desc.case));  
end

%%
