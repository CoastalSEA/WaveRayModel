function wrm_transport_plots(obj,mobj,option)
%                                                                          **************************
%-------function help------------------------------------------------------Currently under development
% NAME
%   wrm_transport_plots.m
% PURPOSE
%   Use sediment transport results for a set of points along the coast to 
%   examine drift rates, the divergence of drift and the Peclet number
%   (indicates balance of advection and diffusion). 
%   The mean and summary plot options can also be used to plot multi-point 
%   wave data and wave roses.
% USAGE
%   wrm_transport_plotsl(obj,mobj,option)
% INPUT
%   obj - WRM_SedimentTransport or WRM_WaveModel class instance
%   mobj - WRM modelUI class instance (only used for constants)
%   option - selected plot options. Currently the options include:
%            'Annual Mean Drift' - at each point the mean and standard deviation 
%               of the selected variable is computed for all years, summer  
%               or winter and the results of the mean and +/- st.dev. 
%               are plotted as a function of point number (non-dimensional 
%               equivalent to distance along the shore).
%            'Monthly Mean Drift' - for a selected point plot the monthly 
%               means in each year of the dataset
%            'Summary Point Drift' - for each point plot the monthly and 
%               annual drift values and the positive and negative contributions 
%               (uses littoraldriftstats from the Derive Output function library).
%            'Summary Shore Drift' - creates surface plots of a selected 
%               variable downsampled to months or years by applying a suitable
%               function (mean, stdec, sum, min, max, etc) and plotted as a 
%               function of point number (non-dimensional equivalent to distance
%               along the shore) and time.
%            'Monthly Peclet Ratio' - plots to examine Peclet ratio along-shore
%               and over time using monthly/annual sampling (see Kahl, et al, 2024)
%            'Cluster Peclet Ratio'- plots to examine Peclet ratio along-shore
%               and over time using cluster sampling (see Kahl, et al, 2024)
%            'Wave-Drift Tables' - not yet implemented                      
% OUTPUT 
%   plot options for drift at multiple points as detailed above
% NOTES
%   called as part of WaveRayModel from WRM_SedimentTranport
%   also referred to in menu options as 'Multi-point Plots'
% SEE ALSO
%   Kahl, et al (2024). Characterizing longshore transport potential and 
%   divergence of drift to inform beach loss trends. Coastal Engineering, 
%   189, 104473. https://doi.org/10.1016/j.coastaleng.2024.104473 
%
% Author: Ian Townend
% CoastalSEA (c)May 2025
%--------------------------------------------------------------------------
%  
    msgtxt = 'WRM_SedimentTransport class required for this option';
    switch option
        case 'Annual Mean Drift'
            ann_mean_drift(obj);
        case 'Binned Mean Drift'
            bin_mean_drift(obj);
        case 'Summary Point Drift'
            summary_point_drift(obj,msgtxt);
        case 'Summary Shore Drift'
            summary_shore_drift(obj);
        case 'Binned Peclet Ratio'
            drift_peclet(obj,msgtxt);
        case 'Cluster Peclet Ratio'
            cluster_peclet(obj,msgtxt);
        case 'Wave Rose Plots'
            multi_rose_plots(obj);
    end
end

%%
function ann_mean_drift(obj)
    %plot the annual mean drift and standard deviation for all years
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end

    %select summmer/winter and +ve/-ve drift
    selection = get_var_sampling([1,1,1,1],false);  %no stats selection
    if isempty(selection), return; end              
    selection(1:2) = {1,1};  %force selection of all years
    [binvar,~] = subsample_variable(dst,varsel,selection);

    meanVar = zeros(1,npnts); upper = meanVar; lower = upper; 
    peclet = upper; uppct = upper; lowpct = upper;
    for i=1:npnts
        Var = binvar{i};
        meanVar(i) = mean(Var,'omitnan');
        stdVar = std(Var,'omitnan');
        upper(i) = meanVar(i)+stdVar;
        lower(i) = meanVar(i)-stdVar;
        peclet(i) = meanVar(i)./stdVar;
        uppct(i) = prctile(Var,95);
        lowpct(i) = prctile(Var,5);
    end
    downcoast = meanVar; upcoast = meanVar;
    downcoast(peclet>-1) = NaN;     %downcoast advection fo Pe<-1
    upcoast(peclet<1) = NaN;        %upcoast advection fo Pe>1

    if selection{4}==2
        lower = zeros(size(lower));
    elseif selection{4}==3
        upper = zeros(size(upper));
    end

    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);
    loc = 1:npnts;
    grey = mcolor('light grey');
    fill(ax,[loc, fliplr(loc)],[upper,fliplr(lower)],grey,...
        'FaceAlpha',0.5,'EdgeColor',grey,'DisplayName',sprintf('Stdev %s',varsel.name));
    hold on
    plot(ax,loc,meanVar,'-k','DisplayName',sprintf('Mean %s',varsel.name));
    if varsel.isdrift && (any(~isnan(downcoast)) || any(~isnan(upcoast)))
        plot(ax,loc,downcoast,'-og','DisplayName','Pe<-1','LineWidth',0.8,'MarkerSize',4);
        plot(ax,loc,upcoast,'-ob','DisplayName','Pe>1','LineWidth',0.8,'MarkerSize',4);
    end
    plot(ax,loc,uppct,'-.b','DisplayName',sprintf('95 percentile %s',varsel.name));
    plot(ax,loc,lowpct,'-.b','DisplayName',sprintf('5 percentile %s',varsel.name));
    hold off
    xlabel('Position along shore')
    ylabel(varsel.labl)
    title(sprintf('Case: %s',dst.(pntnames{1}).Description));    
    seltxt = getSelectionText(selection);
    subtitle(sprintf('Mean and Standard deviation for %s',seltxt))
    legend
end

%%
function bin_mean_drift(obj)
    %plot the selected interval (eg monthly) mean drift for each year 
    %at a selected point
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end

    %select summmer/winter and +ve/-ve drift
    selection = get_var_sampling([4,2,1,1],false);  %no stats selection
    if isempty(selection), return; end              
    [binvar,bintime] = subsample_variable(dst,varsel,selection);
    
    for i=1:npnts                        %loop for each point
        nint = size(binvar,3);           %number of intervals
        nper = size(binvar,2);           %number of periods
        for j=1:nper
            for k=1:nint
                meanVar = mean(binvar{i,j,k},'omitnan');
                stdVar = std(binvar{i,j,k},'omitnan');
                peclet = meanVar./stdVar;
                monthlyMean(i,j,k) = meanVar;
                monthlyPeclet(i,j,k) = peclet;
            end
        end
    end
    downcoast = monthlyMean; upcoast = monthlyMean;
    downcoast(monthlyPeclet>-1) = NaN;     %downcoast advection for Pe<=-1
    upcoast(monthlyPeclet<1) = NaN;        %upcoast advection for Pe>=1
    mnmxMn = minmax(monthlyMean);

    seltxt = getSelectionText(selection);
    bins = {'All','Year', 'Quarter', 'Month', 'Week', 'Dai', 'Hour'};
    bintxt = bins{selection{1}};
    ok = 1;
    while ok>0
        [selpnt,ok] = listdlg('Name','Plot profile', ...
                                 'PromptString','Select variable', ...
                                 'ListSize',[200,300], ...
                                 'SelectionMode','single', ...
                                 'ListString',pntnames);  
        if ok==0, continue; end 

        point.var = squeeze(monthlyMean(selpnt,:,:));
        point.down = squeeze(downcoast(selpnt,:,:));
        point.up = squeeze(upcoast(selpnt,:,:));
        point.pec = squeeze(monthlyPeclet(selpnt,:,:));
        mtime.int = unique(bintime.intervals);
        mtime.per = string(bintime.periods);
             
        plotxt = {dst.(pntnames{selpnt}).Description,pntnames{selpnt},...
                                                            seltxt,bintxt};
        %plot the interval mean drift as a line plot for selected point
        plotPeriodLines(mtime,point,varsel,plotxt);
 
        %plot the interval mean drift as a surface plot for selected point
        plotPeriodSurface(mtime,point,mnmxMn,varsel,plotxt);

        %plot the peclet ratio as a surface for the selected point
        % plotPeriodPeclet(mtime,point,varsel,plotxt);  %mainly just a check
    end
end
%%
function summary_point_drift(obj,msgtxt)
   %summary plot of monthly and annual drift at a point (littoraldriftstats)
   if ~isa(obj,'WRM_SedimentTransport'), getdialog(msgtxt); return; end

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
    %summary plot of selected statistical property and period for all 
    %alongshore points
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end
    mtime = dst.(pntnames{1}).RowNames;
    dt = mode(diff(mtime));

    sample = setDownsampleSettings();
    if isempty(sample), return; end

    if strcmp(sample.func,'pct95')
        func = '@(x) prctile(x,95)';
    elseif strcmp(sample.func,'pct5')
        func = '@(x) prctile(x,5)';
    else
        func = sample.func;
    end

    nfunc = ['nan',func];
    for i=1:npnts
        Var = dst.(pntnames{i}).(varsel.name);
        if strcmp(sample.func,'sum')
            Var = Var*seconds(dt);
        end
        if strcmp(sample.type,'+ve')
            Var(Var<0) = NaN;     %mask all negative values
        elseif strcmp(sample.type,'-ve')
            Var(Var>0) = NaN;     %mask all positive values
        end
       

        [mt,mvar(:,i)] = downsample(Var,mtime,sample.binsize,nfunc); 
    end

    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);
    [X,Y] = meshgrid(1:npnts,datenum(mt)); %#ok<DATNM>
    contourf(ax,X,Y,mvar);
    % Format the axes to display datetime
    if strcmp(sample.binsize,'month')
        datetick('y', 'mmm-yy', 'keepticks'); %#ok<DATIC>
    else
        datetick('y', 'yyyy', 'keepticks'); %#ok<DATIC>
    end
    axis tight
    colormap(cmap_selection(19));
    hc = colorbar;
    hc.Label.String = varsel.labl;
    xlabel('Point number')
    ylabel(sprintf('Time (%s)',sample.binsize))
    title(sprintf('Case: %s',dst.(pntnames{1}).Description));
    subtitle(sprintf('Downsampled %s %s using %s(%s)',sample.type,...
                                      varsel.name,sample.func,sample.binsize));
end

%%
function drift_peclet(obj,msgtxt)
    %plots to examine Peclet ratio using monthly/annual sampling(see Kahl, et al, 2024)
    if ~isa(obj,'WRM_SedimentTransport'), getdialog(msgtxt); return; end

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end
    %select summmer/winter and +ve/-ve drift
    selection = get_var_sampling([4,2,1,1],false);  %no stats selection
    if isempty(selection), return; end              
    [binvar,bintime] = subsample_variable(dst,varsel,selection);
            % 
            % dst = obj.Data;
            % pntnames = fieldnames(dst);
            % npnts = length(pntnames);
            % var = getVariable(dst,pntnames,1);
            % if isempty(var), return; end
        
            % mtime = dst.(pntnames{1}).RowNames;
            % monthlyMean = zeros(npnts,nper*nint); 
            % monthlyStdev = monthlyMean; monthlyPeclet = monthlyMean;
            % annualMean = zerod(npnts,nper); 
            % annualStdev = annualMean; annualPeclet = annualMean;
    for i=1:npnts
                % Var = dst.(pntnames{i}).(varsel.name);        
                % Var(abs(Var)<calms.value) = NaN; %remove near zero values
                % [~,binvar,bintime] = binned_variable(Var,mtime,'week','year');
        nint = size(binvar,3);           %number of intervals
        nper = size(binvar,2);           %number of periods
        nyr = 0;
        for j=1:nper
            for k=1:nint
                meanVar = mean(binvar{i,j,k},'omitnan');
                stdVar = std(binvar{i,j,k},'omitnan');
                    % nint = size(binvar,2);           %number of months
                    % nper = size(binvar,1);           %number of years
                    % nyr = 0;
                    % for j=1:nper
                    %     for k=1:nint
                ptle95 =  prctile(binvar{i,j,k},95);
                ptle05 =  prctile(binvar{i,j,k},5);
                    % meanVar = mean(binvar{j,k},'omitnan');
                    % stdVar = std(binvar{j,k},'omitnan');
                peclet = meanVar/stdVar;
                if isnan(peclet) || isinf(peclet)
                    peclet = 0;
                elseif peclet>-varsel.pecthr && peclet<varsel.pecthr
                    peclet = 0;
                end
                interval95tile(i,nyr+k) = ptle95;
                interval05tile(i,nyr+k) = ptle05;
                intervalMean(i,nyr+k) = meanVar;
                intervalStdev(i,nyr+k) = stdVar;
                intervalPeclet(i,nyr+k) = peclet;
            end            
            annualMean(i,j) = mean(intervalMean(i,(nyr+1:nyr+nint)),'omitnan');
            annualStdev(i,j) = sqrt(sum(intervalStdev(i,(nyr+1:nyr+nint)).^2,'omitnan')/nint^2);
            peclet = annualMean(i,j)/annualStdev(i,j);
            if isnan(peclet) || isinf(peclet)
                peclet = 0;
            elseif peclet>-varsel.pecthr && peclet<varsel.pecthr
                peclet = 0;
            end
            annualPeclet(i,j) = peclet;
            nyr = nyr+nint;
        end
    end

    bins = {'All','Year', 'Quarter', 'Month', 'Week', 'Dai', 'Hour'};
    bintxt = bins{selection{1}};
    ptxt = ': (Peclet ratio: >1 blue o; <1 yellow o)';
    subtxt = @(w,x,y,z) sprintf('%sly %s (Calms <%s m^3/yr) %s',w,x,y,z);
    % startdate = datetime(year(bintime.periods(1)),1,1);   %force full year
    % enddate = datetime(year(bintime.periods(end)),12,31); %to match variable
    bins = datenum(bintime.intstart);  
    [X,Y] = meshgrid(1:npnts,bins);
    pointdown = intervalPeclet; pointup = intervalPeclet;
    pointdown(pointdown>-varsel.pecthr) = NaN; pointup(pointup<varsel.pecthr) = NaN;

    %plot monthly mean peclet ratio as a scatter plot (position,time)
    desctxt = sprintf('Peclet ratio (<-%.1f or >%.1f)',varsel.pecthr,varsel.pecthr);
    desc = struct('case',dst.(pntnames{1}).Description,'var',desctxt);    
    axm = plotPeclet(pointup,bins,desc,0);  
    hold on
    scatter3(axm,X,Y,-pointdown',10,'y','filled','MarkerEdgeColor','k')
    hold off  
    view(2)
    datetick('y', 'yyyy'); %#ok<DATIC>
    subtitle(axm,subtxt(bintxt,'peclet ratio',varsel.calms.text,ptxt))

    %annual mean peclet ratio as a surface plot (position,time)
    annualPeclet(annualPeclet==0) = NaN;
    axa = plotPeclet(annualPeclet,bintime.periods,desc,1);
    subtitle(axa,subtxt('Annual','peclet',varsel.calms.text,ptxt))

    %plot monthly mean as a surface plot (position,time)
    desc = struct('case',dst.(pntnames{1}).Description,'var','Monthly mean');
    axm = plotPeclet(intervalMean,bins,desc,1);  
    hold on
    scatter3(axm,X,Y,-pointdown'*axm.ZLim(2),10,'MarkerEdgeColor','y')
    scatter3(axm,X,Y,pointup'*axm.ZLim(2),10,'MarkerEdgeColor',[0, 0.447, 0.741])    
    hold off
    subtitle(axm,subtxt(bintxt,'mean',varsel.calms.text,''))

    %plot std.dev. as a surface plot (position,time)
    desc = struct('case',dst.(pntnames{1}).Description,'var','Monthly Std.dev.');
    axm = plotPeclet(intervalStdev,bins,desc,1);  
    axm.CLim = [-0.01,0.01];
    subtitle(axm,subtxt(bintxt,'Std.dev.',varsel.calms.text,''))

    %plot monthly 95%tile as a surface plot (position,time)
    desc = struct('case',dst.(pntnames{1}).Description,'var','Monthly 95 & 5 percentile');
    interval95tile(abs(interval95tile)<1e-3) = 0;
    interval05tile(abs(interval05tile)<1e-3) = 0;
    monthlyPrctile = interval95tile+interval05tile; %assumes no overlap
    monthlyPrctile(monthlyPrctile==0) = NaN;
    axm = plotPeclet(monthlyPrctile,bins,desc,1);  
    subtitle(axm,subtxt(bintxt,'percentile',varsel.calms.text,'95%(+ve) and 5%(-ve)'))
end

%%
function cluster_peclet(obj,msgtxt)    
    %plots to Peclet ratio using cluster sampling (see Kahl, et al, 2024)
    if ~isa(obj,'WRM_SedimentTransport'), getdialog(msgtxt); return; end

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames,1); %selects Qs without prompting user
    if isempty(varsel), return; end
    mtime = dst.(pntnames{1}).RowNames;

    ans0 = questdlg('Use absolute values of drift or +/- values?','Clusters',...
                                            'abs(Qs)','+/-(Qs)','abs(Qs)');

    %NB: options defines the variables used in get clusters and includes 
    %additional variables used in mergeSelection for posnegClusters
    options = setClusterOptions(dst.(pntnames{1}).(varsel.name));
    if strcmp(ans0,'abs(Qs)')        
        [cluster,options] = absClusters(options,dst,calms,pecthr);
    else
        [cluster,options] = posnegClusters(options,dst,calms,pecthr);        
    end

    clustpoint = []; clustpec = []; clustints = [];
    for i=1:npnts
        nint = length([cluster.Ints{i,:}]);
        clustpoint = [clustpoint;repmat(i,nint,1)];
        clustints = [clustints;[cluster.Ints{i,:}]'];
        clustpec = [clustpec;[cluster.Peclet{i,:}]']; 
    end
  
    %plot cluster mean peclet ratio as a surface plot (position,time)
    desctxt = sprintf('Peclet ratio (<-%.1f or >%.1f)',pecthr,pecthr);
    desc = struct('case',dst.(pntnames{1}).Description,'var',desctxt);
    bintime = mtime(1):caldays(1):mtime(end);
    x = 1:npnts;
    % Define a grid for interpolation
    [xq, yq] = meshgrid(x, datenum(bintime));    %#ok<DATNM>
    % Interpolate scattered data onto the grid
    zq = griddata(clustpoint,datenum(clustints),clustpec, xq, yq, 'linear'); %#ok<DATNM> % 'linear', 'nearest', or 'cubic'
    axc = plotPeclet(zq',bintime,desc,1);    
    axc.CLim = [-2,2];
    % save('clusterplot',"xq","yq","zq","bintime","desc");
    txt1 = 'Clusters: Downdrift advection (Pe>1); Updrift advection (Pe<-1)';
    if strcmp(ans0,'abs(Qs)') 
        txt2 = sprintf('Absolute - Calms <%s m^3/yr; Threshold: %0.4f; Interval: %dd; Min duration: %dd',...
                    calms.text,options.threshold,...
                    options.clint,options.mincluster);
    else
        txt2 = sprintf('Pos/Neg - Calms <%s m^3/yr; Threshold: %0.4f/%0.4f;\n          Interval: %dd/%dd; Min duration: %dd/%dd',...
                    calms.text,...
                    options.pos.threshold,options.neg.threshold,...
                    options.pos.clint, options.neg.clint,...
                    options.pos.mincluster, options.neg.mincluster);
    end
    subtitle(axc,sprintf('%s\n%s',txt1,txt2))   
end
%%
function multi_rose_plots(obj)
    %plot multiple wave roses in single go
    dst = obj.Data;
    pntnames = fieldnames(dst);
    varname = dst.(pntnames{1}).VariableNames;
    vardesc = dst.(pntnames{1}).VariableDescriptions;

    promptxt = {'To add reference line, enter angle to degTN:',...
                'Variable intensity subdivisions (eg 0 0.5 1 ...)',...
                'Percentage circles to draw (eg 10 20 30)',...
                'Number of direction intervals (default is 36)'};
    rinp = {'','','',''};

    ok = 0;
    while ok<1
        %get the variable and points to use
        varsel = getVariable(dst,pntnames,[],0); %select from all and don't prompt for thresholds
        if isempty(varsel), ok = 1; continue; end
        [selpnt,sok] = listdlg('Name','Plot profile', ...
                            'PromptString','Select variable', ...
                            'ListSize',[200,300], ...
                            'SelectionMode','multiple', ...
                            'ListString',pntnames);
        if sok==0, ok = 1; continue; end

        %get the direction to use

        idvar = find(contains(vardesc,'direction'));
        seldir = 1;
        if numel(idvar)>1  %user needs to select direction
            seldir = listdlg('Name','Directions', ...
            'PromptString','Select direction to use','ListSize',[300,150],... ...
            'SelectionMode','single','ListString',vardesc(idvar));
        end
        if isempty(seldir), seldir = 1; end

        %set the rose plot scaling settings
        rinp = inputdlg(promptxt,'Rose plot',1,rinp);
        rose.theta = parseInput(rinp{1}); %if vector should be same length as numel(selpnt)
        rose.di = parseInput(rinp{2});
        rose.ci = parseInput(rinp{3});
        rose.nd = parseInput(rinp{4});
        %loop to create plots for selected points
        for i=1:numel(selpnt)
            ipnt = selpnt(i);
            hfig = figure('Name','Rose plot','Tag','PlotFig');
            figax = axes(hfig); %#ok<LAXES>
            dir = dst.(pntnames{ipnt}).(varname{idvar(seldir)}); %selected direction variable
            var = dst.(pntnames{ipnt}).(varsel.name);            %selected variable    
            casedesc = dst.(pntnames{ipnt}).Description;
            %title using variable-case-point
            %titletxt = sprintf('%s for %s at %s',varsel.desc,casedesc,pntnames{ipnt});
            %title using case-point-shore_angle
            titletxt = sprintf('%s at %s, theta=%d dTN',casedesc,pntnames{ipnt},rose.theta(i));
            wind_rose(dir,var,'parent',figax,'dtype','meteo',...
                'shore',rose.theta(i),'nd',rose.nd,'di',rose.di,'ci',rose.ci,...
                'labtitle',titletxt,'lablegend',varsel.labl);
        end
    end

    %-nested function------------------------------------------
    function var = parseInput(vartxt)
        if isempty(vartxt)
            var = [];
        else
            var = str2num(vartxt); %#ok<ST2NM> parsing scalar and vector
        end
    end
end

%% ------------------------------------------------------------------------
% Utility functions for data sampling
%--------------------------------------------------------------------------
function [cluster,userops] = absClusters(options,dst,calms,pecthr)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames,1);  %selects Qs without prompting user
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
    
                Var = dst.(pntnames{sel}).(varsel.name); 
                vardst = getDSTable(dst.(pntnames{sel}),'VariableNames',varsel.name);
                vardst.(varsel.name) = abs(vardst.(varsel.name));
                [idcls,userops] = getclusters(vardst,options);
                userops.mincluster = options.mincluster;
                userops.isplot = options.isplot;
                mergeAbsClusters(Var,mtime,idcls,userops);
                medges = mergeAbsClusters(Var,mtime,idcls,userops);
                if isempty(medges)
                    getdialog('No clusters found. Change threshold or minimum duration of cluster')
                else
                    ans1 = questdlg('Use selected options or examine another point?',...
                             'Clusters','Use selected','New point','Use selected');
                    if strcmp(ans1,'Use selected'), ok = 0; end
                end
        end
    
        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');        
    elseif  strcmp(ans2,'No')
        userops = options;
    end

    if strcmp(ans2,'Quit'), return; end

    userops.isplot = false; %supress plots in for loop
    for i=1:npnts        
        Var = dst.(pntnames{i}).(varsel.name);        
        Var(abs(Var)<calms.value) = NaN; %remove near zero values
        vardst = getDSTable(dst.(pntnames{i}),'VariableNames',varsel.name);
        vardst.(varsel.name) = abs(vardst.(varsel.name)); %use absolute values for intervals

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
            elseif peclet>-pecthr && peclet<pecthr %#ok<BDSCI>
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
function [cluster,userops]  = posnegClusters(options,dst,calms,pecthr)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames,1);  %selects Qs without prompting user
    if isempty(var), return; end
    mtime = dst.(pntnames{1}).RowNames;

    %default to use same options for postive and negative drift
    userops.pos = options; userops.neg = options; 
    cluster = [];

    %set a minimum trheold to remove zero values
    ans2 = questdlg('Check settings for selected points?','Clusters','Yes','No','Quit','Yes');
    if strcmp(ans2,'Quit')
        return; 
    elseif strcmp(ans2,'Yes')
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

                Var = dst.(pntnames{sel}).(varsel.name); 
                posdst = getDSTable(dst.(pntnames{sel}),'VariableNames',varsel.name);
                [idpos,userops.pos] = getclusters(posdst,options);
                userops.pos.mincluster = options.mincluster;
                userops.pos.isplot = options.isplot;
                negdst = posdst;
                negdst.(varsel.name) = negdst.(varsel.name)*-1;
                [idneg,userops.neg] = getclusters(negdst,options);
                userops.neg.mincluster = options.mincluster;
                userops.neg.isplot = options.isplot;
                mergePosNegClusters(Var,mtime,idpos,idneg,options);

                ans1 = questdlg('Use selected options or examine another point?',...
                         'Clusters','Use selected','New point','Use selected');
                if strcmp(ans1,'Use selected'), ok = 0; end
        end

        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');
    elseif strcmp(ans2,'No')
        %use same options for postive and negative drift
    end

    if strcmp(ans2,'Quit'), return; end

    options.isplot = false;
    for i=1:npnts        
        Var = dst.(pntnames{i}).(varsel.name);        
        Var(abs(Var)<calms.value) = NaN; %remove near zero values
        posdst = getDSTable(dst.(pntnames{i}),'VariableNames',varsel.name);

        % find clusters based on results from peak selection
        idposcls = getVarClusters(posdst,userops.pos);
        negdst = posdst;
        negdst.(varsel.name) = negdst.(varsel.name)*-1;
        idnegcls = getVarClusters(negdst,userops.neg);

        %merge any overlaps to define intervals to be used
        medges = mergePosNegClusters(Var,mtime,idposcls,idnegcls,options); %only uses mincluster field in options
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
            elseif peclet>-pecthr && peclet<pecthr %#ok<BDSCI>
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
    % find clustecrs based on results from peak selection
    pk_date = mtime(idpks);    %datetime of peak
    pk_vals = var(idpks);  %value of peak
    idcls = clusters(pk_date,pk_vals,days(opts.clint));
end

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

    if ~isempty(medges) && opts.isplot
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

    if ~isempty(medges) && opts.isplot
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

%% ------------------------------------------------------------------------
% Utility functions for variable selection
%--------------------------------------------------------------------------
function [calms,pecthr] = calmsThreshold()
    %set the calms threshold to apply to the data
    calmsthreshold = 100;  %"calms" are drift rates less than threshold
                           % 100m^3/yr ~= 3e-6 m^3/s; 
    promptxt = {'Calms threshold (Hs (m), Qs (m^3/yr, etc):','Peclet plotting threshold'};
    defaults = {num2str(calmsthreshold),'1'};
    answer = inputdlg(promptxt,'Drift',1,defaults);
    if isempty(answer), answer = defaults; end
    calms.value = str2double(answer{1})/31556952;  %value from mobj.Constants.y2s
    calms.text = answer{1};
    pecthr = str2double(answer{2});
end

%%
function varsel = getVariable(dst,pntnames,sel,isthr)
    %select a variable to use in the plot
    varname = dst.(pntnames{1}).VariableNames;
    vardesc = dst.(pntnames{1}).VariableDescriptions; 
    varlabl = dst.(pntnames{1}).VariableLabels; 

    if nargin<3 || isempty(sel)
        [sel,ok] = listdlg('Name','Plot profile', ...
                                     'PromptString','Select variable', ...
                                     'ListSize',[200,80], ...
                                     'SelectionMode','single', ...
                                     'ListString',vardesc);    
        if ok==0, varsel = []; return; end 
    end

    varsel.name = varname{sel};
    varsel.desc = vardesc{sel};
    varsel.labl = varlabl{sel};

    varsel.isdrift = false;
    if contains(varsel.name,'Q'), varsel.isdrift = true; end
    % 
    if nargin<4 || isthr
        [varsel.calms,varsel.pecthr] = calmsThreshold();
    end
end

%%
function sample = setDownsampleSettings()
    %check that the current sit parameters settings are correct
    %modifications used to update RunParams stored with Case
    type = 'all';
    period = 'year';
    func = 'sum';
    promptxt = {'Variable sampling (all, +ve, -ve)',...
                'Bin size (year, month, etc)',...
                'Function (sum, mean, min, max, etc + pct95, pct5)'};
    defaults = {type,period,func};
    data = inputdlg(promptxt,'Drift settings',1,defaults);
    if isempty(data), sample = []; return; end  %no change to default settings
    sample.type = data{1};
    sample.binsize = data{2};
    sample.func = data{3}; 
end

%%
function seltxt = getSelectionText(selection)
    %use selection to construct text for plot titles
    seltxt = 'All';
    if selection{3}>1
        stxt = {'All','Winters','Summers'};
        seltxt = sprintf('%s %s',seltxt,stxt{selection{3}});
    else
        seltxt = [seltxt,' years'];
    end
        
    if selection{4}>1
        stxt = {'All','+ve only','-ve only'};
        seltxt = sprintf('%s, %s',seltxt,stxt{selection{4}});
    end
end

%% ------------------------------------------------------------------------
% Utility functions for plotting 
%--------------------------------------------------------------------------

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
function [ax,hs] = plotPeclet(var,bintime,desc,issurf)
    %monthly or cluster mean peclet ratio as a surface plot (position,time)
    npnts = size(var,1);
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);    
    grid on
    [X,Y] = meshgrid(1:npnts,bintime); 
    if issurf
        hs = surf(ax,X,Y,var');
        shading interp
        ax.Layer = 'top';
        colormap(cmap_selection(19));
        hc = colorbar;
        hc.Label.String = desc.var;        
    else
        hs = scatter3(ax,X,Y,var',10,'b','filled','MarkerEdgeColor','w');
    end
    view(2)
    axis tight
    datetick('y', 'yyyy'); %#ok<DATIC>
    xlabel('Position along shore')
    ylabel('Year')
    title(sprintf('Case: %s',desc.case));  
end

%%
function ax = plotPeriodLines(mtime,point,varsel,plotxt)
    %plot lines of the selected variable and intervals in each year period
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);       
    hold on
    plot(ax,mtime.int,point.var(1,:),'DisplayName','Year','ButtonDownFcn',@godisplay);
    if varsel.isdrift
        plot(ax,mtime.int,point.down(1,:),'-og','DisplayName','Pe<-1',...
            'LineWidth',0.8,'MarkerSize',4,'ButtonDownFcn',@godisplay);
        plot(ax,mtime.int,point.up(1,:),'-ob','DisplayName','Pe>1',...
            'LineWidth',0.8,'MarkerSize',4,'ButtonDownFcn',@godisplay);
    end
    %
    for i=2:size(point.var,1)
        p1 = plot(ax,mtime.int,point.var(i,:),'DisplayName',mtime.per(i),'ButtonDownFcn',@godisplay);
        p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
        if varsel.isdrift
            p1 = plot(ax,mtime.int,point.down(i,:),'-og','LineWidth',0.8,...
            'MarkerSize',4,'DisplayName',mtime.per(i),'ButtonDownFcn',@godisplay);
            p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            p1 = plot(ax,mtime.int,point.up(i,:),'-ob','LineWidth',0.8,...
            'MarkerSize',4,'DisplayName',mtime.per(i),'ButtonDownFcn',@godisplay);
            p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end 
    end
    p1 = plot(ax,ax.XLim,[0,0],'--k');  %zero line
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold off
    xlabel(plotxt{4})
    ylabel(varsel.labl)
    title(sprintf('%s for %s (%s) %s',varsel.desc,plotxt{1:3}));             
    subtitle(sprintf('%sly Means for each year of data set',plotxt{4}))
    legend
end

%%
function ax = plotPeriodSurface(mtime,point,mnmxMn,varsel,plotxt)
    %plot a surface of the selected variable and any peclet events    
    hf = figure('Name','SedTrans','Tag','PlotFig');
    sax = axes(hf);
    [m,n] = size(point.var);
    [X,Y] = meshgrid(1:n,1:m);
    surf(sax,X,Y,point.var)
    shading interp
    view(2)
    axis tight
    ax.Layer = 'top';
    colormap(cmap_selection(19));
    hc = colorbar;
    hc.Label.String = varsel.labl;
    xlabel(plotxt{4})
    ylabel('Year')
    idx = str2double(sax.YTickLabel);
    sax.YTickLabel = mtime.per(idx);
    sax.CLim = mnmxMn;
    title(sprintf('%s for %s (%s) %s',varsel.desc,plotxt{1:3})); 
    subtitle(sprintf('%sly Means (Peclet ratio: >1 blue o; <1 yellow o)',plotxt{4}))
    hold on
    %add points where peclet exceeds above the surface
    scatter3(sax,X,Y,-point.down*2,20,'y','filled','MarkerEdgeColor','k')
    scatter3(sax,X,Y,point.up*2,20,'b','filled','MarkerEdgeColor','w')
    hold off
end

%%
function ax = plotPeriodPeclet(mtime,point,varsel,plotxt)
    %plot a surface of the period peclet value and the peclet events
    hf = figure('Name','SedTrans','Tag','PlotFig');
    sax = axes(hf);
    [m,n] = size(point.pec);
    [X,Y] = meshgrid(1:n,1:m);
    surf(sax,X,Y,point.pec)
    shading interp
    view(2)
    axis tight
    ax.Layer = 'top';
    colormap(cmap_selection(19));
    hc = colorbar;
    hc.Label.String = varsel.labl;
    xlabel(plotxt{4})
    ylabel('Year')
    idx = str2double(sax.YTickLabel);
    sax.YTickLabel = mtime.per(idx);
    sax.CLim = [-2,2];
    title(sprintf('%s for %s (%s) %s',varsel.desc,plotxt{1:3}));
    subtitle(sprintf('%sly peclet ratio (Peclet ratio: >1 filled red o; <1 red o)',plotxt{4}))
    hold on
    %add points where peclet excceeds above the surface

    scatter3(sax,X,Y,point.down./point.down*2,20,'r')
    scatter3(sax,X,Y,point.up./point.up*2,20,'r','filled')
    hold off
end