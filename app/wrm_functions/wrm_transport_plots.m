function wrm_transport_plots(obj,mobj,option) %#ok<INUSD>
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
        %list as per WRM_SedimentTransport.transportPlots 'listxt' variable
        case 'Annual Mean Drift'
            ann_mean_drift(obj);
        case 'Binned Mean Drift'
            bin_mean_drift(obj);
        case 'Summary Point Drift'
            summary_point_drift(obj,msgtxt);
        case 'Summary Shore Drift'
            summary_shore_drift(obj);
        case 'Duration Exceedance'
            duration_exceedance(obj);
        case 'Binned Statistics'
            monthly_drift_peclet(obj,msgtxt);
        case 'Absolute Cluster Statistics'
            cluster_peclet(obj,msgtxt,0);
        case 'Pos/Neg Cluster Statistics'
            cluster_peclet(obj,msgtxt,1);
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
    ok = 0;
    while ok<1
        selection = get_var_sampling([1,1,1,1],false);  %no stats selection
        if isempty(selection), return; end              
        selection(1:2) = {1,1};  %force selection of all years
        [binvar,~] = subsample_variable(dst,varsel,selection);
        if ~isempty(binvar), ok = 1; end
    end

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
    title(sprintf('Case: %s',varsel.case));    
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
    ok = 0;
    while ok<1
        selection = get_var_sampling([4,2,1,1],false);  %no stats selection
        if isempty(selection), return; end              
        [binvar,bintime] = subsample_variable(dst,varsel,selection);
        if ~isempty(binvar), ok = 1; end
    end

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
             
        plotxt = {varsel.case,pntnames{selpnt},seltxt,bintxt};                                                            
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
    grid on
    
    [X,Y] = meshgrid(1:npnts,datenum(mt)); %#ok<DATNM>
    surf(ax,X,Y,mvar);
    shading interp
    view(2)
     
    % Format the axes to display datetime
    ax.YTick = datenum(mtime(1):calyears(5):mtime(end)); %#ok<DATNM>
    datetick('y', 'yyyy', 'keepticks'); %#ok<DATIC>
    % if strcmp(sample.binsize,'month')
    %     datetick('y', 'mmm-yy', 'keepticks'); %#ok<DATIC>
    % else
    %     datetick('y', 'yyyy', 'keepticks'); %#ok<DATIC>
    % end
    axis tight
    ax.Layer = 'top'; %moves grid above surface  
    
    zmap = struct('Z',mvar,'zeroLevel',0);
    colormap(cmap_selection(23,zmap));
    hc = colorbar;
    zlabel = varsel.labl;
    if strcmp(sample.func,'sum')
        zlabel = sprintf('Total transport/%s (m^3)',sample.binsize); 
    end
    hc.Label.String = zlabel;
    xlabel('Point number')
    ylabel(sprintf('Time (%sly bins)',sample.binsize))
    title(sprintf('Case: %s',varsel.case));
    varsel.subtitle = sprintf('Downsampled %s %s using %s(%s)',sample.type,...
                                      varsel.name,sample.func,sample.binsize);
    subtitle(varsel.subtitle);
    %summary of drift by reach
    reachPlot(npnts,mt,mvar,varsel,sample);
end

%%
function duration_exceedance(obj)
    %plot the duration of exceedance for selected thresholds
    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end
    mtime = dst.(pntnames{1}).RowNames;

    tdur = years(1);       %annual
    tstep = years(0.08);   %monthly

    txt = sprintf('Set threshold and sampling durations\nDefine threshold:');
    promptxt = {txt,'Positive or Negative (1/0)','Sampling period (y,d,h,m,s):',...
                'Time step interval (y,d,h,m,s)','Method (mean, sum,... xx prctile}'};
    defaults = ['0','1',cellstr(tdur),cellstr(tstep),'mean'];  %'95 prctile'
    answer = inputdlg(promptxt,'MovingTime',1,defaults);
    if isempty(answer), return; end
    threshold = str2double(answer{1});
    ispos = logical(str2double(answer{2}));
    tdur = str2duration(answer{3});
    tstep = str2duration(answer{4});
    method = answer{5};

    for i=1:npnts
        var = dst.(pntnames{i}).(varsel.name);
        if ~ispos % for peak negative values invert variable
            var = -var;
            %invert threshold if specified as a negative value
            if threshold<0; threshold = -threshold; end
        end

        [stid,edid] = zero_crossing(var,threshold);
        if isempty(stid)
            hw = warndlg('No zero-crossings found'); waitfor(hw); continue;
        elseif numel(stid)<=1
            hw = warndlg('Only one zero-crossing found. Try changing threshold');
            waitfor(hw); continue;
        end

        if stid(1)>edid(1)   %correct order for exceedances above threshold
            stid = stid(1:end-1);
            edid = edid(2:end);
        end
        vardur = mtime(edid)-mtime(stid);
        vardur.Format = 'h';
        vardur = hours(vardur);
        [tmi,vmi] = movingtime(vardur,mtime(stid),tdur,tstep,method,0); %no prompt required
        if i==1 || numel(vmi)==size(vm,1)            
            tm = tmi;
        else
            vmi = interp1(tmi,vmi,tm,'linear');
        end
        vm(:,i) = vmi;
    end

    vartxt = sprintf('%s %s duration (hrs) in each period',varsel.desc,method);
    desc = struct('case',varsel.case,'var',vartxt);  
    vm(vm==0) = NaN;
    [ax,~] = plot2Dvariable(vm',tm,desc,1);
    if ispos, sgntxt = 'Positive'; else, sgntxt = 'Negative'; end
    ax.Subtitle.String = sprintf('%s %s; Sampling period %s; Time step %s',...
                                sgntxt,varsel.name,char(tdur),char(tstep));
    ax.Title.String = sprintf('%s: %.3g threshold',ax.Title.String,threshold);
end

%%
function monthly_drift_peclet(obj,msgtxt)
    %plots to examine Peclet ratio using monthly/annual sampling(see Kahl, et al, 2024)
    if ~isa(obj,'WRM_SedimentTransport'), getdialog(msgtxt); return; end

    dst = obj.Data;    
    pntnames = fieldnames(dst);
    %mtime = dst.(pntnames{1}).RowNames;
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames);
    if isempty(varsel), return; end

    %select summmer/winter and +ve/-ve drift
    ok = 0;
    while ok<1
        selection = get_var_sampling([4,2,1,1],false);  %no stats selection
        if isempty(selection), return; end              
        [binvar,bintime] = subsample_variable(dst,varsel,selection);
        if ~isempty(binvar), ok = 1; end
    end

    hw = waitbar(0,'Processing point 0');
    for i=1:npnts
        nint = size(binvar,3);           %number of intervals
        nper = size(binvar,2);           %number of periods
        nyr = 0;
        for j=1:nper
            annualData = [];
            for k=1:nint
                intstats(i,nyr+k) = wrm_interval_statistics(binvar{i,j,k});
                annualData = [annualData;binvar{i,j,k}];
            end
            nyr = nyr+nint;
            annstats(i,j) = wrm_interval_statistics(annualData);
        end

        % nyr = 0;
        % for j=1:nper
        %     annualData = [];
        %     for k=1:nint
        %         meanVar = mean(binvar{i,j,k},'omitnan');
        %         stdVar = std(binvar{i,j,k},'omitnan');
        %         ptle95 =  prctile(binvar{i,j,k},95);
        %         ptle05 =  prctile(binvar{i,j,k},5);
        %         peclet = meanVar/stdVar;
        %         peclet = checkPecletLimits(peclet,varsel,NaN,peclet);
        % 
        % 
        %         interval95tile(i,nyr+k) = ptle95;
        %         interval05tile(i,nyr+k) = ptle05;
        %         intervalMean(i,nyr+k) = meanVar;
        %         intervalStdev(i,nyr+k) = stdVar;
        %         intervalPeclet(i,nyr+k) = peclet;
        %         annualData = [annualData;binvar{i,j,k}];
        %     end   
        %     annualMean(i,j) = mean(annualData,'omitnan');
        %     annualStdev(i,j) = std(annualData,'omitnan');
        %     annpeclet = annualMean(i,j)/annualStdev(i,j);
        %     annpeclet = checkPecletLimits(annpeclet,varsel,0,0);
        %     annualPeclet(i,j) = annpeclet;
        %     nyr = nyr+nint; 
        % end
        waitbar(i/npnts,hw,sprintf('Processing point %d',i));
    end
    delete(hw)

    % intsz = num2cell(size(intstats));
    % interval95tile = reshape([intstats(:).ptle05],intsz{:});
    % interval05tile = reshape([intstats(:).ptle05],intsz{:});
    % intervalMean = reshape([intstats(:).mean],intsz{:});
    % intervalStdev = reshape([intstats(:).stdev],intsz{:});
    % intervalPeclet = reshape([intstats(:).Pe],intsz{:});
    % intervalIndex = reshape([intstats(:).index],intsz{:});
    % annsz = num2cell(size(annstats));
    % annualMean = reshape([annstats(:).mean],annsz{:});
    % annualPeclet = reshape([annstats(:).Pe],annsz{:});
    % annualIndex = reshape([annstats(:).index],annsz{:});

    %meta-data for bin selection
    bins = {'All','Year', 'Quarter', 'Month', 'Week', 'Day', 'Hour'};
    Season = {'All','Winters','Summers'};
    Direction = {'All','+ve only','-ve only'};
    bintxt = {bins{selection{1}}, bins{selection{2}}};
    seltxt = sprintf('%s for %s directions',Season{selection{3}},Direction{selection{4}});

    %meta-data for limit variables
    varsubtxt = @(w,x,y,z) sprintf('%sly %s (Calms <%s m^3/yr) %s',w,x,y,z);
    %indsubtxt = @(v,w,x,y,z) sprintf('%sly %s (Calms <%s m^3/yr) %s\n%s',v,w,x,y,z);
    % 
    % bintimes = datenum(bintime.intstart);  %#ok<DATNM>
    % [X,Y] = meshgrid(1:npnts,bintimes);       
    % zq = intervalPeclet';
    % zq(zq>2) = 2;
    % zq(zq<-2) = -2;
    % zq(zq==0) = NaN;
    % pointdown = zq; pointup = zq;
    % pointdown(pointdown>-varsel.pecthr) = NaN; pointup(pointup<varsel.pecthr) = NaN;
    % clup = annualPeclet'; cldn = -annualPeclet';
    % clup(clup<varsel.pecthr) = NaN; cldn(cldn<varsel.pecthr) = NaN; 
    % clup = clup./clup*2; cldn = cldn./cldn*2;

    %----------------------------------------------------------------------
    % plots selected options
    %----------------------------------------------------------------------
    plotoptions = {'Interval surface','Period surface',...
                   'Interval points','Period points',...
                   'Pos/Neg points','Reach summary'};
    ok = 0;
    while ok<1
        sel = listdlg('Name','Plot options', ...
            'PromptString','Select plot type:','ListSize',[200,150],... ...
            'SelectionMode','single','ListString',plotoptions);
        if isempty(sel) || strcmp(sel,'Quit'), ok = 1; continue; end

        numtimes = datenum(bintime.intstart);  %#ok<DATNM>
        varsel.bintxt = bintxt{1};
        varsel.iscluster = 0;
        if contains(plotoptions{sel},'Interval')
            [plotvar,meta] = getPlotVariable(intstats,varsel);               
        elseif contains(plotoptions{sel},'Period')
            varsel.bintxt = bintxt{2};
            [plotvar,meta] = getPlotVariable(annstats,varsel);
            numtimes = datenum(bintime.periods);   %#ok<DATNM>            
        else
            [plotvar,meta] = getPlotVariable(intstats,varsel); %???????
        end
        if isempty(plotvar), ok = 1; continue; end

        switch plotoptions{sel}
            case {'Interval surface','Period surface'}
                axm = plot2Dvariable(plotvar,numtimes,meta,1); 
            case {'Interval points','Period points'}
                axm = plot2Dvariable(plotvar,numtimes,meta,2);
            case 'Pos/Neg points'
                axm = plot2Dvariable(plotvar,numtimes,meta,3);
             
                % [X,Y] = meshgrid(1:npnts,bintimes);
                % hold on
                % scatter3(axc,X,Y,-pointdown*axc.ZLim(2),6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
                % scatter3(axc,X,Y,pointup*axc.ZLim(2),6,'b','filled','Marker','square','MarkerEdgeColor','none')  
                % hold off
        end
        subtitle(axm,varsubtxt(meta.bintxt,meta.name,meta.calms.text,seltxt))
        setAxisLimits(axm,npnts,'1980','2025');  %bespoke ********************

        % switch plotoptions{selection}
        %     case 'Peclet points'
        %     plot interval peclet ratio as a scatter plot (position,time)
        %     desctxt = sprintf('Peclet ratio (<-%.1f or >%.1f)',varsel.pecthr,varsel.pecthr);
        %     desc = struct('case',varsel.case,'var',desctxt);    
        %     axm = plotPeclet(pointup',bins,desc,0);  
        %     hold on
        %     scatter3(axm,X,Y,-pointdown,6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        %     hold off  
        %     view(2)
        %     datetick('y', 'yyyy'); %#ok<DATIC>
        %     subtitle(axm,pecsubtxt(bintxt{1},'peclet ratio',varsel.calms.text,seltxt,ptxt))
        %     setAxisLimits(axm,npnts,'1980','2025'); %bespoke ********************
        %     case 'Peclet surface'
        %     %plot peclet ratio as a color map surface (position,time)
        %     hfig = figure('Tag','PlotFig');
        %     axp = axes(hfig); %#ok<LAXES>
        %     zqnan = checkPecletLimits(zq,varsel,NaN,NaN);
        %     pcolor(axp,X,Y,zqnan,'EdgeColor','none','FaceColor','flat');
        %     cmap = colormap('turbo');
        %     colormap(flipud(cmap));
        %     datetick('y', 'yyyy'); %#ok<DATIC>
        %     xlabel('Position along shore')
        %     ylabel('Year')
        %     pectxt = 'Peclet ratio: >1 blue; <1 red';
        %     title(sprintf('Case: %s',varsel.case));  
        %     subtitle(axp,indsubtxt(bintxt{1},'peclet ratio',varsel.calms.text,seltxt,pectxt))
        %     setAxisLimits(axp,npnts,'1980','2025'); %bespoke ********************
        % 
        %     case 'Annual mean'
        %     %annual mean of variable as a surface plot + peclet ratio points (position,time)
        %     desc = struct('case',varsel.case,'var',...
        %                             ['Annual mean of ',lower(varsel.labl)]);
        %     annbins = datenum(bintime.periods); %#ok<DATNM>
        %     [xq, yq] = meshgrid(1:npnts,annbins); 
        %     axa = plot2Dvariable(annualMean,annbins,desc,1); 
        %     hold(axa,'on')
        %     scatter3(axa,xq,yq,clup*axa.ZLim(2),6,'b','filled','Marker','square','MarkerEdgeColor','none')
        %     scatter3(axa,xq,yq,cldn*axa.ZLim(2),6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        %     hold(axa,'off') 
        %     subtitle(axa,indsubtxt('Year','mean',varsel.calms.text,seltxt,ptxt))
        %     % axa.CLim = [-2,2];
        %     setAxisLimits(axa,npnts,'1980','2025'); %bespoke ********************
        % 
        %     case 'Interval mean'
        %     %plot interval mean as a surface plot (position,time) with scatter points
        %     desc = struct('case',varsel.case,'var',...
        %                      [bintxt{1},'ly mean of ',lower(varsel.labl)]);
        %     axc = plot2Dvariable(intervalMean,bintimes,desc,1);  
        %     hold on
        %     scatter3(axc,X,Y,-pointdown*axc.ZLim(2),6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        %     scatter3(axc,X,Y,pointup*axc.ZLim(2),6,'b','filled','Marker','square','MarkerEdgeColor','none')  
        %     hold off
        %     subtitle(axc,indsubtxt(bintxt{1},'mean',varsel.calms.text,seltxt,ptxt))
        %     setAxisLimits(axc,npnts,'1980','2025'); %bespoke ******************
        % 
        %     case  'Interval St.dev'
        %     %plot std.dev. as a surface plot (position,time)
        %     desc = struct('case',varsel.case,'var',...
        %                    [bintxt{1},'ly Std.dev. of',lower(varsel.labl)]);
        %     axs = plot2Dvariable(intervalStdev,bintimes,desc,1);  
        %     axs.CLim = [-0.01,0.01];
        %     subtitle(axs,varsubtxt(bintxt{1},'Std.dev.',varsel.calms.text,seltxt))
        %     setAxisLimits(axs,npnts,'1980','2025'); %bespoke ******************
        % 
        %     case 'Interval percentiles'
        %     %plot interval 95%tile as a surface plot (position,time)
        %     desc = struct('case',varsel.case,'var',...
        %           [bintxt{1},'ly 95 & 5 percentile of ',lower(varsel.labl)]);
        %     interval95tile(abs(interval95tile)<1e-3) = 0;
        %     interval05tile(abs(interval05tile)<1e-3) = 0;
        %     monthlyPrctile = interval95tile+interval05tile; %assumes no overlap
        %     monthlyPrctile(monthlyPrctile==0) = NaN;
        %     axp = plot2Dvariable(monthlyPrctile,bintimes,desc,1);  
        %     subtitle(axp,indsubtxt(bintxt{1},'percentile',varsel.calms.text,seltxt,'95%(+ve) and 5%(-ve)'))
        %     setAxisLimits(axp,npnts,'1980','2025'); %bespoke ******************
        % 
        %     case 'Reach summary'
        %     varsel.subtitle = varsubtxt(bintxt{1},'mean',varsel.calms.text,seltxt);
        %     reachPecletPlot(npnts,bintime.intstart,intervalPeclet',varsel);
        % end
    end
end

%%
function cluster_peclet(obj,msgtxt,isposneg)    
    %plots to Peclet ratio using cluster sampling (see Kahl, et al, 2024)
    if ~isa(obj,'WRM_SedimentTransport'), getdialog(msgtxt); return; end

    dst = obj.Data;
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    varsel = getVariable(dst,pntnames,1); %selects Qs without prompting user
    if isempty(varsel), return; end
    mtime = dst.(pntnames{1}).RowNames;

    % ans0 = questdlg('Use absolute values of drift or +/- values?','Clusters',...
    %                                         'abs(Qs)','+/-(Qs)','abs(Qs)');

    %NB: options defines the variables used in get clusters and includes 
    %additional variables used in mergeSelection for posnegClusters
    options = setClusterOptions(dst.(pntnames{1}).(varsel.name));
    % if strcmp(ans0,'abs(Qs)')    
    if isposneg
        [cluster,clustints,options] = posnegClusters(options,dst,varsel);  
    else
        [cluster,clustints,options] = absClusters(options,dst,varsel);        
    end
    if isempty(cluster), return; end

    % clustpoint = []; clustpec = []; clustints = []; clustmean = [];
    % for i=1:npnts
    %     nint = length([cluster.Ints{i,:}]);
    %     clustpoint = [clustpoint;repmat(i,nint,1)];
    %     clustints = [clustints;[cluster.Ints{i,:}]'];
    %     clustpec = [clustpec;[cluster.Peclet{i,:}]'];
    %     clustmean = [clustmean;[cluster.Mean{i,:}]'];
    % end
    % for i=1:npnts
    %     Clustints{i} = [cluster.Ints{i,:}];
    %     Clustpec{i} = [cluster.Peclet{i,:}];
    % end  
    % Clustints = clustints;
    % for i=1:npnts
    %     Clustpec{i} = [cluster(i,:).Pe];
    % end

    %plot cluster mean peclet ratio as a surface plot (position,time)
    % desctxt = sprintf('Peclet ratio (<-%.1f or >%.1f)',varsel.pecthr,varsel.pecthr);
    % desc = struct('case',varsel.case,'var',desctxt);
    bintime = mtime(1):caldays(1):mtime(end); 
    numtimes = datenum(bintime); %#ok<DATNM>
    x = 1:npnts;

    % y = datenum(bintimes); %#ok<DATNM>
    % [xq, yq] = meshgrid(x,y); 

    %meta-data for bin selection
    subtxt = @(x,y,z) sprintf('%s (Calms <%s m^3/yr) %s',x,y,z);

    %----------------------------------------------------------------------
    % plots selected options
    %----------------------------------------------------------------------
    plotoptions = {'Cluster surface','Cluster points','Reach summary'};

    % Define a grid for interpolation
    varsel.bintxt = '';
    varsel.iscluster = 1;
    
 
    
    % zq = zq';
    % %get up and down drift peclets
    % clup = zq; cldn = -zq;
    % clup(clup<varsel.pecthr) = NaN; cldn(cldn<varsel.pecthr) = NaN;   
    % clup = clup./clup*2; cldn = cldn./cldn*2;

    ok = 0;
    while ok<1
        selection = listdlg('Name','Plot options', ...
            'PromptString','Select plot type:','ListSize',[200,150],... ...
            'SelectionMode','single','ListString',plotoptions);
        if isempty(selection) || strcmp(selection,'Quit'), ok = 1; continue; end
        [cellvar,meta] = getPlotVariable(cluster,varsel);
        % meta = varsel;
        % Interpolate scattered data onto the grid

        plotvar = interval_interpolant(x,clustints,cellvar, x, bintime);
        % Clustvar = interval_interpolant(x,Clustints,Clustpec, x, bintime);
        % axm = plot2Dvariable(Clustvar,numtimes,meta,1); 
        %------------------------------------------------------------------
        % plots as selection options
        %------------------------------------------------------------------
        switch plotoptions{selection}
            case 'Cluster surface'
                axm = plot2Dvariable(plotvar,numtimes,meta,1); 
            case 'Cluster points' 
                axm = plot2Dvariable(plotvar,numtimes,meta,2); 
            case 'Reach summary'
        end

        %subtitle(axm,subtxt(meta.bintxt,meta.name,varsel.calms.text));      
        % Summary output text
        txt1 = subtxt(meta.bintxt,meta.name,varsel.calms.text);
        if isposneg
            txt2 = sprintf('Pos/Neg clusters; Calms <%s m^3/yr; Threshold: %0.4f/%0.4f;\n          Interval: %dd/%dd; Min duration: %dd/%dd',...
                        varsel.calms.text,...
                        options.pos.threshold,options.neg.threshold,...
                        options.pos.clint, options.neg.clint,...
                        options.pos.mincluster, options.neg.mincluster);
        else
            txt2 = sprintf('Absolute clusters; Calms <%s m^3/yr; Threshold: %0.4f; Interval: %dd; Min duration: %dd',...
                        varsel.calms.text,options.threshold,...
                        options.clint,options.mincluster);                    
        end
        subtitle(axm,sprintf('%s\n%s',txt1,txt2))  
        setAxisLimits(axm,npnts,'1980','2025'); %bespoke **********


        % % switch plotoptions{selection}
        %     case 'Peclet points' 
        %         %plot just the  peclet points when >1
        %         axm = plot2Dvariable(clup',y,desc,0); 
        %         hold(axm,'on')
        %         scatter3(axm,xq,yq,cldn,6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        %         hold(axm,'off') 
        %         view(2)
        %         datetick('y', 'yyyy'); %#ok<DATIC>
        %         subtitle(axm,subtxt('Peclet ratio',varsel.calms.text,ptxt))
        %         setAxisLimits(axm,npnts,'1980','2025'); %bespoke ******************
        % 
        %     case 'Peclet surface'
        %         %plot peclet ratio as a color map surface (position,time)
        %         hfig = figure('Tag','PlotFig');
        %         axp = axes(hfig); %#ok<LAXES>
        %         zqnan = checkPecletLimits(zq,varsel,NaN,NaN);  
        %         pcolor(axp,xq,yq,zqnan,'EdgeColor','none','FaceColor','flat');
        %         cmap = colormap('turbo');
        %         colormap(flipud(cmap));
        %         datetick('y', 'yyyy'); %#ok<DATIC>
        %         xlabel('Position along shore')
        %         ylabel('Year')
        %         pectxt = ': (Peclet ratio: >1 blue; <1 red)';
        %         title(sprintf('Case: %s',desc.case));  
        %         subtitle(axp,subtxt('Clusters peclet ratio',varsel.calms.text,pectxt))
        %         axp.CLim = [-2,2];
        %         setAxisLimits(axp,npnts,'1980','2025'); %bespoke **********
        % 
        %     case 'Cluster mean'
        %         %plot surface of interval values with peclet points when >1
        %         desc.var = varsel.labl;                
        %         F = scatteredInterpolant(clustpoint,datenum(clustints),clustmean,'nearest'); %#ok<DATNM> %faster than griddata
        %         zqc = F(xq,yq);
        %         axc = plot2Dvariable(zqc',y,desc,1);  
        %         %axc.CLim = [-options.threshold*3,options.threshold*3];
        %         hold(axc,'on')
        %         scatter3(axc,xq,yq,clup,6,'b','filled','Marker','square','MarkerEdgeColor','none')
        %         scatter3(axc,xq,yq,cldn,6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        %         hold(axc,'off') 
        %         setAxisLimits(axc,npnts,'1980','2025'); %bespoke ******************
        % 
        %         % Summary output text
        %         txt1 = 'Clusters: Downdrift advection (Pe>1, blue); Updrift advection (Pe<-1, yellow)';
        %         if strcmp(ans0,'abs(Qs)') 
        %             txt2 = sprintf('Absolute; Calms <%s m^3/yr; Threshold: %0.4f; Interval: %dd; Min duration: %dd',...
        %                         varsel.calms.text,options.threshold,...
        %                         options.clint,options.mincluster);
        %         else
        %             txt2 = sprintf('Pos/Neg; Calms <%s m^3/yr; Threshold: %0.4f/%0.4f;\n          Interval: %dd/%dd; Min duration: %dd/%dd',...
        %                         varsel.calms.text,...
        %                         options.pos.threshold,options.neg.threshold,...
        %                         options.pos.clint, options.neg.clint,...
        %                         options.pos.mincluster, options.neg.mincluster);
        %         end
        %         subtitle(axc,sprintf('%s\n%s',txt1,txt2))   
        %     case 'Reach summary'            
        %         varsel.subtitle = subtxt('Peclet ratio cluster mean',varsel.calms.text,'');
        %         reachPecletPlot(npnts,bintime,zq,varsel);
        % end
    end
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
            if isempty(rose.theta)
                %title using variable-case-point
                titletxt = sprintf('%s for %s at %s',varsel.desc,varsel.case,pntnames{ipnt});
                theta = [];
            else
                %title using case-point-shore_angle
                titletxt = sprintf('%s at %s, theta=%d dTN',varsek.case,pntnames{ipnt},rose.theta(i));
                theta = rose.theta(i);
            end
            wind_rose(dir,var,'parent',figax,'dtype','meteo',...
                'shore',theta,'nd',rose.nd,'di',rose.di,'ci',rose.ci,...
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
function [cluster,ints,userops] = absClusters(options,dst,varsel)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    % varsel = getVariable(dst,pntnames,1);  %selects Qs without prompting user
    % if isempty(varsel), return; end
    mtime = dst.(pntnames{1}).RowNames;

    ans2 = questdlg('Check settings for selected points?','Clusters','Yes','No','Quit','Yes');
    if strcmp(ans2,'Yes')
        promptxt = {'Accept figures are used to adjust the cluster selection';...
                    'The plot shows absolute drift values (all +ve)';...
                    'Use Threshold and Time between clusters to adjust number of clusters'};
        getdialog(promptxt,[],3);
    
        ok = 1;
        while ok>0
            [vardst,varsel] = getPointData(dst,varsel);
               % [sel,ok] = listdlg('Name','Plot profile', ...
               %                       'PromptString','Select variable', ...
               %                       'ListSize',[200,300], ...
               %                       'SelectionMode','single', ...
               %                       'ListString',pntnames);  
               %  if ok==0, continue; end 
               % 
               %  Var = dst.(pntnames{sel}).(varsel.name); 
               %  varsel.point = sel;    %add point to variable selection struct
               %  vardst = getDSTable(dst.(pntnames{sel}),'VariableNames',varsel.name);

                vardst.(varsel.name) = abs(vardst.(varsel.name));                
                promptxt = sprintf('Cluster definition for %s at Point %d',...
                                       vardst.VariableNames{1},varsel.point);
                [idcls,userops] = getclusters(vardst,options,promptxt);
                userops.mincluster = options.mincluster;
                userops.isplot = options.isplot;                     
                [medges,isok] = mergeAbsClusters(vardst.(varsel.name),mtime,idcls,userops,varsel);
                if isempty(medges)
                    getdialog('No clusters found. Change threshold or minimum duration of cluster')
                elseif isok
                    ok = 0;
                end
        end
    
        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');        
    elseif  strcmp(ans2,'No')
        userops = options;
    end

    if strcmp(ans2,'Quit'), cluster = []; userops = options; return; end

    userops.isplot = false; %supress plots in for loop

    % cluster(npnts,1) = struct('mean',[],'stdev',[],'ptle95',[],'ptle05',[],...
    %                'Pe',[],'dirBias',[],'netQ',[],'skewness',[],...
    %                'classification','','index',[],'nPoints',[]);
    hw = waitbar(0,'Processing point 0');
    for i=1:npnts        
        Var = dst.(pntnames{i}).(varsel.name);        
        Var(abs(Var)<varsel.calms.value) = NaN; %remove near zero values
        vardst = getDSTable(dst.(pntnames{i}),'VariableNames',varsel.name);
        vardst.(varsel.name) = abs(vardst.(varsel.name)); %use absolute values for intervals

        % find clusters based on results from peak selection
        idcls = getVarClusters(vardst,userops);        
        %merge any overlaps to define intervals to be used
        medges = mergeAbsClusters(Var,mtime,idcls,userops);
        if isempty(medges)
            medges = [mtime(1),mtime(end)];
        end
        
        %find the indices of the variable within each interval
        [intervals,intstart] = discretize(mtime,medges);
        nint = length(intstart);

        %compute the statistics for the point over each interval 
        ints{i} = intstart;
        for j=1:nint
            idint = intervals==j;
            cluster(i,j) = wrm_interval_statistics(Var(idint)); %#ok<AGROW>            
        end

                % for j=1:nint
                %     idint = intervals==j;
            	% 	binvar = Var(idint);
                %     if isempty(binvar)
                %         meanVar = NaN; stdVar = NaN; peclet = NaN;
                %     else
                %         meanVar = mean(binvar,'omitnan');
                %         stdVar = std(binvar,'omitnan');
                %         peclet = meanVar./stdVar;
                %         peclet = checkPecletLimits(peclet,varsel,NaN,peclet);
                %     end
                % 
                %     cluster.Ints{i,j} = intstart(j); %#ok<*AGROW>
                %     cluster.Mean{i,j} = meanVar;
                %     cluster.Stdev{i,j} = stdVar;
                %     cluster.Peclet{i,j} = peclet;
                % end
        waitbar(i/npnts,hw,sprintf('Processing point %d',i));
    end
    delete(hw)
end

%%
function [cluster,ints,userops]  = posnegClusters(options,dst,varsel)
    %select varaiable and get time data
    pntnames = fieldnames(dst);
    npnts = length(pntnames);
    mtime = dst.(pntnames{1}).RowNames;

    %default to use same options for postive and negative drift
    userops.pos = options; userops.neg = options;     

    %set a minimum trheold to remove zero values
    ans2 = questdlg('Check settings for selected points?','Clusters','Yes','No','Quit','Yes');
    if strcmp(ans2,'Quit')
        cluster = []; ints = []; return; 
    elseif strcmp(ans2,'Yes')
        promptxt = {'Accept figures are used to adjust the threshold selection';...
                    'The first plot sets the positive threshold';...
                    'The second plot sets the negative threshold';...
                    '(NB: variable is inverted in the 2nd plot)'};
        getdialog(promptxt,[],3);

        ok = 1;
        while ok>0
            [posdst,varsel] = getPointData(dst,varsel);
            % [sel,ok] = listdlg('Name','Plot profile', ...
            %                       'PromptString','Select variable', ...
            %                       'ListSize',[200,300], ...
            %                       'SelectionMode','single', ...
            %                       'ListString',pntnames);
            %  if ok==0, continue; end
            %
            %  Var = dst.(pntnames{sel}).(varsel.name);
            %  varsel.point = sel;    %add point to variable selection struct
            % posdst = getDSTable(dst.(pntnames{sel}),'VariableNames',varsel.name);
            promptxt = sprintf('Positive %s cluster definition for Point %d',...
                posdst.VariableNames{1},varsel.point);
            [idpos,userops.pos] = getclusters(posdst,options,promptxt);
            userops.pos.mincluster = options.mincluster;
            userops.pos.isplot = options.isplot;

            negdst = copy(posdst);
            negdst.(varsel.name) = negdst.(varsel.name)*-1;
            promptxt = sprintf('Negative %s cluster definition for Point %d',...
                negdst.VariableNames{1},varsel.point);
            [idneg,userops.neg] = getclusters(negdst,options,promptxt);
            userops.neg.mincluster = options.mincluster;
            userops.neg.isplot = options.isplot;
            [medges,isok] = mergePosNegClusters(posdst.(varsel.name),mtime,idpos,idneg,userops,varsel);
            if isempty(medges)
                getdialog('No clusters found. Change threshold or minimum duration of cluster')
            elseif isok
                ok = 0;
            end
        end

        ans2 = questdlg('Proceed with analysis of all points using last set of options?',...
                                         'Clusters','Proceed','Quit','Proceed');
    elseif strcmp(ans2,'No')
        %use same options for postive and negative drift
    end

    if strcmp(ans2,'Quit'), cluster = []; ints = []; return; end

    userops.pos.isplot = false; userops.neg.isplot = false;
    hw = waitbar(0,'Processing point 0');
    for i=1:npnts        
        Var = dst.(pntnames{i}).(varsel.name);        
        Var(abs(Var)<varsel.calms.value) = NaN; %remove near zero values
        posdst = getDSTable(dst.(pntnames{i}),'VariableNames',varsel.name);

        % find clusters based on results from peak selection
        idposcls = getVarClusters(posdst,userops.pos);
        negdst = posdst;
        negdst.(varsel.name) = negdst.(varsel.name)*-1;
        idnegcls = getVarClusters(negdst,userops.neg);

        %merge any overlaps to define intervals to be used
        userops.pos.isplot = false;
        medges = mergePosNegClusters(Var,mtime,idposcls,idnegcls,userops,varsel); %only uses mincluster field in options
        if isempty(medges)
            medges = [mtime(1),mtime(end)];       
        end
        
        %find the indices of the variable within each interval
        [intervals,intstart] = discretize(mtime,medges);
        nint = length(intstart);
        ints{i} = intstart; %#ok<AGROW>
        %compute the statistics for the point over each interval 
        for j=1:nint
            idint = intervals==j;
            cluster(i,j) = wrm_interval_statistics(Var(idint));
        end
        waitbar(i/npnts,hw,sprintf('Processing point %d',i));
                % for j=1:nint
                %     idint = intervals==j;
                %     binvar = Var(idint);
                %     meanVar = mean(binvar,'omitnan');
                %     stdVar = std(binvar,'omitnan');
                %     peclet = meanVar./stdVar;
                % 
                %     %set diffusion values to 0
                %     peclet = checkPecletLimits(peclet,varsel,NaN,peclet);
                %     %peclet = checkPecletLimits(peclet,varsel,0,0);
                % 
                %     cluster.Ints{i,j} = intstart(j); %#ok<*AGROW>
                %     cluster.Mean{i,j} = meanVar;
                %     cluster.Stdev{i,j} = stdVar;
                %     cluster.Peclet{i,j} = peclet;
                % end
    end
    delete(hw)
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
    pk_vals = var(idpks);      %value of peak
    idcls = clusters(pk_date,pk_vals,days(opts.clint));
end

%%
function [medges,isok] = mergeAbsClusters(var,mtime,idpos,opts,varsel)
    %merge  absolute cluster selections to a single set of edges
    mincluster = opts.mincluster*24;      %min length of a cluster (h)
    dt = mode(diff(mtime));
    mincls = floor(mincluster/hours(dt));  
    func = @(x) length(x)<mincls;
    postimes = {idpos(:).date};
    posshort = cellfun(func,postimes,"UniformOutput",false);
    postimes([posshort{:}]) = [];

    mdates.posstart = cellfun(@(x) x(1),postimes); %first date in each cell
    mdates.posend = cellfun(@(x) x(end),postimes); %last date in each cell
    medges = sort(unique([mdates.posstart,mdates.posend])); %sorted edges

    if ~isempty(medges) && opts.isplot
        %plotEdges(var,mtime,medges,opts,'Intervals to used for statistics');
        varsel.title = 'Intervals to use for statistics';
        isok = get_plotEdges(var,mtime,medges,opts,varsel);
    else
        isok = true;
    end
end

%%
function [medges,isok] = mergePosNegClusters(var,mtime,idpos,idneg,opts,varsel)
    %merge positive and negative cluster selections to a single set of edges
    dt = mode(diff(mtime));  
    mincluster = opts.pos.mincluster*24;      %min length of a cluster (h)
    mincls = floor(mincluster/hours(dt));  
    func = @(x) length(x)<mincls;
    postimes = {idpos(:).date};
    posshort = cellfun(func,postimes,"UniformOutput",false);
    postimes([posshort{:}]) = [];

    mincluster = opts.neg.mincluster*24;      %min length of a cluster (h)
    mincls = floor(mincluster/hours(dt));  
    func = @(x) length(x)<mincls;    
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

    if ~isempty(medges) && opts.pos.isplot
        varsel.title = 'Intervals to use for statistics';
        isok = get_plotEdges(var,mtime,medges,opts,varsel);
    else
        isok = true;
    end
end

%%
function mdates = mergeOverlaps(var,mtime,mdates,opts)
    %find any overlaps and merge any that are short
    pos.start = mdates.posstart;
    pos.end = mdates.posend;
    neg.start = mdates.negstart;
    neg.end = mdates.negend;

    overlaps = findOverlaps(pos,neg);
    if any(overlaps,'all')        
        [row, col] = find(overlaps);
        %plotMergedVar(var,mtime,posstart(row),posend(row),negstart(col),negend(col),'Unmerged cluster overlaps');
        poso.start = pos.start(row); poso.end = pos.end(row);
        nego.start = neg.start(col); nego.end = neg.end(col);
        if opts.pos.isplot
            plotMergedVar(var,mtime,poso,nego,opts,'Unmerged cluster overlaps');
        end

        for k = 1:numel(row)
            pn = [pos.start(row(k)), neg.start(col(k)), pos.end(row(k)), neg.end(col(k))];
            pos.start(row(k)) = min(pn); %assign to pos and remove neg
            neg.start(col(k)) = NaT;      %makes no difference as pos&neg going
            pos.end(row(k)) = max(pn);   %be merged
            neg.end(col(k)) = NaT;
        end
        neg.start(isnat(neg.start)) = [];
        neg.end(isnat(neg.end)) = [];
        % for k = 1:length(row)
        %     overlap_start = max(pos.start(row(k)), neg.start(col(k)));
        %     overlap_end = min(pos.end(row(k)), neg.end(col(k)));
        %     overlap_length = hours(overlap_end-overlap_start);
        %     mincluster = (opts.pos.mincluster+opts.neg.mincluster)/2;
        %     if overlap_length<mincluster*24    %split between the two
        %         if pos.start(row(k))<neg.start(col(k)) && pos.end(row(k))<neg.end(col(k))
        %             pos.end(row(k)) = pos.end(row(k))-hours(overlap_length/2+0.0);
        %             neg.start(col(k)) = neg.start(col(k))+hours(overlap_length/2+0.0);
        %         elseif pos.start(row(k))>neg.start(col(k)) && pos.end(row(k))>neg.end(col(k))
        %             pos.start(row(k)) = pos.start(row(k))+hours(overlap_length/2+0.0);
        %             neg.end(col(k)) = neg.end(col(k))-hours(overlap_length/2+0.0);
        %         end
        %     else
        %         overlap_start = min(pos.start(row(k)), neg.start(col(k)));
        %         overlap_end = max(pos.end(row(k)), neg.end(col(k)));
        %         pos = removeOverlap(pos,overlap_start,overlap_end,row(k));
        %         neg = removeOverlap(neg,overlap_start,overlap_end,col(k));
        %         fprintf('Overlaps from %d: %s to %s merged\n', k, ...
        %                                        overlap_start, overlap_end);
        %     end            
        % end
        fprintf('%d overlaps have been merged\n', numel(row));
        %check that all have been removed
        % overlaps = findOverlaps(pos,neg);
        % [row, col] = find(overlaps);
        % if ~isempty(row)
        %     if opts.pos.isplot
        %         poso.start = pos.start(row); poso.end = pos.end(row);
        %         nego.start = neg.start(col); nego.end = neg.end(col);
        %         plotMergedVar(var,mtime,poso,nego,opts,'Cluster overlaps to be subdivided');
        %     end            
        % end
        %update struct with merged intervals
        mdates = struct('posstart',pos.start,'posend',pos.end,'negstart',neg.start,'negend',neg.end);
    end

    %-nested function------------------------------------------------------
    function overlaps = findOverlaps(pos,neg)
        % Initialize a logical matrix to store overlaps
        pstart = numel(pos.start);
        nstart = numel(neg.start);
        overlaps = false(pstart,nstart);
        
        % Check for overlaps between intervals
        for i = 1:pstart
            for j = 1:nstart
                overlaps(i, j) = (pos.start(i) <= neg.end(j)) && (neg.start(j) <= pos.end(i));
            end
        end
    end

    % %-nested function------------------------------------------------------
    % function var = removeOverlap(var,lapstart,lapend,idv)
    %     if isbetween(var.start(idv),lapstart,lapend,'open')
    %         var.start(idv) = [];
    %     end
    %     if isbetween(var.end(idv),lapstart,lapend,'open')
    %         var.end(idv) = [];
    %     end
    % end
end

%% ------------------------------------------------------------------------
% Utility functions for variable selection
%--------------------------------------------------------------------------
function [calms,pecthr] = calmsThreshold(varunits)
    %set the calms threshold to apply to the data
    calmsthreshold = 10;  %"calms" are drift rates less than threshold
                          % 10m^3/yr ~= 3e-7 m^3/s;
                          % or wave height or runup less than threshold
    promptxt = {'Calms threshold (Hs(m), Qs(m^3/yr), etc):','Peclet plotting threshold'};
    defaults = {num2str(calmsthreshold),'1'};
    answer = inputdlg(promptxt,'Drift',1,defaults);
    if isempty(answer), answer = defaults; end
    if contains(varunits,'m^3/s')   %input is in m^3/yr
        calms.value = str2double(answer{1})/31556952;  %value from mobj.Constants.y2s
    end
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
    
    varsel.case = dst.(pntnames{1}).Description;
    varsel.name = varname{sel};
    varsel.desc = vardesc{sel};
    varsel.labl = varlabl{sel};

    varsel.isdrift = false;
    if contains(varsel.name,'Q'), varsel.isdrift = true; end
    % 
    if nargin<4 || isthr
        varunits = dst.(pntnames{1}).VariableUnits{sel};
        [varsel.calms,varsel.pecthr] = calmsThreshold(varunits);
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
function options = setClusterOptions(data,opts)
    %define the options used in a peaks and cluster data selection
    if nargin<2 || isempty(opts)
        threshold = round(mean(data,'omitnan'),1,"significant");
        default = {num2str(threshold),'1','0','15','5'};                   
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
function [vardst,varsel] = getPointData(dst,varsel)
    %
    pntnames = fieldnames(dst);
   [sel,ok] = listdlg('Name','Plot profile', ...
                         'PromptString','Select variable', ...
                         'ListSize',[200,300], ...
                         'SelectionMode','single', ...
                         'ListString',pntnames);  
    if ok==0, return; end 

    varsel.point = sel;    %add point to variable selection struct
    vardst = getDSTable(dst.(pntnames{sel}),'VariableNames',varsel.name);
end

%%
function [nrch,stpnts] = setReachPoints(ndpnt)
    %UI to set definition of reaches for summary peclet plot    
    promptxt = {'Number of reaches',sprintf('ID of start points (<%d)',ndpnt)};
    defaults = {'1','1'};
     ok = 0;
     while ok<1   
        inp = inputdlg(promptxt,'Reaches',1,defaults);
        if isempty(inp), return; end
        nrch = str2double(inp{1});
        stpnts = str2num(inp{2}); %#ok<ST2NM> vector input
        if numel(stpnts)==nrch
            ok = 1;
        else
            warndlg('Number of start points must match number of reaches')
        end
    end
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

%%
function pecvar = checkPecletLimits(pecvar,varsel,nullvalue,diffvalue)
    %check peclet value and set to nullvalue if test 1 fails and to
    %diffvalue if test 2 fails
    [m,n] = size(pecvar);
    for i=1:m
        for j=1:n
            peclet = pecvar(i,j);
            if isnan(peclet) || isinf(peclet)
                pecvar(i,j) = nullvalue;
            elseif peclet>-varsel.pecthr && peclet<varsel.pecthr
                if isequal(size(pecvar), size(diffvalue))
                    pecvar(i,j) = diffvalue(i,j);
                elseif isscalar(diffvalue)
                    pecvar(i,j) = diffvalue;
                else
                    errordlg('diffvalue must be scalar or same size as pecvar in checkPecletLimtis')
                end
            end
        end
    end
end

%%
function [plotVar,meta] = getPlotVariable(stats,varsel)
    %select variable to plot and extract from stats struct
    %varoptions are based on the varnames used in the stats struct
    varnames = fieldnames(stats);
    varoptions = {'Mean','Standard deviation','95th Percentile',...
                  '5th Percentile','Peclet ratio','Directional bias index',...
                  'Net transport index','Skewness','Classification',...
                  'Classification index','No. of points'};  
    varlabels = {'mean','st.dev','5/95 percentiles','5/95 percentiles',...
                 'Peclet ratio','bias index','net transport index',...
                 'skewness','classification','class index','no. points'};

    selection = listdlg('Name','Variable options', ...
                        'PromptString','Select variable:','ListSize',[200,150],... ...
                        'SelectionMode','single','ListString',varoptions);
    if isempty(selection), plotVar = []; meta = ''; return; end

    %force use of index for plotting Classification
    if selection==9, selection = 10; end

    var2use = varnames{selection};

    if varsel.iscluster
        %variable is a different length for each point. pad with NaNs to
        %allow variable specific changes
        sz = num2cell(size(stats));
        plotVar = NaN(sz{:}); varlength = zeros(1,sz{1});
        if strcmp(var2use,'ptle95') ||  strcmp(var2use,'ptle05')
            addVar = plotVar;
            for i=1:sz{1}
                varlength(i) = numel([stats(i,:).ptle05]);
                plotVar(i,1:varlength(i)) = [stats(i,:).ptle05];
                addVar(i,1:varlength(i)) = [stats(i,:).ptle95];
            end
        else
            for i=1:sz{1}
                varlength(i) = numel([stats(i,:).(var2use)]);
                plotVar(i,1:varlength(i)) = [stats(i,:).(var2use)];
            end
        end
    else
        sz = num2cell(size(stats));        
        if strcmp(var2use,'ptle95') || strcmp(var2use,'ptle05')
            plotVar = reshape([stats(:).ptle05],sz{:});
            addVar = reshape([stats(:).ptle95],sz{:});
        else
            plotVar = reshape([stats(:).(var2use)],sz{:});
        end
    end

    if strcmp(var2use,'Pe')
         plotVar = checkPecletLimits(plotVar,varsel,NaN,plotVar);
    elseif strcmp(var2use,'ptle95') || strcmp(var2use,'ptle05') 
        %remove small values and merge 95 (positive) and 05 (negative) estimates
        min0 = mean(addVar,'all','omitnan')-std(plotVar,0,'all','omitnan');
        addVar(abs(addVar)<min0) = 0;  %remove very small values
        min0 = mean(plotVar,'all','omitnan')-std(plotVar,0,'all','omitnan');
        plotVar(abs(plotVar)<min0) = 0;           
        plotVar = plotVar+addVar; %assumes no overlap
        plotVar(plotVar==0) = NaN;  
    end

    meta = varsel; 
    meta.statname = varnames{selection};
    meta.statdesc = varoptions{selection};
    meta.statlabl = varlabels{selection};    

    if varsel.iscluster
        %repack variable as a variable length vector in a cell for each point
        temp = plotVar; 
        plotVar = cell(1,sz{1});
        for i=1:sz{1}
            plotVar{i} = [temp(i,1:varlength(i))];
        end
        meta.statxt = sprintf('%s %s',meta.statlabl,lower(meta.labl));
    else
        meta.statxt = sprintf('%sly %s %s',meta.bintxt,meta.statlabl,lower(meta.labl));
    end
end
%% ------------------------------------------------------------------------
% Utility functions for plotting 
%--------------------------------------------------------------------------

function [ax,hs] = plot2Dvariable(var,bintime,desc,options)
    %plot a 2D variable as a surface plot or scatter points (position,time)
    %options: 1 - surface plot; 2 - points with colormap; 3 - points pos/neg
    npnts = size(var,1);
    hf = figure('Name','SedTrans','Tag','PlotFig');
    ax = axes(hf);    
    grid on
    [X,Y] = meshgrid(1:npnts,bintime); 
    if options==1
        hs = surf(ax,X,Y,var');
        shading interp
        ax.Layer = 'top';
        setColormap(var(:),0)
        hc = colorbar;
        hc.Label.String = desc.statxt;  

    elseif options==2
        Xv = X(:);
        Yv = Y(:);    
        Z = var'; Zv = Z(:);
        hs = scatter3(ax,Xv,Yv,Zv,6,Zv,'filled','Marker','square','MarkerEdgeColor','none');                
        setColormap(Zv,1);
        hc = colorbar;
        hc.Label.String = desc.statxt;   

    else
        Xv = X(:);
        Yv = Y(:);    
        Z = var'; Zv = Z(:);
        posZ = Zv>=0; negZ = Zv<0;
        blue = [0,0.447,0.741];
        hs = scatter3(ax,Xv(posZ),Yv(posZ),Zv(posZ),6,blue,'filled','Marker','square','MarkerEdgeColor','none');                
        hold on
        yellow = [0.929,0.694,0.125];
        scatter3(ax,Xv(negZ),Yv(negZ),Zv(negZ),6,yellow,'filled','Marker','square','MarkerEdgeColor','none')
        hold off
    end
    
    view(2)
    axis tight
    datetick('y', 'yyyy'); %#ok<DATIC>
    xlabel('Position along shore')
    ylabel('Year')
    title(sprintf('Case: %s',desc.case));  
end

%%
function plotMergedVar(var,mtime,pos,neg,opts,titxt)
    %plot the merged selection
    hf = figure('Name','SedTrans','Tag','StatFig');
    ax = axes(hf); 
    plot(ax,mtime,var,'Color',[0.75,0.75,0.75],'LineWidth',0.2)
    yy = ylim;
    posy = [0,yy(2)];
    negy = [0,yy(1)];
    hold on
    plot(ax,ax.XLim,[1,1]*opts.pos.threshold,'Color',[0.7,0.7,0.7])
    plot(ax,ax.XLim,[-1,-1]*opts.neg.threshold,'Color',[0.7,0.7,0.7])
    if numel(pos.start)==2 
         %needed if there are only 2 points to avoid plotting diagonal
        for i=1:2
            plot([pos.start(i), pos.start(i)],posy,'-','Color',"#77AC30",'LineWidth',1);%#7E2F8E
            plot([pos.end(i),pos.end(i)],posy,'--','Color',"#77AC30",'LineWidth',1);
            plot([neg.start(i),neg.start(i)],negy,'-','Color','#A2142F','LineWidth',1);
            plot([neg.end(i),neg.end(i)],negy,'--','Color','#A2142F','LineWidth',1); 
        end
    else
        plot([pos.start', pos.start'],posy,'-','Color',"#77AC30",'LineWidth',1);%#7E2F8E
        plot([pos.end',pos.end'],posy,'--','Color',"#77AC30",'LineWidth',1);
        plot([neg.start',neg.start'],negy,'-','Color','#A2142F','LineWidth',1);
        plot([neg.end',neg.end'],negy,'--','Color','#A2142F','LineWidth',1);  
    end
    hold off
    xlabel('Time')
    ylabel('Selected drift variable')
    title(titxt)
end

%%
function isaccept = get_plotEdges(var,mtime,medges,opts,varsel)
    %plot figure of edges with yes/no accept buttons
    promptxt = 'Accept edges definition';
    [h_plt,h_but] = acceptfigure(varsel.title,promptxt,'StatFig');%StatFig is Tag to allow group delete
    h_ax = axes(h_plt);
    [~,ncls] = plotEdges(var,mtime,medges,opts,h_ax);

    waitfor(h_but,'Tag');
    if ~ishandle(h_but)   %this handles the user deleting figure window
        isaccept = []; 
    elseif strcmp(h_but.Tag,'Yes')
        isaccept = true;
        panelText(h_but,ncls,varsel,opts);
    else
        isaccept = false;
    end  

    function panelText(h_but,ncls,varsel,opts)
        %set text to display selected parameters
        hbyn = findobj(h_but,'Tag','YesNo');
        delete(hbyn)
        if isfield(opts,'pos')
            % npoint = opts.pos.point;
            threshold = sprintf('%.5f / %.5f',opts.pos.threshold,opts.neg.threshold);
            method = sprintf('%d / %d',opts.pos.method,opts.neg.method);
            tint = sprintf('%g / %g',opts.pos.tint,opts.neg.tint);
            clint = sprintf('%g / %g',opts.pos.clint,opts.neg.clint);
            mdur = sprintf('%g / %g',opts.pos.mincluster,opts.neg.mincluster);
        else
            % npoint = opts.point;
            threshold = sprintf('%.5f',opts.threshold);
            method = sprintf('%g',opts.method);
            tint = sprintf('%g',opts.tint);
            clint = sprintf('%g',opts.clint);
            mdur = sprintf('%g',opts.mincluster);
        end
        h_but.Title = sprintf('Selected parameters for Point %d: %s',varsel.point,varsel.case);
        txt1 = sprintf('Threshold for peaks = %s',threshold);
        txt2 = sprintf('Selection method = %s',method);
        txt3 = sprintf('Peak separation = %s hrs',tint);
        txt4 = sprintf('Cluster separation = %s days',clint);
        txt5 = sprintf('Min. duration = %s days',mdur);
        txt6 = sprintf('             Clusters: No. = %d;  No./year = %.2f;  Av.dur. = %.1f days;  Percent obs period = %.1f %%',...
                       ncls.count,ncls.rate,ncls.avdur,ncls.pcntdur);
        selparams = sprintf('%s; %s; %s; %s; %s\n%s',txt1,txt2,txt3,txt4,txt5,txt6);
        uicontrol('Parent',h_but,'Tag','YesNo',...
            'Style', 'text', 'String', selparams,...
            'Units','normalized', ...
            'Position', [0.01 0.01 0.99 0.8]);
        %write output to command window for collation
        fprintf('Case, Point, Threshold, Separation, Min.Dur, No., No/yr, Av.dur., %% time\n')
        fprintf('%s, %d, %s, %s, %s, %d, %.2f, %.1f, %.1f\n',varsel.case,...
                            varsel.point,threshold,clint,mdur,...
                            ncls.count,ncls.rate,ncls.avdur,ncls.pcntdur)
    end
end

%%
function [ax,numcls] = plotEdges(var,mtime,medges,opts,ax)
    %plot the merged edges to be used to compute the statitics
    %ie the start or end of each cluster
    if nargin<5   
        hf = figure('Name','SedTrans','Tag','StatFig');
        ax = axes(hf); 
    end
    mvar = max(abs(var),[],'omitnan')/4;
    if isfield(opts,'pos')
        posthreshold = opts.pos.threshold;
        negthreshold = opts.neg.threshold;
    else
        posthreshold = opts.threshold;
        negthreshold = opts.threshold;
    end
    posvar = var; posvar(var<posthreshold) = NaN;
    negvar = var; negvar(var>-negthreshold) = NaN;
    hold(ax,'on')
        plot(ax,mtime,posvar,'Color',[0.75,0.75,0.75],'LineWidth',0.2)     %positive cluster variable
        plot(ax,mtime,negvar,'Color',[0.75,0.75,0.75],'LineWidth',0.2)     %negative cluster variable
        plot(ax,[medges(1:2:end)',medges(1:2:end)'],[-mvar,mvar],'-','Color',"#77AC30")  %cluster edges
        plot(ax,[medges(2:2:end)',medges(2:2:end)'],[-mvar,mvar],'-','Color',"#D95319")  %cluster edges
        plot(ax,medges,0,'.','Color',"#0072BD",'MarkerSize',4)
    hold(ax,'off')
    xlabel('Time')
    ylabel('Selected drift variable')
    
    obsduration = mtime(end)-mtime(1);
    cluster_durations = medges(2:2:end)-medges(1:2:end);
    numcls.count = numel(medges(1:2:end));
    numcls.rate = numcls.count/years(obsduration);
    numcls.avdur = days(mean(cluster_durations,'omitnan'));
    numcls.pcntdur = sum(cluster_durations,'omitnan')/obsduration*100;
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
    yellow = [0.929,0.694,0.125];
    scatter3(sax,X,Y,-point.down*2,25,yellow,'LineWidth',1)
    scatter3(sax,X,Y,point.up*2,25,'b','LineWidth',1)
    hold off
end

%%
function reachPlot(npnts,bintime,var,varsel,sample)
    %plot reach mean values of variable for selected reaches 
    [nrch,stpnts] = setReachPoints(npnts);
    ndpnts = [stpnts(2:end)-1,npnts];  
    ntime = numel(bintime);
    rchmean = zeros(ntime,nrch);
    for i=1:ntime
        for j=1:nrch
            rchmean(i,j) = mean(var(i,stpnts(j):ndpnts(j)),'omitnan');
            rchstd(i,j) = std(var(i,stpnts(j):ndpnts(j)),'omitnan');
        end
    end
    %plot results
    hfig = figure('Tag','PlotFig');
    axs = axes(hfig);
    colororder(axs,'gem12') %"gem12" is extended version of default palette
    colorlist = axs.ColorOrder;
    axs.ColorOrder = colorlist(1:nrch,:);
    hold(axs,'on')
    for j=1:nrch
        rmean = mean(rchmean(:,j),'omitnan');
        rmeanstd = mean(rchstd(:,j),'omitnan');
        rname = sprintf('Reach %d, mean=%.0f',j,rmean);
        hp = plot(axs,bintime,rchmean(:,j),'DisplayName',rname);
        hp.SeriesIndex = j;
        %add mean and std lines
        hp = plot(axs,axs.XLim,[1,1]*rmean,'--');
        hp.SeriesIndex = j;
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

        hp = plot(axs,axs.XLim,[1,1]*(rmean+rmeanstd),'-.');
        hp.SeriesIndex = j;
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

        hp = plot(axs,axs.XLim,[1,1]*(rmean-rmeanstd),'-.');
        hp.SeriesIndex = j;
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    hold(axs,'off')
    xlabel('Time')
    ylbl = varsel.labl;
    if strcmp(sample.func,'sum')
        ylbl = sprintf('Total transport/%s (m^3)',sample.binsize); 
    end
    ylabel(ylbl)
    title(sprintf('Case: %s',varsel.case))
    subtitle(varsel.subtitle)
    legend
end

%%
function reachPecletPlot(npnts,bintime,intervalPeclet,varsel)
    %plot summary mean values of Peclets>1 for selected reaches 
    [nrch,stpnts] = setReachPoints(npnts);
    ndpnts = [stpnts(2:end)-1,npnts];  
    ntime = numel(bintime);
    rchmean = zeros(ntime,nrch); rchpecp = NaN(ntime,nrch); rchpecn = rchpecp;
    for i=1:ntime
        for j=1:nrch
            rchmean(i,j) = mean(intervalPeclet(i,stpnts(j):ndpnts(j)),'omitnan');
            if any(intervalPeclet(i,stpnts(j):ndpnts(j))>varsel.pecthr)
                rchpecp(i,j) = nrch-j+1; 
            end
            if any(intervalPeclet(i,stpnts(j):ndpnts(j))<-varsel.pecthr)
                rchpecn(i,j) = -j;%-(nrch-j+1); 
            end
        end
    end
    %plot results
    hfig = figure('Tag','PlotFig');
    axs = axes(hfig);
    hold(axs,'on')
    for j=1:nrch
        plot(axs,bintime,rchmean(:,j),'DisplayName',sprintf('Reach %d',j))
    end
    ypec = [1,1]*varsel.pecthr;
    hp = plot(axs,axs.XLim,ypec,'--','Color',[0.75,0.75,0.75]);
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp = plot(axs,axs.XLim,-ypec,'--','Color',[0.75,0.75,0.75]);
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hp = plot(axs,axs.XLim,[0,0],'--','Color',[0.5,0.5,0.5]);
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold(axs,'off')
    xlabel('Time')
    ylabel('Peclet ratio')
    title(sprintf('Case: %s',varsel.case))
    subtitle(varsel.subtitle)
    legend

    hfig = figure('Tag','PlotFig');
    axb = axes(hfig);
    colorlist = axb.ColorOrder;
    axb.ColorOrder = colorlist(1:nrch,:);
    axb.YTickLabel = string([-fliplr(1:nrch),0,fliplr(1:nrch)]');
    axb.YTick = -nrch:1:nrch;

    hold(axb,'on')
    for j=1:nrch
        hs1 = stem(axb,bintime,rchpecp(:,j),'Linewidth',1,...
                               'Marker','.','DisplayName',sprintf('Reach %d',j));
        hs1.SeriesIndex = j;
        jj = nrch-j+1; %reverse plotting order
        hs2 = stem(axb,bintime,rchpecn(:,jj),'LineWidth',1,'Marker','.');
        %hs2.SeriesIndex =  hs1.SeriesIndex;  %force the same color
        hs2.SeriesIndex = jj;
        hs2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    hold(axb,'off')
    xlabel('Time')
    ylabel('Reach number & Peclet event -1>P>1')
    title(sprintf('Case: %s',varsel.case))
    varsel.subtitle = replace(varsel.subtitle,'mean','events');
    subtitle(varsel.subtitle)
    legend
end

%%
function setAxisLimits(ax,npnts,std,nnd)
    %format the axes limits to fixed ranges 
    ax.XLim = [1,npnts];
    ax.YLim = [datenum(['01-01-',std]),datenum(['12-31-',nnd])]; %#ok<DATNM>
end

%%
function setColormap(Z,islimits)
    %select the color map to use. default is 19-BuGnYl, for fixed intervals
    %use 25 and set clim to +/-max(minmax(Z))
    mnmx = minmax(Z);
    limits = max(mnmx);
    promptxt = {'Colormap selection (25 for intervals)','Limits'};
    defaults = {'19',num2str([-limits,limits])};
    inp = inputdlg(promptxt,'Colormap',1,defaults);
    if isempty(inp)
        mapidx = 19; 
    else
        mapidx = str2double(inp{1});
        limits = str2num(inp{2}); %#ok<ST2NM>
    end

    if mapidx==24 || mapidx==23
        zoptions.Z = Z;
        zoptions.zeroLevel = 0;
    elseif mapidx==25
        zoptions = limits(2);
    else
        zoptions = [];
    end

    colormap(cmap_selection(mapidx,zoptions));
    if islimits
        clim(limits);
    end    
end