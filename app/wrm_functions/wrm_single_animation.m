function wrm_single_animation(mobj,sobj,plttxt,ptype)
%
%-------function help------------------------------------------------------
% NAME
%   wrm_single_animation.m
% PURPOSE
%   animation of model spectra timeseries
% USAGE
%   wrm_single_animation(mobj,tsdst,plttxt)
% INPUTS
%   mobj - instance of Model class (eg WaveRayModel, CoastalTools)
%   sobj - array of ctWaveSpectra class
%   plttxt - cell array of labels for variable, axes and title of plot
%   ptype - XY or Polar plot options
% OUTPUT
%   animation figure
% SEE ALSO
%   code derived from muiPlots.newAnimation and implemented for animation
% NOTES
%   obj is a data class instance and pobj is a plot class muiPlots instance
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2025
%--------------------------------------------------------------------------
%
    hfig = figure('Name','Animation', 'Units','normalized', ...
                    'Resize','on','HandleVisibility','on','Visible','on',...
                    'Tag','PlotFig'); %,'Position',[0.38,0.42,0.30,0.42]);
    %create an instance of muiPlots and populate the properties that are   
    %needed for the newAnimation method
    if isa(mobj.mUI.Plots,'muiPlots')
        pobj = mobj.mUI.Plots;           %get existing instance of muiPlot object       
        clearPreviousPlotData(pobj);
    else
        pobj = muiPlots.get_muiPlots();  %create new instance of muiPlot object  
    end       
    pobj.Plot.CurrentFig = hfig;
    pobj.Plot.FigNum = hfig.Number;
    pobj.ModelMovie = [];
    pobj.Title = plttxt{4};

    spectra = [sobj(:).Spectrum];        %unpack timeseries of spectra
    params = vertcat(sobj(:).Params);    %unpack integral properties
    %extract the timeseries data and dimensions for plot
    pobj.Data.X = 1./spectra(1).freq;
    pobj.Data.Y = spectra(1).dir;
    pobj.Data.T = [spectra(:).date];
    %pobj.Data.Z = [spectra(:).SG];
    nrec = length(pobj.Data.T);
    [n,m] = size(spectra(1).SG);
    Z = zeros(nrec,n,m);
    parfor i = 1:nrec
        Z(i, :, :) = reshape(spectra(i).SG, 1, n, m);
    end  
    pobj.Data.Z = Z;  clear Z
    pobj.Data.Waves = [params.Hs,params.Tp,params.Dir];

    ax = setupAnimation(sobj,pobj,plttxt,ptype);
    if ~isvalid(pobj.Plot.CurrentFig), return; end

    getAnimation(pobj,ax,hfig);
    ax.UserData = pobj.Data;  %store data set in UserData to
                              %allow user to switch between plots
    %add replay and slider
    hm = setControlPanel(pobj,hfig,length(pobj.Data.T),string(pobj.Data.T(1)));
    %assign bespoke run_Movie handle
    %hm(1).Callback = @(src,evt)run_Movie(pobj,src,evt);
    hm(3).Callback = @(src,evt)run_Movie(pobj,src,evt);

    %assign muiPlots instance to handle
    mobj.mUI.Plots = pobj;
end

%%
function ax = setupAnimation(sobj,pobj,plttxt,ptype)
    %initialise 3Dplot and setup animation variables
    hfig = pobj.Plot.CurrentFig;
    ax = axes(hfig);      
    ax = getPlot(sobj(1),plttxt,ptype,ax);
    t = pobj.Data.T;  
    w = pobj.Data.Waves;
    ax.Subtitle.String = sprintf('Time = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g',...
                       string(t(1)),w(1,1),w(1,2),w(1,3));
    if ~isvalid(hfig), return; end
    hfig.Visible = 'on';

    %assign axes properties                
    ax.ZLimMode = 'manual'; %fix limits of z-axis
    zlim = minmax([pobj.Data.Z(:)]); 
    ax.ZLim = [abs(zlim(1)),zlim(2)];
    ax.NextPlot = 'replaceChildren';
    ax.Tag = 'PlotFigAxes'; 
    ax.CLim = ax.ZLim;

    %assign data sources
    hp1 = findobj(ax.Children,'Type','surface');   
    hp1.CDataSource = 'var1'; 
    hd = findobj(ax,'Tag','DirPk');
    hd.YDataSource = 'var2';
    ht = findobj(ax,'Tag','TpPk');
    ht.XDataSource = 'var3';

    %adjust position of plots and add title
    if strcmp(ptype,'XY')                    %XY plot
        ax.Position([2,4]) = [0.20,0.68];      %make space for slider bar
    else                                     %polar plot
        % hfig.Position(3) = 0.6;
        So = pobj.Data.Z; 
        nrec = size(So,1);
        ff = pobj.Data.X; dd = pobj.Data.Y;
        parfor i=1:nrec
            Po(i,:,:) = reshapePolarGrid(ff,dd,So{i});
        end
        pobj.Data.Z = {Po}; 
    end
end

%%
function getAnimation(pobj,ax,hfig)
    %generate an animation for user selection.
    t = pobj.Data.T;  
    w = pobj.Data.Waves;
    nrec = length(t);
    Mframes(nrec) = struct('cdata',[],'colormap',[]);
    Mframes(1) = getframe(gcf); %NB print function allows more control of 
    hp1 = findobj(ax.Children,'Type','surface');
    hd = findobj(ax,'Tag','DirPk');
    ht = findobj(ax,'Tag','TpPk');
    for i=2:nrec
        var1 = squeeze(pobj.Data.Z(i,:,:)); %#ok<NASGU>        
        refreshdata(hp1,'caller')   
        var2 = [1,1]*w(i,3);  %#ok<NASGU>
        refreshdata(hd,'caller')  
        var3 = [1,1]*w(i,2);  %#ok<NASGU>
        refreshdata(ht,'caller')  
        ax.Subtitle.String = sprintf('Time = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g',...
                       string(t(i)),w(i,1),w(i,2),w(i,3));
        drawnow;                 
        Mframes(i) = getframe(gcf); 
        %NB print function allows more control of resolution 
    end
    idm = size(pobj.ModelMovie,1);            
    pobj.ModelMovie{idm+1,1} = hfig.Number;
    pobj.ModelMovie{idm+1,2} = Mframes;   %save movie to class property
end 

%%
function P = reshapePolarGrid(Period,Dir,S)
    %format spectral array to be in format required by polarplot3d
    wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
    %interpolate var(phi,T) onto plot domain defined by tints,rints
    %intervals match those set to initialise plot in SpectralTransfer.polar_plot
    tints = linspace(0,2*pi,360); 
    rints = linspace(1,25,25);
    [Tq,Rq] = meshgrid(tints,rints); 
    warning('off',wid)
    P = griddata(deg2rad(Dir),Period,S',Tq,Rq);
    P(isnan(P)) = 0;  %fill blank sector so that it plots the period labels
    warning('on',wid)
end

%%
function run_Movie(pobj,src,~)
    %callback function for animation figure buttons and slider
    %runMovie in muiPlot is used for re-running and saving movie
    %bespoke slider update required to update lines and surface
    hfig = src.Parent;
    t = pobj.Data.T;  
    w = pobj.Data.Waves;        
    val = ceil(src.Value);          %slider value
    %get figure axis, extract variable and refresh plot
    ax = findobj(hfig,'Tag','PlotFigAxes');
    hp = findobj(ax.Children,'Type','surface');
    hd = findobj(ax,'Tag','DirPk');
    ht = findobj(ax,'Tag','TpPk');

    time = ax.UserData.T(val);   %time slice selected
    var1 = ax.UserData.Z;

    var1 = squeeze(var1(val,:,:)); %#ok<NASGU>
    refreshdata(hp,'caller')
    var2 = [1,1]*w(val,3);  %#ok<NASGU>
    refreshdata(hd,'caller')  
    var3 = [1,1]*w(val,2);  %#ok<NASGU>
    refreshdata(ht,'caller')  
    ax.Subtitle.String = sprintf('Time = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g',...
                   string(t(val)),w(val,1),w(val,2),w(val,3));
    drawnow;
    %update slider selection text
    stxt = findobj(hfig,'Tag','FrameTime');
    stxt.String = string(time);
end