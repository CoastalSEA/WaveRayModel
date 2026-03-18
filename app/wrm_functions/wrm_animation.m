function wrm_animation(mobj,sptobj,tsdst,offobj,inobj)
%
%-------function help------------------------------------------------------
% NAME
%   wrm_animation.m
% PURPOSE
%   animation of model spectra timeseries
% USAGE
%   wrm_animation(obj,sptobj,tsdst,offobj,inobj)
% INPUTS
%   mobj - instance of WaveRayModel class
%   sptobj - instance of SpectralTransfer class
%   tsdst - dstable with selected offshore wave timeseries data
%   offobj - array of offshore wave spectra (ctWaveSpectrum class)
%   inobj - array of inshore wave spectra (ctWaveSpectrum class)
% OUTPUT
%   animation figure
% SEE ALSO
%   code derived from muiPlots.newAnimation and implemented to show 2
%   subplots in an animation
% NOTES
%   obj is a data class instance and pobj is a plot class muiPlots instance
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%--------------------------------------------------------------------------
%
    hfig = figure('Name','Animation', 'Units','normalized', ...
                    'Resize','on','HandleVisibility','on','Visible','off',...
                    'Position',[0.38,0.42,0.30,0.42],'Tag','PlotFig');
    %create an instance of muiPlots and populate the properties that are   
    %needed for the newAnimation method
    if isa(mobj.mUI.Plots,'muiPlots')
        pobj = mobj.mUI.Plots;    %get existing instance          
        clearPreviousPlotData(pobj);
    else
        pobj = muiPlots.get_muiPlots();                   %create new instance
    end       
    pobj.Plot.CurrentFig = hfig;
    pobj.Plot.FigNum = hfig.Number;
    pobj.ModelMovie = [];
    pobj.Title = sprintf('Case: %s',tsdst(1).Description);

    %extract the timeseries data and dimensions for plot
    offsp = [offobj(:).Spectrum];
    insp = [inobj(:).Spectrum];
    pobj.Data.X = 1./offsp(1).freq;
    pobj.Data.Y = offsp(1).dir;
    pobj.Data.T = tsdst.RowNames; 
    nd = numel(offsp(1).dir); nf = numel(offsp(1).freq); nt = numel(tsdst.RowNames);
    off = zeros(nt,nd,nf); in = off; xyz.off = zeros(nt,3); xyz.in = xyz.off;
    %reshape as an array of dimensions [nt,nd,nf]
    for i=1:nt
        off(i,:,:) = offsp(i).SG;
        in(i,:,:) = insp(i).SG;
        [x,y,z] = markerPosition(pobj.Data.X,pobj.Data.Y,offsp(i).SG);
        xyz.off(i,:) = [x,y,z];
        [x,y,z] = markerPosition(pobj.Data.X,pobj.Data.Y,insp(i).SG);
        xyz.in(i,:) = [x,y,z];
    end
    pobj.Data.Z = struct('off',off,'in',in);
    
    % pobj.Data.Z = {reshape([offsp.SG],[],nd,nf), reshape([insp.SG],[],nd,nf)};
       
    if any(strcmp(tsdst.VariableNames,'Tp'))
        pobj.Data.Waves = [tsdst.Hs,tsdst.Tp,tsdst.Dir];
    else
        Hs = tsdst.Hs;
        [~,idx] = max(tsdst.S,[],2); 
        Tp = 1./tsdst.Dimensions.freq(idx);
        pkDir = tsdst.Dir(idx);
        pobj.Data.Waves = [Hs,Tp,pkDir];
    end

    answer = questdlg('What type of plot','OI spectrum','XY','Polar','XY');
    if strcmp(answer,'XY')
        pobj.MetaData = true;          %Cartesian dir-freq plot
    else
        pobj.MetaData = false;         %Polar dir-freq plot  
    end

    isfixed = true;
    answer = questdlg('Allow z-axis scale to vary?','Animation','Yes','No','Yes');
    if strcmp(answer,'Yes'), isfixed = false; end

    [s1,s2] = setupAnimation(sptobj,pobj,isfixed);
    if ~isvalid(pobj.Plot.CurrentFig), return; end

    getAnimation(pobj,s1,s2,xyz,hfig);
    s1.UserData = pobj.Data;  %store data set in UserData to
                              %allow user to switch between plots
    s1.UserData.xyz = xyz;
    %add replay and slider
    wrmControlPanel(pobj,hfig,length(pobj.Data.T),string(pobj.Data.T(1)));

    %assign muiPlots instance to handle
    mobj.mUI.Plots = pobj;
end

%%
function [s1,s2] = setupAnimation(sptobj,pobj,isfixed)
    %initialise 3Dplot and setup animation variables
    hfig = pobj.Plot.CurrentFig;
    figax = axes(hfig); 
    var1 = squeeze(pobj.Data.Z.off(1,:,:)); 
    var2 = squeeze(pobj.Data.Z.in(1,:,:)); 
    hfig.Visible = 'on';
    [s1,s2] = off_in_plot(sptobj,pobj.Data.X,pobj.Data.Y,var1,var2,...
                                                    figax,pobj.MetaData);
    if ~isvalid(hfig), return; end

    zMax = max([pobj.Data.Z.off],[],'all')/2;
    %assign axes properties                
    % s1.ZLimMode = 'manual'; %fix limits of z-axis
    % s1.ZLim = minmax(pobj.Data.Z.off);   
    s1.NextPlot = 'replaceChildren';
    s1.Tag = 'PlotFigAxes1'; 
    if isfixed
     s1.CLim(2) = zMax; 
    end
    hp1 = findobj(s1.Children,'Type','surface');    
    hp1.CDataSource = 'var1';    %assign data source
    if pobj.MetaData             %XY plot - add markers
        hp3 = findobj(s1.Children,'Tag','DirPk');
        hp3.XDataSource = 'var3'; 
        hp3.YDataSource = 'var4';
    end

    %s2.ZLimMode = 'manual'; %fix limits of z-axis
    %s2.ZLim = minmax(pobj.Data.Z.in); 
    s2.NextPlot = 'replaceChildren';
    s2.Tag = 'PlotFigAxes2'; 
    if isfixed
        s2.CLim(2) = zMax;
    end
    hp2 = findobj(s2.Children,'Type','surface');     
    hp2.CDataSource = 'var2';    %assign data source   
    if pobj.MetaData             %XY plot - add markers
        hp4 = findobj(s2.Children,'Tag','DirPk');
        hp4.XDataSource = 'var5'; 
        hp4.YDataSource = 'var6';
    end

    %adjust position of plots and add title
    if pobj.MetaData                         %holds logical for isXY plot
        s1.Position = [0.13,0.58,0.70,0.34]; %make space for slider bar
        s2.Position = [0.13,0.15,0.70,0.34]; %make space for slider bar
    else                                     %polar plot
        hfig.Position(3) = 0.6;
        So = pobj.Data.Z.off; Si = pobj.Data.Z.in;
        nrec = size(So,1);
        ff = pobj.Data.X; dd = pobj.Data.Y;
        spobj = ctWaveSpectrum;
        parfor i=1:nrec                      %parfor loop
            Po(i,:,:) = reshapePolarGrid(ff,dd,squeeze(So(i,:,:)),spobj);
            Pi(i,:,:) = reshapePolarGrid(ff,dd,squeeze(Si(i,:,:)),spobj);
        end
        pobj.Data.Z = struct('off',Po,'in',Pi);
    end

    w = pobj.Data.Waves;
    sg = sgtitle(sprintf('%s \nTime = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g\n',pobj.Title,...
                     string(pobj.Data.T(1)),w(1,1),w(1,2),w(1,3)));
    sg.FontSize = 10;
    sg.Margin = 1;
    sg.Tag = 'PlotFigTitle';
end

%%
function getAnimation(pobj,s1,s2,xyz,hfig)
    %generate an animation for user selection.
    t = pobj.Data.T;  
    var = pobj.Data.Z;
    w = pobj.Data.Waves;
    nrec = length(t);
    Mframes(nrec) = struct('cdata',[],'colormap',[]);
    Mframes(1) = getframe(gcf); %NB print function allows more control of 
    hp1 = findobj(s1.Children,'Type','surface');  
    hp2 = findobj(s2.Children,'Type','surface');
    if pobj.MetaData            %XY plot - add markers
        hp3 = findobj(s1.Children,'Tag','DirPk');
        hp4 = findobj(s2.Children,'Tag','DirPk');
    end

    for i=2:nrec
        var1 = squeeze(var.off(i,:,:)); %#ok<NASGU> 
        refreshdata(hp1,'caller')
        var2 = squeeze(var.in(i,:,:));  %#ok<NASGU> 
        refreshdata(hp2,'caller')   
        if pobj.MetaData                %XY plot - add markers
            var3 = xyz.off(i,1);        %#ok<NASGU>
            var4 = xyz.off(i,2);        %#ok<NASGU>
            refreshdata(hp3,'caller')            
            var5 = xyz.in(i,1);         %#ok<NASGU>
            var6 = xyz.in(i,2);         %#ok<NASGU>
            refreshdata(hp4,'caller')  
        end
        sg = findobj(s1.Parent.Children,'Tag','PlotFigTitle');
        sg.String = sprintf('%s \nTime = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g\n',...
                       pobj.Title,string(t(i)),w(i,1),w(i,2),w(i,3));
        drawnow;                 
        Mframes(i) = getframe(gcf); 
        %NB print function allows more control of resolution 
    end
    idm = size(pobj.ModelMovie,1);            
    pobj.ModelMovie{idm+1,1} = hfig.Number;
    pobj.ModelMovie{idm+1,2} = Mframes;   %save movie to class property
end 

%%
function [s1,s2] = off_in_plot(obj,T,dir,var0,vari,ax,isXY)
    %plot offshore and inshore spectra
    if nargin<6
        hf = figure('Name','SpecTrans','Tag','PlotFig');
        ax = axes(hf);
    end
    labeli = 'Inshore direction (degTN)';
    label0 = 'Offshore direction (degTN)';
    labelx = 'Wave period (s)';
    shorenorm = obj.Data.Inshore.UserData.ShoreAngle+90;         
    grey = mcolor('light grey');
    % Position of peak marker - not used because needs updating at each timestep

    if isXY
        s1 = subplot(2,1,1,ax);
        surf_plot(obj,T,dir,var0,{'Spectral Energy (m^2s)',labelx,label0,'Off'},s1)
        hold on
            [x,y,z] = markerPosition(T,dir,var0);
            plot3(s1,x,y,z,'+y','MarkerSize',12,'Tag','DirPk')
        hold off        
        s2 = subplot(2,1,2);
        surf_plot(obj,T,dir,vari,{'Spectral Energy (m^2s)',labelx,labeli,'In'},s2)
        hold on
            [x,y,z] = markerPosition(T,dir,vari);
            plot3(s2,x,y,z,'+y','MarkerSize',12,'Tag','DirPk')                             
            xn = minmax(T); yn = [shorenorm,shorenorm]; 
            plot3(s2,xn,yn,z*[1,1],'Color',grey,'LineStyle','--','LineWidth',1);
        hold off
        s2.YLim = s1.YLim;
    else
        s1 = subplot(1,2,1,ax);
        polar_plot(obj,T,dir,var0,{'Spectral Energy (m^2s)',labelx,label0,'Off'},s1)
        s2 = subplot(1,2,2);
        polar_plot(obj,T,dir,vari,{'Spectral Energy (m^2s)',labelx,labeli,'In'},s2)
        hold on
            ang = compass2trig(shorenorm);
            xn = [0,max(T)*cos(ang)]; yn = [0,max(T)*sin(ang)];
            zn = max(vari,[],'all');
            plot3(s2,xn,yn,zn*[1,1],'Color',grey,'LineStyle','--','LineWidth',1);
        hold off                
        s2.YLim = s1.YLim;
    end
end 

%%
function hm = wrmControlPanel(pobj,hfig,nrec,t0)
    %intialise button to re-run animation and slider to scroll through
    %add playback button
    buttxt = 'Run';   
    butpos = [0.05,0.01,0.05,0.05];
    butcall = @(src,evt)wrm_runmovie(pobj,src,evt);
    buttag = 'runMovie';
    buttip = 'Press to rerun animation';         
    hm(1) = setactionbutton(hfig,buttxt,butpos,butcall,buttag,buttip);
    buttxt = 'Save';   
    butpos = [0.12,0.01,0.05,0.05];
    butcall = @(src,evt)wrm_runmovie(pobj,src,evt);
    buttag = 'saveMovie';
    buttip = 'Press to save animation';  
    hm(2) = setactionbutton(hfig,buttxt,butpos,butcall,buttag,buttip);
    %add slide control            
    hm(2) = uicontrol('Parent',hfig,...
            'Style','slider','Value',1,... 
            'Min',1,'Max',nrec,'sliderstep',[1 1]/nrec,...
            'Callback', @(src,evt)wrm_runmovie(pobj,src,evt),...
            'HorizontalAlignment', 'center',...
            'Units','normalized', 'Position', [0.2,0.015,0.6,0.04],...
            'Tag','stepMovie');
    hm(3) = uicontrol('Parent',hfig,...
            'Style','text','String',t0,'Units','normalized',... 
            'Position',[0.8,0.015,0.15,0.03],'Tag','FrameTime');
end

%%
function P = reshapePolarGrid(Period,Dir,S,sobj)
    %format spectral array to be in format required by polarplot3d
    wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
    %interpolate var(phi,T) onto plot domain defined by tints,rints
    %intervals match those set to initialise plot in SpectralTransfer.polar_plot
    tints = linspace(sobj(1).Interp.tlim{:}); 
    rints = linspace(sobj(1).Interp.rlim{:});
    [Tq,Rq] = meshgrid(tints,rints); 
    warning('off',wid)
    P = griddata(deg2rad(Dir),Period,S',Tq,Rq);
    P(isnan(P)) = 0;  %fill blank sector so that it plots the period labels
    warning('on',wid)
end

%%
function [x,y,z] = markerPosition(T,dir,var)
    [sz] = size(var);
    [z,ind] = max(var,[],'all');
    [idy,idx] = ind2sub(sz,ind);
    x = T(idx);
    y = dir(idy);
end
    