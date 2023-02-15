classdef RayTracks < muiDataSet
%-------class help---------------------------------------------------------
% NAME
%  RayTracks.m
% PURPOSE
%   Class description - Class for constructing array of wave ray tracks 
%   as a function of wave direction, wave period and water level, working 
%   in forward or backward ray tring mode. Used in the WaveRayodel app
% NOTES
%   Method is based on Abernethy C L and Gilbert G, 1975, Refraction of 
%   wave spectra, Report No: INT 117,pp. 1-166, Wallingford, UK.
% SEE ALSO
%   WaveRayModel.m, RayTrack.m, Ray.m, arc_ray
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
%     
    properties
        Tracks     %array of Ray objects for Nr rays, Nt periods and Nz water levels
        tol        %struct: tol.angle tolerance to test for angles that are 
                   %multiples of pi/2 and tol.dist for distances from axis
    end
    
    methods (Access = private)
        function obj = RayTracks()                
            %class constructor
            obj.tol.angle = 0.01;    %tolerance to test for angles that are multiples of pi/2 (0.06 deg) 
            obj.tol.dist = 0.001;    %tolerance to test for distances from axis in local coordinates
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,src)
            %function to run a simple 2D diffusion model
            obj = RayTracks;                          

            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in WaveRayModel
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;                    
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
            runobj = mobj.Inputs.WRM_RunParams;
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData'};
            grdobj = selectCaseObj(muicat,[],gridclasses,promptxt);
            grdrec = caseRec(muicat,grdobj.CaseIndex);
            if isempty(grdobj), return; end
            %assign the run parameters to the model instance
            switch src.Text
                case 'Forward Rays'
                    mobj.ModelInputs.RayTracks = {'WRM_RunParams','WRM_FT_Params'};
                case 'Backward Rays'
                    mobj.ModelInputs.RayTracks = {'WRM_RunParams','WRM_BT_Params'};
            end
            setRunParam(obj,mobj,grdrec); %input caserecs passed as varargin 

            grid = getGrid(grdobj,1);
            if isempty(grid.z), return; end

            [X,Y] = meshgrid(grid.x,grid.y);            
            delta = grid.x(2)-grid.x(1);
            obj.tol.dist = 1/1000/delta;
            cgrid = struct('X',X,'Y',Y,'z',grid.z);

            %arrays of waver periods and water levels, 
            %frequency at log spaced intervals, 
            T = runobj.PeriodRange;
            if length(T)>1
                frng = num2cell(1./runobj.PeriodRange); %
                f = log10(logspace(frng{:},runobj.nPeriod))';
                T = round(1./f,1); %round to one decimal place
            end
            %water levels at linear intervals
            zwl = runobj.WaterLevelRange;
            if length(zwl)>1
                WLrng = num2cell(runobj.WaterLevelRange);
                zwl = linspace(WLrng{:},runobj.nWaterLevel)';  
                zwl = round(zwl,1);
            end   
            %cutoff depth for wave ray tracing
            hlimit = runobj.hCutOff;

            %get the celerity, group celerity and celeirty gradient grids
            cgrid = celerity_grid(cgrid,T,zwl,delta);
            
            switch src.Text
                case 'Forward Rays'
                    [rays,rownames] = forwardTrack(obj,cgrid,T,zwl,hlimit);
                    modeltype = 'forward_model';
                case 'Backward Rays'
                    [rays,rownames] = backwardTrack(obj,cgrid,T,zwl,hlimit);
                    modeltype = 'backward_model';
            end
            % hf.Visible = 'on';
            raytable = setTable(obj,rays);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dsp = modelDSproperties(obj,modeltype);
            dst = dstable(raytable,'RowNames',rownames','DSproperties',dsp);
            dst.Dimensions.Period = T;           %wave period range
            dst.Dimensions.WaterLevel = zwl;     %water level range
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = modeltype;
            setDataSetRecord(obj,muicat,dst,modeltype);
            getdialog('Run complete');
        end
%%
        function checkStart(mobj)
            %create a plot of the start points on the selected bathymetry
            %grid with arrows showing the defined wave direction
            runobj = mobj.Inputs.WRM_FT_Params;
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData'};
            grdobj = selectCaseObj(mobj.Cases,[],gridclasses,promptxt);
            if isempty(grdobj), return; end
            grid = getGrid(grdobj,1);
            if isempty(grid.z), return; end

            %line vector of equal spaced start points
            [x_start,y_start] = RayTracks.getStartPoints(runobj.leftXY,...
                                                 runobj.rightXY, runobj.nRay);
            %transform wave direction to grid (trigonometric) direction
            [~,alpha] = compass2trig(runobj.dir0TN);

            hf = figure;  src = axes(hf);
            tabPlot(grdobj,src)
            %add start line points to figure
            ax = findobj(hf.Children,'Type','axes');
            hold on
            plot (ax,x_start,y_start,'-ok')
            hold off
            %add start direction to figure 
            RayTracks.plotWaveDirection(ax,x_start,y_start,alpha);   
        end
    end
%%
    methods
        function tabPlot(obj,src,mobj) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %need base grid based on ID save in RunParams           
            %if deleted warn user
            %select wave priod, water level and direction (for forward ray)
            %plot start points (first row of ray table)

            muicat = mobj.Cases;
            %retrieve bathymetry grid and plot
            caserec = caseRec(muicat,obj.RunParam.WRM_Bathy.caseid); 
            gobj = getCase(muicat,caserec);
            %set >Figure button and create axes
            tabcb  = @(src,evdat)tabPlot(obj,src,mobj); 
            ax = tabfigureplot(obj,src,tabcb,false);
            ax.NextPlot = 'add';
            %use tabPlot in WRM_Bathy to plot the bathymetry grid
            tabPlot(gobj,ax);

            %add the start arrows
            [~,~,delta] = getGridDimensions(gobj.RunParam.GD_GridProps);
            if strcmp(obj.Data.Dataset.MetaData,'forward_model')
                plotArrow(obj,ax,delta);
            else
                plotStartPoint(obj,ax)
            end
            %add the rays
            [~,np,nq] = size(obj.Data.Dataset.xr);
            if np>1 || nq>1
                options = get_selection(obj);
                if isempty(options), return; end    %user cancelled selection
            else
                options = [1,1,2];
            end
            plotRay(obj,ax,options);
            %update title
            title(obj.Data.Dataset.Description);
        end
    end 
%%    
    methods (Access = private)        
        function agrid = subSampleGrid(~,cgrid,j,k)
            %select grids for give wave period (j) and water level (k)
            agrid.X = cgrid.X;
            agrid.Y = cgrid.Y;
            agrid.h = cgrid.h(:,:,k);
            agrid.c = cgrid.c(:,:,j,k);
            agrid.cg = cgrid.cg(:,:,j,k);
            agrid.dcx = cgrid.dcx(:,:,j,k);
            agrid.dcy = cgrid.dcy(:,:,j,k);
        end
%%
        function [rays,rownames] = forwardTrack(obj,cgrid,T,zwl,hlimit)
            %construct set of ray tracks for given wave direction and a set
            %of start points using forward wave ray tracing
            
            ftrobj = obj.RunParam.WRM_FT_Params;
            %line vector of equal spaced start points
            [x_start,y_start] = RayTracks.getStartPoints(ftrobj.leftXY,...
                                              ftrobj.rightXY, ftrobj.nRay);
            %transform direction 'from' to direction ray is travelling 'to'
            %and convert direction degTN to grid (trigonometric) direction
            [~,alpha] = compass2trig(ftrobj.dir0TN);

%             hf = figure('Name','Search','Tag','PlotFig');
%             ax = axes(hf);
%             set(ax,'xgrid','on')
%             set(ax,'ygrid','on')

            nr = ftrobj.nRay;
            np = length(T);
            nq = length(zwl);         
            rays{nr,np,nq} = Ray;
            rownames = 1:nr;
            parfor i=rownames           %ray number
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,cgrid,j,k);
                        xys = [x_start(i),y_start(i)];
                        rayobj = Ray.setRay(agrid,xys,alpha,hlimit,obj.tol);
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
%                         hold on
%                         xr = rayobj.Track.xr;
%                         yr = rayobj.Track.yr;
%                         plot(ax,xr,yr,'-k')
%                         hold off
                        rays{i,j,k} = rayobj;
                    end
                end
            end
        end
%%
        function [rays,rownames] = backwardTrack(obj,cgrid,T,zwl,hlimit)
            %construct set of ray tracks from a point and a set inshore wave
            %directions using backward wave ray tracing

            btrobj = obj.RunParam.WRM_BT_Params;
            xys = btrobj.StartPoint;
            delta = cgrid.X(1,2)-cgrid.X(1,1);
            if all(rem(xys,delta)==0)
                %start point is on a node so move by a small increment - 1m
                xys = xys+1;
            end

            %get the array of inshore wave directions
            phi = btrobj.DirectionRange;
            if length(phi)>1
                dirng = num2cell(btrobj.DirectionRange);
                phi = linspace(dirng{:},btrobj.nDirections)';  
            end  
            rownames = round(phi',1);

            %transform wave direction to grid (trigonometric) direction
            alpha = mod(compass2trig(phi),2*pi);

            nd = btrobj.nDirections;
            np = length(T);
            nq = length(zwl);         
            rays{nd,np,nq} = Ray;
            parfor i=1:nd               %ray direction
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,cgrid,j,k);
                        rayobj = Ray.setRay(agrid,xys,alpha(i),hlimit,obj.tol);
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
                        rays{i,j,k} = rayobj;
                    end
                end
            end

        end
%%
        function plotArrow(obj,ax,delta)
            %plot an arrow at the start of each ray
            nrays = height(obj.Data.Dataset);
            xs = zeros(size(obj.Data.Dataset)); ys = xs; alpha = xs;
            for i=1:nrays
                xs(i) = obj.Data.Dataset.xr{i,1,1}(1);
                ys(i) = obj.Data.Dataset.yr{i,1,1}(1);
                alpha(i) = obj.Data.Dataset.alpha{i,1,1}(1);
            end
            hold(ax,"on")
            [ua,va] = pol2cart(alpha,delta);
            hq = quiver(ax,xs-ua,ys-va,ua,va,'AutoScale',0);
            hq.Color = 'r';
            hq.MaxHeadSize = 2.0;
            hq.LineWidth = 1.0;
            hold(ax,"off")
        end
%%
        function plotStartPoint(obj,ax)
            %plot point at start of backwardt tracking rays
            xs = obj.Data.Dataset.xr{1,1,1}(1);
            ys = obj.Data.Dataset.yr{1,1,1}(1);
            hold(ax,"on")
            plot(ax,xs,ys,'or')
            plot(ax,xs,ys,'.g')
            hold(ax,'off')
        end
%%
        function plotRay(obj,ax,opt)
            %plot each ray on the selected forward track
            %each ray can be a different length so plot iteratively
            cvar = {'celerity','cgroup'};
            nrays = height(obj.Data.Dataset);
            hold(ax,"on")
            if opt(3)>1
                for i=1:nrays
                    xr = obj.Data.Dataset.xr{i,opt(1),opt(2)};
                    yr = obj.Data.Dataset.yr{i,opt(1),opt(2)};
                    plot(ax,xr,yr,'-k');
    %                 if opt(3)==1
    %                     plot(ax,xr,yr,'-k');
    %                 else
    %                     var = obj.Data.Dataset.(cvar{opt(3)-1}){i,opt(1),opt(2)}; 
                            %need to gather point and plot as a set***
                            %see https://uk.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure-in-r2014a-and-earlier
    %                     %create axes for the ray surfaces axes
%                 hsax = axes;
%                 %set visibility for axes to 'off' so it appears transparent
%                 axis(hsax,'off')
%                 colormap(hsax,cool);
%                 s = scatter(hsax,xr,yr,[],var,'fill');
%                 cbCM = colorbar(hsax,'Location','east');
%                 %link the two overlaying axes so they match at all times to remain accurate
%                 linkaxes([ax,hsax]);

    %                     cmap = interpolate_cbrewer(cool,'linear',length(xr));
    % 
    %                     s.CData = cmap;
    %                 end
                end
            else
                nperiods = length(obj.Data.Dataset.Dimensions.Period);
                nwls = length(obj.Data.Dataset.Dimensions.WaterLevel);
                for i=1:nrays
                    for j=1:nperiods
                        for k=1:nwls
                            xr = obj.Data.Dataset.xr{i,j,k};
                            yr = obj.Data.Dataset.yr{i,j,k};
                            plot(ax,xr,yr,'-k')
                        end
                    end
                end
            end

            hold(ax,"off")
        end
%%       
function options = get_selection(obj)
            %get index of period, water level and variable to use in plots or model
            %   Defined using varargin for the following fields
            %    FigureTitle     - title for the UI figure
            %    PromptText      - text to guide user on selection to make
            %    InputFields     - text prompt for input fields to be displayed
            %    Style           - uicontrols for each input field (same no. as input fields)
            %    ControlButtons  - text for buttons to edit or update selection 
            %    DefaultInputs   - default text or selection lists
            %    UserData        - data assigned to UserData of uicontrol
            %    DataObject      - data object to use for selection
            %    SelectedVar     - index vector to define case,dataset,variable selection  
            %    ActionButtons   - text for buttons to take action based on selection
            %    Position        - poosition and size of figure (normalized units)
            T = obj.Data.Dataset.Dimensions.Period;
            zwl = obj.Data.Dataset.Dimensions.WaterLevel;
            var = {'all lines','selected line','celerity','group celerity'};
            selection = inputgui('FigureTitle','Celerity',...
                                 'InputFields',{'Wave Period','Water level'...
                                                               'Variable'},...
                                 'Style',{'popupmenu','popupmenu','popupmenu'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{string(T),string(zwl),var},...
                                 'PromptText','Select values to use');
            if isempty(selection)
                options = []; 
            else
                ki = selection{1};
                li = selection{2};
                mi = selection{3};
                options = [ki,li,mi];
            end  
        end
%%
        function newtable = setTable(~,rays)
            %make the table contents of a ray a single row
            [nr,np,nq] = size(rays);
            newtable = [];
            varnames = {'xr','yr','alpha','depth','celerity','cgroup'};         
            for i=1:nr                  %ray number
                perzwl = cell(1,np,nq);
                varcells = repmat({perzwl},1,6);
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        atable = rays{i,j,k}.Track.DataTable;
                        for l=1:6
                            varcells{1,l}(1,j,k) = {atable.(varnames{l})};
                        end
                    end
                end
                raytable = table(varcells{:},'VariableNames',varnames);
                if isempty(newtable)
                    newtable = raytable;
                else
                    newtable = [newtable;raytable]; %#ok<AGROW> 
                end
            end
        end       
%%
        function dsp = modelDSproperties(~,modeltype) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors                
            dsp.Variables = struct(...  
                'Name',{'xr','yr','alpha','depth','celerity','cgroup'},...
                'Description',{'X-Position','Y-Position','Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'Unit',{'m','m','rad','m','m/s','m/s'},...
                'Label',{'X-Position','Y-Position','Direction (rad)',...
                         'Water depth (m)','Celerity (m/s)','Group Celerity (m/s)'},...
                'QCflag',repmat({'model'},1,6));              
            dsp.Dimensions = struct(...    
                'Name',{'Period','WaterLevel'},...
                'Description',{'Wave Period','Water Level'},...
                'Unit',{'m','mOD'},...
                'Label',{'Wave Period','Water Level'},...
                'Format',{'-','-'});  
            switch modeltype
                case 'forward_model'
                     dsp.Row = struct(...                        
                        'Name',{'Ray'},...
                        'Description',{'Forward ray'},...
                        'Unit',{'-'},...
                        'Label',{'Ray number'},...
                        'Format',{'-'});  
                case 'backward_model'
                    dsp.Row = struct(...
                        'Name',{'InDir'},...
                        'Description',{'Inshore Direction'},...
                        'Unit',{'degTN'},...
                        'Label',{'Direction (degTN)'},...
                        'Format',{'-'});
            end
        end
    end   
%%
    methods (Static, Access=private)
         function [x,y] = getStartPoints(xy1,xy2,npnts)
            %define co-ordinates of equi-spaced points along line
            %Matlab Forum: 63233-interpolating-the-2d-line-to-make-the-new-coordinates-equi-distant?s_tid=ta_ans_results
            pathXY = [xy1;xy2];
            stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
            stepLengths = [0; stepLengths] ;% add the starting point
            cumulativeLen = cumsum(stepLengths);
            finalStepLocs = linspace(0,cumulativeLen(end), npnts);
            finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
            x = finalPathXY(:,1);
            y = finalPathXY(:,2);
        end
%%
        function plotWaveDirection(ax,xnd,ynd,dir)
            %plot arrows to show the initial wave direction
            %xnd,ynd define position of arrow head and dir the direction 'to'
            if isscalar(xnd)
                del = 200;                %default length of 200m
            else
                del = xnd(2)-xnd(1);
            end
            dx = del*cos(dir);
            dy = del*sin(dir);
            xst = xnd-dx;
            yst = ynd-dy;
            uni = ones(size(xst));
            %plot array of arrows
            hold on
            quiver(ax,xst,yst,uni*dx,uni*dy,'AutoScale',0)
            hold off
        end   
    end
end