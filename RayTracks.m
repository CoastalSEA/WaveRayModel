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
    end
    
    methods (Access = private)
        function obj = RayTracks()                
            %class constructor
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
            dsp = modelDSproperties(obj);
            
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
            wavegrid = struct('X',X,'Y',Y,'z',grid.z);

            %arrays of waver periods and water levels, 
            %periods at log spaced intervals, 
            T = runobj.PeriodRange;
            if length(T)>1
                Trng = num2cell(runobj.PeriodRange);
                T = log10(logspace(Trng{:},runobj.nPeriod))';
                T = round(T,1); %round to one decimal place
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
            wavegrid = celerity_grid(wavegrid,T,zwl,delta);
            
            % hf = figure('Name','Search','Tag','PlotFig','Visible','off');
            % ax = axes(hf);
            % set(ax,'xgrid','on')
            % set(ax,'ygrid','on')

%             %start of Forward tracking loop
% 
%             %line vector of equal spaced start points
%             [x_start,y_start] = RayTracks.getStartPoints(runobj.leftXY,...
%                                               runobj.rightXY, runobj.nRay);
% 
%             %transform wave direction to grid (trigonometric) direction
%             [~,alpha] = compass2trig(runobj.dir0TN);
% 
%             nr = runobj.nRay;
%             np = length(T);
%             nq = length(zwl);
%             hlimit = 0.1;            
%             rays{nr,np,nq} = Ray;
%             parfor i=1:nr               %ray number
%                 for j=1:np              %wave period
%                     for k=1:nq          %water level
%                         agrid = subSampleGrid(obj,wavegrid,j,k);
%                         xys = [x_start(i),y_start(i)];
%                         rayobj = Ray.setRay(agrid,xys,alpha,hlimit);
%                         if isempty(rayobj)
%                             rays{i,j,k}.Track.DataTable = emptytable;
%                             continue; 
%                         end
%                         % hold on
%                         % xr = rayobj.Track.xr;
%                         % yr = rayobj.Track.yr;
%                         % plot(ax,xr,yr,'-k')
%                         % hold off
%                         rays{i,j,k} = rayobj;
%                     end
%                 end
%             end
            switch src.Text
                case 'Forward Rays'
                    [rays,rownames] = forwardTrack(obj,wavegrid,T,zwl,hlimit);
                    modeltype = 'forward_model';
                case 'Backward Rays'
                    [rays,rownames] = backwardTrack(obj,wavegrid,T,zwl,hlimit);
                    modeltype = 'backward_model';
            end
            % hf.Visible = 'on';
            raytable = setTable(obj,rays);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(raytable,'RowNames',rownames','DSproperties',dsp);
            dst.Dimensions.Period = T;           %wave period range
            dst.Dimensions.WaterLevel = zwl;     %water level range
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = 'Any additional information to be saved';
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
            plotArrow(obj,ax,delta);
            %add the rays
            [~,np,nq] = size(obj.Data.Dataset.xr);
            if np>1 || nq>1
                [idp,idq] = get_selection(obj);
            else
                idp = 1; idq = 1;
            end
            plotRay(obj,ax,idp,idq);
            %update title
            title(obj.Data.Dataset.Description);
        end
%%
        function agrid = subSampleGrid(~,wavegrid,j,k)
            %select grids for give wave period (j) and water level (k)
            agrid.X = wavegrid.X;
            agrid.Y = wavegrid.Y;
            agrid.h = wavegrid.h(:,:,k);
            agrid.c = wavegrid.c(:,:,j,k);
            agrid.cg = wavegrid.cg(:,:,j,k);
            agrid.dcx = wavegrid.dcx(:,:,j,k);
            agrid.dcy = wavegrid.dcy(:,:,j,k);
        end
    end 
%%    
    methods (Access = private)
        function [rays,rownames] = forwardTrack(obj,wavegrid,T,zwl,hlimit)
            %construct set of ray tracks for given wave direction and a set
            %of start points using forward wave ray tracing
            
            ftrobj = obj.RunParam.WRM_FT_Params;
            %line vector of equal spaced start points
            [x_start,y_start] = RayTracks.getStartPoints(ftrobj.leftXY,...
                                              ftrobj.rightXY, ftrobj.nRay);

            %transform wave direction to grid (trigonometric) direction
            [~,alpha] = compass2trig(ftrobj.dir0TN);

            % hf = figure('Name','Search','Tag','PlotFig','Visible','off');
            % ax = axes(hf);
            % set(ax,'xgrid','on')
            % set(ax,'ygrid','on')

            nr = ftrobj.nRay;
            np = length(T);
            nq = length(zwl);         
            rays{nr,np,nq} = Ray;
            rownames = 1:nr;
            for i=rownames           %ray number
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,wavegrid,j,k);
                        xys = [x_start(i),y_start(i)];
                        rayobj = Ray.setRay(agrid,xys,alpha,hlimit);
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
                        % hold on
                        % xr = rayobj.Track.xr;
                        % yr = rayobj.Track.yr;
                        % plot(ax,xr,yr,'-k')
                        % hold off
                        rays{i,j,k} = rayobj;
                    end
                end
            end
        end
%%
        function [rays,rownames] = backwardTrack(obj,wavegrid,T,zwl,hlimit)
            %construct set of ray tracks from a point and a set inshore wave
            %directions using backward wave ray tracing

            btrobj = obj.RunParam.WRM_BT_Params;
            xys = btrobj.StartPoint;
            delta = wavegrid.X(1,2)-wavegrid.X(1,1);
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
            rownames = round(phi,1);

            %transform wave direction to grid (trigonometric) direction
            alpha = mod(compass2trig(phi),2*pi);

            nd = btrobj.nDirections;
            np = length(T);
            nq = length(zwl);         
            rays{nd,np,nq} = Ray;
            for i=1:nd               %ray direction
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,wavegrid,j,k);
                        rayobj = Ray.setRay(agrid,xys,alpha(i),hlimit);
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
                        % hold on
                        % xr = rayobj.Track.xr;
                        % yr = rayobj.Track.yr;
                        % plot(ax,xr,yr,'-k')
                        % hold off
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
        function plotRay(obj,ax,idp,idq)
            %plot each ray on the selected forward track
            %each ray can be a different length so plot iteratively
            nrays = height(obj.Data.Dataset);
            hold(ax,"on")
            for i=1:nrays
                xr = obj.Data.Dataset.xr{i,idp,idq};
                yr = obj.Data.Dataset.yr{i,idp,idq};
                plot(ax,xr,yr,'-k');
            end
            hold(ax,"off")
        end
%%       
        function [ki,li] = get_selection(obj)
            %get index of period and water level to use in plots or model
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
            selection = inputgui('FigureTitle','Celerity',...
                                 'InputFields',{'Wave Period','Water level'},...
                                 'Style',{'popupmenu','popupmenu'},...
                                 'ActionButtons', {'Select','Close'},...
                                 'DefaultInputs',{string(T),string(zwl)},...
                                 'PromptText','Select values to use');
            if isempty(selection)
                ki = []; li = []; 
            else
                ki = selection{1};
                li = selection{2};
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
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
%             dsp.Variables = struct(...                     
%                 'Name',{'xr','yr','alpha','k','quad','edge'},...
%                 'Description',{'X-Position','Y-Position','Direction',...
%                            'Grid index','Grid quadrant','Element edge'},...
%                 'Unit',{'m','m','rad','-','-','-'},...
%                 'Label',{'X-Position','Y-Position','Direction',...
%                            'Grid index','Grid quadrant','Element edge'},...
%                 'QCflag',repmat({'model'},1,6)); 
            dsp.Variables = struct(...  
                'Name',{'xr','yr','alpha','depth','celerity','cgroup'},...
                'Description',{'X-Position','Y-Position','Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'Unit',{'m','m','rad','m','m/s','m/s'},...
                'Label',{'X-Position','Y-Position','Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'QCflag',repmat({'model'},1,6)); 
            dsp.Row = struct(...
                'Name',{'Ray'},...
                'Description',{'Forward ray'},...
                'Unit',{'-'},...
                'Label',{'Ray number'},...
                'Format',{'-'});   
            dsp.Dimensions = struct(...    
                'Name',{'Period','WaterLevel'},...
                'Description',{'Wave Period','Water Level'},...
                'Unit',{'m','mOD'},...
                'Label',{'Wave Period','Water Level'},...
                'Format',{'-','-'});  
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