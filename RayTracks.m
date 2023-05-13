classdef RayTracks < muiDataSet
%-------class help---------------------------------------------------------
% NAME
%  RayTracks.m
% PURPOSE
%   Class description - Class for constructing array of wave ray tracks 
%   as a function of wave direction, wave period and water level, working 
%   in forward or backward ray tracking mode. Used in the WaveRayodel app
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
        tol          %distance tolerance (m) scaled in setTolerances to 
                     %create struct for angle and distance tolerances 
                     %dependent on the model grid size
    end
    
    methods (Access = private)
        function obj = RayTracks()                
            %class constructor 
            obj.tol.dist = 1e-9;
            obj.tol.angle = atan(obj.tol.dist);
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,src)
            %function to run a simple forwared and backward ray tracing
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
            gridclasses = {'WRM_Bathy','GD_ImportData','WRM_Mesh'};
            grdobj = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(grdobj), return; end
            grdrec = caseRec(muicat,grdobj.CaseIndex);            
            %assign the run parameters to the model instance
            switch src.Text
                case 'Forward Rays'
                    mobj.ModelInputs.RayTracks = {'WRM_RunParams','WRM_FT_Params'};
                case 'Backward Rays'
                    mobj.ModelInputs.RayTracks = {'WRM_RunParams','WRM_BT_Params'};
            end
            setRunParam(obj,mobj,grdrec);      %input caserecs passed as varargin 

            answer = questdlg('Write log to file?','Rays','Yes','No','No');
            islog = false;
            if strcmp(answer,'Yes'), islog = true; end
            
            if isprop(grdobj,'kind')
                cgrid.Tri = grdobj.Data.Dataset.Tri{1};
                cgrid.z = grdobj.Data.Dataset.zlevel{1};
            else
                grid = getGrid(grdobj,1);
                if isempty(grid.z), return; end
    
                [X,Y] = meshgrid(grid.x,grid.y);            
                delta = grid.x(2)-grid.x(1);
                cgrid = struct('X',X,'Y',Y,'z',grid.z);
            end

            %arrays of waver periods and water levels            
            T = runobj.PeriodRange;
            if length(T)>1
                %frequency at log spaced intervals, 
                frng = num2cell(1./runobj.PeriodRange); %
                f = log10(logspace(frng{:},runobj.nPeriod))';
                T = round(1./f,1); %round to one decimal place
            end
            
            zwl = runobj.WaterLevelRange;
            if length(zwl)>1
                %water levels at linear intervals
                WLrng = num2cell(runobj.WaterLevelRange);
                zwl = linspace(WLrng{:},runobj.nWaterLevel)';  
                zwl = round(zwl,1);
            end   

            %cutoff depth for wave ray tracing
            hlimit = runobj.hCutOff;

            %get the celerity, group celerity and celeirty gradient grids
            if isprop(grdobj,'kind')
                cgrid = celerity_mesh(cgrid,T,zwl);
            else
                cgrid = celerity_grid(cgrid,T,zwl,delta);
            end
            
            switch src.Text
                case 'Forward Rays'
                    [rays,rownames] = forwardTrack(obj,cgrid,T,zwl,hlimit,islog);
                    modeltype = 'forward_model';
                    modelprop = obj.RunParam.WRM_FT_Params.dir0TN;
                case 'Backward Rays'
                    [rays,rownames,depths] = backwardTrack(obj,cgrid,T,zwl,hlimit,islog);                    
                    modeltype = 'backward_model';
                    modelprop = mean(depths,'omitnan');
            end
            
            [raytable,outflag] = setTable(obj,rays);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable is an array in the 'results' cell array
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
            dst.UserData.prop = modelprop; %direction or av. depth of point
            dst.UserData.flag = outflag;   %flag to define ray termination
            setDataSetRecord(obj,muicat,dst,modeltype);
            getdialog('Run complete');
        end
%--------------------------------------------------------------------------
% End of runModel
%--------------------------------------------------------------------------
%%
        function checkStart(mobj)
            %create a plot of the start points on the selected bathymetry
            %grid with arrows showing the defined wave direction
            runobj = mobj.Inputs.WRM_FT_Params;
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData','WRM_Mesh'};
            grdobj = selectCaseObj(mobj.Cases,[],gridclasses,promptxt);
            if isempty(grdobj), return; end

            %line vector of equal spaced start points
            [x_start,y_start] = RayTracks.getStartPoints(runobj.leftXY,...
                                                 runobj.rightXY, runobj.nRay);
            %transform wave direction to grid (trigonometric) direction
            [~,alpha] = compass2trig(runobj.dir0TN);

            hf = figure('Name','Start','Tag','PlotFig');  
            src = axes(hf);
            tabPlot(grdobj,src)
            %add start line points to figure
            ax = findobj(hf.Children,'Type','axes');
            axis tight
            hold on
            plot (ax,x_start,y_start,'-ok')
            hold off            
            %add start direction to figure 
            RayTracks.plotWaveDirection(ax,x_start,y_start,alpha);   
        end
%%
        function startDepth(mobj)
            %check the depths of the bacakward ray start point for the range of
            %water levels specified
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData','WRM_Mesh'};
            grdobj = selectCaseObj(mobj.Cases,[],gridclasses,promptxt);
            if isempty(grdobj), return; end

            spnt = mobj.Inputs.WRM_BT_Params.StartPoint;
            if isa(grdobj,'WRM_Mesh')
                method = 'natural';
                mesh = grdobj.Data.Dataset.Tri{1};
                zlevel = grdobj.Data.Dataset.zlevel{1};
                X = mesh.Points(:,1); Y = mesh.Points(:,2);
                bed = griddata(X,Y,zlevel,spnt(1),spnt(2),method);    
            else
                grid = getGrid(grdobj,1);
                if isempty(grid.z), return; end                
                [X,Y] = meshgrid(grid.x,grid.y);
                bed = interp2(X,Y,grid.z',spnt(1),spnt(2),'linear');
            end

            wls =  mobj.Inputs.WRM_RunParams.WaterLevelRange;
            if isscalar(wls)
                msgtxt = sprintf('Water depth at start point: %.2f m',wls-bed); 
            else
                msgtxt = sprintf('Water depths at start point: min=%.2f m; max=%0.2f m',...
                                               wls(1)-bed,wls(2)-bed);
            end
        
            hf = figure('Name','Start','Tag','PlotFig');  
            src = axes(hf);
            tabPlot(grdobj,src)
            %add start point to figure
            ax = findobj(hf.Children,'Type','axes');
            axis tight
            hold on
            plot (ax,spnt(1),spnt(2),'xr','LineWidth',1,'MarkerSize',8)
            plot (ax,spnt(1),spnt(2),'ok','LineWidth',1,'MarkerSize',8)            
            xtxt = ax.XLim(2)*0.05;
            ytxt = ax.YLim(2)*0.95;
            text(ax,xtxt,ytxt,msgtxt);
            hold off              
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
            if isfield(obj.RunParam,'WRM_Bathy')
                caserec = caseRec(muicat,obj.RunParam.WRM_Bathy.caseid); 
            elseif isfield(obj.RunParam,'GD_ImportData')
                caserec = caseRec(muicat,obj.RunParam.GD_ImportData.caseid); 
            elseif isfield(obj.RunParam,'WRM_Mesh')
                caserec = caseRec(muicat,obj.RunParam.WRM_Mesh.caseid); 
            else
                warndlg('Unknown grid/mesh input type'); return
            end
            gobj = getCase(muicat,caserec);

            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src,mobj);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end
            %use tabPlot in WRM_Bathy to plot the bathymetry grid
            tabPlot(gobj,ax);
            %use the YlGnBu colormap generated in cbrewer. This then needs
            %to be interpolated to get a smooth surface plot
            ax = gca;
            cmap = cmap_selection(20);       
            [interpcmap]=interpolate_cbrewer(cmap, 'spline', 200);            
            colormap(ax,interpcmap)
            shading interp
            
            %add the start arrows
            if isprop(gobj,'kind')
                delta = gobj.minshore;
            else
                grid = getGrid(gobj);
                delta = abs(grid.x(2)-grid.x(1));
            end

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
            ax = plotRay(obj,ax,options);
            %update title
            if options(3)~=2
                T = obj.Data.Dataset.Dimensions.Period(options(1));
                zwl = obj.Data.Dataset.Dimensions.WaterLevel(options(2));
                ttxt = sprintf('Rays for %s, period %.3g s; swl %.3g mOD',...
                                   obj.Data.Dataset.Description,T,zwl);
            else
                ttxt = sprintf('Rays for %s',obj.Data.Dataset.Description);
            end
            title(ax,ttxt);
        end
    end 
%%    
    methods (Access = private)   
%         function setTolerances(obj,delta)
%             %set the distance, radius and angle tolerances based on grid
%             %size
%             atol = obj.tol;  %constructor dimension setting of tolerance
%             obj.tol = []; 
%             obj.tol.dist = atol;
%             obj.tol.angle = atan(obj.tol.dist);
%             %obj.tol.radius = 0.2865/obj.tol.angle;            
%         end
%%
        function agrid = subSampleGrid(~,cgrid,j,k,ismesh)
            %select grids for give wave period (j) and water level (k)
            if ismesh
                agrid.Tri = cgrid.Tri;
                agrid.h = cgrid.h(:,k);
                agrid.c = cgrid.c(:,j,k);
                agrid.cg = cgrid.cg(:,j,k);
                agrid.dcx = cgrid.dcx(:,j,k);
                agrid.dcy = cgrid.dcy(:,j,k);
            else
                agrid.X = cgrid.X;
                agrid.Y = cgrid.Y;
                agrid.h = cgrid.h(:,:,k);
                agrid.c = cgrid.c(:,:,j,k);
                agrid.cg = cgrid.cg(:,:,j,k);
                agrid.dcx = cgrid.dcx(:,:,j,k);
                agrid.dcy = cgrid.dcy(:,:,j,k);
            end
        end
%%
        function [rays,rownames] = forwardTrack(obj,cgrid,T,zwl,hlimit,islog)
            %construct set of ray tracks for given wave direction and a set
            %of start points using forward wave ray tracing
            
            ftrobj = obj.RunParam.WRM_FT_Params;
            %line vector of equal spaced start points
            [x_start,y_start] = RayTracks.getStartPoints(ftrobj.leftXY,...
                                              ftrobj.rightXY, ftrobj.nRay);
            %transform direction 'from' to direction ray is travelling 'to'
            %and convert direction degTN to grid (trigonometric) direction
            [~,alpha] = compass2trig(ftrobj.dir0TN);

            filename = sprintf('Forwardtrack_log_%s.txt',char(datetime,"ddMMMyy_HH-mm"));       
            %check plot - comment out when using parfor in loop
            % hf = figure('Name','Search','Tag','PlotFig');
            % ax = axes(hf);
            % set(ax,'xgrid','on')
            % set(ax,'ygrid','on')
            %--------------------------------------------------------------
            if isfield(cgrid,'Tri'),ismesh=true; else, ismesh=false; end

            nr = ftrobj.nRay;
            np = length(T);
            nq = length(zwl);         
            rays{nr,np,nq} = Ray;
            rownames = 1:nr;
            nrec = nr*np*nq;
            hpw = PoolWaitbar(nrec, 'Processing Rays');
            parfor i=rownames           %ray number
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,cgrid,j,k,ismesh);
                        xys = [x_start(i),y_start(i)];
                        rayobj = Ray.setRay(agrid,xys,alpha,hlimit,obj.tol);
                        increment(hpw);
                        npt = height(rayobj.Track.DataTable);
                        if islog
                            lines = sprintf('Ray no: %d, period %.1f, level %.1f, points %d, outflag %d',...
                                                 i,T(j),zwl(k),npt,rayobj.outFlag); %#ok<PFBNS> 
                            writelines(lines,filename,WriteMode="append")
                        end
                        
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
                        %check plot - comment out when using parfor in loop
                        % hold on
                        % xr = rayobj.Track.xr;
                        % yr = rayobj.Track.yr;
                        % plot(ax,xr,yr,'-k')
                        % hold off
                        %--------------------------------------------------
                        rays{i,j,k} = rayobj;                        
                    end
                end
            end
            delete(hpw)
        end
%%
        function [rays,rownames,depths] = backwardTrack(obj,cgrid,T,zwl,hlimit,islog)
            %construct set of ray tracks from a point and a set inshore wave
            %directions using backward wave ray tracing
            btrobj = obj.RunParam.WRM_BT_Params;
            xys = btrobj.StartPoint;

            %get the array of inshore wave directions
            phi = btrobj.DirectionRange;
            if length(phi)>1
                dirng = num2cell(btrobj.DirectionRange);
                phi = linspace(dirng{:},btrobj.nDirections)';  
            end  
            rownames = round(phi',1);

            %transform wave direction to grid (trigonometric) direction
           %[~,alpha]  = compass2trig(phi); %change input 'from' to ray 'to'           
            alpha = mod(compass2trig(phi),2*pi); %ray directions 'to'

            filename = sprintf('Backtrack_log_%s.txt',char(datetime,"ddMMMyy_HH-mm"));                       
            %check plot - comment out when using parfor in loop
            % hf = figure('Name','Search','Tag','PlotFig');
            % ax = axes(hf);
            % set(ax,'xgrid','on')
            % set(ax,'ygrid','on')
            %--------------------------------------------------------------
            if isfield(cgrid,'Tri'),ismesh=true; else, ismesh=false; end

            nd = btrobj.nDirections;
            np = length(T);
            nq = length(zwl);         
            rays{nd,np,nq} = Ray;
            depths = zeros(nd,np,nq);
            nrec = nd*np*nq;
            hpw = PoolWaitbar(nrec, 'Processing Rays');
            parfor i=1:nd               %ray direction
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        agrid = subSampleGrid(obj,cgrid,j,k,ismesh);
                        rayobj = Ray.setRay(agrid,xys,alpha(i),hlimit,obj.tol);
                        increment(hpw);
                        npt = height(rayobj.Track.DataTable);
                        depths(i,j,k) = rayobj.Track.DataTable.depth(1);                        
                        if islog
                            lines = sprintf('Ray dir %.3f, period %.1f, level %.1f, points %d, outflag %d',...
                                                 phi(i),T(j),zwl(k),npt,rayobj.outFlag); %#ok<PFBNS> 
                            writelines(lines,filename,WriteMode="append")
                        end
                                                
                        if isempty(rayobj)
                            rays{i,j,k}.Track.DataTable = [];
                            continue; 
                        end
                        %check plot - comment out when using parfor in loop
                        % hold on
                        % xr = rayobj.Track.xr;
                        % yr = rayobj.Track.yr;
                        % plot(ax,xr,yr,'-k')
                        % hold off   
                        %--------------------------------------------------
                        rays{i,j,k} = rayobj;
                    end
                end
            end
            delete(hpw)
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
            hq.LineWidth = 1;         
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
        function ax = plotRay(obj,ax,opt)
            %plot each ray on the selected forward track
            %each ray can be a different length so plot iteratively
            % Note: the order of setting axis off and positioning colorbar
            % is important for plot overlay with different color maps to work.
            cvar = {'celerity','cgroup','depth'};
            nrays = height(obj.Data.Dataset);
            
            if opt(3)>2    
                %see https://uk.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure-in-r2014a-and-earlier
                %create axes for the ray surfaces
                hsax = axes(ax.Parent); %assign new axes to same figure or tab
                hsax.XLim = ax.XLim;
                hsax.YLim = ax.YLim;                
                colormap(hsax,cool);

                xr = []; yr = []; var = [];
                for i=1:nrays
                    xr = [xr;obj.Data.Dataset.xr{i,opt(1),opt(2)}]; %#ok<AGROW> 
                    yr = [yr;obj.Data.Dataset.yr{i,opt(1),opt(2)}]; %#ok<AGROW> 
                    var = [var;obj.Data.Dataset.(cvar{opt(3)-2}){i,opt(1),opt(2)}]; %#ok<AGROW> 
                end
                stnorm = scalevariable(var,'Normalised');
                sze = ceil(stnorm-2*min(stnorm));
                scatter(hsax,xr,yr,sze,var,'fill');
      
                %set visibility for axes to 'off' so it appears transparent
                axis(hsax,'off')
                linkaxes([ax,hsax]);      %link the two overlaying axes 
                hsax.Position = ax.Position;
                %axis tight %does not work with mesh plotted using 

                %colormap(hsax,flipud(colormap(hsax)))
                cb = colorbar(hsax,'Color',[1,1,1]);
                cb.Position = [0.20,0.15,0.03,0.5];
                cb.Label.String = cvar{opt(3)-2};      
                cb.FontWeight = 'bold';
                clim(hsax,[min(var), max(var)]*1.05);
            elseif opt(3)==1
                hold(ax,"on")
                zmx = ax.ZLim(2);
                for i=1:nrays                    
                    xr = obj.Data.Dataset.xr{i,opt(1),opt(2)};
                    yr = obj.Data.Dataset.yr{i,opt(1),opt(2)};
                    plot(ax,xr,yr,'-k');
                    if obj.Data.Dataset.UserData.flag(i,opt(1),opt(2))==-3
                        %flag=-3; intersection of arc and element not found
                        plot3(ax,xr(end),yr(end),zmx,'*r','MarkerSize',4)
                    elseif obj.Data.Dataset.UserData.flag(i,opt(1),opt(2))==0
                        %flag=-2; radius too small???
                        plot3(ax,xr(end),yr(end),zmx,'vr','MarkerSize',4)
                    elseif obj.Data.Dataset.UserData.flag(i,opt(1),opt(2))==-1
                        %flag=-1; ray returns to shore and stops in shallow water
                        plot3(ax,xr(end),yr(end),zmx,'xr','MarkerSize',4)
                    end
                end
                hold(ax,"off")
            else
                nperiods = length(obj.Data.Dataset.Dimensions.Period);
                nwls = length(obj.Data.Dataset.Dimensions.WaterLevel);
                hold(ax,"on")
                for i=1:nrays
                    for j=1:nperiods
                        for k=1:nwls
                            xr = obj.Data.Dataset.xr{i,j,k};
                            yr = obj.Data.Dataset.yr{i,j,k};
                            plot(ax,xr,yr,'-k')
                        end
                    end
                end
                hold(ax,"off")
            end
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
            var = {'selected case','all cases','celerity','group celerity','depth'};
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
        function [newtable,outflag] = setTable(~,rays)
            %make the table contents of a ray a single row
            [nr,np,nq] = size(rays);
            newtable = [];
            outflag = zeros(nr,np,nq);
            varnames = {'xr','yr','alpha','depth','celerity','cgroup'};         
            for i=1:nr                  %ray number or direction
                perzwl = cell(1,np,nq);
                varcells = repmat({perzwl},1,6);
                for j=1:np              %wave period
                    for k=1:nq          %water level
                        atable = rays{i,j,k}.Track.DataTable;
                        for l=1:6
                            varcells{1,l}(1,j,k) = {atable.(varnames{l})};
                        end
                        outflag(i,j,k) = rays{i,j,k}.outFlag;
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