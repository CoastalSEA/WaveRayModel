classdef WRM_Bathy < muiPropertyUI & GDinterface       
%
%-------class help---------------------------------------------------------
% NAME
%   WRM_RunParams.m
% PURPOSE
%   Generate idealised bathymetries using a linear slope or Dean's profile
%   on a linear or crenulate shoreline. Also option for a mound on a linear
%   slope
% USAGE
%   obj = WRMrunparams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI and GDintereface
%   uses the following external functions:
%   interpolate_cbrewer, cmap_selection, crenulate_bays, deanbeachprofile
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Type: 1=linear; 2=Dean; 3=mound',...
                          'Beach crest or HW level (mOD)',...
                          'Upper beach slope (1:s - enter value for s)',...
                          'Bed level 1km out from SWL (mOD) [or y,z]',...
                          'Linear shore (0) or Crenulate Bay (1)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        BathyType = 2          %bathymetry type: 1=linear; 2=Dean; 3=mound
        BeachCrestLevel = 5    %beach crest level (mOD)
        UpperBeachSlope = 20   %bed slope (1:bs)
        BedLevelat1km = -6     %bed level 1km from shore (mOD) [OR y,z coordinates]
        isbay = true           %logical true if crenulate bay to be used
    end    
%%   
    methods (Access=protected)
        function obj = WRM_Bathy(mobj)             
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            obj = WRM_Bathy(mobj);    
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs','WRM_Bathy',obj);
        end   
%%
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model                         

            obj = WRM_Bathy.setInput(mobj);
            if isempty(obj.BedLevelat1km), return; end
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in WaveRayModel
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model grid parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            %of any existing Case datasets used as inputs to the model. The
            %classes loaded are defined in the main model UI class:
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            %if fully defined RunParam property should contain all the
            %input parameters and dataset links needed to run the model.
            setRunParam(obj,mobj); 
            
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
           grid = getBathy(obj);
           if isempty(grid), obj = []; return; end
           [grid,rotate] = orientGrid(obj,grid); %option to flip or rotate grid 
           newgrid(1,:,:) = grid.z; 
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign X and Y dimensions 
            dims = struct('x',grid.x,'y',grid.y);                                               
            %assign metadata about data source and save grid
            meta.source = metaclass(obj).Name;
            meta.data = sprintf('Rotate option = %d',rotate);
            obj = setGrid(obj,{newgrid},dims,meta);
            
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setCase(muicat,obj,'grid_data');
            getdialog('Run complete');
        end
    end 
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            dst = obj.Data.Grid;
            grid = getGrid(obj,1);
            if isempty(grid.z), return; end
            
            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end
    
            %plot form as a contour plot
            pcolor(ax,grid.x,grid.y,grid.z');
            ax = gd_ax_dir(ax,grid.x,grid.y);
            %use the YlGnBu colormap generated in cbrewer. This then needs
            %to be interpolated to get a smooth surface plot
            cmap = cmap_selection(20);
            [interpcmap]=interpolate_cbrewer(cmap, 'spline', 200);            
            colormap(ax,interpcmap)
            shading interp
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = 'Elevation (mAD)';
            xlabel('Length (m)'); 
            ylabel('Width (m)'); 
            title(dst.Description);
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
        end
    end
%%
    methods (Access=private)
        function grid = getBathy(obj,grid)
            %generate an idealised bathymetry using a linear slope or the
            %Dean profile (default is Dean profile)
            
            gobj = obj.RunParam.GD_GridProps;           %gid dimensions
            [grid.x,grid.y,delx,~] = getGridDimensions(gobj); 
            
            zBC = obj.BeachCrestLevel;    %beach crest level (mOD)
            ubs = obj.UpperBeachSlope;    %bed slope (1:bs)
            z1km = obj.BedLevelat1km;     %bed level 1km out from SWL (mOD)

            %set crosshore profile or mound based on selection
            switch obj.BathyType
                case 1                                 %linear slope
                    bedslope = z1km/1000;
                    zp = [zBC;grid.y(1:end-1)*bedslope];
                case 2                                 %Dean's profile
                    [~,zp,~] = deanbeachprofile(grid.y',zBC,z1km,ubs,false);
                    %point for upper beach above msl added, so remove outer point
                    zp = zp(1:end-1); 
                case 3                                 %linear slope and mound
                    bedslope = z1km/1000;
                    zp = [zBC;grid.y(1:end-1)*bedslope];
                    promptxt = {'Height of mound','Radius of mound'};
                    defaults = {'2',num2str(4*delx)};
                    answer = inputdlg(promptxt,'Mound',1,defaults);
                    if isempty(answer), grid = []; return; end
                    hmound = str2double(answer{1});
                    hrad = str2double(answer{2});
                    xm = grid.x(1)+(grid.x(end)-grid.x(1))/2;
                    ym = grid.y(1)+(grid.y(end)-grid.y(1))/2;
                    Circ = get_circle(obj,hrad,xm,ym); 
                    [X,Y] = meshgrid(grid.x,grid.y);
                    XY = [reshape(X,[],1),reshape(Y,[],1)]; %x,y vectors
                    TFin = isinterior(Circ,XY);
                    grid.z = fliplr(repmat(zp,1,gobj.Xint+1)');
                    idx = reshape(TFin,size(X))';
                    grid.z(idx) = grid.z(idx)+hmound;
                    return;
                otherwise
                    return;
            end
    
            %add bay if required
            if obj.isbay
                inp = setBayProperties(obj,gobj.XaxisLimits(2));
                [E,N,~,~,~] = crenulate_bays(inp);
                delete(inp.ax.Parent)            %delete the crenulate bay figure
                offset = max(N(E>=0))-N(E>=0);   %offset from maximum extent of bay from control line
                offset = interp1(E(E>=0),offset,grid.x,'linear','extrap');
                np = length(grid.y);
                xyz = [];
                for i=1:length(grid.x)
                    xi = repmat(grid.x(i),np,1);
                    xyz = [xyz;[xi,grid.y+offset(i),zp]];  %#ok<AGROW> %created scattered xyz tuples.
                end
                F = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),'linear');
                [X,Y] = meshgrid(grid.x,grid.y);
                h = F(X,Y)';
                h(h>zBC) = zBC;
                grid.z = fliplr(h);               %assume relative to zero datum
            else
                grid.z = fliplr(repmat(zp,1,gobj.Xint+1)');
            end
        end
%%
        function inp = setBayProperties(~,xlim)
            %default data
            inp.Ro = xlim;        %length of control line
            inp.beta = 20;        %Angle between control line and wave crest
            inp.Eo = 0;           %Easting of control point
            inp.No = 0;           %Northing of control point
            inp.alpha = 290;      %Angle of wave crest to TN
            inp.p = 1;            %Spiral rotation (1=clockwise outwards and 0=counter-clockwise)
            hf = figure; inp.ax = axes(hf);
        end
%%
        function Circ = get_circle(~,radius,uc,vc)   
            %define a circle in local coordinates
            n = 50;                                  %number of points in circle
            phi = (0:n-1)*(2*pi/n);                  %angle increments
            u = uc + radius*cos(phi);                %coordinates of circle
            v = vc + radius*sin(phi);
            Circ = polyshape(u,v);
        end
    end
end