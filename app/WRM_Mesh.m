classdef WRM_Mesh < muiPropertyUI & muiDataSet & matlab.mixin.Copyable 
%
%-------class help---------------------------------------------------------
% NAME
%   WRM_Mesh.m
% PURPOSE
%   Generate triangular mesh of nearshore area from a Cartesian grid
% USAGE
%   obj = WRM_Mesh.setInput(mobj); %mobj is a handle to Main UI
%   obj = WRM_Mesh.runModel(mobj);
% SEE ALSO
%   inherits muiPropertyUI and muiDataSet
%
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
% 
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Upper shore level',...
                          'Shore mesh size',...
                          'Boundary mesh size',...
                          'Offshore boundary (1=bottom,2=left,3=top,4=right)',...
                          'Maximum radius-edge ratio',...
                          'Gradient limit',...
                          'Edge-Element length threshold',...
                          'Tria-Element length threshold',...
                          'Include mesh optimisation (1/0)',...
                          };
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        % Tri          %triangulation object
        % zlevels      %bed levels on triangulated mesh
        
    end

    properties (Transient)
        kind = 'DELFRONT';  %'DELFRONT', 'DELAUNAY' type of refinement employed.
        disp = inf;         %refinement verbosity. Set to INF for quiet execution.
    end

    properties
        zhw = 5             %zhw is highest bed level to define shoreline
        minshore = 50       %mesh size along shoreline
        minbound = 200      %mesh size along boundaries        
        offshore = 1        %offshore boundary 1=bottom, 2=left, 3=top, 4=right
        rho2 = 1.025        %the maximum allowable radius-edge ratio
        dhdx = 0.2          %scalar gradient-limit used in limhfn2
        siz1 = 1.333        %normalised rel. length threshold for edge-elements
        siz2 = 1.300        %normalised rel. length threshold for tria-elements  
        isopt = true        %include optimisation if true (ie mesh smooting)
    end   
%%   
    methods (Access=protected)
        function obj = WRM_Mesh(mobj)             
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
            obj = WRM_Mesh(mobj);    
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs','WRM_Mesh',obj);
        end   
%%
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            muicat = mobj.Cases;                      
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
            errtxt = 'No grid data available for selected case';
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData'};
            grdobj = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(grdobj), warndlg(errtxt); return; end

            grid = getGrid(grdobj,1);                        %source grid
            if isempty(grid.z), warndlg(errtxt); return; end

            obj = WRM_Mesh.setInput(mobj);
            if isempty(obj.zhw) 
                warndlg('Some mesh input values not defined (e.g. Upper shore level)')
                return;
            end
            setRunParam(obj,mobj);   

            answer = questdlg('Load or create mesh boundary?','Mesh',...
                              'Load boundary','Load shore','Create','Create');
            if strcmp(answer,'Load boundary')
                %load mesh boundary from a file. Boundary should include
                %islands and boundary has points spaced as required for grid
                [node,edge,part] = WRM_Mesh.loadMeshBoundary;
            else
                if strcmp(answer,'Load shore')
                    shore_xy = WRM_Mesh.loadMeshBoundary;        %x,y column vectors
                    if isempty(shore_xy), return; end
                    gdims = gd_dimensions(grid);
                    shore_xy = clipShore(obj,shore_xy,gdims);    %ensure within grid limits
                    isextract = false;
                else
                    %extract shoreline
                    paneltxt = 'Mesh shoreline selection';
                    shore_xy = PL_Boundary.Figure(grid,paneltxt,1,true);
                    if isempty(shore_xy), return; end            %x,y column vectors
                    isextract = true;
                end
                [node,edge,obj] = getMeshBoundary(obj,shore_xy,grid,isextract); %mesh boundary polygon
                if isempty(node), return; end
                %add internal refinement
                % [node,edge,part] = add internal points
                part{1} = 1:length(edge);
            end

            %now generate grid using mesh2d
            hw = waitbar(0,'Processing mesh');
            [vert,tria] = get_mesh2d(obj,node,edge,part);
            waitbar(0.5,hw)
            %convert to a triangulation object
            Tri = triangulation(tria,vert);
            %find elevations on triangulated mesh
            [X,Y] = meshgrid(grid.x,grid.y);
            waitbar(1,hw)
            zlevels = interp2(X,Y,grid.z',vert(:,1),vert(:,2),'makima');
            delete(hw)

            isok = checkplot(obj,Tri,zlevels);
            if ~isok, return; end
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------   
            dsp = setDSproperties(obj);
            dst = dstable({Tri},{zlevels},'DSproperties',dsp);
            %assign metadata about model            
            dst.Source =  sprintf('%s, using %s',metaclass(obj).Name,...
                                            grdobj.Data.Grid.Description);
            dst.MetaData = '';   
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setDataSetRecord(obj,muicat,dst,'mesh');
            getdialog('Run complete');         
        end
%%
        function [node,edge,part] = loadMeshBoundary()
            %load a mesh boundary from a text file with a two column vector
            %listing of nodes
            node = []; edge = []; part = [];
            [filename,path,~] = getfiles('MultiSelect','off','PromptText','Select file:');
            if filename==0, return; end  %user cancelled
            dataSpec = '%f %f';
            nhead = 1;
            [data,~] = readinputfile([path,filename],nhead,dataSpec); %see muifunctions
            node = [data{1},data{2}];
            range = 1:1:size(node,1);
            links = [2:1:size(node,1),1];
            edge = [range',links'];
            %add internal refinement
            % [node,edge,part] = add internal points
            part{1} = 1:length(edge);
        end
    end 
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            dst = obj.Data.Dataset;
            verts = dst.Tri{1}.Points;
            tria = dst.Tri{1}.ConnectivityList;
            zlevel = dst.zlevel{1};
            
            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end
            axis equal tight
    
            %plot form as a contour plot
            trisurf(tria,verts(:,1),verts(:,2),zlevel);
            view(2)
            %use the YlGnBu colormap generated in cbrewer. This then needs
            %to be interpolated to get a smooth surface plot
            cmap = cmap_selection(20);
            [interpcmap]=interpolate_cbrewer(cmap, 'spline', 200);            
            colormap(ax,interpcmap)
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
        function [node,edge,obj] = getMeshBoundary(obj,shore_xy,grid,isextract)
            %construct defintion of nodes and edges for mesh
            %return obj to capture any change to the class properties
            shore_xy(any(isnan(shore_xy), 2), :) = []; %remove trailing NaNs
            node = setBoundary(obj,shore_xy,grid);            
            % node = [[0,0];newshore_xy;[max(x),0]];

            %edge is an E-by-2 array of edge indexing, where each row in edge
            %represents an edge of the polygon
            range = 1:1:size(node,1);
            links = [2:1:size(node,1),1];
            edge = [range',links'];

            %create accept button figure
            figtitle = 'Mesh boundary selection';
            promptxt = 'Accept mesh boundary?';
            tag = 'MeshFig'; %used for collective deletes of a group
            butnames = {'Yes','Adjust','No'};
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames);
            ax = axes(h_plt);
            hp = plot(ax,node(:,1),node(:,2),'.-');
            %add offset to axis limits around boundary
            nint = obj.minbound; 
            ax.XLim(1) = ax.XLim(1)-nint;  ax.XLim(2) = ax.XLim(2)+nint; 
            ax.YLim(1) = ax.YLim(1)-nint;  ax.YLim(2) = ax.YLim(2)+nint; 

            ok = 0;
            while ok<1
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window
                    node = []; return;
                elseif strcmp(h_but.Tag,'No')
                    node = []; ok = 1;
                elseif strcmp(h_but.Tag,'Adjust')
                    obj = editProperties(obj);
                    if isextract
                        %extract shoreline
                        paneltxt = 'Mesh shoreline extraction';
                        shore_xy = PL_Boundary.Figure(grid,paneltxt,1,true);
                        if isempty(shore_xy)
                            node = []; 
                            delete(h_plt.Parent); return; 
                        end
                        shore_xy(any(isnan(shore_xy), 2), :) = []; %remove trailing NaNs
                    end
                    node = setBoundary(obj,shore_xy,grid);
                    range = 1:1:size(node,1);
                    links = [2:1:size(node,1),1];
                    edge = [range',links'];
                    hp.XData = node(:,1);
                    hp.YData = node(:,2);
                else
                    ok = 1;
                end
            end
            delete(h_plt.Parent)
        end
%%
    function shore_xy = clipShore(obj,shore_xy,gd)
        %ensure that the shoreline data points do not go outside of grid
        if ~istable(gd)
            gd = table(gd{:},'VariableNames',{'xmin','xmax','ymin','ymax'});
        end
        idx = shore_xy(:,1)<gd.xmin | shore_xy(:,1)>gd.xmax;
        idy = shore_xy(:,2)<gd.ymin | shore_xy(:,2)>gd.ymax;
        idd = idx | idy;
        shore_xy(idd,:) = [];    
        icol = 1;                    %dominant axis is x-axis
        if obj.offshore==2 || obj.offshore==4
            icol = 2;                %dominant axis is y-axis
        end
        [~,ic] = unique(shore_xy(:,icol),'sorted'); %sort over dominant axis
        shore_xy = shore_xy(ic,:);

    end
%%
        function node = setBoundary(obj,shore_xy,grid)
            %create the boundary with specified limits
            shorelength = sum(vecnorm(diff(shore_xy),2,2));
            %interpolate shoreline to user specified mesh size
            shorepts = round(shorelength/obj.minshore);
            newshore_xy = curvspace(shore_xy,shorepts);
            
            gd = gd_dimensions(grid);  %table of dimensions for grid
        
            %construct boundaries of mesh using user defined mesh sizes
            %convention for shore vector:
            %increases from xmin for offshore = 1 or 3 and 
            %increases from ymin for offshore = 2 or 4
            %boundargy convention is them:
            %start at xmin for offshore = 1 or 3 and 
            %start at ymin for offshore = 2 or 4
            %vector proceeds clockwise for 1 and 4 and anticlockwise when
            %offshore is 2 or 3           
            if obj.offshore==1 || obj.offshore==3
                if newshore_xy(1,1)>newshore_xy(end,1)  %test x-order
                    newshore_xy = flipud(newshore_xy);
                end
            else
                if newshore_xy(1,2)>newshore_xy(end,2)  %test y-order
                    newshore_xy = flipud(newshore_xy);
                end
            end
            xs = newshore_xy(:,1);
            ys = newshore_xy(:,2);

            if obj.offshore==1 || obj.offshore==3
                if obj.offshore==1
                    ystart = min(gd.ymin,min(ys)-obj.minbound);
                    sy1 = ystart:obj.minbound:ys(1);
                    sy3 = fliplr(gd.ymin:obj.minbound:ys(end));                    
                else
                    ystart = max(gd.ymax,max(ys)+obj.minbound);
                    sy1 = fliplr(ys(1):obj.minbound:ystart);
                    sy3 = ys(end):obj.minbound:ystart;
                end
                sx1 = ones(size(sy1))*gd.xmin;
                sx3 = ones(size(sy3))*gd.xmax;
                sx4 = fliplr(gd.xmin:obj.minbound:gd.xmax);
                sy4 = ones(size(sx4))*ystart;
            else  %offshore is 2 or 4
                if obj.offshore==2
                    xstart = min(gd.xmin,min(xs)-obj.minbound);
                    sx1 = xstart:obj.minbound:xs(1);
                    sx3 = fliplr(gd.xmin:obj.minbound:xs(end));                    
                else
                    xstart = max(gd.xmax,max(xs)+obj.minbound);
                    sx1 = fliplr(xs(1):obj.minbound:xstart);
                    sx3 = xs(end):obj.minbound:xstart;
                end
                sy1 = ones(size(sx1))*gd.ymin;
                sy3 = ones(size(sx3))*gd.ymax;
                sy4 = fliplr(gd.ymin:obj.minbound:gd.ymax);
                sx4 = ones(size(sy4))*xstart;
            end

            %initialise inputs for refine2
            %node defines mesh boundary as an N-by-2 array of polygonal vertices
            node = [sx1',sy1'; xs,ys; sx3',sy3'; sx4',sy4'];
            % Find unique rows
            node = unique(node, 'rows','stable');
        end
%%
        function [vert,tria] = get_mesh2d(obj,node,edge,part)
            %use mesh2d to generate the verticies and triangles for a mesh
            options = getPropertiesStruct(obj);
            options.kind = obj.kind; options.disp = obj.disp;
            %the generation of the mesh requires mesh2D
                %check present???
            %compute a size-estimate for a multiply-connected geometry
            [vlfs,tlfs, hlfs] = lfshfn2(node,edge,part,options);
            %spatial-indexing structure for a 2-simplex triangulation embedded in the two-dimensional plane.
            [slfs] = idxtri2(vlfs,tlfs);
            hfun = @trihfn2;            %evaluate a discrete mesh-size function        
            %construct mesh using (Frontal)-Delaunay-refinement for 
            %two-dimensional, polygonal geometries
            [vert,etri,tria,tnum] = refine2(node,edge,part,options,hfun,vlfs,tlfs,slfs,hlfs);
            %optimize the mesh using smooth2
            if options.isopt
                [vert,~,tria,~] = smooth2(vert,etri,tria,tnum,options);
            end
        end
%%
        function isok = checkplot(~,Tri,zlevel)
            %plot to confirm the mesh is to be saved
            %create accept button figure
            verts = Tri.Points;
            tria = Tri.ConnectivityList;

            figtitle = 'Mesh generation';
            promptxt = 'Save mesh?';
            tag = 'MeshFig'; %used for collective deletes of a group
            butnames = {'Yes','No'};
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames);
            ax = axes(h_plt);
            %plot form as a contour plot
            trisurf(tria,verts(:,1),verts(:,2),zlevel);
            view(ax,2)
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = 'Elevation (mAD)';
            xlabel('Length (m)'); 
            ylabel('Width (m)'); 

            isok = false;
            waitfor(h_but,'Tag');
            if ~ishandle(h_but) 
                %this handles the user deleting figure window    
                return;
            elseif strcmp(h_but.Tag,'Yes')
               isok = true;
            end  
            delete(h_plt.Parent)
        end
%%
        function [dspec,dsprop] = setDSproperties(~)
            %define the variables in the dataset
            %define the metadata properties for the demo data set
            dspec = struct('Variables',[],'Row',[],'Dimensions',[]); dsprop = dspec;    
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dspec.Variables = struct(...
                'Name',{'Tri','zlevel'},...
                'Description',{'Triangulation','Bed elevation'},...
                'Unit',{'-','mAD'},...
                'Label',{'Triangulation','Elevation (mAD)'},...
                'QCflag',repmat({'raw'},1,2)); 
            dspec.Row = struct(...
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});        
            dspec.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});   
        end
    end
end