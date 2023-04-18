classdef WRM_Mesh < muiPropertyUI & muiDataSet      
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
                          'Maximum radius-edge ratio',...
                          'Gradient limit',...
                          'Edge-Element length threshold',...
                          'Tria-Element length threshold',...
                          'Include mesh optimisation (1/0)',...
                          };

        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
%         Tri          %triangulation object
%         zlevels      %bed levels on triangulated mesh
        
    end

    properties (Transient)
        kind = 'DELFRONT';  %'DELFRONT', 'DELAUNAY' type of refinement employed.
        disp = inf;         %refinement verbosity. Set to INF for quiet execution.
    end

    properties
        zhw = 5             %zhw is highest bed level to define shoreline
        minshore = 50       %mesh size along shoreline
        minbound = 1000     %mesh size along boundaries        
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
            setClassObj(mobj,'Inputs','WRM_Bathy',obj);
        end   
%%
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------  
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = WRM_Mesh.setInput(mobj);
            if isempty(obj.zhw), return; end
            muicat = mobj.Cases;
            setRunParam(obj,mobj);             
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
            promptxt = 'Select grid to use for wave model'; 
            gridclasses = {'WRM_Bathy','GD_ImportData'};
            grdobj = selectCaseObj(muicat,[],gridclasses,promptxt);
            if isempty(grdobj), return; end

            grid = getGrid(grdobj,1);                    %source grid
            if isempty(grid.z), return; end

            shore_xy = getShoreline(obj,grid);           %extract shoreline
            if isempty(shore_xy), return; end

            [node,edge] = getMeshBoundary(obj,shore_xy,grid); %mesh boundary polygon
            if isempty(node), return; end
            %add internal refinement
            % [node,edge,part} = add internal points
            part{1} = 1:length(edge);

            %now generate grid using mesh2d
            [vert,tria] = get_mesh2d(obj,node,edge,part);

            %convert to a triangulation object
            Tri = triangulation(tria,vert);
            %find elevations on triangulated mesh
            [X,Y] = meshgrid(grid.x,grid.y);
            zlevels = interp2(X,Y,grid.z',vert(:,1),vert(:,2),'makima');

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
    
            %plot form as a contour plot
            trisurf(tria,verts(:,1),verts(:,2),zlevel);
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
        function shore_xy = getShoreline(obj,grid)
            %get contour at zhw to define shoreline     
            %could add option to select smoothing method and window size**
            [X,Y] = meshgrid(grid.x,grid.y);

            hf = figure('Name','Search','Tag','PlotFig','Visible','off');
            ax = axes(hf);
            surf(ax,X,Y,grid.z'); shading interp; colorbar;
            hold on 
            M = contourc(grid.x,grid.y,grid.z',[obj.zhw,obj.zhw]);
            hold off
            shore_xy = M(:,2:end);
%             idx = shore_xy(1,:)==0 & shore_xy(2,:)==0;
%             shore_xy(:,idx) = [];
            delete(hf)
            
            %create accept button figure
            figtitle = 'Mesh shoreline selection';
            promptxt = 'Accept shoreline?';
            tag = 'MeshFig'; %used for collective deletes of a group
            butnames = {'Yes','Smooth','Quit'};
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames);
            ax = axes(h_plt);
            plot(ax,shore_xy(1,:),shore_xy(2,:))

            ok = 0;
            while ok<1
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window    
                   shore_xy = []; return;
                elseif strcmp(h_but.Tag,'Quit')                    
                   shore_xy = []; ok = 1;
                elseif strcmp(h_but.Tag,'Smooth')
                   %Do something
                   h_but.Tag = '';
                   shore = smoothdata(shore_xy(2,:),'sgolay');
                   hold on 
                   plot(shore_xy(1,:),shore)
                   hold off
                   shore_xy(2,:) = shore;
                else
                   ok = 1;
                end  
            end
            delete(h_plt.Parent)
        end
%%
        function [node,edge] = getMeshBoundary(obj,shore_xy,grid)
            %construct defintion of nodes and edges for mesh
            node = setBoundary(obj,shore_xy);            
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
            ax.XLim = ax.XLim*1.05 - ax.XLim(2)*0.025;
            ax.YLim = ax.YLim*1.05 - ax.YLim(2)*0.025;

            ok = 0;
            while ok<1
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window
                    node = []; return;
                elseif strcmp(h_but.Tag,'No')
                    node = []; ok = 1;
                elseif strcmp(h_but.Tag,'Adjust')
                    obj = editProperties(obj);
                    shore_xy = getShoreline(obj,grid);           %extract shoreline
                    if isempty(shore_xy)
                        node = []; 
                        delete(h_plt.Parent); return; 
                    end
                    node = setBoundary(obj,shore_xy);
                    hp.XData = node(:,1);
                    hp.YData = node(:,2);
                else
                    ok = 1;
                end
            end
            delete(h_plt.Parent)
        end
%%
        function node = setBoundary(obj,shore_xy)
            %create the boundary with specified limits

            shorelength = sum(vecnorm(diff(shore_xy'),2,2));
            %interpolate shoreline to user specified mesh size
            shorepts = round(shorelength/obj.minshore);
            newshore_xy = curvspace(shore_xy',shorepts);

            %construct seaward boundaries of mesh using user define mesh sizes
            xmax = max(newshore_xy(:,1));
            lbound = 0:obj.minbound:newshore_xy(1,2);
            left_bound = zeros(length(lbound),2);
            left_bound(:,2) = lbound;
            rbound = fliplr(0:obj.minbound:newshore_xy(end,2));
            right_bound = ones(length(rbound),2);
            right_bound(:,2) = rbound;
            right_bound(:,1) = right_bound(:,1)*xmax;
            sbound = fliplr(0:obj.minbound:xmax);
            sea_bound = zeros(length(sbound)-2,2);
            sea_bound(:,1) = sbound(2:end-1);

            %initialise inputs for refine2
            %node defines mesh boundary as an N-by-2 array of polygonal vertices
            node = [left_bound;newshore_xy;right_bound;sea_bound];    
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