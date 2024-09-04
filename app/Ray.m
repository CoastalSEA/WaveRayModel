classdef Ray < handle                        
%
%-------class help---------------------------------------------------------
% NAME
%   Ray.m
% PURPOSE
%   Class description - Class for ray objects in the WaveRayModel app
% NOTES
%   Ray method based on Abernethy C L and Gilbert G, 1975, Refraction of 
%   wave spectra, Report No: INT 117,pp. 1-166, Wallingford, UK.
%   This version computes the route track and properties along track. 
%   The table ray is used to hold ray properties. The last row of the 
%   table is used to call arc_ray. Whilst constructing the ray the table
%   holds xr, yr, alpha, k, quad, edge, hr, cr, cgr. Once the ray is
%   complete k, quad, edge are removed from the table before assigning the
%   reultant ray table to the Track property.
% SEE ALSO
%   WaveRayModel.m, RayTracks.m, and arc_ray.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
%     
    properties
        Track         %table of track properties: xr,yr,alpha,k,quad,edge  
        outFlag = 0   %flag to indicate nature of ray termination
                      % 1 - successful ray;   -1 - ray ends on shore
                                        %-2 - radius too small; -3 - intersection not found
    end
    
    methods
        function obj = Ray()                
            %class constructor
        end
    end      
%%
    methods (Static)                
        function obj = setRay(cgrid,xys,alpha,hlimit,tol)
            %function to compute a single ray track
            % cgrid - contains X,Y,h,c,dcx,dcy
            % xys - start point grid coordinates
            % alpha - start angle (radians)
            % hlimit - cutoff depth limit for ray (optional)
            % tol - tolerance around angles that are multiples of pi/2
            if nargin<4
                hlimit = 0.1;
            end

            obj = Ray;                               %Ray class instance 

            if isfield(cgrid,'Tri')
                ray = meshRay(obj,cgrid,xys,alpha,hlimit,tol);
            else
                ray = gridRay(obj,cgrid,xys,alpha,hlimit,tol);
            end
            if obj.outFlag==0
                %find why????????
                fprintf('%d',height(ray))               
            end
            obj.Track = dstable(ray,'DSproperties',modelDSproperties(obj));
        end
    end
%%
    methods (Access=private)  
        function ray = gridRay(obj,cgrid,xys,alpha,hlimit,tol)
            %get ray path using Carteisan grid
            %get first first element intersection from the start point
            ray = startRay(obj,cgrid,xys,alpha,tol);
            if isempty(ray)
                error('Ray solution not found in start_ray')
            elseif ~isa(ray,'table')
                error('Ray start position outside grid domain')
            end   
            tol.radius = 1e3;
            %check plot for finding ray errors - comment out when not required
            % hf = figure('Name','Search','Tag','PlotFig');
            % ax = axes(hf);
            % plot(ax,ray.xr(1),ray.yr(1),'ob')
            %--------------------------------------------------------------

            %loop to get ray track to edge of grid or depth limit
            hr = ray.hr(1);
            while hr>hlimit
                newray = arc_ray(cgrid,ray(end,:),tol);
                if isempty(newray)                         %ray exits grid
                    obj.outFlag = 1; hr = hlimit; continue; 
                elseif isnumeric(newray) && newray<0       %ray error (eg radius too small)
                    obj.outFlag = newray; hr = hlimit; continue; 
                end
                ray = [ray;newray]; %#ok<AGROW> 
                %check plot for finding ray errors - comment out when not required
                % hold on
                % plot(ax,newray.xr,newray.yr,'+k')
                % hold off   
                %----------------------------------------------------------                
                hr = newray.hr;
                if hr<=hlimit, obj.outFlag = -1; end
            end
            ray(:,4:6) = [];   %remove kr quad and r from the ray table
        end
%%

        function ray = meshRay(obj,cmesh,xys,alpha,hlimit,tol)
            %get ray path using triangular mesh held in cmesh (renamed from
            %cgrid)

            %find nearest node to start point
            kr = nearestNeighbor(cmesh.Tri,xys);

            %wave properties at start point
            tol.radius = 1e6;
            r = inf;                                     %intial radius
            xr = xys(1); yr = xys(2); quad = 0;
            [hr,cr,cgr] = raypoint_properties(obj,cmesh,xr,yr);
            %initialise ray table            
            ray = table(xr,yr,alpha,kr,quad,r,hr,cr,cgr);%grid properties of ray position

            %check plot for finding ray errors - comment out when not required
%             hf = figure('Name','Search','Tag','PlotFig');
%             ax = axes(hf);
%             Tri = cmesh.Tri;
%             save('trigrid_input','Tri')
%             pts = cmesh.Tri.Points;
%             tria = cmesh.Tri.ConnectivityList;
%             trimesh(tria,pts(:,1),pts(:,2),'Color','k')
%             hold on
%             plot(ax,ray.xr(1),ray.yr(1),'ob')
%             hold off
            %--------------------------------------------------------------
            
            while hr>hlimit
                newray = mesh_arc_ray(cmesh,ray(end,:),tol);
                if isempty(newray)                         %ray exits grid
                    obj.outFlag = 1; hr = hlimit; continue; 
                elseif isnumeric(newray) && newray<0       %ray error (eg radius too small)
                    obj.outFlag = newray; hr = hlimit; continue; 
                end
                ray = [ray;newray]; %#ok<AGROW> 
                %check plot for finding ray errors - comment out when not required
%                 hold on
%                 plot(ax,newray.xr,newray.yr,'+r')
%                 hold off   
                %----------------------------------------------------------                
                hr = newray.hr;
                if hr<=hlimit, obj.outFlag = -1; end
            end   
            ray(:,4:6) = [];   %remove kr quad and r from the ray table
        end
%%
        function ray = startRay(obj,cgrid,xyr,alpha,tol)
            %compute the first element intersection from the start point
            % grid - contains X,Y,h,c,dcx,dcy
            % xyr - start point grid coordinates
            % alpha - start angle (radians)
            
            delta = cgrid.X(1,2)-cgrid.X(1,1);                  %grid spacing
            XY = [reshape(cgrid.X,[],1),reshape(cgrid.Y,[],1)]; %x,y vectors

            %find nearest node to start point
            kr = dsearchn(XY,xyr);
            %wave properties at start point
            [hr,cr,cgr] = raypoint_properties(obj,cgrid,xyr(1),xyr(2));
            %initialise ray table
            r = inf;                                     %intial radius
            xr = xyr(1); yr = xyr(2); quad = 0;
            ray = table(xr,yr,alpha,kr,quad,r,hr,cr,cgr);%grid properties of ray position

            %coordinates of local origin
            [row,col] = ind2sub(size(cgrid.X),kr);
            xi = cgrid.X(row,col); yi = cgrid.Y(row,col); 
            %get vector to start point in local grid coordinates
            uvr = [(xr-xi)/delta,(yr-yi)/delta];
            %find quadrant for start point
            ison = is_axis_point(ray,uvr,tol);
            ray.quad = get_quadrant(ray,uvr,ison);
            
            %find first intersection point
            [ue,ve] = pol2cart(alpha,sqrt(2));      %vector from start in direction of alpha
                                                    %sqrt(2) ensures it crosses a boundary
            [u0,v0] = pol2cart(alpha,tol.dist);     %offset from start point
            us = uvr(1); vs = uvr(2);               %start point
            lineseg = [us+u0,vs+v0;us+ue,vs+ve];    %ray vector line segment
            r = inf;                                %intial radius

            % figure; ax = axes;
            % plot(ax,xi,yi,'ob')  %origin
            % hold on
            % plot(ax,xr,yr,'xr')
            % plot(ax,xi+lineseg(:,1)*delta,yi+lineseg(:,2)*delta,'--k')
            % hold off
            %check line does not go out of grid from start point                                   
            xye = [xyr(1)+ue*delta,xyr(2)+ve*delta];
            isbound = checkGridBoundary(obj,cgrid,xye);       
            if isbound, ray = NaN; return; end

            %get new ray point where it exits current element
            uvray = get_intersection(ray,lineseg,uvr,tol);
            xr = xi+uvray(1)*delta; yr = yi+uvray(2)*delta;
            kr = dsearchn(XY,[xr,yr]);            
            %wave properties at new ray point
            [hr,cr,cgr] = raypoint_properties(obj,cgrid,xr,yr);      
            %coordinates of start point
            [row,col] = ind2sub(size(cgrid.X),kr);
            xi = cgrid.X(row,col); yi = cgrid.Y(row,col); 
            %get vector to start point in local grid coordinates
            uvr = [(xr-xi)/delta,(yr-yi)/delta];

            % hold on
            % plot(ax,xr,yr,'+g')
            % plot(ax,xi,yi,'vb')
            % hold off
            %find quadrant for new ray point
            ison = is_axis_point(ray,uvr,tol);
            quad = get_quadrant(ray,uvr,ison);
            %add new ray position to table
            newray = table(xr,yr,alpha,kr,quad,r,hr,cr,cgr);%grid properties of ray position
            ray = [ray;newray];
        end
%%
        function edge = get_edge(~,ison)
            %assign and edge when directed along axis
            if ison(2)==1 || ison(2)==3
                edge = 2;
            else 
                edge = 1;
            end
        end
%%
        function [hr,cr,cgr] = raypoint_properties(~,cgrid,xr,yr)
            %interpolate depth, celerity and group celerity at ray point
            if isfield(cgrid,'Tri')
                X = cgrid.Tri.Points(:,1); Y = cgrid.Tri.Points(:,2); 
                hr = griddata(X,Y,cgrid.h,xr,yr);
                cr = griddata(X,Y,cgrid.c,xr,yr);
                cgr = griddata(X,Y,cgrid.cg,xr,yr);                
            else
                X = cgrid.X; Y = cgrid.Y; 
                hr = interp2(X,Y,cgrid.h',xr,yr);
                cr = interp2(X,Y,cgrid.c',xr,yr);
                cgr = interp2(X,Y,cgrid.cg',xr,yr);
            end
        end
%%
        function isbound = checkGridBoundary(~,cgrid,xye)
            %check if point, xye, is outside grid
            x = cgrid.X(1,:); y = cgrid.Y(:,1);
            limxy = [x(1),y(1);x(end),y(end)];                     %limits of grid
            isbound = xye(1)<limxy(1,1) || xye(1)>limxy(2,1) || ...%check x
                      xye(2)<limxy(1,2) || xye(2)>limxy(2,2);      %check y       
        end
%%
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                       % <<Edit metadata to suit model
                'Name',{'xr','yr','alpha','depth','celerity','cgroup'},...
                'Description',{'X-Position','Y-Position','Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'Unit',{'m','m','rad','m','m/s','m/s'},...
                'Label',{'X-Position','Y-Position','Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'QCflag',repmat({'model'},1,6)); 
            dsp.Row = struct(...
                'Name',{'-'},...
                'Description',{'-'},...
                'Unit',{'-'},...
                'Label',{'-'},...
                'Format',{'-'});        
            dsp.Dimensions = struct(...    
                'Name',{'-'},...
                'Description',{'-'},...
                'Unit',{'-'},...
                'Label',{'-'},...
                'Format',{'-'});  
        end  
    end
end