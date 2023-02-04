classdef Ray < handle                        
%
%-------class help---------------------------------------------------------
% NAME
%   Ray.m
% PURPOSE
%   Class description - Class for ray objects in the WaveRayModel app
% NOTES
%   This version only computes the route track. 
%   The table ray is used to hold ray properties. The last row of the 
%   table is used to call arc_ray. Whilst constructing the ray the table
%   holds xr, yr, alpha, k, quad, edge, hr, cr, cgr. Once the ray is
%   complete k, quad, edge are removed from the table before assigning the
%   reultant ray table to the Track property.
% SEE ALSO
%   WaveRayModel.m, RayTrack.m, and arc_ray.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
%     
    properties
        Track     %table of track properties: xr,yr,alpha,k,quad,edge  
    end
    
    methods
        function obj = Ray()                
            %class constructor
        end
    end      
%%
    methods (Static)                
        function obj = setRay(grid,xys,alpha,hlimit)
            %function to compute a single ray track
            % grid - contains X,Y,h,c,dcx,dcy
            % xys - start point grid coordinates
            % alpha - start angle (radians)
            % hlimit - cutoff depth limit for ray (optional)
            if nargin<4
                hlimit = 0.1;
            end

            obj = Ray;                               %Ray class instance 

            %get first first element intersection from the start point
            ray = startRay(obj,grid,xys,alpha);
            if isempty(ray)
                error('Ray solution not found in start_ray')
            elseif ~isa(ray,'table')
                error('Ray start positiion outside grid domain')
            end   

            %loop to get ray track to edge of grid or depth limit
            hr = interp2(grid.X,grid.Y,grid.h',ray.xr,ray.yr,'linear',0);
            while hr>hlimit
                newray = arc_ray(grid,ray(end,:));
                ray = [ray;newray]; %#ok<AGROW> 
                hr = interp2(grid.X,grid.Y,grid.h',ray.xr,ray.yr,'linear',0);
            end
            ray(:,4:6) = [];   %remove k, quad and edge from the ray table
            obj.Track = dstable(ray,'DSproperties',modelDSproperties(obj));
        end
    end
%%
    methods (Access=private)  
        function ray = startRay(obj,grid,xys,alpha)
            %compute the first element intersection from the start point
            % grid - contains X,Y,h,c,dcx,dcy
            % xys - start point grid coordinates
            % alpha - start angle (radians)
            
            delta = grid.X(1,2)-grid.X(1,1);                  %grid spacing
            XY = [reshape(grid.X,[],1),reshape(grid.Y,[],1)]; %x,y vectors

            %find nearest node to start point
            k = dsearchn(XY,[xys(1),xys(2)]);
            [row,col] = ind2sub(size(grid.X),k);
            xi = grid.X(row,col); yi = grid.Y(row,col);       %coordinates of start point
            
            %get vector to start point in local grid coordinates
            us = (xys(1)-xi)/delta;
            vs = (xys(2)-yi)/delta;
            [theta,~] = cart2pol(us,vs);            %vector to start point
            theta= mod(theta,2*pi);
            
            [ue,ve] = pol2cart(alpha,sqrt(2));      %vector from start in direction of alpha
                                                    %sqrt(2) ensures it crosses a boundary
            lineseg = [us,vs;us+ue,vs+ve];          %ray vector line segment

            %check line does not go out of grid from start point                                   
            xye = [xys(1)+ue*delta,xys(2)+ve*delta];
            isbound = checkGridBoundary(obj,grid,xye);       
            if isbound, ray = NaN; return; end

            %define a triangle polygon for the quadrant of the start point
            [Tri,quad] = get_quadrant(theta);
            if isempty(Tri), ray = []; return; end

            % figure; plot(Tri)
            % hold on
            % plot(us,vs,'og')
            % plot(lineseg(:,1),lineseg(:,2))
            % hold off
            
            %find intersection point of line with triangle
            [inside,outside] = intersect(Tri,lineseg); 
            %inside line segment coordinates returned as a two-column matrix (x,y).
            %find the common point in both vectors
            if isempty(inside)
                tol = pi/1000;
                [Tri,quad] = get_quadrant(mod(theta-tol,2*pi));
                [inside,outside] = intersect(Tri,lineseg); 
            end
            [~,idx] = intersect(inside,outside,'rows');
            %use coordinates of point to identify which edge it lies on
            tol = 1/1000/delta;                     %tolerance equivalent to 1mm
            edge = get_edge(inside,idx,tol);
        
            %transform new ray position from local to grid coordinates
            xr = xi+inside(idx,1)*delta; yr = yi+inside(idx,2)*delta;
            [hr,cr,cgr] = raypoint_properties(obj,grid,xr,yr);
            ray = table(xr,yr,alpha,k,quad,edge,hr,cr,cgr);%grid properties of ray position
            %wave properties at start point
            [hs,cs,cgs] = raypoint_properties(obj,grid,xys(1),xys(2));
            %duplicate row and add start point to first row
            ray = [ray;ray];
            ray{1,1} = xys(1); ray{1,2} = xys(2); 
            ray{1,7} = hs; ray{1,8} = cs; ray{1,9} = cgs; 
        end
%%
        function [hr,cr,cgr] = raypoint_properties(~,grid,xr,yr)
            %interpolate depth, celerity and group celerity at ray point
            X = grid.X; Y = grid.Y; 
            hr = interp2(X,Y,grid.h',xr,yr);
            cr = interp2(X,Y,grid.c',xr,yr);
            cgr = interp2(X,Y,grid.cg',xr,yr);
        end
%%
        function isbound = checkGridBoundary(~,grid,xye)
            %check if point, xye, is outside grid
            x = grid.X(1,:); y = grid.Y(:,1);
            limxy = [x(1),y(1);x(end),y(end)];                     %limits of grid
            isbound = xye(1)<limxy(1,1) || xye(1)>limxy(2,1) || ...%check x
                      xye(2)<limxy(1,2) || xye(2)>limxy(2,2);      %check y       
        end
%%
%        function obj = setTable(obj,table)
%             %make the table contents of a ray a single row
%             dst = table;
%             dst(2:end,:) = [];
%             varnames = table.Properties.VariableNames;
%             for i=1:width(table)
%                 dst.(varnames{i}) = table{:,i}';
%             end
%             dsp = setDSproperties(obj);
%             obj.Track = dstable(dst,'DSproperties',dsp);
%        end
%%
%         function dsp = modelDSproperties(~) 
%             %define a dsproperties struct and add the model metadata
%             dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
%             %define each variable to be included in the data table and any
%             %information about the dimensions. dstable Row and Dimensions can
%             %accept most data types but the values in each vector must be unique
%             
%             %struct entries are cell arrays and can be column or row vectors
%             dsp.Variables = struct(...                       % <<Edit metadata to suit model
%                 'Name',{'xr','yr','alpha','k','quad','edge'},...
%                 'Description',{'X-Position','Y-Position','Direction',...
%                            'Grid index','Grid quadrant','Element edge'},...
%                 'Unit',{'m','m','rad','-','-','-'},...
%                 'Label',{'X-Position','Y-Position','Direction',...
%                            'Grid index','Grid quadrant','Element edge'},...
%                 'QCflag',repmat({'model'},1,6)); 
%             dsp.Row = struct(...
%                 'Name',{'-'},...
%                 'Description',{'-'},...
%                 'Unit',{'-'},...
%                 'Label',{'-'},...
%                 'Format',{'-'});        
%             dsp.Dimensions = struct(...    
%                 'Name',{'-'},...
%                 'Description',{'-'},...
%                 'Unit',{'-'},...
%                 'Label',{'-'},...
%                 'Format',{'-'});  
%         end     
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