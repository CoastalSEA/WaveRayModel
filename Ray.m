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
        Track     %table of track properties: xr,yr,alpha,k,quad,edge  
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

            %get first first element intersection from the start point
            ray = startRay(obj,cgrid,xys,alpha,tol);
            if isempty(ray)
                error('Ray solution not found in start_ray')
            elseif ~isa(ray,'table')
                error('Ray start position outside grid domain')
            end   

            %check plot for finding ray errors - comment out when not required
            % hf = figure('Name','Search','Tag','PlotFig');
            % ax = axes(hf);
            % plot(ax,ray.xr(1),ray.yr(1),'ob')
            %--------------------------------------------------------------

            %loop to get ray track to edge of grid or depth limit
            hr = interp2(cgrid.X,cgrid.Y,cgrid.h',ray.xr,ray.yr,'linear',0);
            while hr>hlimit
                newray = arc_ray(cgrid,ray(end,:),tol);
                if isempty(newray), hr = hlimit; continue; end
                ray = [ray;newray]; %#ok<AGROW> 
                %check plot for finding ray errors - comment out when not required
                % hold on
                % plot(ax,newray.xr,newray.yr,'+k')
                % hold off   
                %----------------------------------------------------------
                hr = interp2(cgrid.X,cgrid.Y,cgrid.h',ray.xr,ray.yr,'linear',0);
            end
            ray(:,4:7) = [];   %remove k, quad and edge from the ray table 888888
            obj.Track = dstable(ray,'DSproperties',modelDSproperties(obj));
        end
    end
%%
    methods (Access=private)  
        function ray = startRay(obj,cgrid,xys,alpha,tol)
            %compute the first element intersection from the start point
            % grid - contains X,Y,h,c,dcx,dcy
            % xys - start point grid coordinates
            % alpha - start angle (radians)
            
            delta = cgrid.X(1,2)-cgrid.X(1,1);                  %grid spacing
            XY = [reshape(cgrid.X,[],1),reshape(cgrid.Y,[],1)]; %x,y vectors

            %find nearest node to start point
            k = dsearchn(XY,xys);
            [row,col] = ind2sub(size(cgrid.X),k);
            xi = cgrid.X(row,col); yi = cgrid.Y(row,col);       %coordinates of start point
            
            %get vector to start point in local grid coordinates
            us = (xys(1)-xi)/delta;
            vs = (xys(2)-yi)/delta;

            [ue,ve] = pol2cart(alpha,2);      %vector from start in direction of alpha
                                                    %sqrt(2) ensures it crosses a boundary
            lineseg = [us-ue,vs-ve;us+ue,vs+ve];    %ray vector line segment
            r = inf;                                %intial radius

            %check line does not go out of grid from start point                                   
            xye = [xys(1)+ue*delta,xys(2)+ve*delta];
            isbound = checkGridBoundary(obj,cgrid,xye);       
            if isbound, ray = NaN; return; end

            %check whether point is on grid axes and determine quad
            ray = table(alpha,k,r);
            [ison,uvi] = is_axis_point(ray,[us,vs],tol);
            quad = get_quadrant(ray,[us,vs],ison);

            if ison(1)<1 || ison(1)>2
                %not on axis (in element or on a diagonal)
                [uvr,edge] = get_intersection(quad,lineseg,ray,[us,vs],tol);
                xr = xi+uvr(1)*delta; yr = yi+uvr(2)*delta;
            else
                %on an axis and directed along axis
                edge = get_edge(obj,ison);
                xr = xi+uvi(1)*delta; yr = yi+uvi(2)*delta;
            end            

            [hr,cr,cgr] = raypoint_properties(obj,cgrid,xr,yr);
            ray = table(xr,yr,alpha,k,quad,edge,r,hr,cr,cgr);%grid properties of ray position
            if xr~=xys(1) || yr~=xys(2)
                %wave properties at start point
                [hs,cs,cgs] = raypoint_properties(obj,cgrid,xys(1),xys(2));
                %duplicate row and add start point to first row
                ray = [ray;ray];
                ray{1,1} = xys(1); ray{1,2} = xys(2); 
                ray{1,8} = hs; ray{1,9} = cs; ray{1,10} = cgs; 
            end
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
            X = cgrid.X; Y = cgrid.Y; 
            hr = interp2(X,Y,cgrid.h',xr,yr);
            cr = interp2(X,Y,cgrid.c',xr,yr);
            cgr = interp2(X,Y,cgrid.cg',xr,yr);
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