function newray = arc_ray(cgrid,ray,tol)
%
%-------function help------------------------------------------------------
% NAME
%   arc_ray.m
% PURPOSE
%   compute the exit position and direction of a ray entering a triangular
%   element at the position and direction defined by the incoming ray.
% USAGE
%   newray = arc_ray(cgrid,ray,tol);
% INPUTS
%   cgrid - struct of grid array properties, including:
%          X, Y - x and y coordinates as meshgrid matrices
%          h - water depths
%          c - celerity grid matrix 
%          dcx, dcy - gradients of celerity in x and y directions
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad
%   tol - tolerance to test for angles that are multiples of pi/2 (5.7 deg)
% OUTPUTS
%   newray - table of outgoing ray position (xr,yr), direction, alpha, local
%            node index, k, quadrant being entered, quad, side of element, edge
%            - returned empty if new point is outside grid domain
%            - returned as -2 if radius is too small (<0.5) 
% NOTES
%   Ray method based on Abernethy C L and Gilbert G, 1975, Refraction of 
%   wave spectra, Report No: INT 117,pp. 1-166, Wallingford, UK.
%   Right angled isosceles triangles derived from uniform grid are used 
%   rather than equalateral triangles. Starting from the nearest grid point
%   to the start point, a local grid is constructed comprising a central
%   node (0,0) and the surrounding nodes (1,0),(1,1),(0,1),(-1,1),(-1,0)
%   (-1,-1),(0,-1),(1,-1). The position of the ray at any given time
%   relative to the grid is defined in grid coordinates as xr,yr and in
%   indices as the nearest node, kr (or i,j), the quadrant, quad, 
%   of the triangle. Where quad is the trigonometric quadrant 1-4 of the 
%   ray point. 
% SEE ALSO
%   get_quadrant, get_element, is_axis_point, get_intersection
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    X = cgrid.X; Y = cgrid.Y;
    delta = X(1,2)-X(1,1);                  %grid spacing
    XY = [reshape(cgrid.X,[],1),reshape(cgrid.Y,[],1)]; %x,y vectors

    %variables used in function
    % alpha - angle of ray direction
    % phi - angle of normal to ray direction
    % theta - angle from origin to ray position
    % kr - index of reference grid node
    % xi,yi - position of reference grid node for local coordinates
    % xr,yr - ray position in grid coordinates
    % ur,vr - ray position in local coordinates
    % uc,vc - centre of arc in local coordinates
    % r - radius of arc in local coordinates (R = r*delta)
    
    %find node of ray point
    [row,col] = ind2sub(size(X),ray.kr);
    xi = X(row,col); yi = Y(row,col);       %coordinates of start point

    %get vector to start point in local grid coordinates
    ur = (ray.xr-xi)/delta;
    vr = (ray.yr-yi)/delta;
    
    %get the centre of the arc that is tangential to the ray at ur,vr
    [phi,r,uc,vc] = arc_properties(cgrid,ur,vr,ray,tol);
    uvArc = get_arc(phi,r,uc,vc,ur,vr,tol);
    if any(isnan(uvArc))
        Tri = get_element(ray.quad);
        plot_element(Tri,uvArc,ur,vr,ray.kr)
        error('Arc segment not found at xr=%.0f; yr=%.0f',ray.xr,ray.yr)
    elseif isempty(uvArc)
        % warndlg(sprintf('Radius, r=%0.2f too small at xr=%.0f; yr=%.0f',r,ray.xr,ray.yr));
        newray = -2; return;
    end
    %find the intersection of the arc segment with quad triangle
    uvray = get_intersection(ray,uvArc,[ur,vr],tol);
    xr = xi+uvray(1)*delta; yr = yi+uvray(2)*delta;
    kr = dsearchn(XY,[xr,yr]);
    %wave properties at new ray point
    [hr,cr,cgr] = raypoint_properties(cgrid,xr,yr); 
    %get exit angle (ur,vr entry point into element, uvray exit point)
    alpha =  exit_angle(phi,r,ur,vr,uvray,tol);  ray.alpha = alpha;

    %coordinates of start point
    [row,col] = ind2sub(size(cgrid.X),kr);
    xi = cgrid.X(row,col); yi = cgrid.Y(row,col); 
    %get vector to start point in local grid coordinates
    uvr = [(xr-xi)/delta,(yr-yi)/delta];
    %find quadrant for new ray point
    ison = is_axis_point(ray,uvr,tol);
    quad = get_quadrant(ray,uvr,ison);

    isoutbound = checkGridBoundary(cgrid,[xr,yr]);
    if isoutbound 
        newray = []; return; 
    else
        newray = table(xr,yr,alpha,kr,quad,r,hr,cr,cgr);
    end
end
%%
function [phi,r,uc,vc] = arc_properties(cgrid,ur,vr,ray,tol)
    %find the radius and centre of the arc that the ray follows in element
    X = cgrid.X; Y = cgrid.Y; c = cgrid.c; dcx = cgrid.dcx; dcy = cgrid.dcy;
    delta = X(1,2)-X(1,1);             %grid spacing
    method = 'linear';
    cr = interp2(X,Y,c',ray.xr,ray.yr,method);
    dcrx = interp2(X,Y,dcx',ray.xr,ray.yr,method,0);
    dcry = interp2(X,Y,dcy',ray.xr,ray.yr,method,0);

    %unit normal to ray (convention is left in direction of ray is positive)
    phi = mod(ray.alpha+pi/2,2*pi);    %angle of normal to ray direction    
    [un,vn] = pol2cart(phi,1);         %normal vector in local coordinates

    %radius of arc
    Ndc = un*dcrx+vn*dcry;
    R = -cr/Ndc;                       %radius in grid coordinates
    r = R/delta;                       %radius in local coordinates
    if abs(r)<tol.radius
        uc = ur+r*un;
        vc = vr+r*vn;
    else
        uc = 0; vc = 0;   %radius is large so use straight line segment
    end
end
%%
function [uv_arc] = get_arc(phi,radius,uc,vc,ur,vr,tol)  
    %calculate the coordinates of an arc either side of the radius vecor
    %from uc,vc to ur,vr.
    N = 20;                                       %number of points in half-Arc
    rt2 = sqrt(2);
    rad = abs(radius);
    alpha = mod(phi-pi/2,2*pi);
    
    if rad>=tol.radius
        %straight line segment will suffice       
        [ue,ve] = pol2cart(alpha,rt2); %vector from ray point in direction of alpha
                                       %use 2 to ensure line crosses a boundary
                                       %max element length is sqrt(2) 
        [us,vs] = pol2cart(phi-pi/2,tol.dist);         
        uv_arc = [ur+us,vr+vs;ur+ue,vr+ve];       %ray vector line segment
        return;
    elseif rad<1.0                                %radius is so small that it may not exit element
        uv_arc = []; return;
    end

    phi = phi+pi;                                  %angle of normal from centre of arc
    arcang = rt2/rad;                              %set arc segment based on radius
    %abstract Circle Function For Angles In Radians
    circr = @(radius,angle)  [radius*cos(angle)+uc,  radius*sin(angle)+vc]; 
    r_angl = linspace(phi-arcang,phi-tol.angle, N);%angles Defining left Arc Segment (radians)
    uv_arcl = circr(radius,r_angl');
    r_angr = linspace(phi+tol.angle,phi+arcang, N);%angles Defining right Arc Segment (radians)
    uv_arcr = circr(radius,r_angr');
    %find arc in direction of ray
    ua1 = uv_arcl(1,1)-ur; va1 = uv_arcl(1,2)-vr;  %vector from ray point to first point on arc
    arcl = mod(cart2pol(ua1,va1),2*pi);
    ua2 = uv_arcr(1,1)-ur; va2 = uv_arcr(1,2)-vr;  %vector from ray point to first point on arc
    arcr = mod(cart2pol(ua2,va2),2*pi);
    angtol = 1.0;                                  %angle limit of +/-1.0 rads (57.3 deg)
    bound = [alpha-angtol,alpha+angtol];
   
    if isangletol(arcl,bound)                      %check if arcl is within bound   
        uv_arc = uv_arcl;
    elseif isangletol(arcr,bound)                  %check if arcr is within bound   
        uv_arc = uv_arcr;
    else
        uv_arc = [];
        % warndlg('Arc not found in get_arc')
    end
end
%%
function [hr,cr,cgr] = raypoint_properties(cgrid,xr,yr)
    %interpolate depth, celerity and group celerity at ray point
    X = cgrid.X; Y = cgrid.Y; 
    method = 'linear';
    hr = interp2(X,Y,cgrid.h',xr,yr,method);
    cr = interp2(X,Y,cgrid.c',xr,yr,method);
    cgr = interp2(X,Y,cgrid.cg',xr,yr,method);
end
%%
function alpha = exit_angle(phi,r,ur,vr,uvray,tol)
    %find the angle that is tangent to the arc at the exit point
    % phi,r - angle of normal and radius of arc
    % uv,vr - ray entry point into element
    % uvray - ray exit point out of element
    phi = mod(phi+pi,2*pi);                    %angle of normal from centre of arc
    L = sqrt((ur-uvray(1))^2+(vr-uvray(2))^2); %arc segment length
    if abs(r)<tol.radius
        theta = 2*asin(L/2/r);                 %angle between entry and exit point
    else
        theta = 0;                             %straight line no change in direction
    end
    alpha = phi+theta+pi/2;                    %angle of tangent to ray at exit point
    alpha = mod(alpha,2*pi);
end
%%
function plot_element(Tri,Arc,ur,vr,k)
    %plot the element and arc
    figure('Name','ArcRay','Tag','PlotFig');
    plot(Tri);
    hold on 
    %plot(uc,vc,'or');
    plot(ur,vr,'xr');
    plot(Arc(:,1),Arc(:,2));
    hold off
    xlim([-10,10])
    ylim([-10,10])
    title(sprintf('Element at node index %d',k))
end
%%
function isoutbound = checkGridBoundary(cgrid,xye)
    %check if point, xye, is outside grid
    x = cgrid.X(1,:); y = cgrid.Y(:,1);
    limxy = [x(1),y(1);x(end),y(end)];                     %limits of grid
    isoutbound = xye(1)<limxy(1,1) || xye(1)>limxy(2,1) || ...%check x
              xye(2)<limxy(1,2) || xye(2)>limxy(2,2);         %check y       
end
