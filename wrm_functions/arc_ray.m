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
%         node index, k, quadrant being entered, quad, side of element, edge
%   tol - tolerance to test for angles that are multiples of pi/2 (5.7 deg)
% OUTPUTS
%   newray - table of outgoing ray position (xr,yr), direction, alpha, local
%            node index, k, quadrant being entered, quad, side of element, edge
% NOTES
%   Ray method based on Abernethy C L and Gilbert G, 1975, Refraction of 
%   wave spectra, Report No: INT 117,pp. 1-166, Wallingford, UK.
%   Right angled isosceles triangles derived from uniform grid are used 
%   rather than equalateral triangles. Starting from the nearest grid point
%   to the start point, a local grid is constructed comprising a central
%   node (0,0) and the surrounding nodes (1,0),(1,1),(0,1),(-1,1),(-1,0)
%   (-1,-1),(0,-1),(1,-1). The position of the ray at any given time
%   relative to the grid is defined in grid coordinates as xr,yr and in
%   indices as the nearest node, k (or i,j), the quadrant, q, and the edge
%   of the triangle, e. Where q is the trigonometric quadrant 1-4 of the 
%   ray point, and e has a value of 1-3: e=1 if yr=yi; e=2 if xr=zi; 
%   else e=3.
% SEE ALSO
%   get_quadrant, get_element and next_element 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    X = cgrid.X; Y = cgrid.Y;
    delta = X(1,2)-X(1,1);                  %grid spacing

    %variables used in function
    % alpha - angle of ray direction
    % phi - angle of normal to ray direction
    % theta - angle from origin to ray position
    % k - index of reference grid node
    % xi,yi - position of reference grid node for local coordinates
    % xr,yr - ray position in grid coordinates
    % ur,vr - ray position in local coordinates
    % uc,vc - centre of arc in local coordinates
    % r - radius of arc in local coordinates (R = r*delta)

    %find node of ray point
    [row,col] = ind2sub(size(X),ray.k);
    xi = X(row,col); yi = Y(row,col);       %coordinates of start point

    %get vector to start point in local grid coordinates
    ur = (ray.xr-xi)/delta;
    vr = (ray.yr-yi)/delta;

    %find which element the ray is entering
    dcx = interp2(cgrid.X,cgrid.Y,cgrid.dcx',ray.xr,ray.yr,'linear',0); %gradients at start point
    dcy = interp2(cgrid.X,cgrid.Y,cgrid.dcy',ray.xr,ray.yr,'linear',0);
    [quad,edge,uvi] = next_element(ray,[ur,vr],tol);
    %transform local coordinates if start point is on edge 3 (hypotenuse)
    if ray.edge==3 || sum(abs(uvi))>0
        xi = xi+uvi(1)*delta;       %translate origin to new local origin
        yi = yi+uvi(2)*delta;
        isoutbound = checkGridBoundary(cgrid,[xi,yi]);
        if isoutbound, newray = []; return; end
        k = sub2ind(size(X),row+uvi(2),col+uvi(1));
        ur = (ray.xr-xi)/delta;     %update ray position to new origin
        vr = (ray.yr-yi)/delta;
    else
        k = ray.k;                  %no change in origin
    end

    if quad>4                       %assign values for ray along an axis
        xr = xi+uvi(1)*delta; yr = yi+uvi(2)*delta;
        alpha = ray.alpha;
        r = inf;
    else
        %get the centre of the arc that is tangential to the ray at ur,vr
        [phi,r,uc,vc] = arc_properties(cgrid,ur,vr,ray);
        %create line vector based on defined arc
        xyArc = get_arc(phi,r,uc,vc,ur,vr);
        if any(isnan(xyArc))
            Tri = get_element(quad);
            plot_element(Tri,xyArc,ur,vr,k)
            error('Arc segment not found at xr=%.0f; yr=%.0f',ray.xr,ray.yr)
        elseif isempty(xyArc)
            warndlg(sprintf('Radius too small at xr=%.0f; yr=%.0f',ray.xr,ray.yr));
            newray = []; return;
        end
        %find the intersection of the arc segment with quad triangle
        [uvray,edge] = get_intersection(quad,xyArc,ray.alpha,[ur,vr],tol);
        xr = xi+uvray(1)*delta; yr = yi+uvray(2)*delta;
        %get exit angle
        alpha =  exit_angle(phi,r,ur,vr,uvray);  
    end

    isoutbound = checkGridBoundary(cgrid,[xr,yr]);
    if isoutbound 
        newray = []; return; 
    end

    [hr,cr,cgr] = raypoint_properties(cgrid,xr,yr);
    try
    newray = table(xr,yr,alpha,k,quad,edge,r,hr,cr,cgr);%grid properties of ray position
    catch
        warndlg()
    end
end
%%
function [xy_arc] = get_arc(phi,radius,uc,vc,ur,vr)  
    %calculate the coordinates of an arc either side of the radius vecor
    %from uc,vc to ur,vr.
    N = 1001;                                    %number of points in Arc
    if abs(radius)>1000 
        %straight line segment will suffice        
        [ue,ve] = pol2cart(phi-pi/2,sqrt(2));    %vector from ray point in direction of alpha
                                                 %sqrt(2) ensures it crosses a boundary
        xy_arc = [ur-ue,vr-ve;ur,vr;ur+ue,vr+ve];%ray vector line segment
        return;
    elseif abs(radius)<0.001
        xy_arc = []; return;
    end
    % arcang = 2*asin(1/radius);                 %angle between entry and exit point
    % arcang = 2*asin(0.5/radius);
    arcang = pi/2;
    phi = phi+pi;                                %angle of normal from centre of arc

    r_angl = linspace(phi+arcang,phi-arcang, N); %angles Defining Arc Segment (radians)
    %abstract Circle Function For Angles In Radians
    circr = @(radius,angle)  [radius*cos(angle)+uc,  radius*sin(angle)+vc];     
    xy_arc = circr(radius,r_angl');              %matrix (Nx2) of (x,y) coordinates

    %Arc segment use
    % Arc = polyshape([xy_arc(:,1)],[xy_arc(:,2)]);
    %For Arc sector use 
    % Arc = polyshape([uc;xy_arc(:,1)],[vc;xy_arc(:,2)]);
end
%%
function [phi,r,uc,vc] = arc_properties(cgrid,ur,vr,ray)
    %find the radius and centre of the arc that the ray follows in element
    X = cgrid.X; Y = cgrid.Y; c = cgrid.c; dcx = cgrid.dcx; dcy = cgrid.dcy;
    delta = X(1,2)-X(1,1);             %grid spacing
    cr = interp2(X,Y,c',ray.xr,ray.yr);
    dcrx = interp2(X,Y,dcx',ray.xr,ray.yr,'linear',0);
    dcry = interp2(X,Y,dcy',ray.xr,ray.yr,'linear',0);

    %unit normal to ray (convention is left in direction of ray is positive)
    phi = mod(ray.alpha+pi/2,2*pi);    %angle of normal to ray direction    
    [un,vn] = pol2cart(phi,1);         %normal vector in local coordinates

    %radius of arc
    Ndc = un*dcrx+vn*dcry;
    R = -cr/Ndc;                       %radius in grid coordinates
    r = R/delta;                       %radius in local coordinates
    if abs(r)<1000
        uc = ur+r*un;
        vc = vr+r*vn;
    else
        uc = 0; vc = 0;   %radius is large so use straight line segment
    end
end
%%
function [hr,cr,cgr] = raypoint_properties(cgrid,xr,yr)
    %interpolate depth, celerity and group celerity at ray point
    X = cgrid.X; Y = cgrid.Y; 
    hr = interp2(X,Y,cgrid.h',xr,yr);
    cr = interp2(X,Y,cgrid.c',xr,yr);
    cgr = interp2(X,Y,cgrid.cg',xr,yr);
end
%%
function alpha = exit_angle(phi,r,ur,vr,uvray)
    %find the angle that is tangent to the arc at the exit point
    % phi,r - angle of normal and radius of arc
    % uv,vr - ray entry point into element
    % uvray - ray exit point out of element
    phi = mod(phi+pi,2*pi);                    %angle of normal from centre of arc
    L = sqrt((ur-uvray(1))^2+(vr-uvray(2))^2); %arc segment length
    if abs(r)<1000
        theta = 2*asin(L/2/r);                     %angle between entry and exit point
    else
        theta = 0;
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
              xye(2)<limxy(1,2) || xye(2)>limxy(2,2);      %check y       
end
%%
% function isintol = isclose(angle,test,tol)
%     %check whether test is within tol of angle
%     isintol = abs(angle-test)<tol.angle;
% end
%%
% function [isdir,npt] = checkDirection(lineseg,ur,vr,alpha)
%     %check the direction of a vector relative to direction given by alpha
%     % npt is nearest point to [ur,vr] and this is in the direction of alpha
%     % if isdir is true
%     anglim = 0.2;                         %angle limit of +/-0.2 rads (11.5 deg)
%                                           %could pass tol.angle
%     npt = dsearchn(lineseg,[ur,vr]);
%     un = lineseg(npt,1)-ur;   
%     vn = lineseg(npt,2)-vr; 
%     [xsi,~] = cart2pol(un,vn);            %vector to next point
%     bound = [alpha-anglim,alpha+anglim];
%     isdir = isangletol(xsi,bound);        %check if xsi is within bound
% end

%%
% function Circ = get_circle(radius,uc,vc)   
%     %define a circle in local coordinates
%     n = 50;                                  %number of points in circle
%     phi = (0:n-1)*(2*pi/n);                  %angle increments
%     u = uc + radius*cos(phi);                %coordinates of circle
%     v = vc + radius*sin(phi);
%     Circ = polyshape(u,v);
% end

