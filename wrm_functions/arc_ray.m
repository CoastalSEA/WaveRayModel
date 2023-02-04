function newray = arc_ray(grid,ray)
%
%-------function help------------------------------------------------------
% NAME
%   arc_ray.m
% PURPOSE
%   compute the exit position and direction of a ray entering a triangular
%   element at the position and direction defined by the incoming ray.
% USAGE
%   newray = arc_ray(X,Y,ray,c,dcx,dcy);
% INPUTS
%   grid - struct of grid array properties, including:
%          X, Y - x and y coordinates as meshgrid matrices
%          h - water depths
%          c - celerity grid matrix 
%          dcx, dcy - gradients of celerity in x and y directions
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad, side of element, edge
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
%   ray point, and e has a value of 1-3: e=1 if yr=yi; e-=2 if xr=zi; 
%   else e=3.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    X = grid.X; Y = grid.Y;
    delta = X(1,2)-X(1,1);                  %grid spacing

    %variables used in function
    % alpha - angle of ray direction
    % phi - angle of normal to ray direction
    % k - index of reference grid node
    % xi,yi - position of reference grid node for local coordinates
    % xr,yr - ray position in grid coordinates
    % ur,vr - ray position in local coordinates
    % uc,vc - centre of arc in local coordinates
    % r - radius of arc in local coordinates (R = r*delta)

    %find node of ray point
    [row,col] = ind2sub(size(X),ray.k);
    xi = X(row,col); yi = Y(row,col);       %coordinates of start point

    %find which element the ray is entering
    [uvi,quad] = next_element(ray.quad,ray.edge);
    %transform local coordinates if start point is on edge 3 (hypotenuse)
    if ray.edge==3
        xi = xi+uvi(1)*delta;  %translate origin to new local origin
        yi = yi+uvi(2)*delta;
        k = sub2ind(size(X),row+uvi(2),col+uvi(1));
    else
        k = ray.k;
    end

    %get vector to start point in local grid coordinates
    ur = (ray.xr-xi)/delta;
    vr = (ray.yr-yi)/delta;
%     Tri = get_quadrant(quad);
    Tri = get_element(quad);
    if isempty(Tri)
        error('Ray quadrant not found, or has changed, in arc_ray')
    end

    %get the centre of the arc that is tangential to the ray at ur,vr
    [phi,r,uc,vc] = arc_properties(grid,ur,vr,ray);
    %create line vector based on defined arc
    xyArc = get_arc(phi,r,uc,vc,ur,vr);
    if any(isnan(xyArc))
        plot_element(Tri,xyArc,ur,vr,k)
        error('Arc segment not found')
    end
    % plot_element(Tri,Arc,ur,vr,k)

    %find intersection point of line with triangle
    [inside,outside] = intersect(Tri,xyArc); 
    if isempty(inside)
        [theta,~] = cart2pol(ur,vr);            %vector to start point
        theta = mod(theta,2*pi);
        tol = pi/1000;
        [Tri,quad] = get_quadrant(mod(theta,2*pi));
        [inside,outside] = intersect(Tri,xyArc); 
    end

    %inside line segment coordinates returned as a two-column matrix (x,y).
    %find the common point in both vectors
    tol = 1/delta;                          %tolerance equivalent to 1m
    idx = ~ismembertol(inside,[ur,vr],tol,'ByRows',true);
    inside = inside(idx,:);                 %inside points excluding entry point
    [~,idx] = intersect(inside,outside,'rows');

    if size(idx,1)>1
        %more than one intersection of element
        %find nearest point in direction of travel
        ok = 0;
        nxtpnt = inside(idx,:);
        while ok<1            
            [isdir,npt] = checkDirection(nxtpnt,ur,vr,ray.alpha);
            if isdir
                idx = idx(npt); ok = 1;
            else
                nxtpnt(npt,:) = [];
                if isempty(nxtpnt)
                    %point not found
                    idx = []; ok = 1;
                end
            end
        end
    end

    if isempty(idx)
        plot_element(Tri,xyArc,ur,vr,k)
        error('Intersection with element not found in arc_ray')
    end
    uvray = inside(idx,:);                  %local coordinates of ray exit point

    %exit angle and edge 
    alpha =  exit_angle(phi,r,ur,vr,uvray);
    tol = 1/1000/delta;                     %tolerance equivalent to 1mm
    edge = get_edge(inside,idx,tol);

    %transform new ray position from local to grid coordinates
    xr = xi+uvray(1)*delta; 
    yr = yi+uvray(2)*delta;
    [hr,cr,cgr] = raypoint_properties(grid,xr,yr);
    newray = table(xr,yr,alpha,k,quad,edge,hr,cr,cgr);%grid properties of ray position
end
%%
function [xy_arc] = get_arc(phi,radius,uc,vc,ur,vr)  
    %calculate the coordinates of an arc either side of the radius vecor
    %from uc,vc to ur,vr.
    N = 101;                                   %number of points in Arc
    if abs(radius)>1000
        %straight line segment will suffice        
        [ue,ve] = pol2cart(phi-pi/2,sqrt(2));  %vector from start in direction of alpha
                                               %sqrt(2) ensures it crosses a boundary
        xy_arc = [ur,vr;ur+ue,vr+ve];          %ray vector line segment
        return;
    end
%     arcang = 2*asin(1/radius);
%     arcang = 2*asin(0.5/radius);               %angle between entry and exit point
    arcang = pi/2;
    phi = phi+pi;                              %angle of normal from centre of arc

    r_angl = linspace(phi+arcang,phi-arcang, N); %angles Defining Arc Segment (radians)
    %abstract Circle Function For Angles In Radians
    circr = @(radius,angle)  [radius*cos(angle)+uc,  radius*sin(angle)+vc];     
    xy_arc = circr(radius,r_angl');            %matrix (Nx2) of (x,y) coordinates

%     [isdir,npt] = checkDirection(xy_arc,ur,vr,phi-pi/2);

%     npt = dsearchn(xy_arc,[ur,vr]);
%     ut = xy_arc(npt+1,1)-ur; vt = xy_arc(npt+1,2)-vr;
% 
%     [xsi,~] = cart2pol(ut,vt);            %vector to next point
%     xsi = mod(xsi,2*pi);
% 
%     isdir = xsi<phi && xsi>(mod(phi-pi,2*pi)); %equivalent to alpha +/-pi/2

%     if isdir
%         idr = 1:npt;
%     else
%         idr = npt:size(xy_arc,1);
%     end
%     xy_arc = xy_arc(idr,:);

    %Arc segment use
    % Arc = polyshape([xy_arc(:,1)],[xy_arc(:,2)]);
    %For Arc sector use 
    % Arc = polyshape([uc;xy_arc(:,1)],[vc;xy_arc(:,2)]);
end
%%
function [phi,r,uc,vc] = arc_properties(grid,ur,vr,ray)
    %find the radius and centre of the arc that the ray follows in element
    X = grid.X; Y = grid.Y; c = grid.c; dcx = grid.dcx; dcy = grid.dcy;
    delta = X(1,2)-X(1,1);             %grid spacing
    cr = interp2(X,Y,c',ray.xr,ray.yr);
    dcrx = interp2(X,Y,dcx',ray.xr,ray.yr,'linear',0);
    dcry = interp2(X,Y,dcy',ray.xr,ray.yr,'linear',0);

    %unit normal to ray (convention is left in direction of ray is positive)
    phi = ray.alpha+pi/2;              %angle of normal to ray direction
    offset = 0.001;                    %offset in radians (~0.06 deg)    
    if any(isclose(ray.alpha,[0,pi,2*pi],offset))        
        phi = phi+sign(dcry)*offset;
    end
    [un,vn] = pol2cart(phi,1);         %normal vector in local coordinates

    %radius of arc
    Ndc = un*dcrx+vn*dcry;
    R = -cr/Ndc;                       %radius in grid coordinates
    r = R/delta;                       %radius in local coordinates
    uc = ur+r*un;
    vc = vr+r*vn;
    [ucheck,vcheck] = pol2cart(ray.alpha-pi/2,r);
    tol = eps(R);
    if abs(ur-(uc+ucheck))>tol || abs(vr-(vc+vcheck))>tol
        fprintf('Arc centre coordinates are not within %.2g in arc_properties\n',tol)
    end
end
%%
function [hr,cr,cgr] = raypoint_properties(grid,xr,yr)
    %interpolate depth, celerity and group celerity at ray point
    X = grid.X; Y = grid.Y; 
    hr = interp2(X,Y,grid.h',xr,yr);
    cr = interp2(X,Y,grid.c',xr,yr);
    cgr = interp2(X,Y,grid.cg',xr,yr);
end
%%
function Tri = get_element(quad)
    %define new triangular polyshape based on quadrant being entered
    if quad==1                              %first quadrant
        Tri = polyshape([0,0,1],[0,1,0]);
    elseif quad==2                          %second quadrant
        Tri = polyshape([0,0,-1],[0,1,0]);
    elseif quad==3                          %third quadrant
        Tri = polyshape([0,0,-1],[0,-1,0]);
    elseif quad==4                          %fourth quadrant
        Tri = polyshape([0,0,1],[0,-1,0]);
    else  
        %quadrant not found
        Tri = [];
        return;
    end
end
%%
function [uvi,newquad] = next_element(quad,edge)
    %find the grid definition for the element that the ray is entering
    %in local coordinates (ui,vi), quadrant (edge does not change)
    pi_2 = pi/2;
    uvi = [0,0];
    %sign and magnitude of quadrant change for each case
    if edge==3
        if quad==1 ||quad==2
            sgn = +2;
        else
            sgn = -2;
        end
    elseif (edge==1 && (quad==1 || quad==3)) || (edge==2 && (quad==2 || quad==4))
        sgn = -1;                
    else
        sgn = +1;
    end
    %change in radians to handle 1 to 4 cross-over
    phi = mod(quad*pi_2+sgn*pi_2,2*pi);
    if phi==0, phi = 2*pi; end
    newquad = int8(2*phi/pi);

    %update central node coordinates if point is on edge 3
    if edge==3
        if quad==1
            uvi = [+1,+1];
        elseif quad==2
            uvi = [-1,+1;];
        elseif quad==3
            uvi = [-1,-1];
        elseif quad==4
            uvi = [+1,-1];
        else
            error('New quadrant not found in next_element')
        end
    end
end
%%
function alpha = exit_angle(phi,r,ur,vr,uvray)
    %find the angle that is tangent to the arc at the exit point
    phi = mod(phi+pi,2*pi);                   %angle of normal from centre of arc
    L = sqrt((ur-uvray(1))^2+(vr-uvray(2))^2); %arc segment length
    theta = 2*asin(L/2/r);                     %angle between entry and exit point
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
function isintol = isclose(angle,test,tol)
    %check whether test is within tol of angle
    isintol = abs(angle-test)<tol;
end
%%
function [isdir,npt] = checkDirection(lineseg,ur,vr,alpha)
    %check the direction of a vector relative to direction given by alpha
    % npt is nearest point to [ur,vr] and this is in the direction of alpha
    % if isdir is true
    anglim = 0.2;                         %angle limit of +/-0.2 rads (11.5 deg)
    npt = dsearchn(lineseg,[ur,vr]);
    un = lineseg(npt,1)-ur;   
    vn = lineseg(npt,2)-vr; 
    [xsi,~] = cart2pol(un,vn);            %vector to next point
    xsi = mod(xsi,2*pi); %N
    ub = mod(alpha+anglim,2*pi); %b  
    lb = mod(alpha-anglim,2*pi); %a
    if lb<ub                              %test accounts for wrap at 0-2pi
        isdir = lb<=xsi && xsi<=ub;
    else 
        isdir = lb<=xsi || xsi <=ub;      
    end
end

%%
% function Circ = get_circle(radius,uc,vc)   
%     %define a circle in local coordinates
%     n = 50;                                  %number of points in circle
%     phi = (0:n-1)*(2*pi/n);                  %angle increments
%     u = uc + radius*cos(phi);                %coordinates of circle
%     v = vc + radius*sin(phi);
%     Circ = polyshape(u,v);
% end

