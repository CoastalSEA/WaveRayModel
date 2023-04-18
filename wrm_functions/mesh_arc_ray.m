function newray = mesh_arc_ray(cmesh,ray,tol)
%
%-------function help------------------------------------------------------
% NAME
%   mesh_arc_ray.m
% PURPOSE
%   compute the exit position and direction of a ray entering a triangular
%   mesh element at the position and direction defined by the incoming ray.
% USAGE
%   newray = mesh_arc_ray(cgrid,ray,tol);
% INPUTS
%   cgrid - struct of grid array properties, including:
%          Tri - triangulation object of mesh
%          h - water depths
%          c - celerity grid matrix 
%          dcx, dcy - gradients of celerity in x and y directions
%   ray - table of incoming ray position (xr,yr), direction, alpha, and 
%         local node index, k
%   tol - tolerance to test for angles that are multiples of pi/2 (5.7 deg)
% OUTPUTS
%   newray - table of outgoing ray position (xr,yr), direction, alpha, and
%            local node index, k
%            - returned empty if new point is outside grid domain
%            - returned as -2 if radius is too small (<0.5) 
% NOTES
%   Ray method based on Abernethy C L and Gilbert G, 1975, Refraction of 
%   wave spectra, Report No: INT 117,pp. 1-166, Wallingford, UK.
%   
% SEE ALSO
%   
%
% Author: Ian Townend
% CoastalSEA (c) Apr 2023
%----------------------------------------------------------------------
%

    %variables used in function
    % alpha - angle of ray direction
    % phi - angle of normal to ray direction
    % theta - angle from origin to ray position
    % kr - index of reference grid node (nearest point to xr,yr)
    % xi,yi - position of reference grid node for local coordinates
    % xr,yr - ray position in grid coordinates
    % r - radius of arc in local coordinates (R = r*delta)  

    %find the largest edge connected to kr    
    tol.length = maxEdge(cmesh,ray);

    %get the centre of the arc that is tangential to the ray at xr,yr
    [phi,r,xc,yc] = arc_properties(cmesh,ray,tol);
    xyArc = get_arc(phi,r,xc,yc,ray.xr,ray.yr,tol);
    Tri = get_element(cmesh,ray,tol);
    if isempty(Tri)
        newray = []; return; 
    elseif any(isnan(xyArc))
        plot_element(Tri,xyArc,[ray.xr,ray.yr],ray)
        error('Arc segment not found at xr=%.0f; yr=%.0f',ray.xr,ray.yr)
    elseif isempty(xyArc)
        % warndlg(sprintf('Radius, r=%0.2f too small at xr=%.0f; yr=%.0f',r,ray.xr,ray.yr));
        newray = -2; return;
    end
    xyray = get_intersection(Tri,xyArc,[ray.xr,ray.yr],ray,tol);
    kr = nearestNeighbor(cmesh.Tri,xyray);
    %wave properties at new ray point
    xr = xyray(1); yr = xyray(2); quad = 0;
    [hr,cr,cgr] = raypoint_properties(cmesh,xr,yr); 
    %get exit angle (ur,vr entry point into element, uvray exit point)
    alpha =  exit_angle(phi,r,ray.xr,ray.yr,xyray,tol);  

    isoutbound = checkGridBoundary(cmesh,[xr,yr]);
    if isoutbound 
        newray = []; return; 
    else
        newray = table(xr,yr,alpha,kr,quad,r,hr,cr,cgr);
    end  
end
%%
function [phi,R,xc,yc] = arc_properties(cmesh,ray,tol)
    %find the radius and centre of the arc that the ray follows in element
    pts = cmesh.Tri.Points;
    X = pts(:,1); Y = pts(:,2); c = cmesh.c; dcx = cmesh.dcx; dcy = cmesh.dcy;
    
    method = 'natural';
    cr = griddata(X,Y,c,ray.xr,ray.yr,method);
    dcrx = griddata(X,Y,dcx,ray.xr,ray.yr,method);
    dcry = griddata(X,Y,dcy,ray.xr,ray.yr,method);

    %unit normal to ray (convention is left in direction of ray is positive)
    phi = mod(ray.alpha+pi/2,2*pi);    %angle of normal to ray direction    
    [un,vn] = pol2cart(phi,1);         %normal unit vector

    %radius of arc
    Ndc = un*dcrx+vn*dcry;
    R = -cr/Ndc;                       %radius in grid coordinates
    if abs(R)<tol.radius
        xc = ray.xr+R*un;
        yc = ray.yr+R*vn;
    else
        xc = 0; yc = 0;   %radius is large so use straight line segment
    end
end
%%
function [xy_arc] = get_arc(phi,radius,xc,yc,xr,yr,tol)  
    %calculate the coordinates of an arc either side of the radius vecor
    %from uc,vc to ur,vr.
    N = 20;                                       %number of points in half-Arc
    rt2 = tol.length;                             %maximum edge length
    rad = abs(radius);
    alpha = mod(phi-pi/2,2*pi);
    
    if rad>=tol.radius
        %straight line segment will suffice       
        [xe,ye] = pol2cart(alpha,rt2); %vector from ray point in direction of alpha
                                       %use 2 to ensure line crosses a boundary
                                       %max element length is sqrt(2) 
        [xs,ys] = pol2cart(phi-pi/2,tol.dist);    %small offset from start point      
        xy_arc = [xr+xs,yr+ys;xr+xe,yr+ye];       %ray vector line segment
        return;
    elseif rad<1.0                                %radius is so small that it may not exit element
        xy_arc = []; return;
    end

    phi = phi+pi;                                  %angle of normal from centre of arc
    arcang = rt2/rad;                              %set arc segment based on radius
    %abstract Circle Function For Angles In Radians
    circr = @(radius,angle)  [radius*cos(angle)+xc,  radius*sin(angle)+yc]; 
    r_angl = linspace(phi-arcang,phi-tol.angle, N);%angles Defining left Arc Segment (radians)
    xy_arcl = circr(radius,r_angl');
    r_angr = linspace(phi+tol.angle,phi+arcang, N);%angles Defining right Arc Segment (radians)
    xy_arcr = circr(radius,r_angr');
    %find arc in direction of ray
    ua1 = xy_arcl(1,1)-xr; va1 = xy_arcl(1,2)-yr;  %vector from ray point to first point on arc
    arcl = mod(cart2pol(ua1,va1),2*pi);
    ua2 = xy_arcr(1,1)-xr; va2 = xy_arcr(1,2)-yr;  %vector from ray point to first point on arc
    arcr = mod(cart2pol(ua2,va2),2*pi);
    angtol = 1.0;                                  %angle limit of +/-1.0 rads (57.3 deg)
    bound = [alpha-angtol,alpha+angtol];
   
    if isangletol(arcl,bound)                      %check if arcl is within bound   
        xy_arc = xy_arcl;
    elseif isangletol(arcr,bound)                  %check if arcr is within bound   
        xy_arc = xy_arcr;
    else
        xy_arc = [];
        % warndlg('Arc not found in get_arc')
    end
end
%%
function len = maxEdge(cmesh,ray)
    %find the maximum length of the edges connected to nearest node to
    %current ray position
    pts = cmesh.Tri.Points;

    edgeID = edges(cmesh.Tri);                   %node ids for all edges
    eids = [edgeID(edgeID(:,1)==ray.kr,:);edgeID(edgeID(:,2)==ray.kr,:)];
    eids = eids(eids~=ray.kr);                   %node ids linked to kr

    npts = length(eids);
    a = repmat(pts(ray.kr,:),npts,1);
    b = pts(eids,:);
    len = max(vecnorm(a-b,2,2))*1.1;   %increase by 10% to ensure length
                                       %exits element in cases where 3rd 
                                       %side is longer

%     figure; trisurf(cgrid.Tri.ConnectivityList,cgrid.Tri.Points(:,1),...
%                                             cgrid.Tri.Points(:,2),-cgrid.h);
%     view(2)
%     hold on
%     plot(a(1,1),a(1,2),'xr')
%     plot(b(:,1),b(:,2),'or')
%     hold off
end
%%
function Tri = get_element(cmesh,ray,tol)
    %
    [xt,yt] = pol2cart(ray.alpha,tol.dist*1000); %normal vector in local coordinates *****
    xp = ray.xr+xt;
    yp = ray.yr+yt;
    triID = pointLocation(cmesh.Tri,xp,yp);
    if isnan(triID), Tri = []; return; end
    %probably need to check if along edge

    tri_pt_ids = cmesh.Tri.ConnectivityList(triID,:);
    tri_pts = cmesh.Tri.Points(tri_pt_ids,:);
    Tri = polyshape(tri_pts);
end
%%
function xyray = get_intersection(Tri,lineseg,xyr,ray,tol)
    %find intersection point of line with triangle
    tri = [Tri.Vertices;Tri.Vertices(1,:)];
    xyray = InterX(tri',lineseg')';
    if isempty(xyray) 
        plot_element(Tri,lineseg,xyr,ray)
        errtxt = sprintf('element %d uv [%.3g %.3g] with direction %.3g',ray.kr,xyr(1),xyr(2),ray.alpha);
        error('Intersection with element not found in get_intersection 1\n%s\n',errtxt)
    end    
    entrypoint = ismembertol(xyray,xyr,tol.dist, 'DataScale', 1,'ByRows',true);  %tolerance 
    xyray(entrypoint,:) = [];

    %if more than one intersection of element find nearest point
    if size(xyray,1)>1
        npt = dsearchn(xyray,xyr);
        xyray = xyray(npt,:);
    end

    if isempty(xyray)
        plot_element(Tri,lineseg,xyr,ray)
        errtxt = sprintf('quad %d uv [%.3g %.3g] with direction %.3g',ray.quad,xyr(1),xyr(2),ray.alpha);
        error('Intersection with element not found in get_intersection 2\n%s\n',errtxt)
    end      
end
%%
function [hr,cr,cgr] = raypoint_properties(cmesh,xr,yr)
    %interpolate depth, celerity and group celerity at ray point  
    method = 'natural';
    X = cmesh.Tri.Points(:,1); Y = cmesh.Tri.Points(:,2); 
    hr = griddata(X,Y,cmesh.h,xr,yr,method);
    cr = griddata(X,Y,cmesh.c,xr,yr,method);
    cgr = griddata(X,Y,cmesh.cg,xr,yr,method);               
end
%%
function alpha = exit_angle(phi,r,xr,yr,xyray,tol)
    %find the angle that is tangent to the arc at the exit point
    % phi,r - angle of normal and radius of arc
    % uv,vr - ray entry point into element
    % uvray - ray exit point out of element
    phi = mod(phi+pi,2*pi);                    %angle of normal from centre of arc
    L = sqrt((xr-xyray(1))^2+(yr-xyray(2))^2); %arc segment length
    if abs(r)<tol.radius
        theta = 2*asin(L/2/r);                 %angle between entry and exit point
    else
        theta = 0;                             %straight line no change in direction
    end
    alpha = phi+theta+pi/2;                    %angle of tangent to ray at exit point
    alpha = mod(alpha,2*pi);
end
%%
function isoutbound = checkGridBoundary(cmesh,xye)
    %check if point, xye, is outside grid
    x = cmesh.Tri.Points(:,1); y = cmesh.Tri.Points(:,2);
    xminmax = minmax(x); yminmax = minmax(y);
    limxy = [xminmax',yminmax'];                              %limits of grid               
    isoutbound = xye(1)<limxy(1,1) || xye(1)>limxy(2,1) || ...%check x
              xye(2)<limxy(1,2) || xye(2)>limxy(2,2);         %check y       
end
%%
function plot_element(Tri,Arc,xyr,ray)
    %plot the element and arc
    figure('Name','ArcRay','Tag','PlotFig');
    plot(Tri);
    hold on 
    %plot(uc,vc,'or');
    plot(xyr(1),xyr(2),'xr');
    plot(Arc(:,1),Arc(:,2));
    hold off
    title(sprintf('Element at node index %d',ray.kr))
end
