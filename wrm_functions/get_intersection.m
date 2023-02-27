    
function [uvray,edge] = get_intersection(quad,lineseg,alpha,uvr,tol)
%
%-------function help------------------------------------------------------
% NAME
%   get_intersection.m
% PURPOSE
%    find intersection of a triangle element and a line segment that can be
%    a straight line or an arc segment
% USAGE
%    [uvray,edge] = get_intersection(quad,lineseg,alpha,uvr,tol)
% INPUTS
%   quad - trigonometric quadrant 1-4 of the right angled triangle, or
%          adjacent combination of triangles, 12,23,34,41.
%   lineseg - line or arc segment to be intersected
%   alpha - angle tangential to ray direction
%   uvr - location of ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   uvray - local coordinates of ray exit point
%   edge - side of triangle that the ray intersects
% SEE ALSO
%   get_quadrant, next_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    Tri = get_element(quad);
    if isempty(Tri), uvray = []; return; end

    %find intersection point of line with triangle
    [inside,outside] = intersect(Tri,lineseg); 
    %inside line segment coordinates returned as a two-column matrix (x,y).
    %find the common point in both vectors
    idx = ~ismembertol(inside,uvr,tol.dist,'ByRows',true);  %tolerance equivalent to 1m
    inside = inside(idx,:);                 %inside points excluding entry point
    [~,idx] = intersect(inside,outside,'rows');

    if size(idx,1)>1
        %more than one intersection of element
        %find nearest point in direction of travel
        ok = 0;
        nxtpnt = inside(idx,:);
        while ok<1            
            [isdir,npt] = checkDirection(nxtpnt,uvr,alpha);
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
        plot_element(Tri,lineseg,uvr,quad,alpha)
        error('Intersection with element not found in get_intersection')
    end

    uvray = inside(idx,:);                  %local coordinates of ray exit point
    %use coordinates of point to identify which edge it lies on
    if abs(uvray(2))<=tol.dist              %y<tol ie appox 0
        edge = 1;                               %x-directed edge
    elseif abs(uvray(1))<=tol.dist          %x<tol ie appox 0
        edge = 2;                               %y-directed edge
    else                                    %x~=0 & y~=0
        edge = 3;                               %hypotenuse
    end
end
%%
function [isdir,npt] = checkDirection(lineseg,uvr,alpha)
    %check the direction of a vector relative to direction given by alpha
    % npt is nearest point to [ur,vr] and this is in the direction of alpha
    % if isdir is true
    angtol = 1.0;                         %angle limit of +/-1.0 rads (57.3 deg)
    npt = dsearchn(lineseg,uvr);
    un = lineseg(npt,1)-uvr(1);   
    vn = lineseg(npt,2)-uvr(2); 
    [xsi,~] = cart2pol(un,vn);            %vector to next point
    bound = [alpha-angtol,alpha+angtol];
    isdir = isangletol(xsi,bound);        %check if xsi is within bound
end
%%
function plot_element(Tri,lineseg,uvr,quad,alpha)
    %plot the element and arc
    figure('Name','ArcRay','Tag','PlotFig');
    plot(Tri);
    hold on 
    %plot(uc,vc,'or');
    plot(uvr(1),uvr(2),'xr');
    plot(lineseg(:,1),lineseg(:,2));
    hold off
    xlim([-10,10])
    ylim([-10,10])
    title(sprintf('Element in quad %d with direction %.3g',quad,alpha))
end