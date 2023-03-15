    
function uvray = get_intersection(ray,lineseg,uvr,tol)
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
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad
%   lineseg - line or arc segment to be intersected
%   uvr - location of ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   uvray - local coordinates of ray exit point
% SEE ALSO
%   get_quadrant, get_element, is_axis_point, arc_ray
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    Tri = get_element(ray.quad);
    if isempty(Tri), uvray = []; return; end

    %find intersection point of line with triangle
    tri = [Tri.Vertices;Tri.Vertices(1,:)];
    uvray = InterX(tri',lineseg')';
    if isempty(uvray) 
        plot_element(Tri,lineseg,uvr,ray)
        errtxt = sprintf('quad %d uv [%.3g %.3g] with direction %.3g',ray.quad,uvr(1),uvr(2),ray.alpha);
        error('Intersection with element not found in get_intersection 1\n%s\n',errtxt)
    end    
    entrypoint = ismembertol(uvray,uvr,tol.dist, 'DataScale', 1,'ByRows',true);  %tolerance 
    uvray(entrypoint,:) = [];

    %if more than one intersection of element find nearest point
    if size(uvray,1)>1
        npt = dsearchn(uvray,uvr);
        uvray = uvray(npt,:);
    end

    if isempty(uvray)
        plot_element(Tri,lineseg,uvr,ray)
        errtxt = sprintf('quad %d uv [%.3g %.3g] with direction %.3g',ray.quad,uvr(1),uvr(2),ray.alpha);
        error('Intersection with element not found in get_intersection 2\n%s\n',errtxt)
    end      
end
%%
function plot_element(Tri,lineseg,uvr,ray)
    %plot the element and arc
    [ualp,valp] = pol2cart(ray.alpha,1);    %arrow vector

    figure('Name','ArcRay','Tag','PlotFig');
    plot(Tri);
    hold on 
    %plot(uc,vc,'or');
    plot(uvr(1),uvr(2),'xr');
    plot(lineseg(:,1),lineseg(:,2));
    quiver(uvr(1),uvr(2),ualp,valp);
    hold off
    xlim([-10,10])
    ylim([-10,10])
    title(sprintf('Element in quad %d with direction %.3g',ray.quad,ray.alpha))
end