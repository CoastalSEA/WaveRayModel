function edge = get_edge(invector,idx,tol)
%
%-------function help------------------------------------------------------
% NAME
%   get_edge.m
% PURPOSE
%    use coordinates of intesection point to identify which side or edge
%    of a triangle the point lies on
% USAGE
%    edge = get_edge(invector,idx,tol)
% INPUTS
%   invector - vector of points (x,y) inside the element
%   idx - index of point that intersects a side of the element
%   tol - distance tolerance when testing proximity to x=0 or y=0
% OUTPUTS
%   edge - side of triangle that is intersected by ray
% SEE ALSO
%   get_element, get_quadrant, next_element and arc_ray. 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if idx>size(invector,1) || size(invector,2)~=2
        warndlg('invector must be an [Nx2] with N>=idx')
        edge = []; return;
    end

    if abs(invector(idx,2))<=tol            %y<tol ie appox 0
        edge = 1;                           %x-directed edge
    elseif abs(invector(idx,1))<=tol        %x<tol ie appox 0
        edge = 2;                           %y-directed edge
    else                                    %x~=0 & y~=0
        edge = 3;                           %hypotenuse
    end
end
