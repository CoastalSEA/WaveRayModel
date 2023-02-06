function Tri = get_element(quad)
%
%-------function help------------------------------------------------------
% NAME
%   get_element.m
% PURPOSE
%    define triangular polyshape based on quadrant being entered
% USAGE
%    Tri = get_element(quad);
% INPUTS
%   quad - trigonometric quadrant 1-4 of the right angled triangle, or
%          adjacent combination of triangles, 12,23,34,41.
% OUTPUTS
%   Tri - triangular polyshape for the required quadrant(s)
% SEE ALSO
%   get_edge, get_quadrant, next_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if length(quad)>1
        warndlg('quad must be a scalar integer or double')
        Tri = []; return;
    end
    
    switch quad
        case 1                                  %first quadrant
            Tri = polyshape([0,1,0],[0,0,1]);
        case 2                                  %second quadrant
            Tri = polyshape([0,0,-1],[0,1,0]);
        case 3                                  %third quadrant
            Tri = polyshape([0,-1,0],[0,0,-1]);
        case 4                                  %fourth quadrant
            Tri = polyshape([0,0,1],[0,-1,0]);
        case 12                                 %1st & 2nd quadrants
            Tri = polyshape([1,0,-1],[0,1,0]);
        case 23                                 %2nd & 3rd quadrants
            Tri = polyshape([0,-1,0],[1,0,-1]);
        case 34                                 %3rd & 4th quadrants
            Tri = polyshape([-1,0,1],[0,-1,0]);
        case 41                                 %4th & 1st quadrants
            Tri = polyshape([0,1,0],[-1,0,1]);
        otherwise
            %quadrant not found
            Tri = [];
        return;
    end

end