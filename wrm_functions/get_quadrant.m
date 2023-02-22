function quad = get_quadrant(theta)
%
%-------function help------------------------------------------------------
% NAME
%   get_quadrant.m
% PURPOSE
%    find the quadrant that the start point lies in or on. Once first
%    intersection has been found subsequent quadrants are defined in
%    next_element, which calls get_quadrant if ray direction is aligned to axis 
% USAGE
%    quad = get_quadrant(theta,alpha,uvs,dcs,tol);
% INPUTS
%   theta - angle in radians from local origin to start point (0-2pi)
%           set to ray direction, alpha, if start point is at a node
% OUTPUTS
%   quad - trigonometric quadrant 1-4 of the right angled triangle, or
%          adjacent combination of double triangles, 12,23,34,41.
% NOTES
%   uses tol to test for proximity to multiple of pi/2. if within tol then
%   a double triangle is used rather than a single qudrant.
% SEE ALSO
%   get_element, next_element and arc_ray. 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if length(theta)>1
        warndlg('theta must be a scalar integer or double')
        quad = []; return;
    end

    theta = mod(theta,2*pi);              %ensure theta between 0-2pi
    %use theta to determine required quadrant
    if theta>=0 && theta<pi/2               %first quadrant
        quad = int8(1);
    elseif theta>=pi/2 && theta<pi          %second quadrant
        quad = int8(2);
    elseif theta>=pi && theta<3*pi/2        %third quadrant);
        quad = int8(3);
    elseif theta>=3*pi/2 && theta<2*pi      %fourth quadrant
        quad = int8(4);
    else  
        %quadrant not found
        quad = [];
        return;
    end
end