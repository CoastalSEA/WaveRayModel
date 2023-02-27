function quad = get_quadrant(ray,uvr,ison)
%
%-------function help------------------------------------------------------
% NAME
%   get_quadrant.m
% PURPOSE
%    find the quadrant that the start point lies in or on. Once first
%    intersection has been found subsequent quadrants are defined in
%    next_element, which calls get_quadrant if ray direction is aligned to axis 
% USAGE
%    quad = get_quadrant(ray,uvr,ison)
% INPUTS
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad, side of element, edge
%   uvr - location of ray point [u,v] in local coordinates
% ison(1) 0=not on axis; 1=on x-axis, 2=on y-axis,
%         3=on forward axis-diagonal, 4=on backward axis-diagonal;
% ison(2) direction of ray if on axis: 1=pi/2,2=pi,3=3pi/2,4=0|2pi  
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
    [theta,rs] = cart2pol(uvr(1),uvr(2));           %vector to start point
    if rs==0, theta = ray.alpha; end    
    theta = mod(theta,2*pi);              %ensure theta between 0-2pi

    
    if ison(1)==0                          %not on axis
        quad = getQuad(theta);
    elseif ison(1)>2                       %along diagonal
        quad = int8(ison(2));
    else
        quad = getAxisQuad(ison);          %along x or y axis
    end
end
%%
function quad = getQuad(theta)
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
%%
function quad = getAxisQuad(ison)
    %use isaxispoint to determine which double quadrant to use
    if ison(2)==1
        quad = int8(12);                    %theta close to pi/2
    elseif ison(2)==2
        quad = int8(23);                    %theta close to pi
    elseif ison(2)==3
        quad = int8(34);                    %theta close to 3pi/2
    elseif ison(2)==4
        quad = int8(41);                    %theta close to 0 or 2pi
    end
end


