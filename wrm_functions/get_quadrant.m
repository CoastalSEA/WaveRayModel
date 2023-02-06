function quad = get_quadrant(theta,tol)
%
%-------function help------------------------------------------------------
% NAME
%   get_quadrant.m
% PURPOSE
%    find the quadrant that the start point lies in. Once first
%    intersection has been found subsequent quadrants are defined in
%    next_element
% USAGE
%    quad = get_quadrant(theta);
% INPUTS
%   theta - angle in radians from local origin to start point (0-2pi)
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   quad - trigonometric quadrant 1-4 of the right angled triangle, or
%          adjacent combination of double triangles, 12,23,34,41.
% NOTES
%   uses tol to test for proximity to multiple of pi/2. if within tol then
%   a double triangle is used rather than a single qudrant.
% SEE ALSO
%   get_edge, get_element, next_element and arc_ray. 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if length(theta)>1
        warndlg('theta must be a scalar integer or double')
        quad = []; return;
    end  

    theta = mod(theta,2*pi);                   %ensure theta between 0-2pi

    pi2 = pi/2; isdir = false(1,5);
    for i=0:1:4                                 %check whether close to n.pi/2
        bound =[i*pi2-tol,i*pi2+tol];
        isdir(i+1) = isangletol(theta,bound);
    end
    
    if any(isdir)
        %use theta to determine which double quadrant to use
        if isdir(1)                            
            quad = int8(41);                    %theta close to 0
        elseif isdir(2)
            quad = int8(12);                    %theta close to pi/2
        elseif isdir(3)
            quad = int8(23);                    %theta close to pi
        elseif isdir(4)
            quad = int8(34);                    %theta close to 3pi/2
        elseif isdir(5)
            quad = int8(41);                    %theta close to 2pi
        end
    else
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
end