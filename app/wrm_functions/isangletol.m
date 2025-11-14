function isdir = isangletol(theta,bound)
%
%-------function help------------------------------------------------------
% NAME
%   isangletol.m
% PURPOSE
%    boolean check of whether an angle lies between upper and lower bounds
%    defined as specific angles, or a tolerance
% USAGE
%    isdir = isangletol(theta,bound)
% INPUTS
%   theta - array of angles to be checked in radians
%   bound - scalar or [1x2] element vector. scalar defines a +/- tolerance 
%           relative to theta and vector specifies the angles to be used 
%           for the bounds in radians. 
% OUTPUTS
%   isdir - logical, true if angle lies between the defined bounds
% NOTES
%   If bound is defined as 2 angles relative to a mean direction, theta0, 
%   the function can be used to check whether theta is within the limits 
%   relative to the mean direction.
% SEE ALSO
%   get_quadrant, next_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%
%--------------------------------------------------------------------------
%
    if length(bound)>2
        warndlg('Bound must be scalar or a [1x2] array')
        isdir = [];
        return;       
    end

    theta = mod(theta,2*pi);         %use mod to check that all angles are 
                                     %between 0 and 2pi
    if length(bound)>1               %use specified angles as range
        lb = mod(bound(1),2*pi);     %lower bound
        ub = mod(bound(2),2*pi);     %upper bound        
    else                             %use specified tolerance as range
        lb = mod(theta-bound,2*pi);  %lower bound
        ub = mod(theta+bound,2*pi);  %upper bound
    end

    %test whether theta lies between lower and upper bound
    if lb<ub                              %test accounts for wrap at 0-2pi
        isdir = lb<=theta & theta<=ub;
    else 
        isdir = lb<=theta | theta<=ub;      
    end
end