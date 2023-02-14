function [isaxispoint,theta] = is_axis_point(theta,alpha,uvr,dcs,tol)
%
%-------function help------------------------------------------------------
% NAME
%   is_axis_point.m
% PURPOSE
%    find whether point lies on an axis and is travelling in the direction
%    of that axis
% USAGE
%    [isaxispoint,theta] = is_axis_point(theta,alpha,uvr,dcs,tol);
% INPUTS
%   theta - angle in radians from local origin to start point (0-2pi)
%           set to ray direction, alpha, if start point is at a node
%   alpha - angle tangential to ray direction
%   uvr - location of ray point [u,v] in local coordinates
%   dcs - celerity gradient at start point [dcx,dcy]
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   isaxispoint - logical, true if point is on axis and ray is travellingin
%                 direction of axis
%   theta - angle in radians from local origin to start point (0-2pi) with
%           offset if gradient defines direction of change
% SEE ALSO
%   get_quadrant and next_element
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%----------------------------------------------------------------------
%
    gradtol = 1e-6;                       %tolerance used to chack gradient
    distol = tol;                         %use angle tolerance for distance tolerance 

    pi2 = pi/2; isaxispoint = false(1,5);
    for i=0:1:4                           %check whether theta is close to n.pi/2
        bound =[i*pi2-tol,i*pi2+tol];     %and point lies on axis in that direction        
        isaxisdir = isangletol(theta,bound);
        israydir1 = isangletol(alpha,bound); %ray direction is also along axis
        israydir2 = isangletol(alpha+pi,bound); %ray direction is also along axis
        if i==0 || i==2 || i==4
            isonaxis = uvr(2)<distol;     %0, pi, 2pi and y~0
        else
            isonaxis = uvr(1)<distol;     %pi/2, 3pi/2 and x~0
        end        
        isaxispoint(i+1) = isaxisdir && isonaxis && (israydir1 || israydir2);
    end

    if any(isaxispoint)
        %try to resolve ray direction with celerity gradient 
        phi = mod(alpha+pi/2,2*pi);        %angle of normal to ray direction  
        [un,vn] = pol2cart(phi,1);         %normal unit vector in local coordinates
        cgrad = un*dcs(1)+vn*dcs(2);       %celerity gradient dot product with normal unit vector
        %if gradient is greater than gradtol, adjust theta to the 
        %quadrant based on the element gradient
        if abs(cgrad)>gradtol
            %check if theta in same direction as ray to determine gradient sign
            israydir1 = isangletol(theta,[alpha-tol,alpha+tol]);   %same direction
            israydir2 = isangletol(theta+pi,[alpha-tol,alpha+tol]);%opposite directions
            if israydir1
                theta = theta+sign(cgrad)*2*tol;  %alpha~=theta
            elseif israydir2
                theta = theta-sign(cgrad)*2*tol;  %alpha~=-theta
            end
            theta = mod(theta,2*pi);       %ensure theta between 0-2pi
            isaxispoint = false;
        end
    end
end