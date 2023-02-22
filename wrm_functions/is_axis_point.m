function [quad,edge,uvi] = is_axis_point(alpha,uvr,dcs,tol)
%
%-------function help------------------------------------------------------
% NAME
%   is_axis_point.m
% PURPOSE
%    find whether point lies on an axis and is travelling in the direction
%    of that axis
% USAGE
%    [quad,edge,uvi] = is_axis_point(alpha,uvr,dcs,tol)
% INPUTS
%   alpha - angle tangential to ray direction
%   uvr - location of ray point [u,v] in local coordinates
%   dcs - celerity gradient at start point [dcx,dcy]
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   quad - quadrant ray is passing through relative to the local origin
%   edge - edge is assigned if point is on a grid axis and ray is directed
%          along the axis
%   uvi - central node coordinates relative to current [0,0] local origin
% SEE ALSO
%   get_quadrant and next_element
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%----------------------------------------------------------------------
%
    gradtol = 1e-6;                       %tolerance used to chack gradient

    pi2 = pi/2; 
    atol = tol.angle;
    iquad = int8([41,12,23,34,41]);
    quad = []; edge = []; uvi = [0,0];    %vector to next origin
    isaxispoint = int8(zeros(1,5));
    for i=0:1:4                           
        %check whether raypoint is on an axis
        if i==0 || i==2 || i==4
            isonxaxis = uvr(2)<tol.dist;      %0, pi, 2pi and y~0
            isonyaxis = false;
        else
            isonyaxis = uvr(1)<tol.dist;      %pi/2, 3pi/2 and x~0
            isonxaxis = false;        
        end  
        %check whether alpha is close to n.pi/2
        bound =[i*pi2-atol,i*pi2+atol];
        israydir = isangletol(alpha,bound); %ray direction is also along axis
        %assign quadrant and edge if on an axis and determine direction to
        %next point based on position on axis relative to loal origin
        if israydir && isonxaxis
            isaxispoint(i+1) = 1;           %axis point on x-axis
            quad = iquad(i+1);
            edge = 1;
            if uvr(1)>=0 && i==0 ||i==4      %point is +ve and in direction of x-axis
                uvi = [1,0];
            elseif uvr(1)<=0 && i==2         %point is -ve and in direction of -x-axis
                uvi = [-1,0];
            end
        elseif israydir && isonyaxis
            isaxispoint(i+1) = -1;           %axis point on y-axis
            quad = iquad(i+1);
            edge = 2;
            if uvr(2)>=0 && i==1              %point is +ve and in direction of y-axis
                uvi = [0,1];
            elseif uvr(2)<=0 && i==3          %point is -ve and in direction of -yaxis
                uvi = [0,-1];
            end             
        end
    end

    [theta,rs] = cart2pol(uvr(1),uvr(2));           %vector to start point
    if rs==0, theta = alpha; end

    if any(isaxispoint)
        %try to resolve ray direction with celerity gradient 
        phi = mod(alpha+pi/2,2*pi);         %angle of normal to ray direction  
        [un,vn] = pol2cart(phi,1);          %normal unit vector in local coordinates
        cgrad = un*dcs(1)+vn*dcs(2);        %celerity gradient dot product with normal unit vector
        %if gradient is greater than gradtol, adjust theta to the 
        %quadrant based on the element gradient
        if abs(cgrad)>gradtol
            %check if theta in same direction as ray to determine gradient sign     
            israydir1 = isangletol(theta,[alpha-atol,alpha+atol]);   %same direction
            israydir2 = isangletol(theta+pi,[alpha-atol,alpha+atol]);%opposite directions
            if israydir1
                theta = theta-sign(cgrad)*2*atol;  %alpha~=theta
            elseif israydir2
                theta = theta+sign(cgrad)*2*atol;  %alpha~=-theta
            end
            quad = get_quadrant(theta); %update quadrant based on gradient direction
        end
    else
        quad = get_quadrant(theta); %update quadrant based on gradient direction
    end    
end