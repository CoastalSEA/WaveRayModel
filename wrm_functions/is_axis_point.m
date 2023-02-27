function [ison,uvi] = is_axis_point(ray,uvr,tol)
%-------function help------------------------------------------------------
% NAME
%   is_axis_point.m
% PURPOSE
%    find whether point lies on an axis and is travelling in the direction
%    of that axis
% USAGE
%    [ison,uvi] = is_axis_point(ray,uvr,tol)
% INPUTS
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad, side of element, edge
%   uvr - location of ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   ison - [1x2] vector: ison(1) 0=not on axis; 1=on x-axis, 2=on y-axis,
%          3=on forward axis-diagonal, 4=on bakward axis-diagonal;
%          ison(2) direction of ray if on axis: 1=pi/2,2=pi,3=3pi/2,4=0|2pi
%   uvi - central node coordinates relative to current [0,0] local origin
% SEE ALSO
%   get_quadrant and next_element
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%----------------------------------------------------------------------
%
    pi2 = pi/2; pi4 = pi/4;
    ison = zeros(1,2);
    %check whether raypoint is on an axis
    isonxaxis = IsPointOnLine([-1,0],[1,0],uvr);
    isonyaxis = IsPointOnLine([0,-1],[0,1],uvr);
    isorigin = isonxaxis && isonyaxis;
    isonforediag = IsPointOnLine([-1,-1],[1,1],uvr);
    isonbackdiag = IsPointOnLine([-1,1],[1,-1],uvr);
    
    for i=0:1:4                           
        %check whether alpha is close to n.pi/2
        bound =[i*pi2-tol.angle,i*pi2+tol.angle];
        israydir2 = isangletol(ray.alpha,bound);     %ray direction is also along axis
        %for case when ray is on a node and running along grid diagonal
        bound =[(2*i-1)*pi4-tol.angle,(2*i-1)*pi4+tol.angle];
        israydir4 = isangletol(ray.alpha,bound);     %ray direction is also along axis
        if isorigin && israydir2
            if i==0 || i==2 || i==4
                ison(1) = 1;                         %along x-axis from origin
            else
                ison(1) = 2;                         %along y-axis from origin
            end
        elseif isonxaxis && israydir2 
            ison(1) = 1;                             %axis point & along x-axis
        elseif isonyaxis && israydir2
            ison(1) = 2;                             %axis point & along y-axis    
        elseif isorigin && israydir4
            if i==1 || i==3
                ison(1) = 3;                         %along forward diagonal from origin
            else
                ison(1) = 4;                         %along backward diagonal from origin
            end  
        elseif isonforediag && israydir4
            ison(1) = 3;                             %diagonal point along forward diagonal
        elseif isonbackdiag && israydir4    
            ison(1) = 4;                             %diagonal point along backward diagonal
        end
        %
        if ison(1)>0
            %identifier of axis direction or ray: 1=pi/2,2=pi,3=3pi/2,4=0 or 2pi
            ison(2) = (i>0)*i+(i==0)*4;              %4 if i=0 otherwise i
            uvi = get_uvi(uvr,ison,i);
            break
        else
            uvi = [0,0];                             %point in element
        end
    end
end
%%
function uvi = get_uvi(uvr,ison,i)
    %update central node coordinates if point is on a grid axis
    %in both cases only move origin if moving away from origin
    uvi = [0,0];
    if ison(1)>0 && ison(1)<3   %point is on x or y-axis and directed along axis        
        if uvr(1)>=0 && i==0 || i==4     %point is +ve and in direction of x-axis
            uvi = [1,0];
        elseif uvr(1)<=0 && i==2         %point is -ve and in direction of -x-axis
            uvi = [-1,0];
        elseif uvr(2)>=0 && i==1         %point is +ve and in direction of y-axis
            uvi = [0,1];
        elseif uvr(2)<=0 && i==3         %point is -ve and in direction of -yaxis
            uvi = [0,-1];
        end

    elseif ison(1)>2            %point is on diagonal and directed along diagonal
        switch ison(2)
            case 1       
                if uvr(1)>0 && uvr(2)>0
                    uvi = [1,1];         %diagonal across first quadrant
                end
            case 2
                if uvr(1)<0 && uvr(2)>0
                    uvi = [-1,1];         %diagonal across second quadrant
                end
            case 3
                if uvr(1)<0 && uvr(2)<0
                    uvi = [-1,-1];      %diagonal across third quadrant
                end
            case 4
                if uvr(1)>0 && uvr(2)<0
                    uvi = [1,-1];       %diagonal across fourth quadrant            
                end 
        end
    end
end
%%
function ison = IsPointOnLine(xy1,xy2,xy3)
    % Line equation: y = m*x + b;  by Jan on Matlab Forum
    Limit = 100 * eps(max(abs([xy1,xy2,xy3])));
    if xy1(1) ~= xy2(1)
      m = (xy2(2)-xy1(2))/(xy2(1)-xy1(1));
      yy3 = m*xy3(1) + xy1(2) - m*xy1(1);
      ison   = (abs(xy3(2) - yy3) < 100 * Limit);
    else
      ison   = (xy3(1) < Limit);
    end
end
