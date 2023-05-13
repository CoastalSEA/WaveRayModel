function ison = is_axis_point(ray,uvr,tol)
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
%         node index, k, quadrant being entered, quad
%   uvr - location of ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   ison - [1x2] vector: ison(1) 0=not along axis; 1=on x-axis, 2=on y-axis,
%          3=on forward axis-diagonal, 4=on backward axis-diagonal,
%          -1=origin but not along an axis;
%          if ison(1)>0
%          ison(2) direction of ray if on axis: 1=pi/2,2=pi,3=3pi/2,4=0|2pi
%          if ison(1)=0
%          ison(2) direction of ray: north,west,south,east relative to axis
% SEE ALSO
%   get_quadrant, get_element, arc_ray, get_intersection
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
        %check each direction for solution and exit loop if found
        %check whether alpha is close to n.pi/2
        bound =[i*pi2-tol.angle,i*pi2+tol.angle];
        israydir2 = isangletol(ray.alpha,bound);  %ray direction is also along axis
        %for case when ray is on a node and running along grid diagonal
        bound =[(2*i-1)*pi4-tol.angle,(2*i-1)*pi4+tol.angle];
        israydir4 = isangletol(ray.alpha,bound);  %ray direction is also along axis
        if isorigin && israydir2
            if i==0 || i==2 || i==4
                ison(1) = 1;                      %along x-axis from origin
            else
                ison(1) = 2;                      %along y-axis from origin
            end
        elseif isonxaxis && israydir2 && (i==0 || i==2 || i==4)
            ison(1) = 1;                          %axis point & along x-axis
        elseif isonyaxis && israydir2 && (i==1 || i==3)
            ison(1) = 2;                          %axis point & along y-axis    
        elseif isorigin && israydir4
            if i==1 || i==3
                ison(1) = 3;                      %along forward diagonal from origin
            else
                ison(1) = 4;                      %along backward diagonal from origin
            end 
        elseif isonforediag && israydir4
            ison(1) = 3;                          %diagonal point along forward diagonal
        elseif isonbackdiag && israydir4    
            ison(1) = 4;                          %diagonal point along backward diagonal        
        end
        %
        if ison(1)>0
            %identifier of axis direction or ray: 1=pi/2,2=pi,3=3pi/2,4=0 or 2pi
            ison(2) = (i>0)*i+(i==0)*4;           %4 if i=0 otherwise i
            break                                 %solution found exit loop        
        end
    end

    %no solution for ray along axis but point is on an axis
    if ison(1)==0
        if isorigin
            ison(2) = -1;                             %origin point
        elseif isonxaxis                          
            if ray.alpha>0 && ray.alpha<pi            %on x-axis but not along it
                ison(2) = 1;                          %north going
            else
                ison(2) = 3;                          %south going
            end
        elseif isonyaxis                     
            if ray.alpha>pi/2 && ray.alpha<3*pi/2     %on y-axis but not along it
                ison(2) = 2;                          %west going
            else
                ison(2) = 4;                          %east going
            end
        end
    end                               
end
%%
function ison = IsPointOnLine(xy1,xy2,xy3)
    %test whether point xy3 is on the line defined by xy1 and xy2
    % Line equation: y = m*x + b;  by Jan on Matlab Forum
    Limit = 100 * eps(max(abs([xy1,xy2,xy3])));
    if xy1(1) ~= xy2(1)
      m = (xy2(2)-xy1(2))/(xy2(1)-xy1(1));
      yy3 = m*xy3(1) + xy1(2) - m*xy1(1);
      ison   = (abs(xy3(2) - yy3) < 100 * Limit);
    else
      ison   = abs(xy3(1)) < Limit;
    end
end