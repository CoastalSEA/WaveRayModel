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
%         node index, k, quadrant being entered, quad
%   uvr - location of ray point [u,v] in local coordinates
%   ison(1) 0=not on axis; 1=on x-axis, 2=on y-axis,
%           3=on forward axis-diagonal, 4=on backward axis-diagonal;
%   ison(1)>0: ison(2) direction of ray if on axis: 1=pi/2,2=pi,3=3pi/2,4=0|2pi
%   ison(1)=0: ison(2) direction of ray: north,west,south,east relative to axis 
% OUTPUTS
%   quad - trigonometric quadrant 1-4 of the right angled triangle, or
%          adjacent combination of double triangles, 12,23,34,41.
% NOTES
%   uses tol to test for proximity to multiples of pi/2. if within tol then
%   a double triangle is used rather than a single qudrant.
% SEE ALSO
%   get_element, is_axis_point, arc_ray, get_intersection
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if sum(ison)==0                           %point not on either axis  
        [theta,rs] = cart2pol(uvr(1),uvr(2)); %vector to start point
        if rs==0, theta = ray.alpha; end    
        theta = mod(theta,2*pi);              %ensure theta between 0-2pi        
        quad = getQuad(theta);   
    elseif ison(1)==0 && ison(2)<0            %point at origin but not along axis
        quad = getQuad(ray.alpha);   
    elseif ison(1)==0 && ison(2)>0            %point on axis but not along it
        quad = getAxisQuad(ison,uvr);
%     elseif ison(1)>2                          %along diagonal
%         quad = int8(ison(2));
    else                                      %point on x or y axis and along it
        quad = getAlongAxisQuad(ison);        
    end        
end
%%
function quad = getQuad(theta)
    %use theta to determine required quadrant
    if theta>=0 && theta<pi/2               %first quadrant
        quad = int8(1);
    elseif theta>=pi/2 && theta<pi          %second quadrant
        quad = int8(2);
    elseif theta>=pi && theta<3*pi/2        %third quadrant
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
    function quad = getAxisQuad(ison,uvr)
        %point on an axis but direction is not aligned along axis
        if ison(2)==1 || ison(2)==3         %on x-axis
            if uvr(1)>0                     %x>0
                if ison(2)==1               %north going
                    quad = int8(1);
                elseif ison(2)==3           %south going
                    quad = int8(4);                    
                end                    
            else
                if ison(2)==1               %north going
                    quad = int8(2);
                elseif ison(2)==3           %south going
                    quad = int8(3);                    
                end
            end
        else                                %on y-axis
            if uvr(2)>0                     %y>0
                if ison(2)==2               %west going
                    quad = int8(2);
                elseif ison(2)==4           %east going
                    quad = int8(1);                    
                end                    
            else
                if ison(2)==2               %west going
                    quad = int8(3);
                elseif ison(2)==4           %east going
                    quad = int8(4);                    
                end
            end
        end
    end
%%
    function quad = getAlongAxisQuad(ison)
    %unable to resolve quadrant becuase point on axis and alinged along it
    %assign a double quadrant to progress ray to next node unless the ray
    %curvature moves it off axis (tested in next_element)
    if ison(2)==1
        quad = int8(12);                    %theta close to pi/2
    elseif ison(2)==2
        quad = int8(23);                    %theta close to pi
    elseif ison(2)==3
        quad = int8(34);                    %theta close to 3pi/2
    elseif ison(2)==4
        quad = int8(41);                    %theta close to 0 or 2pi
    else
        %quadrant not found
        quad = [];
        return;
    end
end


