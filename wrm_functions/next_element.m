function [newquad,edge,uvi] = next_element(ray,uvr,tol)
%
%-------function help------------------------------------------------------
% NAME
%   next_element.m
% PURPOSE
%    find the grid definition for the element that the ray is entering
%    in local coordinates (ui,vi) and quadrant (edge does not change)
% USAGE
%    [newquad,edge,uvi] = next_element(ray,uvr,tol)
% INPUTS
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad, side of element, edge
%         quad - trigonometric quadrant that the ray is leaving defined 
%                as 1-4 of the right angled triangle, or
%                adjacent combination of triangles, 12,23,34,41.
%         edge - the side of the triangle that the ray intersects and has a 
%                value of 1-3: edge=1 if yr=yi; edge=2 if xr=zi; else edge=3.
%   uvr - location of ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   newquad - quadrant ray is entering relative to the new local origin
%   edge - side of triangle that ray is moving in to
%   uvi - central node coordinates relative to current [0,0] local origin
% SEE ALSO
%   get_quadrant, get_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if ~isa(ray,'table')
        warndlg('Input must be a single row Ray table')
    end
    uvi = [0,0];

    quad = ray.quad; edge = ray.edge; %element being exited

    %resolve double quad cases
    if quad>4
        [newquad,edge,uvi] = resolve_quad(ray,uvr,tol);  
    elseif edge==0
        newquad = find_diag_quad(ray);
        diagquad = int8(mod(quad+1,4)+1);
        if newquad~=diagquad
            fprintf('Direction quad %d and Diagonal quad %d\n',newquad,diagquad)
        end
    else
        pi_2 = pi/2;    
        %sign and magnitude of quadrant change for each case
        if edge==3
            if quad==1 || quad==2
                sgn = +2;
            else
                sgn = -2;
            end
        elseif (edge==1 && (quad==1 || quad==3)) || (edge==2 && (quad==2 || quad==4))
            sgn = -1;                
        else
            sgn = +1;
        end
        %change in radians to handle 1 to 4 cross-over
        quad = double(quad);                %convert from int8 for arithmetic
        phi = mod(quad*pi_2+sgn*pi_2,2*pi);
        if phi==0, phi = 2*pi; end
        newquad = int8(2*phi/pi);
    end

    %update central node coordinates if point is on edge 3
    if edge==3
        switch quad
            case 1 
                uvi = [+1,+1];
            case 2
                uvi = [-1,+1;];
            case 3
                uvi = [-1,-1];
            case 4
                uvi = [+1,-1];
        otherwise
            error('New quadrant not found in next_element')
        end
    end
end
%%
function [quad,edge,uvi] = resolve_quad(ray,uvr,tol)
    %ray is along axis try using radius to resolve
        
        %quadrant based on the celerity gradient radius
        if abs(ray.r)<1000               %matches tol.angle tolerance 
            %check if theta in same direction as ray to determine gradient sign     
            ray.alpha = ray.alpha+sign(ray.r)*2*tol.angle;
            [ison,uvi] = is_axis_point(ray,uvr,tol);
            quad = get_quadrant(ray,uvr,ison); %update quadrant based on gradient direction
            if quad>4
                edge = ray.edge;
            else
                edge = get_edge(ison);
            end
        else
            %directed along axis with radius~=straight line, no change
            quad = ray.quad; edge = ray.edge; uvi = [0,0];
        end
end
%%
function edge = get_edge(ison)
    %assign and edge when directed along axis
    if ison(2)==1 || ison(2)==3
        edge = 2;
    else 
        edge = 1;
    end
end
%%
function quad = find_diag_quad(ray)
    %use the ray direction to find the quad being enetered as ray passes
    %through or close to the origin
    alpha = ray.alpha;
    if alpha>=0 && alpha<pi/2               %first quadrant
        quad = int8(1);
    elseif alpha>=pi/2 && alpha<pi          %second quadrant
        quad = int8(2);
    elseif alpha>=pi && alpha<3*pi/2        %third quadrant
        quad = int8(3);
    elseif alpha>=3*pi/2 && alpha<2*pi      %fourth quadrant
        quad = int8(4);
    else  
        %quadrant not found
        quad = [];
        return;
    end
end