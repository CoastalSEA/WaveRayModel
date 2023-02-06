function [uvi,newquad] = next_element(ray,uvr,tol)
%
%-------function help------------------------------------------------------
% NAME
%   next_element.m
% PURPOSE
%    find the grid definition for the element that the ray is entering
%    in local coordinates (ui,vi) and quadrant (edge does not change)
% USAGE
%    [uvi,newquad] = next_element(ray)
% INPUTS
%   ray - table of incoming ray position (xr,yr), direction, alpha, local
%         node index, k, quadrant being entered, quad, side of element, edge
%         quad - trigonometric quadrant that the ray is leaving defined 
%                as 1-4 of the right angled triangle, or
%                adjacent combination of triangles, 12,23,34,41.
%         edge - the side of the triangle that the ray intersects and has a 
%                value of 1-3: edge=1 if yr=yi; edge=2 if xr=zi; else edge=3.
%   uvr - location of exit ray point [u,v] in local coordinates
%   tol - tolerance around angles that are multiples of pi/2
% OUTPUTS
%   uvi - central node coordinates relative to current [0,0] local origin
%   newquad - quadrant ray is entering relative to the new origin
% SEE ALSO
%   get_edge, get_quadrant, get_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if ~isa(ray,'table')
        warndlg('Input must be a single row Ray table')
    end
    uvi = [0,0];

    quad = ray.quad; edge = ray.edge; alpha = ray.alpha;  %element being exited

    %resolve double quad cases
    if quad>4
        [ua,va] = pol2cart(alpha,1);  %unit vector in ray direction
        [theta,rs] = cart2pol(uvr(1)+ua,uvr(2)+va);
        %[theta,rs] = cart2pol(uvr(1),uvr(2));

        newquad = get_quadrant(theta,tol);
        if newquad>4
%             isposdir = 1-abs(theta-alpha)*2/pi<tol;
%             if isposdir && rs>0.5
% %                 if rs>0.5
%                     switch newquad
%                         case 12
%                             uvi = [sign(uvr(1)),0];
%                         case 23
%                             uvi = [0,sign(uvr(2))];
%                         case 34
%                             uvi = [sign(uvr(1)),0];
%                         case 41
%                             uvi = [0,sign(uvr(2))];
%                     end
%             else
                    switch newquad
                        case 12
                            uvi = [0,1];
                        case 23
                            uvi = [-1,0];
                        case 34
                            uvi = [0,-1];
                        case 41
                            uvi = [1,0];
                    end
% %                 end
%             end
            return;
        end            
        

% quad = get_quadrant(alpha,tol);

        
            %travelling along grid line with no arc (ie normal to bed contours)
%             newquad = ray.quad;
%             raydir = get_quadrant(alpha,tol);
%             if raydir>4
%                 if newquad==12 && raydir==12
%                     uvi = [0,1];
%                 elseif newquad==23 && raydir==23
%                     uvi = [-1,0];
%                 elseif newquad==34 && raydir==34
%                     uvi = [0,-1];
%                 elseif newquad==41 && raydir==41
%                     uvi = [1,0];
%                 end


%             end

    end

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