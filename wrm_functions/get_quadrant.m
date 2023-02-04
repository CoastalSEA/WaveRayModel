function [Tri,quad] = get_quadrant(theta)

    %define a triangle polygon for the quadrant of the start point
    if isa(theta,'double')                          %find quad
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
            Tri = []; quad = [];
            return;
        end
    else
        quad = theta;                            %quad is known 
    end

    %define new triangular polyshape based on quadrant being entered
    if quad==1                              %first quadrant
        Tri = polyshape([0,0,1],[0,1,0]);
    elseif quad==2                          %second quadrant
        Tri = polyshape([0,0,-1],[0,1,0]);
    elseif quad==3                          %third quadrant
        Tri = polyshape([0,0,-1],[0,-1,0]);
    elseif quad==4                          %fourth quadrant
        Tri = polyshape([0,0,1],[0,-1,0]);
    else  
        %quadrant not found
        Tri = [];
        return;
    end
end