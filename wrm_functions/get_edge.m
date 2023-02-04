function edge = get_edge(invector,idx,tol)
    %use coordinates of intesection point to identify which edge it lies on
    if abs(invector(idx,2))<=tol            %y<tol ie appox 0
        edge = 1;                           %x-directed edge
    elseif abs(invector(idx,1))<=tol        %x<tol ie appox 0
        edge = 2;                           %y-directed edge
    else                                    %x~=0 & y~=0
        edge = 3;                           %hypotenuse
    end
end
