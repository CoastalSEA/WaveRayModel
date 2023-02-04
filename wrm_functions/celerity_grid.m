function grid = celerity_grid(grid,T,zwl,delta)
%
%-------function help------------------------------------------------------
% NAME
%   celerity_grid.m
% PURPOSE
%    calculate the celerity and group celerity over a bathymetry grid for
%    a range of wave periods and a range of water levels
% USAGE
%    grid = celerity_grid(grid,T,zwl,delta)
% INPUTS
%   grid.z - [m,n] array of bed elevations relative to zero datum (m)
%   T - [1,p] array of wave periods (s)
%   zwl - [1,q] water level elevations relative to zero datum (m)
%   delta - grid scpaing interval of uniform grid (delx=dely)
% OUTPUTS
%   Adds the following variables to the 'grid' struct
%   h - [m,n,q] array of water depths
%   c  - [m,n,p,q] array of wave celerity using linear wave theory
%   cg - [m,n,p,q] array of wave group celerity using linear wave theory
%   dcx - [m,n,p,q] array of wave celerity gradient in the x-direction
%   dcy - [m,n,p,q] array of wave celerity gradient in the y-direction
% NOTES
%   1) For method see Hunt,ASCE,WW4,1974,p457-459
%   2) Only returns gradients for a grid when delta is specified. However,
%      function can be called using an array, or vector of depths (ie h[m,1]),
%      to get the wave celerity and group celerity
% SEE ALSO
%   celerity.m, refraction.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    [nm,nn] = size(grid.z);
    np = length(T);
    nq = length(zwl);
    %loop over water levels and wave periods for each column (x) of the grid

    c = zeros(nm,nn,np,nq); cg = c; dcx = c; dcy = c; depth = zeros(nm,nn,nq);
    for l=1:nq
        h = zwl(l)-grid.z;
        h(h<0) = 0;
        depth(:,:,l) = h;
        for k=1:np
            Tper = T(k);
            parfor j=1:nn
                c(:,j,k,l) = celerity(Tper,depth(:,j,l));            
                cg_fact = 4*pi()*depth(:,j,l)./(c(:,j,k,l)*Tper);
                cg(:,j,k,l) = c(:,j,k,l)/2.*(1+cg_fact./sinh(cg_fact));            
            end
            if nargin>3
                [dcy(:,:,k,l),dcx(:,:,k,l)] = gradient(c(:,:,k,l),delta);
            end
        end
    end
    grid.c = c; grid.cg = cg; 
    grid.dcx = dcx; grid.dcy = dcy; 
    grid.h = depth;
end