function cmesh = celerity_mesh(cmesh,T,zwl)
%
%-------function help------------------------------------------------------
% NAME
%   celerity_mesh.m
% PURPOSE
%    calculate the celerity and group celerity over a bathymetry mesh for
%    a range of wave periods and a range of water levels
% USAGE
%    cmesh = celerity_mesh(cmesh,T,zwl,delta)
% INPUTS
%   cmesh.z - [m,n] array of bed elevations relative to zero datum (m)
%   T - [1,p] array of wave periods (s)
%   zwl - [1,q] water level elevations relative to zero datum (m)
%   delta - grid spacing interval of uniform grid (delx=dely)
% OUTPUTS
%   Adds the following variables to the 'cmesh' struct
%   h - [m,q] array of water depths
%   c  - [m,p,q] array of wave celerity using linear wave theory
%   cg - [m,p,q] array of wave group celerity using linear wave theory
%   dcx - [m,p,q] array of wave celerity gradient in the x-direction
%   dcy - [m,p,q] array of wave celerity gradient in the y-direction
% NOTES
%   1) For method see Hunt,ASCE,WW4,1974,p457-459
%   2) Only returns gradients for a grid when delta is specified. However,
%      function can be called using an array, or vector of depths (ie h[m,1]),
%      to get the wave celerity and group celerity
% SEE ALSO
%   uses celerity.m, cf refraction.m. see test_wavemodel.m for examples of
%   the two uses (grid and vector).
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    nm = length(cmesh.z);
    np = length(T);
    nq = length(zwl);
    %loop over water levels and wave periods for each column (x) of the grid
    pts = cmesh.Tri.Points;
    tria = cmesh.Tri.ConnectivityList;
    c = zeros(nm,np,nq); cg = c; dcx = c; dcy = c; depth = zeros(nm,nq);
    for l=1:nq                       %water levels
        h = zwl(l)-cmesh.z;
        h(h<0) = 0;
        depth(:,l) = h;
        for k=1:np                   %wave periods
            Tper = T(k);
            c(:,k,l) = celerity(Tper,h);            
            cg_fact = 4*pi()*h./(c(:,k,l)*Tper);
            cg(:,k,l) = c(:,k,l)/2.*(1+cg_fact./sinh(cg_fact));            
            [dcx(:,k,l),dcy(:,k,l)] = trigradient(pts(:,1),pts(:,2),...
                                                   c(:,k,l), tria);
        end
    end
    cmesh.c = c;           cmesh.h = depth;
    cg(isnan(cg)) = 0;     dcx(isnan(dcx)) = 0;    dcy(isnan(dcy)) = 0; 
    cmesh.cg = cg;         cmesh.dcx = dcx;        cmesh.dcy = dcy; 
end