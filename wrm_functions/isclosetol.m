function [isclose,dist] = isclosetol(pnts1,pnts2,tol)
%
%-------function help------------------------------------------------------
% NAME
%   isclosetol.m
% PURPOSE
%    boolean check of whether an x,y point lies close to another x,y point
% USAGE
%    isdir = isclosetol(theta,bound)
% INPUTS
%   pnts1 - single x,y point, or [nx2] vectors of x,y coordinates of points
%   pnts2 - [nx2] vectors of x,y coordinates of points to be checked
%   tol - distance to be used as tolerance to define close
% OUTPUTS
%   isclose - logical, true if points are within specified tolerance
%   dist = distance between points
% SEE ALSO
%   get_edge, get_quadrant, next_element and arc_ray.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if (size(pnts1,2)~=2 || size(pnts2,2)~=2) || ...
       (size(pnts1,1)>1 && any(size(pnts1)~=size(pnts2)))
        warndlg('Points must be [Nx2] arrays')
        isclose = []; dist = [];
        return;       
    end

    dist = sqrt((pnts2(:,1)-pnts1(:,1)).^2 + (pnts2(:,2)-pnts1(:,2)).^2);
    isclose = dist<=tol;
end