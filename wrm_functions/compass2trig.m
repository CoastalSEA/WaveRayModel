function [dirfm,dirto] = compass2trig(theta,israd)
%
%-------function help------------------------------------------------------
% NAME
%   compass2trig.m
% PURPOSE
%    Convert compass directions to trigonometric angles or vice versa
% USAGE
%    [dirfm,dirto] = compass2trig(theta,israd);
% INPUTS
%   theta - direction in degrees True North, or trigonometric angle in radians
%   israd - logical true if input is in radians (optional, default=false)
% OUTPUTS
%   dirfm - angle vector is "from" in radians if input is degrees and vice versa  
%   dirto - angle vector is "to" in radians if input is degrees and vice versa
% NOTES
%   Thanks to Mayra  
%   Matlab Forum: 383063-convert-wind-direction-in-true-north-to-math-convention?s_tid=ta_ans_results
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%----------------------------------------------------------------------
%
    if nargin<2 || ~israd
        theta = deg2rad(theta);
        israd = false;
    end

    [u,v] = pol2cart(theta,1);  %cartesian co-ordinates
    %swapping u,v sets theta to counterclockwise and 0° at right.
    [dirfm,~] = cart2pol(v,u);           %vector direction "from"
    dirto = mod(dirfm+pi,2*pi);          %vector direction "to"

    if israd
        dirfm = rad2deg(dirfm);          %output in degrees if input is in radians
        dirto = rad2deg(dirto);          
    end
end