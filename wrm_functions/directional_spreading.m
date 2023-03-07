function [Gspread,theta] = directional_spreading(dir,nint,nspread,iscos)
% 
%-------function help------------------------------------------------------
% NAME
%   directional_spreading.m
% PURPOSE
%	sample a directional spreading function at selected direction intervals
% USAGE
%   [Gspread,theta] = directional_spreading(dir,nint,nspread,iscos)
% INPUTS
%   dir - mean wave direction (degTN)
%   nint - number of directions to sample from spreading function (degTN) 
%          or vector of directions (degTN) 
%   nspread - direction spreading index (-)
%   iscos - true  uses SPM ie cosine function;
%           false uses Donelan ie secant function;
% OUTPUT
%   Gspread - direction distribution for given mean direction
%   theta = the angles used to define the function (degTN)
% NOTES
%   Donelan, Hamilton and Hui, R.Soc, 1985, A315, 509-562
%   US Army Corps, Shore Protection Manual, 1984
% SEE ALSO
%   wave_spectrum.m. Used in WaveRayModel
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%---------------------------------------------------------------------------
%
    beta = 2.28;                %frequency dependent coefficient (from Donelan et al)
    gamma  = 1;                 %SPM scaling coefficient

    if isscalar(nint)
        if rem(nint,2)>0, nint=nint+1; end
        angles = linspace(0,pi/2,nint/2); 
        angles = [fliplr(-angles),angles(2:end)]; %ensures range centred on 0
    else
        angles = deg2rad(nint-dir);               %user defined intervals        
    end
    
    if iscos        
        gfun = @(ang) cos(gamma*ang).^nspread;    %SMP cosine function
    else
        gfun = @(ang) sech(beta*ang).^nspread;    %Donelan secant function
    end

    G = gfun(angles);
    Gspread = G./trapz(angles,G);                 %normalise by function integral

    theta = mod(rad2deg(angles)+dir,360);

    if ~isscalar(nint)
        bound = [-pi/2,pi/2];
        idx = isangletol(angles,bound);           %only use angles within +/-pi/2 of dir
        Gspread(~idx) = 0;
    end
end