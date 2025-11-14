function [Gspread,dir] = directional_spreading(dir0,nint,nspread,iscos)
% 
%-------function help------------------------------------------------------
% NAME
%   directional_spreading.m
% PURPOSE
%	sample a directional spreading function at selected direction intervals
% USAGE
%   [Gspread,theta] = directional_spreading(dir,nint,nspread,iscos)
% INPUTS
%   dir0 - mean wave direction (degTN)
%   nint - number of directions to sample from spreading function (degTN) 
%          or vector of directions (degTN) 
%   nspread - direction spreading index (-)
%   iscos - true  uses SPM ie cosine function;
%           false uses Donelan ie secant function;
% OUTPUT
%   Gspread - direction distribution for given mean direction
%   dir = the angles used to define the function (degTN)
% NOTES
%   Donelan, Hamilton and Hui, R.Soc, 1985, A315, 509-562
%   US Army Corps, Shore Protection Manual (SPM), 1984
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
        theta = linspace(0,pi/2,nint/2); 
        theta = [fliplr(-theta),theta(2:end)];   %ensures range centred on 0
    else
        theta = deg2rad(nint-dir0);              %user defined intervals        
    end
    
    if iscos        
        gfun = @(ang) cos(gamma*ang).^nspread;   %SPM cosine function
    else
        gfun = @(ang) sech(beta*ang).^nspread;   %Donelan secant function
    end

    G = gfun(theta');
    % Gspread = G./trapz(angles,G);              %normalise by function integral
    if ~isscalar(nint)
        bound = [-pi/2,pi/2];
        idx = isangletol(theta,bound);           %only use angles within +/-pi/2 of dir
        G(~idx) = 0;
    end

    w = trapz_weights_periodic(theta);
    Gspread = G./sum(G.*w);                      %normalise by function integral
    dir = mod(rad2deg(theta)+dir0,360);
    %checks 
    % sum(Gspread.*w)                            %integral should =1
    % figure; plot(dir,Gspread)
end