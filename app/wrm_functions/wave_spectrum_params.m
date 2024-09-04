function params = wave_spectrum_params(SG,freq,dir)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum_params.m
% PURPOSE
%   integrate a 2-D spectra to obtain wave parameters
% USAGE
%   [Hs,Dir] = wave_spectrum_params(SG,f,dir)
% INPUTS
%   SG - array of offshore direction-frequency spectral energy
%   freq - frequency vector
%   dir  - direction vector
% OUTPUT
%   params - table containing:
%               Hs - significant wave height
%               m0 - zero moment
%               Dir0 - mean wave direction
%               Dirpk - Direction of peak energy
%               Tp - period peak energy
%               T2 - mean period
% SEE ALSO
%   get_inshore_spectrum.m and SpectralTransfer.m in WaveRayModel
%   
% Author: Ian Townend
% CoastalSEA (c) March 2023
%--------------------------------------------------------------------------
%
    if isrow(dir), dir = dir'; end                   %force a column vector 
    if iscolumn(freq), freq = freq'; end             %force a row vector 

    raddir = deg2rad(dir);
    m0 = trapz(raddir,abs(trapz(freq,SG,2)));        %zero moment
    Hs = 4*sqrt(m0);                                 %significant wave height
    SGdir = trapz(raddir,abs(trapz(freq,raddir.*SG,2))); %direction weighted moment
    Dir0 = mod(rad2deg(SGdir/m0),360);               %mean direction
    %find peak of spectrum 
    [~,idf] = max(SG,[],'All');                      %index of peak energy
    [idir,ifrq] = ind2sub([length(dir),length(freq)],idf);
    Dirpk = dir(idir);                               %direction at inshore peak
    Tp =1/freq(ifrq);                                %period at inshore peak
    SGf2 = trapz(raddir,abs(trapz(freq,(freq.^2).*SG,2))); %offshore f^2 moment
    T2 = sqrt(m0/SGf2);

    params = table(Hs,m0,Dir0,Dirpk,Tp,T2);  
end