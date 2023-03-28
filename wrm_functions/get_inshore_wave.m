function output = get_inshore_wave(SGo,SGi,Dims,intable,sp)
%
%-------function help------------------------------------------------------
% NAME
%   get_inshore_wave.m
% PURPOSE
%   integrate the 2-D spectra to obtain wave parameters and transfer coefficients  
% USAGE
%   output = get_inshore_wave(SGo,SGi,Dims,inp)
% INPUTS
%   SGo - array of offshore direction-frequency spectral energy
%   SGi - array of inshore direction-frequency spectral energy
%   Dims - dimensions used for the spectral arrays, struct with fields:
%          freq - interpolated frequency vector (1/1-30s)
%          beta - interpolated direction vector (0-360 deg)
%          depi - inshore or minimum depth
%   intable - input table of wave parameters
%   sp - defines type of model to use and model parameters
% OUTPUT
%   output - table of wave statistical parameters and ratios including:
%            Hsi,T2i,Diri,Tpi,Diripk,kw,kt2,ktp,kd,swl,depi                                                            
% SEE ALSO
%   get_inshore_spectrum.m and SpectralTransfer.m in WaveRayModel
%   
% Author: Ian Townend
% CoastalSEA (c) March 2023
%--------------------------------------------------------------------------
%
    if isempty(SGo)
        varnames = {'Hsi','T2i','Diri','Tpi','Diripk','kw','kt2','ktp','kd','swl','depi'};
        nans = num2cell(NaN(1,11));
        output = table(nans{:},'VariableNames',varnames);
        return;
    end   

    %spectrum dimensions and parameter settings
    freq = Dims.freq; dir = Dims.dir; depi = Dims.depi; swl = intable.swl;
    dir_int = abs(median(diff(dir))); %interval used to interpolate directions (deg)
    radint = deg2rad(dir_int);
    raddir = deg2rad(dir);
    m0 = trapz(radint,abs(trapz(freq,SGo,2)));  %integral of offshore spectrum
    m0i = trapz(radint,abs(trapz(freq,SGi,2))); %integral of inshore spectrum
    SGdiro = trapz(raddir,abs(trapz(freq,raddir'.*SGo,2)));%offshore direction moment
    SGdiri = trapz(raddir,abs(trapz(freq,raddir'.*SGi,2)));%inshore direction moment

    %input parameters
    if sp.ismodel
        Hso = intable.Hs;
        Tp = intable.Tp;
        Dir0 = intable.Dir;
    else
        Hso = intable.Hs;     
        Dir0 = rad2deg(SGdiro/m0);
        %[~,idd] = min(abs(intable.Dir-Dir0));
        %mnTp =1/sp.freq(idd);
        %check Hs by integrating spectral energy: 4.sqrt(mo)
        % cfHso = 4*sqrt(m0);
        [~,idx] = max(intable.S); 
        Tp = 1/sp.freq(idx);
        %pkDir0 = intable.Dir(idx);        
    end
    
    %transfer coefficients
    kw = sqrt(m0i/m0);                %wave transfer coefficient
    [~,idf] = max(SGi,[],'All');      %index of peak energy
    [idir,ifrq] = ind2sub([length(dir),length(freq)],idf);
    Diripk = dir(idir);               %direction at inshore peak
    Tpi =1/freq(ifrq);                %period at inshore peak
    ktp = Tpi/Tp;                     %peak period coefficient
    Diri = rad2deg(SGdiri/m0i);       %inshore mean direction
    kd = Diri-Dir0;                   %direction shift

    SGif2 = trapz(radint,abs(trapz(freq,(freq.^2).*SGi,2))); %inshore f^2 moment
    T2i = sqrt(m0i/SGif2);

    SGof2 = trapz(radint,abs(trapz(freq,(freq.^2).*SGi,2))); %offshore f^2 moment
    T2o = sqrt(m0/SGof2);
    kt2 = T2i/T2o;                     %mean period coefficient

    Hsi = Hso*kw;                      %inshore wave height (Hmo)
    % Hmi = 4*sqrt(m0i);               %check of inshore Hs
    % Hmo = 4*sqrt(m0);                %check of offshore Hs
    % disp([Hsi,Hmi,inp.Hso,Hmo])

    output = table(Hsi,T2i,Diri,Tpi,Diripk,kw,kt2,ktp,kd,swl,depi);                                                       
end