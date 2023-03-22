function output = get_inshore_wave(SGo,SGi,Dims,inp,sp)
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
%   Dims - dimensions used for the spectral arrays
%   inp - input parameters wave parameters
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
    fri = Dims.f; xsi = Dims.xsi; depi = Dims.depi; swl = inp.swl;
    if sp.ismodel
        Hso = inp.Hs;
        Tp = inp.Tp;
        Dir0 = inp.Dir;
    else
        [~,idx] = max(inp.dst.Spectra.S);        
        Dir0 = inp.dst.Spectra.Dir(idx);
        Hso = inp.dst.Properties.Hs;
        Tp = 1/inp.dst.Spectra.Dimensions.freq(idx);
    end

    dir_int = abs(median(diff(xsi))); %interval used to interpolate directions (deg)
    radint = deg2rad(dir_int);

    m0 = trapz(radint,abs(trapz(fri,SGo,2))); %integral of offshore spectrum
    m0i = trapz(radint,abs(trapz(fri,SGi,2))); %integral of inshore spectrum
    kw = sqrt(m0i/m0);                %wave transfer coefficient
    [~,idf] = max(SGi,[],'All');     %index of peak energy
    [idir,ifrq] = ind2sub([length(xsi),length(fri)],idf);
    Diripk = xsi(idir);              %direction at inshore peak
    Tpi =1/fri(ifrq);                %period at inshore peak
    ktp = Tpi/Tp;                   %peak period coefficient

    radxso = deg2rad(xsi);
    SGdir = trapz(radxso,abs(trapz(fri,radxso'.*SGi,2)));%direction moment
    Diri = rad2deg(SGdir/m0i);
    kd = Diri-Dir0;

    SGif2 = trapz(radint,abs(trapz(fri,(fri.^2).*SGi,2))); %inshore f^2 moment
    T2i = sqrt(m0i/SGif2);

    SGof2 = trapz(radint,abs(trapz(fri,(fri.^2).*SGi,2)));  %offshore f^2 moment
    T2o = sqrt(m0/SGof2);
    kt2 = T2i/T2o;                  %mean period coefficient

    Hsi = Hso*kw;                %inshore wave height (Hmo)
    % Hmi = 4*sqrt(m0i);              %check of inshore Hs
    % Hmo = 4*sqrt(m0);               %check of offshore Hs
    % disp([Hsi,Hmi,inp.Hso,Hmo])

    output = table(Hsi,T2i,Diri,Tpi,Diripk,kw,kt2,ktp,kd,swl,depi);                                                       
end