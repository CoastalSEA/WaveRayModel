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
    p0 = wave_spectrum_params(SGo,freq,dir); %parameters of offshore spectrum
    pi = wave_spectrum_params(SGi,freq,dir); %parameters of inshore spectrum

    %input parameters
    if sp.ismodel
        Hso = intable.Hs;
        Tpo = intable.Tp;
        Dir0 = intable.Dir;
    else
        Hso = intable.Hs;     
        Tpo = p0.Tp;
        Dir0 = p0.Dir0;
    end
    T2o = p0.T2;

    %transfer coefficients
    kw = sqrt(pi.m0/p0.m0);
    Diri = pi.Dir0;
    Diripk = pi.Dirpk;
    Tpi = pi.Tp;
    T2i = pi.T2;
    ktp = Tpi/Tpo;                    %peak period coefficient
    kd = Diri-Dir0;                   %direction shift
    kt2 = T2i/T2o;                    %mean period coefficient

    Hsi = Hso*kw;                      %inshore wave height (Hmo)
    % Hmi = 4*sqrt(m0i);               %check of inshore Hs
    % Hmo = 4*sqrt(m0);                %check of offshore Hs
    % disp([Hsi,Hmi,inp.Hso,Hmo])

    output = table(Hsi,T2i,Diri,Tpi,Diripk,kw,kt2,ktp,kd,swl,depi);                                                       
end