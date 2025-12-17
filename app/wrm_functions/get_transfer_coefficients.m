function output = get_transfer_coefficients(offobj,inobj)
%
%-------function help------------------------------------------------------
% NAME
%   get_transfer_coefficients.m
% PURPOSE
%   compute transfer coefficients from offshore and inshore parameters
% USAGE
%   output = get_transfer_coefficients(offobj,inobj)
% INPUTS
%   offobj - instance of ctWaveSpectra that defines offshore spectrum
%   inobj - instance of ctWaveSpectra of the associated inshore spectrum
% OUTPUT
%   output - table of wave transfer coefficients for:
%            kw - wave height, kt2 - T2 wave period, 
%            ktp - peak wave period, kd - wave direction                                                          
% SEE ALSO
%   get_inshore_spectrum.m and SpectralTransfer.m in WaveRayModel.
%   Replaces get_inshore_wave.m to use ctWaveSpectra as inputs and only 
%   return the inshore transfer coefficients (inshore wave parameters are
%   held in the inobj.Params property of the ctWaveSpectra class)
%   
% Author: Ian Townend
% CoastalSEA (c) Dec 2025
%--------------------------------------------------------------------------
%
    if isempty(offobj.Spectrum.SG)
        varnames = {'kw','kt2','ktp','kd'};
        nans = num2cell(NaN(1,4));
        output = table(nans{:},'VariableNames',varnames);
        return;
    end   

    %spectrum dimensions and parameter settings
    p0 = offobj.Params;   %parameters of offshore spectrum
    pi = inobj.Params;    %parameters of inshore spectrum
    % inp = offobj.inpData;
    % fprintf('Hs = %.2f/%.2f, Dir = %.1f/%.1f\n',inp.Hs,p0.Hs,inp.Dir,p0.Dir)

    %input parameters
    Tpo = p0.Tp;
    Dir0 = p0.Dir;
    T2o = p0.T2;

    %transfer coefficients
    kw = sqrt(pi.m0/p0.m0);
    Diri = pi.Dir;
    Tpi = pi.Tp;
    T2i = pi.T2;
    ktp = Tpi/Tpo;                    %peak period coefficient
    kd = Diri-Dir0;                   %direction shift
    kt2 = T2i/T2o;                    %mean period coefficient
    output = table(kw,kt2,ktp,kd);   
end