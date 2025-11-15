function stats = setSpectrum(obj,obsfreq)
%
%-------function help------------------------------------------------------
% NAME
%   get_inshore_spectrum.m
% PURPOSE
%   construct the offshore and inshore spectra for given wave conditions or
%   wave buoy spectral data
% USAGE
%   [SGo,SGi,Dims] = get_inshore_spectrum(sptobj,intable,sp)
% INPUTS
%   obj - instance of a ctWaveSpectra class
%   obsfreq - the frequency bins to be used for the ouput
% OUTPUT
%   stats - a table with variables: 'S','Dir','Spr','Skew','Kurt' and a row
%           for each frequency bin
% NOTES
%   see test_setSpectrum for the code to read a spectrum file, create a
%   spectrum and then call setSpectrum to decompose spectrum into
%   statistics across a set of frequency bins (f: S,Dir,Spr,Skew,Kurt).
%   A stand-alone version with alternative methods can be run using
%   test_get_set_spectrum.m
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2025
%--------------------------------------------------------------------------
%
    %reduce a detailed spectrum to the format defined by the Datawell buoy
    %spt file format
    
    if ~isa(obj,'ctWaveSpectra')
        warndlg('Incorrect input to setSpectrum'); 
        stats = []; return;
    end          
    SGin = obj.Spectrum.SG;              %spectrum to be saved
    dir = obj.Spectrum.dir;        
    freq = obj.Spectrum.freq;        
    
    % interpolate spectrum to 64 frequency intervals for spt format
    SG = zeros(size(SGin,1),length(obsfreq));
    for i = 1:size(SGin,1)
        SG(i,:) = interp1(freq, SGin(i,:), obsfreq, 'pchip', 0);                
    end 

    stats = smoothedMomentEstimates(dir,SG,0,'none');
end

%%
function stats = smoothedMomentEstimates(dir, SG, winLen, method)
    %use a smoothed vector-sum method of moments to estimate the statistical
    %properties
    % Smooth m1, m2 along frequency to stabilize mean & spread.
    % method: 'movavg', 'sgolay', or 'adsmooth'
    if nargin < 3 || isempty(winLen), winLen = 7; end
    if nargin < 4 || isempty(method), method = 'movavg'; end

    % Ensure [I x J]
    if size(SG,1) ~= numel(dir), SG = SG.'; end
    theta = deg2rad(mod(dir(:),360));  % [I x 1]
    wp = trapz_weights_periodic(theta);
    Sf = sum(SG.*wp);

    W = sum(SG,1);                          % Total energy per frequency
    Wsafe = max(W, eps);
    P = SG ./ Wsafe;                        % normalized weights per freq

    % Trigonometric moments (vectorized)
    m1 = sum(P .* exp(1i*theta), 1);         % [1 x J]
    m2 = sum(P .* exp(2i*theta), 1);         % [1 x J]

    % Smooth complex series along frequency
    switch lower(method)
        case 'movavg'
            win = ones(1,winLen) / winLen;
            m1s = conv(m1, win, 'same');
            m2s = conv(m2, win, 'same');
        case 'sgolay'
            % SG: odd window, poly order <= window-1
            polyOrder = min(3, winLen-2);
            [~, g] = sgolay(polyOrder, winLen);
            m1s = conv(m1, g(:,1)', 'same');
            m2s = conv(m2, g(:,1)', 'same');
        case 'adsmooth'
            m1s = adapt_smooth_complex(m1, abs(m1), winLen);
            m2s = adapt_smooth_complex(m2, abs(m1), winLen); % tie smoothing to m1 concentration
        otherwise
            m1s = m1; m2s = m2;
    end

    % Mean direction and concentration
    mu   = angle(m1s);             % mean direction (radians)
    rho  = abs(m1s);               % concentration
    stdc = sqrt(-2*log(max(rho,eps)));

    % Spread components via m2 relative to mu
    skew = imag(m2s .* exp(-2i*mu)) ./ max((1-rho).^(3/2), eps);
    kurt = real(m2s .* exp(-2i*mu)) ./ max((1-rho).^2,   eps);

    %store results
    sts.sf = Sf';
    sts.mean = mod(rad2deg(mu(:)),360);
    % sts.rho  = rho(:);
    % sts.std  = stdc(:);          % radians
    sts.spr = mod(rad2deg(stdc(:)),360);
    sts.skew = skew(:);
    sts.kurt = kurt(:);
    stats = struct2table(sts);
    stats.Properties.VariableNames =  {'S','Dir','Spr','Skew','Kurt'};
end

%%
function y = adapt_smooth_complex(x, rho, baseLen)
    % x: [1 x J] complex series, rho: [1 x J], baseLen: min window
    J = numel(x); y = zeros(1,J);
    for j = 1:J
        k = max(baseLen, round(baseLen + 10*(1 - rho(j))));
        h = max(1, floor(k/2));
        i0 = max(1, j-h); i1 = min(J, j+h);
        w = ones(1, i1-i0+1);
        y(j) = sum(x(i0:i1) .* w) / sum(w);
    end
end