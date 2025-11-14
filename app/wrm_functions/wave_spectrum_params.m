function [params,diagnost] = wave_spectrum_params(SG,freq,dir,isdirmodes)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum_params.m
% PURPOSE
%   integrate a 2-D spectra to obtain wave parameters
% USAGE
%   paramtable = wave_spectrum_params(SG,f,dir);
%   [paramtable,diagnost] = wave_spectrum_params(SG,f,dir,true);
%   paramtable = wave_spectrum_params(obj)
%   [paramtable,diagnost]  = wave_spectrum_params(obj,true)
% INPUTS
%   SG - array of direction-frequency spectral energy
%   freq - frequency vector
%   dir - direction vector
%   isdirmodes - diagnostic data output if true (optonal)
%   OR
%   obj - instance of the ctWaveSpectra class
%   isdirmodes - diagnostic data output if true (optional)
% OUTPUT
%   params - table containing:
%               Hs - significant wave height
%               m0 - zero moment
%               Dir - mean wave direction
%               Sp - peak spectral density of omni-directional spectrum, S(f)
%               Tp - period at peak energy of omni-directional spectrum
%               Dp - direction at peak energy of omni-directional spectrum
%               Sfdpk - peask spectral density of frequency-direction spectrum, S(theta,f)
%               Tfdpk - period at peak energy of frequency-direction spectrum 
%               Dfdpk - direction at peak energy of frequency-direction spectrum
%               T2 - mean period
%   diagnost - table containing

% SEE ALSO
%   SpectralTransfer.m in WaveRayModel and ctWaveSpectra in CoastalClasses
%   
% Author: Ian Townend
% CoastalSEA (c) March 2023
%--------------------------------------------------------------------------
%
    if nargin<4, isdirmodes = false; diagnost = []; end

    if isa(SG,'ctWaveSpectra') %unpack ctWaveSprectra object
        if nargin<2
            isdirmodes = false; diagnost = []; 
        else
            isdirmodes = freq;
        end        
        obj = SG;
        freq = obj.Spectrum.freq;
        dir = obj.Spectrum.dir;
        SG = obj.Spectrum.SG;
        clear obj
    end
    theta = deg2rad(dir(:));                         %radians and column vector  

    w = trapz_weights_periodic(theta);
    m0 = trapz(freq,sum(SG.*w,1)); 
    Hs = 4*sqrt(m0);                                 %significant wave height
    [Dir,rho,R,W] = mean_direction_global(dir,SG);  %mean direction

    %find peak of the frequency-direction (f-d) spectrum 
    [Sfdpk,idf] = max(SG,[],'All');                  %f-d peak energy and index
    [idir,ifrq] = ind2sub([length(theta),length(freq)],idf);
    Dfdpk = rad2deg(theta(idir));                    %f-d direction at peak
    Tfdpk =1/freq(ifrq);                             %f-d period at peak
    SGf2 = sum(abs(trapz(freq,SG.*(freq.^2),2)).*w); %f^2 moment
    T2 = sqrt(m0/SGf2);

    %get the omni-directional spectrum and values at the peak
    Sf = sum(SG.*w,1);                               %omni-directional spectrum   
    [Sp,ifpk] = max(Sf);                             %peak energy density
    Tp = 1/freq(ifpk);                               %period at peak
    [~,idpk] = max(SG(:,ifpk));                      %direction at peak frequency
    Dp = rad2deg(theta(idpk));
    params = table(Hs,m0,Dir,Sp,Tp,Dp,Sfdpk,Tfdpk,Dfdpk,T2);

    if isdirmodes
        %global mean diagnostic ouput
        diagnost.rho = rho;   diagnost.R = R;   diagnost.W = W;
        %frequency mean directions
        [dir_f,rho_f,R_f,W_f] = mean_direction_per_freq(dir, SG);
        diagnost.dir_f = dir_f;   diagnost.rho_f = rho_f; 
        diagnost.R_f = R_f;       diagnost.W_f = W_f;
        %modes - bearing, energy, and fraction for kpeaks peaks per frequency.
        kpeaks = 2; winlen = 9;
        diagnost.modes = dominant_direction_modes(dir,SG,kpeaks,winlen); 
    end
end

%%
function [dir0,rho,R,W] = mean_direction_global(dir,SG)
    %compute mean direction (compass convention 0-360 degTN) from a 
    %frequency-drection spectrum, SG, with dir direction bins. Returns:
    % dir0 - global mean direction (degTN)
    % rho - measure of concentration. Low rho means mean is unstable or the
    %       the distribution is broad/multimodal. If rho<0.3 examine dominant 
    %       mode(s) instead of a single mean.
    % R - magnitude of the resultant vector obtained by summing all the 
    %     weighted unit vectors that represent your directional distribution.
    %     the weights are spectral densities and as R is not normalizedm it 
    %     is an energy‑weighted resultant length.
    if size(SG,1) ~= numel(dir), SG = SG.'; end
    dir = mod(dir(:)', 360);
    theta = deg2rad(dir);

    X = sin(theta)';     % [1 x i]
    Y = cos(theta)';     % [1 x i]

    w_dir = sum(SG,2);   % total energy per direction [i x 1]
    C = sum(X.*w_dir);   % resultant x
    S = sum(Y.*w_dir);   % resultant y

    R = hypot(C,S);
    W = sum(SG(:));
    rho = R/max(W,eps);
    dir0 = mod(rad2deg(atan2(C,S)),360);
end

%%
function [theta_mean_deg_f, rho_f, R_f, W_f] = mean_direction_per_freq(dir, SG)
    %similar to mean_direction_global but returns results for each frequency
    % dir: compass bearings in degrees, 0° = North, clockwise to 360°
    % SG: spectral density S(dir,f)

    % Ensure orientation [dir,freq]
    if size(SG,1) ~= numel(dir), SG = SG.'; end

    dir = mod(dir(:)',360);
    theta = deg2rad(dir);

    % Compass mapping: x = sin(b), y = cos(b)
    X = sin(theta)';  
    Y = cos(theta)';  

    % Weighted sums over direction for each frequency j
    Cj = X' * SG;                   
    Sj = Y' * SG;                    
    R_f = hypot(Cj, Sj);            
    W_f = sum(SG, 1);          
    rho_f = R_f ./ max(W_f, eps);    

    % Mean bearing in compass convention
    theta_mean_deg_f = mod(rad2deg(atan2(Cj, Sj)), 360).';
    rho_f = rho_f.'; R_f = R_f.'; W_f = W_f.';
end

%%
function modes = dominant_direction_modes(dir,SG,kpeaks,winLen)
    % Returns bearing, energy, and fraction for kpeaks peaks per frequency.
    % dir: compass bearings in degrees, 0° = North, clockwise to 360°
    % SG: spectral density S(dir,f)    
    % winLen: optional integer window length (odd recommended, e.g., 9)
    if nargin < 4 || isempty(winLen), winLen = 9; end

    % Ensure [I x J] orientation
    if size(SG,1) ~= numel(dir), SG = SG.'; end
    i = numel(dir); 
    j = size(SG,2);
    dir = mod(dir(:), 360);

    % Circular smoothing along direction using Hann
    win = hann_window(winLen);                  %Hann requires Signal Processing Toolbox
    Sdf_s = conv2(SG,win,'same');             % smooth columns (dir axis)
    Sdf_s = (Sdf_s+Sdf_s([2:end,1],:))/2;   % light seam blend

    modes(j).bearing_deg = [];
    for j = 1:j
        s = Sdf_s(:,j);
        isPeak = (s>s([end 1:end-1])) & (s>=s([2:end 1]));
        pkIdx = find(isPeak);
        [~,order] = sort(s(pkIdx), 'descend');
        pkIdx = pkIdx(order);
        pkIdx = pkIdx(1:min(kpeaks,numel(pkIdx)));

        modes(j).bearing_deg = dir(pkIdx).';
        modes(j).energy      = s(pkIdx).';
        modes(j).frac        = s(pkIdx).'/sum(s);
    end
end

%%
function w = hann_window(N)
    % Hann window: w[n] = 0.5 - 0.5*cos(2*pi*n/(N-1)), n=0..N-1
    % Normalized so sum(w) = 1
    if N <= 1
        w = 1;
        return
    end
    n = (0:N-1)';
    w = 0.5 - 0.5 * cos(2*pi*n/(N-1));
    w = w / sum(w);
end
