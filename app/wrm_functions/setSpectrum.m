function stats = setSpectrum(obj,obsfreq,iserrplt)
    %reduce a detailed spectrum to the format defined by the Datawell buoy
    %spt file format
    
    if ~isa(obj,'ctWaveSpectra')
        warndlg('Incorrect input to setSpectrum'); 
        stats = []; diagnost = []; return;
    end          
    SGin = obj.Spectrum.SG;              %spectrum to be saved
    dir = obj.Spectrum.dir;        
    freq = obj.Spectrum.freq;        
    theta = deg2rad(dir(:));               %force column vector
    
    % interpolate spectrum to 64 frequency intervals for spt format
    SG = zeros(size(SGin,1),length(obsfreq));
    for i = 1:size(SGin,1)
        SG(i,:) = interp1(freq, SGin(i,:), obsfreq, 'pchip', 0);                
    end 

    stats = smoothedMomentEstimates(dir,SG,0,'none');
end

% %%
% function [stats,diagnost] = harmonicEstimates(theta,SG,iserrplt)
%     %use the first harmonics of the spectrum to estimate the statistical
%     %properties
% 
%     %tolerances and constants
%     tol_m1 = eps;         %tolerances for edge-case handling
%     tol_sigma = 1e-8;
%     %windowSize = 5;       %window for smoothing in directional dimension
% 
%     theta = theta(:);   %force column vector
%     w = trapz_weights_periodic(theta);
%     fact = 1/pi;  % model normalization factor from your D expression
%     nfreq = size(SG,2);
%     A1=zeros(1,nfreq); B1=A1; A2=A1; B2=A1; theta0=A1; m1=A1; m2=A1; n2=A1; 
%     Sf=A1; Dir=A1; Spr=A1; Skew=A1; Kurt=A1; rms=A1;
% 
%     %extract direction distribution, get harmonics and estimate parameters
%     for f=1:nfreq
%         Gf = SG(:,f);
%         % Gf = movmean(Gf, windowSize, 'omitnan'); %smoothing option
%         % spectral density per frequency is integral over theta of D (if D normalized to 1 across theta, S = 1)
%         Sf(1,f) = sum(Gf.*w);
%         % Gf = Gf/Sf(1,f);
%         % projections (use D as given; if D is already normalized to S you may scale)
%         A1(f) = fact * sum(Gf.*cos(theta).*w);
%         B1(f) = fact * sum(Gf.*sin(theta).*w);
%         A2(f) = fact * sum(Gf.*cos(2*theta).*w);
%         B2(f) = fact * sum(Gf.*sin(2*theta).*w);
% 
%         % first harmonic amplitude and mean direction
%         theta0(f) = mod(atan2(B1(f),A1(f)),2*pi);   % radians
%         m1(f) = sqrt(A1(f)^2+B1(f)^2);
% 
%         % debug option to check Phase Alignment of Directional Spectrum
%         % NB produces a plot for every frequency
%         % checkPhaseAlignment(theta,theta0(f));
% 
%         % rotate second harmonic by -2*theta0 to express about (theta-theta0)
%         A2r =  A2(f)*cos(2*theta0(f)) + B2(f)*sin(2*theta0(f));
%         B2r = -A2(f)*sin(2*theta0(f)) + B2(f)*cos(2*theta0(f));
%         m2(f) = A2r;
%         n2(f) = B2r;
% 
%         % invert to physical parameters with edge-case handling
%         if m1(f)<tol_m1
%             % first harmonic vanished => mean direction undefined, spread fixed by m1
%             Dir(f) = NaN;
%         else
%             Dir(f) = rad2deg(theta0(f));
%         end
% 
%         % sigma from m1: sigma = sqrt(2*(1-m1))
%         sigma = sqrt(max(0,2*(1-m1(f))));
%         Spr(f) = rad2deg(sigma);
%         %fits a smooth unimodal model to the directional spectrum, to
%         %provide a robust estimate of spread even when (m1) is very small
%         % [mu_hat, kappa_hat] = fit_vonmises(theta, Gf, w);
%         % sigma = sqrt(-2*log(besseli(1, kappa_hat) / besseli(0, kappa_hat)));
%         % Spr(f) = mod(rad2deg(sigma),360);
% 
%         if sigma<=tol_sigma
%             Kurt(f) = NaN;   % singular: division by sigma^4
%         else
%             Kurt(f) = (2*m2(f)+ 6-8*m1(f))/(sigma^4);
%         end
% 
%         denom = 0.5*(1-m2(f));
%         if denom<=0
%             Skew(f) = NaN;   % invalid / complex
%         else
%             Skew(f) = - n2(f)/(denom^(3/2));
%         end     
%         %rms error
%         Ghat = (1/pi)*(0.5+m1(f)*cos(theta-theta0(f))+...
%                m2(f)*cos(2*(theta-theta0(f)))+n2(f)*sin(2*(theta-theta0(f))));
%         Ghat_scaled = Ghat*Sf(f);        
%         rms(f) = sqrt(mean((Gf-Ghat_scaled).^2));
%     end
%     varnames = {'S','Dir','Spr','Skew','Kurt'};
%     stats = table(Sf',Dir',Spr',Skew',Kurt','VariableNames',varnames);
%     diagnost = struct('A1',A1,'B1',B1,'A2',A2,'B2',B2,'theta0',theta0,'m1',m1,'m2',m2,'n2',n2,'rms',rms);
% 
%     if iserrplt
%         figure('Tag','PlotFig')
%         plot(1:nfreq,rms,'x')
%         title('RMS Error')
%         xlabel('Frequency index')
%         ylabel('RMS error')
%     end
% end
% 
% %%
% function stats = momentEstimates(dir,SG)
%     %use the vector-sum method of moments to estimate the statistical
%     %properties
%     % dir: [1 x I] compass bearings (0°=North, clockwise)
%     % SG:  [I x J] spectrum (dir x freq)
%     % stats: struct with fields mean, rho, std, skew, kurt, each [J x 1]
% 
%     if size(SG,1) ~= numel(dir)
%         SG = SG.'; %pure transpose without modifying the complex values
%     end
% 
%     theta = deg2rad(mod(dir(:),360)); % [I x 1]
%     w = trapz_weights_periodic(theta);
%     nfreq = size(SG,2);
% 
%     sts.sf = zeros(1,nfreq);
%     sts.mean = zeros(1,nfreq);
%     % m1f  = zeros(1,nfreq);
%     % m2f  = zeros(1,nfreq);
%     sts.spr  = zeros(1,nfreq);
%     sts.skew = zeros(1,nfreq);
%     sts.kurt = zeros(1,nfreq);
% 
%     for f = 1:nfreq
%         Gf = SG(:,f);       
%         Sf = sum(Gf.*w);
% 
%         W = sum(Gf);
%         if W <= 0
%             sts.mean(f) = 0; continue;
%         end
%         p = Gf / W;
% 
%         % moments
%         m1 = sum(p .* exp(1i*theta));
%         m2 = sum(p .* exp(2i*theta));
% 
%         % mean direction and concentration
%         mu = angle(m1);
%         rho = abs(m1);
% 
% % m1f(f) = m1;
% % m2f(f) = m2;
% 
%         % spread
%         circ_std = sqrt(-2*log(max(rho,eps)));
% 
%         % skewness and kurtosis
%         skew = imag(m2 * exp(-2i*mu)) / max((1-rho)^(3/2), eps);
%         kurt = real(m2 * exp(-2i*mu)) / max((1-rho)^2, eps);
% 
%         % store
%         sts.sf(f) = Sf; 
%         sts.mean(f) = mod(rad2deg(mu),360);
%         % sts.rho(f)  = rho;
%         % sts.std(f)  = circ_std;
%         sts.spr(f) = mod(rad2deg(circ_std),360);
%         sts.skew(f) = skew;
%         sts.kurt(f) = kurt;
%     end
%     stats = struct2table(sts);
%     stats.Properties.VariableNames =  {'S','Dir','Spr','Skew','Kurt'};
% 
% 
% %figure('Tag','PlotFig'); plot(1:nfreq,m1f,1:nfreq,m2f);
% 
% end

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

% %%
% function [mu_hat, kappa_hat] = fit_vonmises(theta, D, w)
%     % theta: angular vector in radians
%     % D: directional spectrum values (normalized)
%     % w: integration weights (e.g., trapz_weights_periodic(theta))
% 
%     % Compute first trigonometric moment
%     C = sum(D .* cos(theta) .* w);
%     S = sum(D .* sin(theta) .* w);
%     R = sqrt(C^2 + S^2);
%     mu_hat = atan2(S, C);
% 
%     % Estimate kappa from R using approximation
%     if R < 0.53
%         kappa_hat = 2*R + R^3 + (5*R^5)/6;
%     elseif R >= 0.53 && R < 0.85
%         kappa_hat = -0.4 + 1.39*R + 0.43/(1 - R);
%     else
%         kappa_hat = 1/(R^3 - 4*R^2 + 3*R);
%     end
% 
%     % Optional: refine kappa by numerical optimization (not shown)
% end
% 
% %%
% function checkPhaseAlignment(theta,theta0)
%     %Check Phase Alignment of Directional Spectrum
%     % Rotate angles
%     dtheta = mod(theta-theta0(f)+pi, 2*pi)-pi;
%     % Interpolate spectrum onto rotated angles
%     D_rot = interp1(dtheta, D, theta, 'linear', 0);
%     % Plot original and rotated spectrum
%     figure('Name','Direction spread','Tag','PlotFig');
%     rd = mod(rad2deg(theta),360);
%     plot(rd, D, 'b-', rd, D_rot, 'r--');
%     xlabel('Angle (deg)'); ylabel('Directional Spectrum');
%     title(sprintf('Frequency %d: Original (blue) vs Rotated (red) Spectrum, theta0=%.1f',f,rad2deg(theta0(f))));
%     legend('Original', 'Rotated');
% end