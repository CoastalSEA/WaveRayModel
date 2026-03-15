function [S,gamma] = wave_spectrum(stype,f,inputs)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum.m
% PURPOSE
%   Calculate a f-D wave spectrum using selected model (Bretschneider open ocean, 
%   Pierson-Moskowitz fully developed, JONSWAP fetch limited, and 
%   TMA shallow water).
% USAGE
%   S = wave_spectrum(stype,f,inputs)
% INPUTS
%   stype - type of spectrum. options are 'Bretschneider open ocean', 
%           'Pierson-Moskowitz fully developed', 'JONSWAP fetch limited', 
%           or 'TMA shallow water'
%   f - frequencies at which the spectrum is to be defined <scalar or vector>
%  inputs - a struct that depends on use. When deriving the spectrum using 
%           wind speed and fetch the struct includes 'wind',Uw,zW,Fetch,df
%           When using wave records the fields are 'wave',Hmo,Tp,gamma,
%   For wind input (JONSWAP and TMA only)
%       source - 'Wind'
%       Uw - wind speed (m/s) <scalar or vector>
%       zW - elevation of wind speed measurement (m)
%       Fetch- dominant fetch length (m) <scalar or vector>
%       df - average water depth over fetch (m) (default is deep water)
%            <scalar or vector> but must be same length as Uw and Fetch
%   For wave input
%       source - 'Wave'
%       Hs - significant wave height [4sqrt(mo)] (m) <scalar or vector>
%       Tp - peak wave period (s)  <can be vector - same length as Hs>
%       T2 - zero upcrossing wave period (sqrt(m0/m2)): used to calculate
%            gamma, or gamma can be specified
%       gamma - spectrum shape parameter (defaults to 3.3 if input as 0 
%               unless T2 is also included)
%
%   When using 'TMA shallow water', wind or wave struct includes 
%   ds - for tma wave: depths at site (m) <scalar or vector> 
%        default is deep water if ds is not included, isempty or zero
%        NB: wind or wave vector data must be the same length
% OUTPUT
%   S - spectral energy density at specified frequencies [nrec,nf] (m^2s)
%   gamma - spectrum shape parameter in case set in TMA or Jonswap
% NOTES
%   The spectrum definitions in Carter, 1982, MIAS publication No.4 are
%   used with additional information from Hughes, CERC-84-7 (p10-12); 
%   Bouws et al, 1985; Donelan, Hamilton and Hui, R.Soc, 1985, eqn 6.1) and
%   Hunt,ASCE,WW4,1974,p457-459
%   NB: 1) gamma only used for JONSWAP and TMA spectra
%       2) Pierson-Moskowitz can use just Hmo or just Tp. If specifying Tp
%          enter do not include Hmo, or assign [] and vice versa.
% SEE ALSO
%   uses celerity.m. cf tma_spectrum.m which outputs [Hs,Tp,Tz]
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2023
%--------------------------------------------------------------------------
%   
    S = []; gamma = [];
    %ensure legacy compatability and format used in ctWaveSpectrum
    if isfield(inputs,'input'), inputs.source = inputs.input; end

    if strcmp(inputs.source,'Wind')            %check for valid inputs
        if inputs.Uw==0 || inputs.Fetch==0, return; end 
        if inputs.zW==0, inputs.zW = 10; end
    else
        if inputs.Hs==0 || inputs.Tp==0, return; end    
    end
    
    switch stype
        case 'Bretschneider open ocean'
            S = bretschneider(f,inputs);
        case 'Pierson-Moskowitz fully developed'
            S = pierson_moskowitz(f,inputs);
        case 'JONSWAP fetch limited'
            [S,gamma] = jonswap(f,inputs);
        case 'TMA shallow water'
            [S,gamma] = tma(f,inputs);
        otherwise
            warndlg('Unknown spectrum type')
    end

    %checks
    % m0 = trapz(f,S);     %first moment
    % Hs = 4*sqrt(m0);     %significant wave height
    % [Sp,ifpk] = max(S);  %peak energy density
    % Tp = 1/f(ifpk);      %period at peak
end

%%
function S = bretschneider(f,inp)
    %construct the Bretschneider spectrum
    if ~strcmp(inp.source,'Wave')
        warndlg('Bretschneider option only available with wave type input')
        S = []; return;
    end
    Hmo = inp.Hs;
    Tp = inp.Tp;
    %using Carter eq(11)
    Func = @(f) 0.31*Hmo.^2.*Tp.*(Tp.*f).^-5.*exp(-1.25./(Tp.*f).^4);
    S = zeros(length(Tp),length(f));
    for i=1:length(f)
        S(:,i) = Func(f(i));
    end
end

%%
function S = pierson_moskowitz(f,inp)
    %construct the Pierson-Moskowitz spectrum
    if ~strcmp(inp.source,'Wave')
        warndlg('Pierson-Moskowitz option only available with wave type input')
        S = []; return;
    end
    Hmo = inp.Hs;
    Tp = inp.Tp;

    if isempty(Hmo)
        Func = @(f) 5e-4*(f).^-5.*exp(-1.25./(Tp.*f).^4);     %Carter eq(15)
        S = zeros(length(Tp),length(f));
    else 
        Func = @(f) 5e-4*(f).^-5.*exp(-2e-3./(Hmo.^2)./f.^4); %Carter eq(14)
        S = zeros(length(Hmo),length(f));
    end

    for i=1:length(f)
        S(:,i) = Func(f(i));
    end
end

%%
function [S,gamma] = jonswap(f,inputs)
    %construct the JONSWAP spectrum using Carter eq(16) with alpha based on
    %Hughes for wind and Carter for wave type input
    [fp,alpha,gamma] = get_input(inputs,f,0);
    if isempty(fp), return; end
    sigma = [0.07,0.09];
    S = getSpectrum(f,fp,alpha,gamma,sigma,[],false);
end

%%
function [S,gamma] = tma(f,inp)
    %construct the TMA spectrum using Kitaigorodskii limit (see Bouws et al, 
    %or Hughes for details)
    [fp,alpha,gamma] = get_input(inp,f,1);
    if isempty(fp), return; end 
    sigma = [0.07,0.09];
    S = getSpectrum(f,fp,alpha,gamma,sigma,inp.ds,true);
end

%%
function S = getSpectrum(f,fp,alpha,gamma,sigma,ds,istma)
    %construct the Jonswap or TMA spectrum
    g = 9.81;
    sig = @(f) sigma(1).*(f<=fp) + sigma(2).*(f>fp);   %Hughes eq(5 & 25), or
    q = @(f) exp((-(f-fp).^2)./(2.*sig(f).^2.*fp.^2)); %Carter eq(16)
    const = g^2*(2*pi)^-4;
    Jonswap = @(f) const*alpha.*(f.^-5).*(exp(-1.25*(f./fp).^-4)).*(gamma.^q(f));

    S = zeros(length(fp),length(f));
    if istma                                          %TMA spectrum
        if isempty(ds) || ds<=0
            ds = (g./fp./2./pi)./fp/2;                %use deep water 
        end 
        Func = @(f) kit_limit(f,ds).*Jonswap(f);      %Hughes eq(15)
        for i=1:length(f)
            S(:,i) = Func(f(i));                          
        end

    else                                              %Jonswap spectrum
        for i=1:length(f)
            S(:,i) = Jonswap(f(i));
        end
    end
end

%%
function [fp,alpha,gamma] = get_input(inp,f,istma)
    %unpack inp for the wind and wave cases
    g = 9.81;
    switch inp.source
        case 'Wind'
            % Adjust wind speed to 10m using power law profile
            inp.U = inp.Uw*(10/inp.zW)^(1/7);
            gFU2 = g*inp.Fetch/inp.U.^2;

            [Tp,Lp] = peak_period(inp);    fp = 1./Tp;
            %fp = 3.5*(g/U)*(gFU2)^-0.33; Tp = 1/fp; %Hughes eq(8)

            %define alpha and gamma for TMA or Jonswap
            if istma
                kp = 2*pi*inp.U.^2./g./Lp;         %Hughes eq(24)
                alpha = 0.0078*kp.^0.49;           %Hughes eq(22)
                gamma = 2.47*kp.^0.39;             %Hughes eq(23)
            else
                alpha = 0.076*(gFU2).^-0.22;       %Hughes eq(6)
                gamma = 7.0*(gFU2).^-0.143;        %Hughes eq(8)
            end

        case 'Wave'
            Hmo = inp.Hs;
            [Tp,Lp] = peak_period(inp);    fp = 1./Tp;
            sigma = [0.07,0.09];
            bnds = [1,7];
            gamma0 = 2;
            lambda = 0.1;

            if istma
                alpha = (pi*Hmo/Lp).^2;            %Hughes eq(29)
                if inp.gamma==0
                    gamma = 6614*(Hmo/4/Lp).^1.59; %Hughes eq(28)
                else
                    gamma = inp.gamma;
                end
                if gamma<=bnds(1)|| gamma>=bnds(2),  gamma = 3.3; end  %use default

            else
                if inp.gamma==0 &&  isfield(inp,'T1')
                    %estimate the best fit for gamma. Does not always fall
                    %within realistic region and default of 3.3 is used                   
                    %[gamma,gamma_sensitivity] = gamma_estimation(inp,f,fp,sigma,bnds,istma);
                    %[gamma, gamma_sensitivity] = jonswap_gamma_estimation(inp,bnds);
                    %gamma = fit_gamma_multimoment(inp,bnds,gamma0,lambda);
                    % dg = max((gamma-gamma_sensitivity));
                    % if dg>0.5
                    %     fprintf('gamma sensivity is %.1f %% \n',dg*100)
                    % end
                    % 
                    % I0 = spectral_moment(0,gamma,sigma);       
                    % I1 = spectral_moment(1,gamma,sigma);
                    I0I1 = inp.T1/Tp;                       %Carter eq(23)
                    if I0I1>0.7 && I0I1<0.9
                        gamma = 45.3*(I0I1)^14.6;           %fit to MIAS No.4, Table 1
                    else
                        gamma = 1.0;
                    end

                elseif inp.gamma>0
                    gamma = inp.gamma;
                else
                    gamma = 3.3;
                end
                %if gamma<=bnds(1)|| gamma>=bnds(2),  gamma = 1.0; end  %use default

                I0 = spectral_moment(0,gamma,sigma);            %Carter eq(18)
                alpha = (2*pi).^4*Hmo.^2.*fp.^4./(16*g.^2.*I0); %Carter eq(20)                
            end

        otherwise
            warndlg('Invalid source - should be wind or wave')
            fp = []; alpha = []; gamma = [];
    end
end

%%
function [Tp,Lp] = peak_period(inp)
    %find peak period based on wind or wave input
    g = 9.81;
    if strcmp(inp.source,'Wind')
        Tp = 0.54*g^-0.77*inp.U.^0.54.*inp.Fetch.^0.23; %Donelan, 1985
    else
        Tp = inp.Tp;
    end

    %get the wave length at peak frequency
    Lp = (g*Tp./2./pi).*Tp;                   %use deep water celerity
    if isfield(inp,'df') && ~isempty(inp.df)
        if isscalar(inp.df) && inp.df>0       %scalar depth over fetch
            Lp = celerity(Tp,inp.df).*Tp;     %use Hunt eq.for celerity
        elseif length(inp.df)==length(Tp)     %vector depths over fetch
            for i=1:length(inp.df)
                if inp.df(i)>0
                    Lp(i) = celerity(Tp(i),inp.df(i)).*Tp(i);  %use Hunt eq.for celerity
                end
            end
        elseif length(inp.df)~=length(Tp)
            errordlg('Input vectors to wave_spectrum must be the same length')
        else
            %defaults to deepwater for all other combinations
        end
    end 
end

%%
function In = spectral_moment(nm,gamma,sigma)
    %find the nm spectral moment of the Jonswap spectrum for given gamma
    sig = @(f) sigma(1).*(f<=1) + sigma(2).*(f>1);   
    q = @(f) exp((-(f-1).^2)./(2.*sig(f).^2));
    Func = @(f) (f.^nm).*(f.^-5).*exp(-1.25.*f.^-4).*gamma.^q(f);
    In = integral(Func,0,Inf,'RelTol',1e-3,'AbsTol',1e-3);
end

%%
function phi = kit_limit(f,ds)
    % Calculate the Kitaigorodskii limit to the spectrum at frquency f
    % f - wave frequency (1/s)
    % ds - water depth at site (m)
    % phi - frquency dependent Kitaigorodskii adjustment to be applied to the 
    % JONSWAP spectrum to take accound of depth limiting effects
    % using Thompson and Vincent, 1983 approximation as given in
    % Hughes, 1984, Eq.(13) and (15)   
    % <duplicated in ctWaveSpectrum>
    g = 9.81;
    omega = 2*pi*f.*sqrt(ds/g);
    % if omega<=1 (vectorised) or omega>2
    phi = 0.5*omega.^2.*(omega<=1) + 1.*(omega>2);        
    phi = phi + (1 - 0.5*(2-omega).^2).*(omega>1 & omega<=2);
end

%%
function [gamma,gamma_sensitivity] = gamma_estimation(inp,f,fp,sigma,bnds,istma)
    %estimate gamma using Hs, Tp and the sea state component value of T1
    % Inputs
    Hs = inp.Hs; T1 = inp.T1; ds = inp.ds; 
    lb = bnds(1)-1e-3; ub = bnds(2)+1e-3; 
    gamma_sensitivity = nan(1,3); % store gamma for alt sigma sets

    % fun_I0 = @(gamma,sigma) spectral_moment(0,gamma,sigma);                         %Carter eq(18)
    % fun_alp = @(g,s) (2*pi).^4*Hs.^2.*fp.^4./(16*g.^2.*fun_I0(g,s)); %Carter eq(20)  
    obj1 = @(gamma) abs(computeT1(fp,gamma,sigma) - T1);
    obj2 = @(gamma) abs(computeHs(f,fp,gamma,sigma,Hs,ds,istma) - Hs);
    obj = @(gamma) abs(obj2(gamma)/obj1(gamma)^2 -Hs/T1^2);      %steepness
    % search gamma in [lb,ub]    
    gamma = fminbnd(obj, lb, ub, optimset('TolX',1e-6,'Display','off'));


    % sensitivity: try three sigma variants
    % gamma_sensitivity(1,1) = fminbnd(@(g) abs(computeT1(fp,g,[0.06,0.08],ds,istma)-T1),lb,ub);
    % gamma_sensitivity(1,2) = fminbnd(@(g) abs(computeT1(Hs,f,fp,g,[0.07,0.09],ds,istma)-T1),lb,ub);
    % gamma_sensitivity(1,3) = fminbnd(@(g) abs(computeT1(Hs,f,fp,g,[0.08,0.10],ds,istma)-T1),lb,ub);
end
%%
function T1 = computeT1(fp,gamma,sigma)
    %estimate T1 for given spectrum parameters 
    I0 = spectral_moment(0,gamma,sigma);                %Carter eq(18)
    I1 = spectral_moment(1,gamma,sigma);
    T1 = I0/I1/fp;
end
%%
function Hs = computeHs(f,fp,gamma,sigma,Hs,ds,istma)
    %estimate T1 for given spectrum parameters 
    g= 9.81;
    % const = g^2*(2*pi)^-4;
    I0 = spectral_moment(0,gamma,sigma);           %Carter eq(18)
    alpha = (2*pi).^4*Hs.^2.*fp.^4./(16*g.^2.*I0); %Carter eq(20) 
    % Hs = sqrt(16*alpha*const*(fp^-4)*I0);
    S = getSpectrum(f,fp,alpha,gamma,sigma,ds,istma);
    Hs = sqrt(16*trapz(f,S));
end

%%
function gamma_est = fit_gamma_multimoment(inp,bnds,g0,lambda)%(Hs, Tp, T1_obs, T2_obs, T10_obs)

    % Inputs
    Hs = inp.Hs; Tp = inp.Tp; T1_obs = inp.T1; T2_obs = inp.T2; T10_obs = inp.T10;
    flim = [0.025,0.58]; %frequency limits
    f = [flim(1):0.005:0.1,0.11:0.01:flim(2)];  %observed frequency intervals (spt format)
    % Weights and prior
    w1  = 1; w2 = 1; w10 = 1;
    usePrior = true;
    % gamma0   = 3.0;
    % lambda   = 0.2;   % tune this

    % Objective in gamma
    obj = @(g) objective_gamma(g, Hs, Tp, f, ...
                               T1_obs, T2_obs, T10_obs, ...
                               w1, w2, w10, usePrior, g0, lambda);

    % Bounded search
    gamma_est = fminbnd(obj, bnds(1),bnds(2), optimset('TolX',1e-6,'Display','off'));
end
%%
function J = objective_gamma(gamma, Hs, Tp, f, ...
                             T1_obs, T2_obs, T10_obs, ...
                             w1, w2, w10, usePrior, gamma0, lambda)
    % Standard sigmas
    sigma1 = 0.07; sigma2 = 0.09;
    fp = 1/Tp;
    sigma = sigma1*(f<=fp) + sigma2*(f>fp);

    % JONSWAP shape (alpha chosen to match Hs via m0)
    % Start with unit alpha, then rescale
    Sshape = jonswap_shape_unit(f, fp, gamma, sigma);
    m0_unit = trapz(f, Sshape);
    alpha = (Hs/4)^2 / m0_unit;   % enforce m0 -> Hs

    S = alpha * Sshape;

    % Moments
    m0  = trapz(f, S);
    m1  = trapz(f, f    .* S);
    m2  = trapz(f, f.^2 .* S);
    m_1 = trapz(f, f.^(-1) .* S);

    T1_mod  = m0 / m1;
    T2_mod  = sqrt(m0 / m2);
    T10_mod = m_1 / m0;

    e1  = (T1_mod  - T1_obs ) / T1_obs;
    e2  = (T2_mod  - T2_obs ) / T2_obs;
    e10 = (T10_mod - T10_obs) / T10_obs;

    J = w1*e1.^2 + w2*e2.^2 + w10*e10.^2;

    if usePrior
        J = J + lambda * ((gamma - gamma0)/gamma0).^2;
    end
end
%%
function Sshape = jonswap_shape_unit(f, fp, gamma, sigma)
    % Unit-alpha JONSWAP shape (m^2/Hz per unit alpha)
    r  = (f./fp);
    Sshape = r.^(-5) .* exp(-1.25 * r.^(-4)) .* ...
             gamma.^(exp(- ( (r - 1).^2 ) ./ (2*sigma.^2)));
end



% function [gamma, gamma_sensitivity] = jonswap_gamma_estimation(inp,bnds)
%     %estimate gamma using Hs, Tp and the sea state component value of T1
%     % Defaults
%     sigma1 = 0.07; sigma2 = 0.09;
%     % Inputs
%     Hs = inp.Hs; Tp = inp.Tp; T1 = inp.T1;
%     flim = [0.025,0.58]; %frequency limits
%     f = [flim(1):0.005:0.1,0.11:0.01:flim(2)];  %obbserved frequency intervals (spt format)
% 
%     %f = linspace(1/(5*Tp), 2/Tp, 2000);   % frequency grid (adjust as needed)
%     ncomp = numel(Hs);
%     gamma = nan(size(Hs));
%     gamma_sensitivity = nan(size(Hs,2),3); % store gamma for alt sigma sets
% 
%     for i=1:ncomp
%         Hs = Hs(i);
%         T1_obs = T1(i);
% 
%         % objective: for given gamma compute model T1 and return error
%         obj = @(gamma) abs( compute_T1_from_jonswap(Hs, Tp, gamma, f, sigma1, sigma2) - T1_obs );
% 
%         % search gamma in [lb,ub]
%         lb = bnds(1)-1e-3; ub = bnds(2)+1e-3; 
%         gamma_est = fminbnd(obj,lb,ub,optimset('TolX',1e-6,'Display','off'));
%         gamma(i) = gamma_est;
% 
%         % sensitivity: try three sigma variants
%         gamma_sensitivity(i,1) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.06,0.08)-T1_obs),lb,ub);
%         gamma_sensitivity(i,2) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.07,0.09)-T1_obs),lb,ub);
%         gamma_sensitivity(i,3) = fminbnd(@(g) abs(compute_T1_from_jonswap(Hs,Tp,g,f,0.08,0.10)-T1_obs),lb,ub);
%     end
% end
% 
% %%
% function T1 = compute_T1_from_jonswap(Hs, Tp, gamma, f, sigma1, sigma2)
%     fp = 1/Tp;
%     sigma = sigma1*(f<=fp) + sigma2*(f>fp);
%     % JONSWAP shape (Goda / common parameterisation)
%     A = 0.0624 ./ (0.230 + 0.0336*gamma - 0.185./(1.9+gamma)); % empirical alpha prefactor
%     S = A .* Hs.^2 .* Tp.^(-4) .* f.^(-5) .* exp(-1.25*(Tp.*f).^(-4)) .* ...
%         gamma.^(exp(-((Tp.*f - 1).^2)./(2*sigma.^2)));
%     df = f(2)-f(1);
%     m0 = trapz(f, S);
%     m1 = trapz(f, f .* S);
%     T1 = m0 ./ m1;
% end