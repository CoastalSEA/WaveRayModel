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
%       T2 - mean wave period or zero upcrossing period: used to calculate
%            gamma, or gamma can be specified
%       gamma - spectrum shape parameter (defaults to 3.3 if not specified)
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
    g = 9.81;
    [fp,alpha,gamma] = get_input(inputs,0);
    if isempty(fp), return; end
    
    sigma = @(f) 0.07.*(f<=fp) + 0.09.*(f>fp);           %Hughes eq(5 & 25), or
    q = @(f) exp((-(f-fp).^2)./(2.*sigma(f).^2.*fp.^2)); %Carter eq(16)
    const = g^2*(2*pi)^-4;
    Jonswap = @(f) const*alpha.*(f.^-5).*(exp(-1.25*(f./fp).^-4)).*(gamma.^q(f));

    S = zeros(length(fp),length(f));
    for i=1:length(f)
        S(:,i) = Jonswap(f(i));
    end
end

%%
function [S,gamma] = tma(f,inp)
    %construct the TMA spectrum using Kitaigorodskii limit (see Bouws et al, 
    %or Hughes for details)
    g = 9.81;
    [fp,alpha,gamma] = get_input(inp,1);
    if isempty(fp), return; end 

    sigma = @(f) 0.07.*(f<=fp) + 0.09.*(f>fp);           %Hughes eq(5 & 25), or
    q = @(f) exp((-(f-fp).^2)./(2.*sigma(f).^2.*fp.^2)); %Carter eq(16)
    cn = g^2*(2*pi)^-4;
    Jonswap = @(f) cn*alpha.*(f.^-5).*(exp(-1.25*(f./fp).^-4)).*(gamma.^q(f));

    if isfield(inp,'ds') && ~isempty(inp.ds) && inp.ds>0
        ds = inp.ds;
    else
        ds = (g./fp./2./pi)./fp/2;                %use deep water 
    end 

    Func = @(f) kit_limit(f,ds).*Jonswap(f);      %Hughes eq(15)
    S = zeros(length(fp),length(f));
    for i=1:length(f)
        S(:,i) = Func(f(i));
        % Sf(:,i) = Jonswap(f(i));
    end
    % figure('Tag','PlotFig'); plot(f,S,f,Sf);
    % fprintf('TMA m0=%f, Jonswap m0=%f\n',trapz(f,S),trapz(f,Sf))
end

%%
function [fp,alpha,gamma] = get_input(inp,istma)
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

            if istma
                alpha = (pi*Hmo/Lp).^2;            %Hughes eq(29)
                if inp.gamma==0
                    gamma = 6614*(Hmo/4/Lp).^1.59; %Hughes eq(28)
                else
                    gamma = inp.gamma;
                end
            else
                if inp.gamma==0 &&  isfield(inp,'T2')
                    gamma = 70*(inp.T2/inp.Tp)^12.23;
                elseif inp.gamma>0
                    gamma = inp.gamma;
                else
                    gamma = 3.3;
                end

                Io = spectral_moment(0,gamma);                  %Carter eq(18)
                alpha = (2*pi).^4*Hmo.^2.*fp.^4./(16*g.^2.*Io); %Carter eq(20)
            end

            if gamma<1
                gamma = 0.9999;
            elseif gamma>7.9
                gamma = 7.9999;
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
function In = spectral_moment(nm,gamma)
    %find the nm spectral moment of the Jonswap spectrum for given gamma
    sigma = @(f) 0.07.*(f<=1) + 0.09.*(f>1);   
    q = @(f) exp((-(f-1).^2)./(2.*sigma(f).^2));
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
