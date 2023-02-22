function S = wave_spectrum(stype,f,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   wave_spectrum.m
% PURPOSE
%   Calculate the spectral energy at a number of frequencies using a
%   selection of spectrum definitions
% USAGE
%   S = wave_spectrum(stype,f,varargin)
% INPUTS
%   stype - type of spectrum. options are 'Bretschneider open ocean', 
%           'Pierson-Moskowitz fully developed', 'JONSWAP fetch limited', 
%           or 'TMA shallow water'
%   f - frequencies at which the spectrum is to be defined <scalar or vector>
%   varargin - depends on use. When deriving the spectrum using wind speed
%              and fetch the order is 'wind',Uw,zw,Fch as defined bleow. 
%              When using wave records the order is 'wave',Hmo,Tp,gamma
%   For wind input (JONSWAP and TMA only)
%       Uw - wind speed (m/s) <scalar or vector>
%       zw - elevation of wind speed measurement (m)
%       Fch- dominant fetch length (m) <scalar or vector>
%   For wave input
%       Hmo - significant wave height [4sqrt(mo)] (m) <scalar or vector>
%       Tp peak wave period (s)  <can be vector - same length as Hs>
%       gamma - spectrum shape parameter (defaults to 3.3 if not specified)
%   When using 'TMA shallow water'
%   df - average water depth over fetch (m) (default is deep water)
%   ds - depths at site (m) <scalar or vector>
%           NB: wind or wave vector data must be the same length
% OUTPUT
%   S - spectral energy density at specified frequencies [nrec,nf] (m^2s)
% NOTES
%   The spectrum definitions in Carter, 1982, MIAS publication No.4 are
%   used with additional information from Hughes, CERC-84-7 (p10-12); 
%   Bouws et al, 1985; Donelan, Hamilton and Hui, R.Soc, 1985, eqn 6.1) and
%   Hunt,ASCE,WW4,1974,p457-459
%   NB: 1) gamma only used for JONSWAP and TMA spectra
%       2) Pierson-Moskowitz can use just Hmo or just Tp. If specifying Tp
%          enter Hmo as [] and vice versa.
% SEE ALSO
%   uses celerity.m. cf tma_spectrum.m which outputs [Hs,Tp,Tz]
%
% Author: Ian Townend
% CoastalSEA (c)Feb 2023
%--------------------------------------------------------------------------
%
    switch stype
        case 'Bretschneider open ocean'
            S = bretschneider(f,varargin{:});
        case 'Pierson-Moskowitz fully developed'
            S = pierson_moskowitz(f,varargin{:});
        case 'JONSWAP fetch limited'
            S = jonswap(f,varargin{:});
        case 'TMA shallow water'
             S = tma(f,varargin{:});
        otherwise
            warndlg('Unknown spectrum type')
            S = [];
    end
end
%%
function S = bretschneider(f,varargin)
    %construct the Bretschneider spectrum
    if ~strcmp(varargin{1},'wave')
        warndlg('Bretschneider option only available with wave type input')
        return;
    end
    Hmo = varargin{2};
    Tp = varargin{3};
    %using Carter eq(11)
    Func = @(f) 0.31*Hmo.^2.*Tp.*(Tp.*f).^-5.*exp(-1.25./(Tp.*f).^4);
    S = zeros(length(Tp),length(f));
    for i=1:length(f)
        S(:,i) = Func(f(i));
    end
end
%%
function S = pierson_moskowitz(f,varargin)
    %construct the Pierson-Moskowitz spectrum
    if ~strcmp(varargin{1},'wave')
        warndlg('Pierson-Moskowitz option only available with wave type input')
        return;
    end
    Hmo = varargin{2};
    Tp = varargin{3};

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
function S = jonswap(f,varargin)
    %construct the JONSWAP spectrum using Carter eq(16) with alpha based on
    %Hughes for wind and Carter for wave type input
    g = 9.81;
    [fp,alpha,gamma] = get_input(varargin{:});
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
function S = tma(f,varargin)
    %construct the TMA spectrum using Kitaigorodskii limit (see Bouws et al, 
    %or Hughes for details)
    g = 9.81;
    [fp,alpha,gamma] = get_input(varargin{:});
    if isempty(fp), return; end 
    sigma = @(f) 0.07.*(f<=fp) + 0.09.*(f>fp);           %Hughes eq(5 & 25), or
    q = @(f) exp((-(f-fp).^2)./(2.*sigma(f).^2.*fp.^2)); %Carter eq(16)
    cn = g^2*(2*pi)^-4;
    Jonswap = @(f) cn*alpha.*(f.^-5).*(exp(-1.25*(f./fp).^-4)).*(gamma.^q(f));
    if length(varargin)>5
        ds = varargin{5};
    else
        ds = (g./fp./2./pi)./fp/2;                %use deep water 
    end 

    Func = @(f) kit_limit(f,ds).*Jonswap(f);      %Hughes eq(15)
    S = zeros(length(fp),length(f));
    for i=1:length(f)
        S(:,i) = Func(f(i));
    end
end
%%
function [fp,alpha,gamma] = get_input(varargin)
    %unpack varargin for the wind and wave cases
    g = 9.81;
    switch varargin{1}
        case 'Wind'
            Uw = varargin{2};
            zw = varargin{3};
            Fch = varargin{4};
            % Adjust wind speed to 10m using power law profile
            U = Uw*(10/zw)^(1/7);

            Tp = 0.54*g^-0.77*U.^0.54.*Fch.^0.23;
            fp = 1./Tp;
            if length(varargin)>4
                df = varargin{5};
                Lp = celerity(Tp,df).*Tp;         %use Hunt eq.for celerity
            else
                Lp = (g*Tp./2./pi).*Tp;           %use deep water celerity
            end            
            kp = 2*pi*U.^2./g./Lp;                %Hughes eq(24)
            alpha = 0.0078*kp.^0.49;              %Hughes eq(22)
            gamma = 2.47*kp.^0.39;                %Hughes eq(23)
        case 'Wave'
            Hmo = varargin{2};
            Tp = varargin{3};
            gamma = varargin{4};

            fp = 1./Tp;
            Io = spectral_moment(0,gamma);                 %Carter eq(18)
            alpha = (2*pi).^4*Hmo.^2.*fp.^4./(16*g.^2.*Io); %Carter eq(20)
        otherwise
            warndlg('Invalid source - should be wind or wave')
            fp = []; alpha = []; gamma = [];
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
    g = 9.81;
    omega = 2*pi*f.*sqrt(ds/g);
    % if omega<=1 (vectorised) or omega>2
    phi = 0.5*omega.^2.*(omega<=1) + 1.*(omega>2);        
    phi = phi + (1 - 0.5*(2-omega).^2).*(omega>1 & omega<=2);
end
