function [SGo,Dims] = get_offshore_spectrum(intable,sp,Dims)
%
%-------function help------------------------------------------------------
% NAME
%   get_offshore_spectrum.m
% PURPOSE
%   construct the offshore using wave buoy spectral data or wave properties 
% USAGE
%   [SGo,SGi,Dims] = get_inshore_spectrum(sptobj,intable,sp)
% INPUTS
                %   sptobj - SpectralTransfer class object
                %            Data - inshore and offshore transfer tables
                %            interp - direction and period interpolation intervals
                        %no longer passed
                        %   transtable - spectral transfer table created by SpectralTransfer class        
                        %   ShoreAngle - angle of shoreline (degTN), use NaN to exclude this limit
                        %   hlimit - limiting depth (m)
%   intable - input parameters - table of timeseries of wave data or spectra
%   sp - defines type of model to use and model parameters
% OUTPUT
%   SGo - array of offshore direction-frequency spectral energy
%   SGc - array of constructed direction-frequency spectral energy
%   Dims - dimensions used for the spectral arrays
% NOTES
%   
% Author: Ian Townend
% CoastalSEA (c) Oct 2025
%--------------------------------------------------------------------------
%
    SGo = []; 
    dir_int = 0.5;  %interval used to interpolate directions (deg)
    per_int = 0.5;  %interval used to interpolate periods (s)
    per_range = 30; %upper bound of period range

    %frequencies at 1/per_int period intervals
    freq = 1./(1:per_int:per_range);         
    %offshore direction at dir_int degree intervals
    beta = 0:dir_int:359.5;

    %variable definitions:
    % direction ranges for construction of offshore directional spreading
    % diro - offshore model direction range (degTN)
    % dspectra - measured spectra offshore direction range (degTN)

    %setup input definition using model or measured definitions
    if sp.ismodel
        sp_inputs = setInputParams(intable,sp);
        if isempty(sp_inputs), return; end

        %directional spreading factor for selected function  
        dir0 = intable.Dir0;                       %mean wave direction 
        diro = dir0-90:dir_int:dir0+90;
        G = directional_spreading(dir0,diro,sp.nspread,sp.spread);
        G = interp1(diro,G,beta,'linear',0);
        % figure; plot(xso,G);

        if strcmp(sp.form,'TMA shallow water')
            [Dims,sp_inputs] = saturation_depths(offdst,G,Dims,...
                                sp_inputs,InDir,beta,freq,swl,idx,radint);
        end

        %spectral energy for selected wave spectrum
        S = wave_spectrum(sp.form,freq,sp_inputs{:});
        % Hs = 4*sqrt(abs(trapz(fri,S))); %check = input Hsn
        % figure; plot(fri,sf);
        SGo = G'*S;
    else
        fspectra = sp.freq;
        S = intable.S;
        SpecDir = intable.Dir;
        intable = {fspectra,intable.Dir,intable.Spr,intable.Skew,intable.Kurt};
        %use range of directional means to define direction range which can
        %be greater than 180 if bidrectional sea-state
        Dmnmx = minmax(SpecDir,'omitnan');
        dspectra = Dmnmx(1)-90:dir_int:Dmnmx(2)+90;
        G = datawell_directional_spectra(dspectra,false,intable{:}); %isplot=false
        [XoFo,XiFi] = transferDims(dspectra,fspectra,beta,freq); 
        Gint = interpn(XoFo{:},G,XiFi{:},'linear',0);  
        Sint = interp1(fspectra,S,freq,'linear',0);  
        SGo = Gint.*Sint;
    end

    Dims.freq = freq; Dims.dir = beta;  %return dimensions used
end



%%
function spectrum_inputs = setInputParams(input,sp)
    %setup input definition using model definitions
    spectrum_inputs = [];
    %check for invalid conditions when timeseries wave or wind 
    %data used to define conditions
    if istable(input) 
        idvar = ismatch(input.Properties.VariableNames,{'Hs','Dir','Tp','swl'});
        if strcmp(sp.source,'Wave') && any(isnan(input{1,idvar})) %offshore waves  
            return;
        elseif strcmp(sp.source,'Wind') && any(isnan(input{1,[1,2]})) %wind input  
            return;
        end
    end   

    %setup model spectrum inputs
    if strcmp(sp.source,'Wave')
        spectrum_inputs = {sp.source,input.Hs,input.Tp,sp.gamma};
    elseif strcmp(sp.source,'Wind')
        spectrum_inputs = {sp.source,input.AvSpeed,...
                                    input.zW,input.Fetch,sp.gamma};
    end
end

%%
function [PFW,XoFoWo,fw,fowo] = transferDims(InDir,fray,beta,freq)
    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
        [P,F] = ndgrid(InDir,fray);
        PFW = {P,F};  
        fw = {fray};
        %frequency, direction and water level arrays to interpolate to 
        %var(xso,fro)
        [Xso,Fro] = ndgrid(beta,freq);  
        XoFoWo = {Xso,Fro};  
        fowo = {freq};
end