function [SGo,SGi,Dims] = get_inshore_spectrum(transtable,input,sp)
%
%-------function help------------------------------------------------------
% NAME
%   get_inshore_spectrum.m
% PURPOSE
%   construct the offshore and inshore spectra for given wave conditions or
%   wave buoy spectral data
% USAGE
%   [SGo,SGi,Dims] = get_inshore_spectrum(transtable,input,sp)
% INPUTS
%   transtable - spectral transfer table created by SpectralTransfer class
        %   ShoreAngle - angle of shoreline (degTN), use NaN to exclude this limit
        %   hlimit - limiting depth (m)
%   input - input parameters - timeseries of wave data or spectra
%   sp - defines type of model to use and model parameters
% OUTPUT
%   SGo - array of offshore direction-frequency spectral energy
%   SGi - array of inshore direction-frequency spectral energy
%   Dims - dimensions used for the spectral arrays
% NOTES
%   
% Author: Ian Townend
% CoastalSEA (c) March 2023
%--------------------------------------------------------------------------
%
    %
    SGo = []; SGi = [];
    dir_int = 0.5;     %interval used to interpolate directions (deg)
    radint = deg2rad(dir_int);
    per_int = 0.5;     %interval used to interpolate periods (s)
    
    %extract transfer tables and get shoaling coefficients
    indst = transtable.Inshore;             %inshore properties       
    depths = indst.UserData.Depths;         %inshore water depths            
    indst = indst.DataTable;
    offdst = transtable.Offshore;           %offshore properties  
%     T = offdst.Dimensions.Period;           %wave periods used in ray model
    zwl = offdst.Dimensions.WaterLevel;     %water levels used in ray model
    swl = input.swl;                        %required water level
    if isscalar(zwl)
        Dims.depi = depths;                 %depth at inshore point
    else
        Dims.depi = interp1(zwl,depths,swl);%inshore water depth
    end
    InDir = offdst.RowNames;  %inshore directions are held in the offshore table    

    %check limits for valid rays based on depth and shoreline angle  
%     [idx,bound] = checkLimits(offdst,ShoreAngle,hlimit);
    idx = offdst.depth<=0;
    %frequencies at 1/per_int period intervals
    % freq = 1./(min(T):per_int:max(T));         
    %offshore direction at dir_int degree intervals
    % beta = mod(bound(1):dir_int:bound(2),360);  %range determined by shoreline angle

    %variable definitions:
    % InDir - transfer table inshore directions (degTN)
    % theta - transfer table offshore directions (degTN)
    % beta - interpolated directions 0-360 (degTN)

    % fray - ray model and transfer table frequencies (1/T),  (Hz)
    % fspectra - measured spectra frequencies (Hz)    
    % freq - interpolated frequencies 1/1-30s (Hz)   

    % direction ranges for construction of offshore directional spreading
    % diro - offshore model direction range (degTN)
    % dspectra - measured spectra offshore direction range (degTN)

    %frequencies at 1/per_int period intervals
    freq = 1./(1:per_int:30);         
    %offshore direction at dir_int degree intervals
    beta = 0:dir_int:359.5; 

    [shoal,offdir] = getShoal(offdst,indst,InDir,beta,freq,swl,idx);

    %setup input definition using model or measured definitions
    if sp.ismodel
        sp_inputs = setInputParams(input,sp);
        if isempty(sp_inputs), return; end

        %directional spreading factor for selected function  
        dir0 = input.Dir;                       %mean wave direction 
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
        if istable(input)
            S = input.S;
            SpecDir = input.Dir;
            input = {fspectra,input.Dir,input.Spr,input.Skew,input.Kurt};
        else
            S = input.dst.Spectra.S;
            SpecDir = input.dst.Spectra.Dir;
            input = {input.dst};
        end
        %use range of directional means to define direction range which can
        %be greater than 180 if bidrectional sea-state
        Dmnmx = minmax(SpecDir,'omitnan');
        dspectra = Dmnmx(1)-90:dir_int:Dmnmx(2)+90;
        G = datawell_directional_spectra(dspectra,false,input{:});
        [XoFo,XiFi] = transferDims(dspectra,fspectra,swl,beta,freq,swl); %swl forces scalar use
        Gint = interpn(XoFo{:},G,XiFi{:},'linear',0);  
        Sint = interp1(fspectra,S,freq,'linear',0);  
        SGo = Gint.*Sint;
    end

    SGsh = shoal.*SGo;
    SGsh(isnan(SGsh)) = 0;
    SGi = zeros(size(SGo));
    parfor j=1:length(freq)
        odir = offdir(:,j);
        sgi = interp1(beta,SGsh(:,j),odir);
        sgi(isnan(sgi)) = 0;
        SGi(:,j) = SGi(:,j)+sgi;
    end 

    if strcmp(sp.form,'TMA shallow water') || sp.issat
        phi = kit_limit(freq,Dims.depi); %apply saturation to depth at site
        SGi = phi.*SGi;
    end

    
    Dims.freq = freq; Dims.dir = beta;  %return dimensions used
end
%%
function [shoal,offdir] = getShoal(offdst,indst,InDir,beta,freq,swl,idx)
    %compute the shoaling coeffient and offshore direction arrays
    T = offdst.Dimensions.Period;            %wave periods used in ray model
    fray = 1./T;                                %wave frequencies used in ray model
    zwl = offdst.Dimensions.WaterLevel;      %water levels used in ray model 
    offdst = offdst.DataTable;
    

    %extract direction & celerity data and apply limits
    theta = offdst.theta;  theta(idx) = NaN; %offshore direction    
    c0 = offdst.celerity;  c0(idx) = 0;      %offshore celerity
    cg0 = offdst.cgroup;   cg0(idx) = 0;     %offshore group celerity
    ci = squeeze(indst.celerity);            %inshore celerity
    cgi = squeeze(indst.cgroup);             %inshore group celerity

    %pad the high frequencies with values from the minimum frequency
    if min(T)>3
        addfray = [1,0.5,0.33];
        fray = [addfray';fray];  %pad wave ray frequencies for periods 1-3s
        theta = [repmat(theta(:,1,:),1,3),theta];
        c0 = [repmat(c0(:,1,:),1,3),c0];
        cg0 = [repmat(cg0(:,1,:),1,3),cg0];
        
        addci = repmat(ci(1,:),3,1);         %inshore celerities for highest frequency
        addcgi = repmat(cgi(1,:),3,1);
        cdeepwater = 9.81./addfray/2/pi;     %deepwater celerity = gT/2pi
        cdeep = repmat(cdeepwater',1,size(ci,2));
        addci = min(addci,cdeep);
        addcgi = min(addcgi,cdeep/2);        %deepwater group celerity = c/2

        ci = [addci;ci];                                                    
        cgi = [addcgi;cgi];
    end

    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    [PFW,XoFoWo,fw,fowo] = transferDims(InDir,fray,zwl,beta,freq,swl);

    %interpolate offshore directions and celerity arrays
    offdir = interpn(PFW{:},theta,XoFoWo{:});           %offshore directions array
    c0fdw = interpn(PFW{:},c0,XoFoWo{:},'linear',0);    %offshore celerity frequency,direction,water level array            
    cg0fdw = interpn(PFW{:},cg0,XoFoWo{:},'linear',0);  %offshore group celerity frequency,direction,water level array
    cifw = interpn(fw{:},ci,fowo{:});                   %inshore celerity frequency,water level array
    cgifw = interpn(fw{:},cgi,fowo{:});                 %inshore group celerity frequency,water level array

    %calculate the shoaling coefficient ks(xso,fri)
    ishoal = cifw.*cgifw;                   %inshore c.cg
    oshoal = c0fdw.*cg0fdw;                 %offshore c.cg
    shoal = oshoal./ishoal';                %shoaling factor array            
    % check_plot(1./fro,xso,shoal,{'shoal','Wave period (s)','Offshore direction (degTN)'}); 
end
%%
function [Dims,sp_inputs] = saturation_depths(offdst,G,Dims,sp_inputs,...
                                            InDir,beta,freq,swl,idx,radint)   
    %determine the depths to use if the TMA spectrum or depth saturation is
    %being used
    T = offdst.Dimensions.Period;          %wave periods used in ray model
    fray = 1./T;                           %wave frequencies used in ray model
    zwl = offdst.Dimensions.WaterLevel;    %water levels used in ray model
    offdst = offdst.DataTable;
    hof = offdst.depth;    hof(idx) = 0;   %offshore depth of ray
    hav = offdst.avdepth;  hav(idx) = 0;   %average depth along ray
    hmn = offdst.mindepth; hmn(idx) = 0;   %minimum depth along ray      

    if min(T)>3
        addfray = [1,0.5,0.33];
        fray = [addfray';fray];  %pad wave ray frequencies for periods 1-3s
        hof = [repmat(hof(:,1,:),1,3),hof];
        hav = [repmat(hav(:,1,:),1,3),hav];
        hmn = [repmat(hmn(:,1,:),1,3),hmn];
    end

    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    [PFW,XoFoWo] = transferDims(InDir,fray,zwl,beta,freq,swl);

    hGof = interpn(PFW{:},hof,XoFoWo{:},'linear',0);
    hGav = interpn(PFW{:},hav,XoFoWo{:},'linear',0);
    hGmn = interpn(PFW{:},hmn,XoFoWo{:},'linear',0);  

    %use direction spread to find average depth on rays
    %based on the proportion they contribute to the spectrum
    hGof = trapz(radint,abs(trapz(freq,G'.*hGof,2)));%offshore depth
    hGav = trapz(radint,abs(trapz(freq,G'.*hGav,2)));%average depth
    %for wind inputs offshore spectrum uses average depth 
    %over dominant rays to compute Lp and offshore depth 
    %for saturation limit of offshore spectrum
    sp_inputs = [sp_inputs,{hGav,hGof}]; 
    %for saturation limit use site or dominant rays min depth
    hGmn = trapz(radint,abs(trapz(freq,G'.*hGmn,2)));%minimum depth
    Dims.depi = min([Dims.depi,min(hGmn,[],'All','omitnan')]);
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
end
%%
function [PFW,XoFoWo,fw,fowo] = transferDims(InDir,fray,zwl,beta,freq,swl)
    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    if isscalar(zwl)
        [P,F] = ndgrid(InDir,fray);
        PFW = {P,F};  
        fw = {fray};
        %frequency, direction and water level arrays to interpolate to 
        %var(xso,fro)
        [Xso,Fro] = ndgrid(beta,freq);  
        XoFoWo = {Xso,Fro};  
        fowo = {freq};
    else
        [P,F,W]  = ndgrid(InDir,fray,zwl);
        PFW = {P,F,W}; 
        fw = {fray,zwl};
        %frequency, direction and water level arrays to interpolate to 
        %var(xso,fro) for selected water level swl
        [Xso,Fro,Wln] = ndgrid(beta,freq,swl);  
        XoFoWo = {Xso,Fro,Wln}; 
        fowo = {freq,swl};
    end
end
%%
% function [idx,bound] = checkLimits(offdst,ShoreAngle,hlimit)
%     %find indices of nodes with depths that are too shallow or
%     %offshore directions that are travelling shoreward
%     offdst = offdst.DataTable;
%     id1 = offdst.depth<=hlimit;                    %depth limit for offshore end of rays
%     if ~isnan(ShoreAngle)  %to exclude this limit use NaN when running ray model
%         shoreang = mod(compass2trig(ShoreAngle),2*pi);
%         bound = [shoreang,mod(shoreang+pi,2*pi)];  %shoreline limits
%         dir = compass2trig(offdst.theta);          %convert to trigonometric radians
%         id2 = isangletol(dir,bound);               %index of angles within bound
%     else
%         id2 = false(size(offdst.theta));           %no constraint applied 
%     end
%     idx = logical(id1+id2);                        %combined limits
%     bound = mod(compass2trig(bound,true),360);     %convert back to compass degrees
% end
%%
function spectrum_inputs = setInputParams(input,sp)
    %setup input definition using model definitions
    spectrum_inputs = [];
    %check for invalid conditions when timeseries wave or wind 
    %data used to define conditions
    if istable(input) 
        if strcmp(sp.source,'Wave') && any(isnan(input{1,[1,3,5,9]})) %offshore waves  
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
function check_plot(period,dir,var,varname,ax) %#ok<DEFNU> 
    %check plot for data selection
    if nargin<6
        hf = figure('Name','SpecTrans','Tag','PlotFig');
        ax = axes(hf);
    end
    surf(ax,period,dir,var,'Tag','PlotFigSurface');
    view(2);
    shading interp
    axis tight
    %add the colorbar and labels
    cb = colorbar;
    cb.Label.String = varname{1};
    xlabel(varname{2}); 
    ylabel(varname{3}); 
    if length(varname)>3
        cb.Tag = varname{4};
    else
        cb.Tag = varname{1};
    end
end