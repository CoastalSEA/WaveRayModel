function [SGo,SGi,Dims] = get_inshore_spectrum(transtable,ShoreAngle,hlimit,input,sp)
%
%-------function help------------------------------------------------------
% NAME
%   get_inshore_spectrum.m
% PURPOSE
%   construct the offshore and inshore spectra for given wave conditions or
%   wave buoy spectral data
% USAGE
%   [SGo,SGi,Dims] = get_inshore_spectrum(transtable,rayobj,input,sp)
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
    T = offdst.Dimensions.Period;           %wave periods used in ray model
    zwl = offdst.Dimensions.WaterLevel;     %water levels used in ray model
    swl = input.swl;                        %required water level
    if isscalar(zwl)
        Dims.depi = depths;                 %depth at inshore point
    else
        Dims.depi = interp1(zwl,depths,swl);%inshore water depth
    end
    phi = offdst.RowNames;  %inshore directions are held in the offshore table    

    %check limits for valid rays based on depth and shoreline angle  
    [idx,bound] = checkLimits(offdst,ShoreAngle,hlimit);

    %frequencies at 1/per_int period intervals
    fri = 1./(min(T):per_int:max(T));         
    %offshore direction at dir_int degree intervals
    xsi = mod(bound(1):dir_int:bound(2),360);  %range determined by shoreline angle
    
    [shoal,offdir,fri] = getShoal(offdst,indst,phi,xsi,fri,swl,idx);

    %setup input definition using model or measured definitions
    if sp.ismodel
        [Dims,sp_inputs] = setInputParams(input,sp,Dims);
        if isempty(sp_inputs), return; end

        %directional spreading factor for selected function  
        dir0 = input.Dir;                       %mean wave direction 
        xso = dir0-90:dir_int:dir0+90;
        G = directional_spreading(input.Dir,xso,sp.nspread,sp.spread);
        G = interp1(xso,G,xsi,'linear',0);
        % figure; plot(xso,G);

        if strcmp(sp.form,'TMA shallow water')
            fr0 = 1./(min(T):per_int:max(T));   %unpadded frequency range
            [Dims,sp_inputs] = saturation_depths(offdst,G,Dims,...
                                sp_inputs,phi,xso,fr0,swl,idx,radint);
        end

        %spectral energy for selected wave spectrum
        S = wave_spectrum(sp.form,fri,sp_inputs{:});
        % Hs = 4*sqrt(abs(trapz(fri,S))); %check = input Hsn
        % figure; plot(fri,sf);
        SGo = G'*S;
    else
        fso = input.dst.Spectra.Dimensions.freq;
        S = input.dst.Spectra.S;
        %use range of directional means to define direction range which can
        %be greater than 180 if bidrectional sea-state
        Dmnmx = minmax(input.dst.Spectra.Dir,'omitnan');
        xso = Dmnmx(1)-90:dir_int:Dmnmx(2)+90;
        G = datawell_directional_spectra(xso,false,input.dst);
        [XoFo,XiFi] = transferDims(xso,fso,swl,xsi,fri,swl); %swl forces scalar use
        Gint = interpn(XoFo{:},G,XiFi{:},'linear',0);  
        Sint = interp1(fso,S,fri,'linear',0);  
        SGo = Gint.*Sint;
    end

    SGsh = shoal.*SGo;
    SGsh(isnan(SGsh)) = 0;
    SGi = zeros(size(SGo));
    for j=1:length(fri)
        odir = offdir(:,j);
        sgi = interp1(xsi,SGsh(:,j),odir);
        sgi(isnan(sgi)) = 0;
        SGi(:,j) = SGi(:,j)+sgi;
    end 

    if strcmp(sp.form,'TMA shallow water') || sp.issat
        phi = kit_limit(fri,Dims.depi); %apply saturation to depth at site
        SGi = phi.*SGi;
    end
    Dims.f = fri; Dims.xsi = xsi;           %return dimensions used
end
%%
function [shoal,offdir,fro] = getShoal(offdst,indst,phi,xso,fro,swl,idx)
    %compute the shoaling coeffient and offshore direction arrays
    T = offdst.Dimensions.Period;            %wave periods used in ray model
    f = 1./T;                                %wave frequencies used in ray model
    zwl = offdst.Dimensions.WaterLevel;      %water levels used in ray model 
    offdst = offdst.DataTable;
    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    [PFW,XoFoWo,fw,fowo] = transferDims(phi,f,zwl,xso,fro,swl);

    %extract direction & celerity data and apply limits
    theta = offdst.theta;  theta(idx) = NaN; %offshore direction    
    c0 = offdst.celerity;  c0(idx) = 0;      %offshore celerity
    cg0 = offdst.cgroup;   cg0(idx) = 0;     %offshore group celerity
    ci = squeeze(indst.celerity);            %inshore celerity
    cgi = squeeze(indst.cgroup);             %inshore group celerity

    offdir = interpn(PFW{:},theta,XoFoWo{:});           %offshore directions array
    c0fdw = interpn(PFW{:},c0,XoFoWo{:},'linear',0);    %offshore celerity frequency,direction,water level array            
    cg0fdw = interpn(PFW{:},cg0,XoFoWo{:},'linear',0);  %offshore group celerity frequency,direction,water level array
    cifw = interpn(fw{:},ci,fowo{:});                   %inshore celerity frequency,water level array
    cgifw = interpn(fw{:},cgi,fowo{:});                 %inshore group celerity frequency,water level array

    %pad the high frequencies with values from the minimum frequency
    if min(T)>3
        addfro = [1,0.5,0.33];
        fro = [addfro,fro];      %pad frequencies for periods 1-3s.
        offdir = [repmat(offdir(:,1,:),1,3),offdir];
        c0fdw = [repmat(c0fdw(:,1,:),1,3),c0fdw];
        cg0fdw = [repmat(cg0fdw(:,1,:),1,3),cg0fdw];

        cdeepwater = 9.81./addfro/pi;
        addci = min([repmat(cifw(1),1,3);cdeepwater],[],1);
        addcgi = min([repmat(cgifw(1),1,3);cdeepwater/2],[],1);
        cifw = [addci,cifw'];
        cgifw = [addcgi,cgifw'];
    end

    %calculate the shoaling coefficient ks(xso,fri)
    ishoal = cifw.*cgifw;                   %inshore c.cg
    oshoal = c0fdw.*cg0fdw;                 %offshore c.cg
    shoal = oshoal./ishoal;                 %shoaling factor array            
    % check_plot(1./fro,xso,shoal,{'shoal','Wave period (s)','Offshore direction (degTN)'}); 
end
%%
function [Dims,sp_inputs] = saturation_depths(offdst,G,Dims,sp_inputs,...
                                            phi,xso,fro,swl,idx,radint)                                 
    %determine the depths to use if the TMA spectrum or depth saturation is
    %being used
    T = offdst.Dimensions.Period;            %wave periods used in ray model
    f = 1./T;                                %wave frequencies used in ray model
    zwl = offdst.Dimensions.WaterLevel;    %water levels used in ray model
    offdst = offdst.DataTable;
    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    [PFW,XoFoWo] = transferDims(phi,f,zwl,xso,fro,swl);

    hof = offdst.depth;    hof(idx) = 0;   %offshore depth of ray
    hav = offdst.avdepth;  hav(idx) = 0;   %average depth along ray
    hmn = offdst.mindepth; hmn(idx) = 0;   %minimum depth along ray    
    hGof = interpn(PFW{:},hof,XoFoWo{:},'linear',0);
    hGav = interpn(PFW{:},hav,XoFoWo{:},'linear',0);
    hGmn = interpn(PFW{:},hmn,XoFoWo{:},'linear',0);  

    if min(T)>3
        fro = [[1,0.5,0.33],fro];      %pad frequencies for periods 1-3s.
        hGof = [repmat(hGof(:,1,:),1,3),hGof];
        hGav = [repmat(hGav(:,1,:),1,3),hGav];
        hGmn = [repmat(hGmn(:,1,:),1,3),hGmn];
    end
    %use direction spread to find average depth on rays
    %based on the proportion they contribute to the spectrum
    hGof = trapz(radint,abs(trapz(fro,G'.*hGof,2)));%offshore depth
    hGav = trapz(radint,abs(trapz(fro,G'.*hGav,2)));%average depth
    %for wind inputs offshore spectrum uses average depth 
    %over dominant rays to compute Lp and offshore depth 
    %for saturation limit of offshore spectrum
    sp_inputs = [sp_inputs,{hGav,hGof}]; 
    %for saturation limit use site or dominant rays min depth
    hGmn = trapz(radint,abs(trapz(fro,G'.*hGmn,2)));%minimum depth
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
function [PFW,XoFoWo,fw,fowo] = transferDims(phi,f,zwl,xso,fro,swl)
    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    if isscalar(zwl)
        [P,F] = ndgrid(phi,f);
        PFW = {P,F};  
        fw = {f};
        %frequency, direction and water level arrays to interpolate to 
        %var(xso,fro)
        [Xso,Fro] = ndgrid(xso,fro);  
        XoFoWo = {Xso,Fro};  
        fowo = {fro};
    else
        [P,F,W]  = ndgrid(phi,f,zwl);
        PFW = {P,F,W}; 
        fw = {f,zwl};
        %frequency, direction and water level arrays to interpolate to 
        %var(xso,fro) for selected water level swl
        [Xso,Fro,Wln] = ndgrid(xso,fro,swl);  
        XoFoWo = {Xso,Fro,Wln}; 
        fowo = {fro,swl};
    end
end
%%
function [idx,bound] = checkLimits(offdst,ShoreAngle,hlimit)
    %find indices of nodes with depths that are too shallow or
    %offshore directions that are travelling shoreward
    offdst = offdst.DataTable;
    id1 = offdst.depth<=hlimit;                    %depth limit for offshore end of rays
    if ~isnan(ShoreAngle)  %to exclude this limit use NaN
        shoreang = mod(compass2trig(ShoreAngle),2*pi);
        bound = [shoreang,mod(shoreang+pi,2*pi)];  %shoreline limits
        dir = compass2trig(offdst.theta);          %convert to trigonometric radians
        id2 = isangletol(dir,bound);               %index of angles within bound
    else
        id2 = false(size(offdst.theta));           %no constraint applied 
    end
    idx = logical(id1+id2);                        %combined limits
    bound = mod(compass2trig(bound,true),360);     %convert back to compass degrees
end
%%
function [Dims,spectrum_inputs] = setInputParams(input,sp,Dims)
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
    else
        spectrum_inputs = {sp.source,input.AvSpeed,...
                                    input.zW,input.Fetch,sp.gamma};
    end
end
%%
function check_plot(T,phi,var,varname,ax) %#ok<DEFNU> 
    %check plot for data selection
    if nargin<6
        hf = figure('Name','SpecTrans','Tag','PlotFig');
        ax = axes(hf);
    end
    surf(ax,T,phi,var,'Tag','PlotFigSurface');
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