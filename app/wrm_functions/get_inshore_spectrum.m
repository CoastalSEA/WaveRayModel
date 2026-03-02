function inobj = get_inshore_spectrum(obj,offobj)
%
%-------function help------------------------------------------------------
% NAME
%   get_inshore_spectrum.m
% PURPOSE
%   construct the inshore spectra for given offshore spectra derived from 
%   wave conditions or wave buoy spectral data
% USAGE
%   inobj = get_inshore_spectrum(offobj,sptobj)
% INPUTS
%   obj - SpectralTransfer class object
%   offobj - instance of ctWaveSpectrum that defines offshore spectrum
% OUTPUT
%   inobj - instance of ctWaveSpectrum of the associated inshore spectrum
% NOTES
%   used in SpectralTransfer and WRM_WaveModel as part of WaveRayModel.
%   Modified from original version to use ctWaveSpectrum to manipulate
%   spectra and only handle inshore spectrum, rather that computing
%   offshore and inshore spectra
%
% Author: Ian Townend
% CoastalSEA (c) Dec 2025
%--------------------------------------------------------------------------
%
    inobj = ctWaveSpectrum;
    [dir,freq] = spectrumDimensions(inobj); %default frequency-direction vectors
    %extract transfer tables and get shoaling coefficients
    indst = obj.Data.Inshore;               %inshore properties       
    depths = indst.UserData.Depths;         %inshore water depths            
    %indst = indst.DataTable;
    offdst = obj.Data.Offshore;             %offshore properties  
    zwl = offdst.Dimensions.WaterLevel;     %water levels used in ray model

    if isfield(offobj.inpData,'swl')
        swl = offobj.inpData.swl;            %required water level
    else
        errdlg('swl not set in input to get_inshore_params')
        inobj = []; return
    end

    if isscalar(zwl)
        depi = depths;                      %depth at inshore point
    else
        depi = interp1(zwl,depths,swl);     %inshore water depth
    end

    %check limits for valid rays based on depth and shoreline angle  
    %parfor not finding dynamic property, so use table
    idx = offdst.DataTable.depth<=0 | isnan(offdst.DataTable.depth);   

    [shoal,offdir] = getShoal(obj,offobj,swl,idx);

    %construct inshore spectrum
    SGo = offobj.Spectrum.SG;

    SGsh = shoal.*SGo;
    SGsh(isnan(SGsh)) = 0;
    SGi = zeros(size(SGo));
    for j=1:length(freq)
        odir = offdir(:,j);
        sgi = interp1(dir,SGsh(:,j),odir);
        sgi(isnan(sgi)) = 0;
        SGi(:,j) = SGi(:,j)+sgi;
    end 

    if isfield(offobj.spModel,'form') && strcmp(offobj.spModel.form,'TMA shallow water')
        %apply saturation limit if TMA spectrum (depth determined from rays)
        if size(SGi,1) ~= numel(dir), SGi = SGi.'; end
        theta = deg2rad(mod(dir(:),360));  % [I x 1]
        wp = trapz_weights_periodic(theta);
        Sf = sum(SGi.*wp); 
        Gi = SGi./Sf;  %checksum = sum(Gi.*wp); %check that non-dimensional     
        depi = inshoreDepths(obj,inobj,depi,Gi,swl);
        %apply saturation to depth at site
        phi = kit_limit(offobj,freq,depi); 
        SGi = phi.*SGi;
    end
    inobj.Spectrum.SG = SGi;        %inshore spectrum
    inobj.Spectrum.freq = freq;
    inobj.Spectrum.dir = dir;
    inobj.Spectrum.date = offobj.Spectrum.date;
    inobj.Spectrum.depth = depi;
end

%%
function [shoal,offdir] = getShoal(obj,inobj,swl,idx) %...getShoal(sptobj,dir,freq,swl,idx)     
    %compute the shoaling coeffient and offshore direction arrays   
    offdst = obj.Data.Offshore;              %offshore properties  
    indst = obj.Data.Inshore;                %inshore properties    
    T = offdst.Dimensions.Period;            %wave periods used in ray model
    fray = 1./T;                             %wave frequencies used in ray model
    %zwl = offdst.Dimensions.WaterLevel;      %water levels used in ray model 
    offdst = offdst.DataTable;
    indst = indst.DataTable;

    %extract direction & celerity data and apply limits
    theta = offdst.theta;  theta(idx) = NaN; %offshore direction    
    c0 = offdst.celerity;  c0(idx) = 0;      %offshore celerity
    cg0 = offdst.cgroup;   cg0(idx) = 0;     %offshore group celerity
    ci = squeeze(indst.celerity);            %inshore celerity
    cgi = squeeze(indst.cgroup);             %inshore group celerity

    
    Tdef = [2,3,4,5,6,7,8,9,10,11.8,13.3,15.4,18.2,22.2,28.6,40];
    %pad the high frequencies with values from the maximum frequency
    idx = find(Tdef<min(T)); 
    if ~isempty(idx)
        addfray = 1./Tdef(idx);
        nadd = numel(addfray);
        fray = [addfray';fray];  %pad wave ray frequencies to maximum (1/2s)
        theta = [repmat(theta(:,1,:),1,nadd),theta];
        c0 = [repmat(c0(:,1,:),1,nadd),c0];
        cg0 = [repmat(cg0(:,1,:),1,nadd),cg0];

        if isrow(ci), ci = ci'; end          %force a column vector 
        if isrow(cgi), cgi = cgi'; end       %force a column vector 
        addci = repmat(ci(1,:),nadd,1);         %inshore celerities for highest frequency
        addcgi = repmat(cgi(1,:),nadd,1);
        cdeepwater = 9.81./addfray/2/pi;     %deepwater celerity = gT/2pi
        cdeep = repmat(cdeepwater',1,size(ci,2));
        addci = min(addci,cdeep);
        addcgi = min(addcgi,cdeep/2);        %deepwater group celerity = c/2

        ci = [addci;ci];                                                    
        cgi = [addcgi;cgi];
    end

    %pad the high frequencies with values from the minimum frequency
    idx = find(Tdef>max(T)); 
    if ~isempty(idx)
        addfray = 1./Tdef(idx);
        nadd = numel(addfray);
        fray = [fray;addfray'];  %pad wave ray frequencies to minimum (1/40s)
        theta = [repmat(theta(:,end,:),1,nadd),theta];
        c0 = [c0,repmat(c0(:,end,:),1,nadd)];
        cg0 = [cg0,repmat(cg0(:,end,:),1,nadd)];

        if isrow(ci), ci = ci'; end          %force a column vector 
        if isrow(cgi), cgi = cgi'; end       %force a column vector 
        addci = repmat(ci(end,:),nadd,1);      %inshore celerities for highest frequency
        addcgi = repmat(cgi(end,:),nadd,1);
        cdeepwater = 9.81./addfray/2/pi;     %deepwater celerity = gT/2pi
        cdeep = repmat(cdeepwater',1,size(ci,2));
        addci = min(addci,cdeep);
        addcgi = min(addcgi,cdeep/2);        %deepwater group celerity = c/2

        ci = [ci;addci];                                                    
        cgi = [cgi;addcgi];
    end

    %replicate grid vectors to produce grids for interpn using
    %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
    [PFW,XoFoWo,fw,fowo] = transferDims(obj,inobj,fray,swl);
    %interpolate offshore directions and celerity arrays
    offdir = interpn(PFW{:},theta,XoFoWo{:});           %offshore directions array
    c0fdw = interpn(PFW{:},c0,XoFoWo{:},'linear',0);    %offshore celerity frequency,direction,water level array            
    cg0fdw = interpn(PFW{:},cg0,XoFoWo{:},'linear',0);  %offshore group celerity frequency,direction,water level array
    cifw = interpn(fw{:},ci,fowo{:});                   %inshore celerity frequency,water level array
    cgifw = interpn(fw{:},cgi,fowo{:});                 %inshore group celerity frequency,water level array

    %calculate the shoaling coefficient 
    ishoal = cifw.*cgifw;                               %inshore c.cg
    if iscolumn(ishoal), ishoal = ishoal'; end          %force a column vector 
    oshoal = c0fdw.*cg0fdw;                             %offshore c.cg
    shoal = oshoal./ishoal;                             %shoaling factor array 
            %celerity tables are in ascending period and the spectrum is defined in
            %terms of ascending frequency. Flip array to give ascending frequencies
            %shoal = fliplr(shoal); 
            
    %check plot of surf (make sure commented out when running 
    %coefficientsPlot as the loop will produce thousands of plots!!)
    surfobj = ctWaveSpectraPlots;
    surfobj.Spectrum.SG = shoal; 
    surfobj.Spectrum.freq = inobj.Spectrum.freq; 
    surfobj.Spectrum.dir = inobj.Spectrum.dir;
    surfobj.Plotxt.ttxt = 'Test shoal function';
    surfPlot(surfobj,'on');
end
