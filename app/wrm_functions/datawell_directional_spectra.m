function D = datawell_directional_spectra(dirs,isplot,varargin)
%
%-------function help------------------------------------------------------
% NAME
%   datawell_directional_spectra.m
% PURPOSE
%   Estimates the directional distribution of a wave spectrum 
%   for directions, dirs, given the mean, spread, skewness and kurtosis 
%   parameters as output by datawell buoys SPT file format.
% USAGE
%   D = datawell_directional_spectra(dirs,isplot,m,s,sk,ku)
%   or
%   D = datawell_directional_spectra(dirs,isplot,dst)
% INPUTS
%   dirs - directions to define directional spread (degTN)
%   isplot - option to plot the distribution (default is false)
%   varargin can be:
%     dst - a dstable of an SPT format file loaded using wave_cco_spectra.m
%   or:
%     fr - frequency of statistical directional parameters
%     mn - mean direction (degTN)
%     sp - directional spread (deg)
%     sk - skew of distribution (-)
%     ku - kurtosis of distribution (-)
% OUTPUTS
%   D - directional distribution for given x values and the defined
%   statistical parameters.
% NOTES
%   Uses equations 14.2.24 and 14.2.27 in the Datawell Waves5 Reference 
%   Manual, 2023, p224-5.
% SEE ALSO
%   A. Kuik, G. P. van Vledder and L. Holthuijsen, “A Method for the 
%   Routine Analysis of Pitch-and-Roll Buoy Wave Data,” Journal of Physical 
%   Oceanography, vol. 18, no. 7, p1020–1034, 1 July 1988. 
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2023
%----------------------------------------------------------------------
%
    if isempty(isplot)
        isplot = false;
    end
    dirs = dirs(:);   %force column vecotr

    if isscalar(varargin)
        dst = varargin{1};
        varargin = [{dst.Spectra.Dimensions.freq},table2cell(dst.Spectra.DataTable)];
        % fr = dst.Spectra.Dimensions.freq;
        % mn = dst.Spectra.Dir;
        % sp = dst.Spectra.Spr;
        % sk = dst.Spectra.Skew;
        % ku = dst.Spectra.Kurt;
    elseif length(varargin)~=5
        warndlg('Incorrect input. Requires a dstable or 5 variables')
    end
    fr = varargin{1};
    mn = varargin{2};
    sp = varargin{3};
    sk = varargin{4};
    ku = varargin{5};   

    nfreq = length(fr);
    nmean = length(mn);
    nspr = length(sp);
    nskew = length(sk);
    nkurt = length(ku);

    %check that statistical parameters are scalar or vectors of the same length
    check = @(var,nvar) (isscalar(var) || nvar==nfreq);
    assert(check(mn,nmean) && check(sp,nspr) && check(sk,nskew) && check(ku,nkurt),...           
           'Statistical parameters must be scalar or vectors of the same length');

    %compute the centred Fourier coefficients
    theta0 = mn*pi/180;               %convert mean direction to radians
    sig = sp*pi/180;                  %convert spread to radians
    m1 = 1-0.5*sig.^2;
    m2 = 0.5*(ku.*sig.^4-6+8*m1);
    n2 = -sk.*real((0.5*(1-m2)).^(3/2));
    % denom = 0.5*(1-m2);
    % n2 = zeros(size(denom)); %zeros is as n2 using real above, NaN blanks invalid f
    % valid = denom>0;
    % n2(valid) = -sk(valid).*(denom(valid).^(3/2));

    %compute the directional components for direction and frequency range
    theta = dirs*pi/180;
    ndir = length(theta);
    D = zeros(ndir,nfreq);

    for f=1:nfreq
        if isnan(theta0(f))
            D(:,f) = 0;
        else
            dtheta = mod(theta-theta0(f)+pi,2*pi)-pi; %wrap trap
            dd = (1/pi)*(0.5+m1(f)*cos(dtheta)+m2(f)*cos(2*dtheta)+...                           
                                              n2(f)*sin(2*dtheta)); 
            % idx = isangletol(theat,[theta0(f)-pi,theta0(f)+pi]); %limit to mean direction x0(f)+/-pi
            % dd = dd.*idx;  
            dd(dd<0) = 0;                      %exclude negative values
            %the integral of D(:,f)=1 the following uses wrapped weights for
            %to integrate D rather than trapz
            w = trapz_weights_periodic(theta);
            D(:,f) = dd./sum(dd.*w);  %renormalise because of dd<0 clipping
        end                           %hence need to restore D(:,j)=1
    end    

    %plot directional distribution surface if required
    if isplot         
        getPlot(fr,dirs,D);
    end
end

%%
function getPlot(fr,dirs,D)
     %plot directional distribution surface if required
    hf = figure('Name','Direction spread','Tag','PlotFig');
    ax = axes(hf);
    [X,Y] = meshgrid(fr,dirs);
    surf(ax,X,Y,D);
    view(2);
    shading interp
    axis tight
    ylabel('Direction')
    xlabel('Frequency (Hz)')
    colorbar
end