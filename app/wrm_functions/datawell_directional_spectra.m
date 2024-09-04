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
%     f - frequency of statistical directional parameters
%     m - mean direction (degTN)
%     s - directional spread (deg)
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

    if length(varargin)==1
        dst = varargin{1};
        f = dst.Spectra.Dimensions.freq;
        m = dst.Spectra.Dir;
        s = dst.Spectra.Spr;
        sk = dst.Spectra.Skew;
        ku = dst.Spectra.Kurt;
    elseif length(varargin)==5
        f = varargin{1};
        m = varargin{2};
        s = varargin{3};
        sk = varargin{4};
        ku = varargin{5};
    else
        warndlg('Incorrect input. Requires a dstable or 5 variables')
    end

    nfreq = length(f);
    nmean = length(m);
    nspr = length(s);
    nskew = length(sk);
    nkurt = length(ku);

    %check that statistical parameters are scalar or vectors of the same length
    check = @(var,nvar) (isscalar(var) || nvar==nfreq);
    assert(check(s,nmean) && check(s,nspr) && check(sk,nskew) && check(ku,nkurt),...
        'Statistical parameters must be scalar or vectors of the same length');

    %compute the centred Fourier coefficients
    x0 = m*pi/180;               %convert mean direction to radians
    sig = s*pi/180;              %convert spread to radians
    m1 = 1-0.5*sig.^2;
    m2 = 0.5*(ku.*sig.^4-6+8*m1);
    n2 = -sk.*((1-m2)/2).^1.5;

    %compute the directional components for direction and frequency range
    x = mod(dirs*pi/180,2*pi);
    ndir = length(x);
    D = zeros(ndir,nfreq);
    for i=1:ndir
        for j=1:nfreq
            if x(i)<=x0(j)+pi && x(i)>=x0(j)-pi
                D(i,j) = (0.5+m1(j)*cos(x(i)-x0(j))+...
                          m2(j)*cos(2*(x(i)-x0(j)))+...
                          n2(j)*sin(2*(x(i)-x0(j))))/pi;
            end
        end
    end
    D(D<0) = 0;

    %convert x back to degrees
    x = mod(x*180/pi,360);
    if isplot         
        hf = figure('Name','Direction spread','Tag','PlotFig');
        ax = axes(hf);
        [X,Y] = meshgrid(f,x);
        surf(ax,X,Y,D);
        view(2);
        shading interp
        axis tight
        ylabel('Direction')
        xlabel('Frequency (Hz)')
        colorbar
    end
end