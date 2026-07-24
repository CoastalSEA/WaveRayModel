function wrm_adjust_scatter_points(hfig,pntsz,edgecolor)
%
%-------function help------------------------------------------------------
% NAME
%   wrm_adjust_peclet_points.m
% PURPOSE
%   find multiple axes in a tiled plot and change the marker size of 
%   scatter points
% USAGE
%   >> hfig = gcf;
%   >> wrm_adjust_scatter_points(hfig,pntsz,edgcolor)
% INPUT
%   hfig - handle to figure
%   pntsz - adjust scatter property SizeData to pntsz 
%   edgecolor - adjustscatter property MarkerEdgeColor to edgecolor
% OUTPUT
%   updated tiles in plot
% SEE ALSO
%   compile_tiled_figure, clear_figure_tile, replace_tile, set_tile_legend
%
% Author: Ian Townend & Copilot
% CoastalSEA (c)July 2026
%----------------------------------------------------------------------
%    
    if isa(hfig.Children,'matlab.graphics.layout.TiledChartLayout')
        axs = hfig.Children.Children;
    else
        axs = hfig.Children;
    end
    axs = findobj(axs,'Type','axes');
    nax = numel(axs);
    for i=1:nax
        hp = findobj(axs(i),'Type','scatter');
        nrec = numel(hp);
        sz = num2cell(ones(1,nrec)*pntsz);
        [hp(:).SizeData] = sz{:};
        if nargin==3
            %change edgecolor of negative points
            yellow = [0.929,0.694,0.125];
            for j=1:numel(hp)
                hpecneg = findobj(hp(j),'CData',yellow);
                if isempty(hpecneg), continue; end
                hpecneg.MarkerEdgeColor = edgecolor;
            end
        end
    end
end