function wrm_adjust_peclet_points(hfig,pntsz,edgcolor)
%find multiple axes in a tiled plot and change the marker size of scatter
%points
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
        % yellow = [0.929,0.694,0.125];
        % for j=1:numel(hp)
        %     hpecneg = findobj(hp(j),'CData',yellow);
        %     if isempty(hpecneg), continue; end
        %     hpecneg.SizeData = edgcolor;
        % end
    end
end