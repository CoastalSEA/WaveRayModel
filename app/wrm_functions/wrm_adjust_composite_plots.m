function wrm_adjust_composite_plots(hfig,option)
%find tiles in a tiled plot and adjust (edit as needed)
%   'option' can be any of the following:
%   'Clear off diagonal tile' - in a 3x3 tile array removes tiles 2,4,6,8
%   'Clear colorbar' - delete ALL colorbars
%   'Adjust colormap' - Change the figure colormap
%   'Disable toolbar' - disable hover data tips and floating axes toolbar for all tiles
%
    if isa(hfig.Children,'matlab.graphics.layout.TiledChartLayout')
        htile = hfig.Children.Children;
    else
        htile = hfig.Children;
    end

    switch option
        case 'Clear off diagonal tile'
            if isa(hfig.Children,'matlab.graphics.layout.TiledChartLayout')
                tiledLayoutObj = hfig.Children;
                %remove tiles 
                rowcol = [1,2;2,1;2,3;3,2]; 
                for i=1:size(rowcol,1)
                    clear_figure_tile(tiledLayoutObj, rowcol(i,1), ...
                                            rowcol(i,2), 'remove',true);
                end
            end

        case 'Clear colorbar'
            %delete all colorbars
            hcb = findobj(htile,'Type','colorbar');
            delete(hcb)

        case 'Adjust colormap' 
            %correct colormap
            cmap = cmap_selection();
            colormap(cmap);

        case 'Disable toolbar'
            %disable hover data tips and floating axes toolbar
            axs = findobj(htile,'Type','axes');
            for i=1:numel(axs)
                 ax = axs(i);
                 ax.Toolbar.Visible = 'off';
                 ax.Interactions = [];
            end
    end
end




