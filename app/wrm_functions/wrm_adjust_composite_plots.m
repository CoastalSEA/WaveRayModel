function wrm_adjust_composite_plots(hfig)
%find tiles in a tiled plot and adjust (edit as needed)

    hobj = hfig.Children.Children;
    axs = findobj(hobj,'Type','axes');
    

    if isa(hfig.Children,'matlab.graphics.layout.TiledChartLayout')
        tiledLayoutObj = hfig.Children;
        %remove tiles 
        rowcol = [1,2;2,1;2,3;3,2]; 
        for i=1:size(rowcol,1)
            clear_figure_tile(tiledLayoutObj, rowcol(i,1), rowcol(i,2),...
                                                         'remove',true);
        end

        hobj = hfig.Children.Children;
        hcb = findobj(hobj,'Type','colorbar');
        delete(hcb)

        %correct colormap
        colormap(cmap_selection(19));
    end



end




