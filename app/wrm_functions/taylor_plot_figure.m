function ax = taylor_plot_figure(rLim)
%
%-------function help------------------------------------------------------
% NAME
%   taylor_plot_figure.m
% PURPOSE
%   Create base plot for a Taylor diagram 
% USAGE
%   ax = taylor_plot_figure(rlim);
% INPUTS
%   rlim - the radial limit for the plot
% OUTPUT
%   ax - axes to plot of Taylor diagram
% NOTES
%   Taylor, K, 2001, Summarizing multiple aspects of model performance 
%   in a single diagram, JGR-Atmospheres, V106, D7. 
%   MSE contours based on Exchange Forum code by Guillaume Maze.
% SEE ALSO
%   Function plot_spectrum_model_skill.m in WaveRayModel
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2025
%--------------------------------------------------------------------------
%
    h_fig = figure('Name','Taylor Diagram', ...
        'Units','normalized', ...
        'Resize','on','HandleVisibility','on', ...
        'Tag','StatFig');
    h_fig.Position(1) = 1-h_fig.Position(3)-0.01;  %top right
    h_fig.Position(2) = 1-h_fig.Position(4)-0.12;
    polaraxes(h_fig);
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
    ax.ThetaLim = [0,90];
    ints = [0,0.2,0.4,0.6,0.8,0.9,0.92,0.94,0.96,0.98,0.99,1.0];
    ax.ThetaTick = (asin(ints)*180/pi);
    ax.ThetaTickLabelMode = 'manual';
    ax.ThetaTickLabel = ints;
    ax.RLim = [0,rLim];
    ax.NextPlot = 'add';
    hold on
    %plot normalised RMSE circles at 0.5 intervals
    R = 0.25:0.25:1.0;
    tau = 0:0.05:pi;
    for ii=1:4
        theta = atan((R(ii).*sin(tau))./(1+R(ii).*cos(tau)));
        tct = 2*cos(theta);
        fac = 4*(1-R(ii)^2);
        rad1 = (tct+sqrt(tct.^2-fac))/2;
        rad2 = (tct-sqrt(tct.^2-fac))/2;                
        [~,idrad] = min(rad1);
        rad = [rad1(1:idrad),rad2(idrad+1:end)];
        hp = polarplot(ax,(pi/2-theta),rad,'--',...
            'LineWidth',0.4,'Color',[0.8 0.8 0.8],...
            'DisplayName','RMS error','Tag','RMSgrid');
        if ii>1
        set(get(get(hp,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend 
        end
        text((pi/2-theta(10)),rad(10)-0.05,num2str(R(ii)),...
            'Color',[0.82 0.82 0.82]);
    end            
    hold off
    %Add axis labels (not part of polar plot)
    uiNameValue('Normalized Std. Dev',[0.38 -0.08],0); 
    uiNameValue('Normalized Std. Dev.',[-0.05 0.38],90); 
    uiNameValue('Correlation Coefficient',[0.75 0.85],-50); 
    uicontrol('Style', 'pushbutton', 'String', 'Case list',...
            'Units','normalized','Position', [0.8 0.06 0.15 0.05],...
            'Callback', @uiCaseList);  
    %-nested function------------------------------------------------------
    function uiNameValue(uitext,uipos,uirot)
        %function to generate text label at sppecified position
        %and angle
        h_ui = text(1,1,uitext);
        h_ui.LineStyle = 'none';
        h_ui.String = uitext;
        h_ui.HorizontalAlignment = 'left';
        h_ui.FontUnits = 'normalized';
        h_ui.FontSize = 0.04;
        h_ui.Rotation = uirot;
        h_ui.Units = 'normalized';
        h_ui.Position = uipos;
    end
    %-nested function------------------------------------------------------
    function uiCaseList(~,~)
        %create figure that lists meta-data for all cases listed in
        %the current Talyor Plot legend
        fig = findobj('Name','Taylor Diagram');
        figax = fig.CurrentAxes;
        figpts = findobj(figax,'Type','Line','-not','Tag','RMSgrid');
        figpts = sortplots(figpts);
        nrec = length(figpts);
        tstrings = cell(nrec,1);
        tstrings{1} = [figpts(1).DisplayName,': ',figpts(1).UserData];                                            
        for i=2:nrec
            tstrings{i} = [figpts(i).DisplayName,'; ',figpts(i).UserData];                             
        end
        if isempty(tstrings), return; end
        hg = figure('Name','Taylor Diagram Summary','Units','normalized',...                
            'Resize','on','HandleVisibility','on','Tag','PlotFig');
        hg.Position(1) = 1-hg.Position(3)-0.01;  %top right
        hg.Position(2) = 1-hg.Position(4)-0.52;
        txtpos = [0.025 0.025 0.95 0.85];
        ht = uicontrol(hg,'style','listbox','Units','normalized',...
        'Position',txtpos,'BackgroundColor',[0.92 0.92 0.92], ...
        'tag','CaseList');
        set(ht,'string',tstrings);
        %Create push button to copy data to clipboard
        uicontrol('Parent',hg,...
                'Style','pushbutton',...
                'String', 'Copy to clipboard',...
                'Units','normalized', ...
                'Position', [0.75 0.915 0.20 0.065], ...
                'UserData',tstrings, ...
                'Callback',@(src,evdata)mat2clip(src.UserData));
    end
end