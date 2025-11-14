function ax = taylor_plot_ts(ax,stats,skill,metatxt)
%
%-------function help------------------------------------------------------
% NAME
%   taylor_plot_ts.m
% PURPOSE
%   Add a timeseries of test points to a Taylor diagram (assumes common
%   reference point)
% USAGE
%   ax = taylor_plot_ts(ax,stats,skill,metatxt);
% INPUTS
%   ax - axes of base plot for a Taylor diagraam (taylor_plot_figure.m)
%   stats - 
%   skill - 
%   metatxt - 
% OUTPUT
%   ax - axes to plot of Taylor diagram
% NOTES
%   Taylor, K, 2001, Summarizing multiple aspects of model performance 
%   in a single diagram, JGR-Atmospheres, V106, D7. 
% SEE ALSO
%   Function plot_spectrum_model_skill.m in WaveRayModel
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2025
%--------------------------------------------------------------------------
%

    hold(ax,'on')
    useref = metatxt{1};
    polarplot(ax,pi/2,1,'+','LineWidth',2.0,'ButtonDownFcn',@godisplay,...
                    'DisplayName','Reference','Tag','0','UserData',useref);
   
    nrec = length(stats);
    for i=1:nrec
        %unpack and normalize data to be plotted
        ndteststd = stats(i).teststd/stats(i).refstd; %normalised std
        ndcrmsd = stats(i).crmsd/stats(i).refstd;     %normalised centred root mean square differences
        acoscor = asin(stats(i).corrcoef(1,2));
        bias = stats(i).testmean-stats(i).refmean;
        corr = stats(i).corrcoef(1,2);
        %check statistcs
        %RMSD-sqrt(testSTD^2+refSTD^2-2*testSTD*refSTD*COR)=0 (NB:refSTD=1)
        check = (ndcrmsd-sqrt(ndteststd^2+1-...
                        2*ndteststd*stats(i).corrcoef(1,2)));
        if check>0.1
            warndlg('Statistics do not agree. Error in plotTaylor');
            continue;
        end    
        datetxt = string(stats(i).date);
        %user data to construct table
        restxt = sprintf('%s: bias= %.3f; corr= %.3f; ndstd= %.3f with skill S.G= %1.3g',...
                              datetxt,bias,corr,ndteststd,stats(i).global);
        if ~isempty(stats(i).local)
            if skill.iter
                itxt = 'for all cells';
            else
                itxt = 'with no overlaps';
            end
            usertxt = sprintf('%s; S.L= %1.3g (Ro=%1.2g, n=%1.1g, W=%d, %s)',...
                     restxt,stats(i).local,skill.Ro,skill.n,skill.W,itxt);
        else
            usertxt = restxt;
        end
        %add point to plot
        hp = polarplot(ax,acoscor,ndteststd,'+','DisplayName',num2str(i),...           
                 'ButtonDownFcn',@godisplay,'Tag',num2str(i),'UserData',usertxt);
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end
    
    
    hold(ax,'off')
    %
legend(ax,'show','Location','northeastoutside');
        % %re-impose suppression of grid lines in the Taylor diagram plot
        % hp = findobj(figax,'Type','Line');
        % hgrd = findobj(hp,'Tag','RMSgrid'); 
        % hp = findobj(hp,'-not','Tag','RMSgrid'); 
        % newhp = vertcat(hp,hgrd(1));
 
end