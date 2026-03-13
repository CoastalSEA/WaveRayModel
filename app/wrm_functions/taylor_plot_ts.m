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
    %unpack and normalize data to be plotted
    ndteststd = [stats(:).teststd]./[stats(:).refstd]; %normalised std
    ndcrmsd = [stats(:).crmsd]./[stats(:).refstd];     %normalised centred root mean square differences
    bias = [stats(:).testmean]-[stats(:).refmean];
    corrcoefs = [stats(:).corrcoef];
    corr = corrcoefs(2,1:2:end);                       %extact [2,1] value from square matrix
    acoscor = asin(corr);
    %check statistcs
    %RMSD-sqrt(testSTD^2+refSTD^2-2*testSTD*refSTD*COR)=0    (NB:refSTD=1)
    check =  (ndcrmsd-sqrt(ndteststd.^2+1-2*ndteststd.*corr));

    nrec = length(stats);
    for i=1:nrec
        %check statistcs are valid 
        if check(i)>0.1
            warndlg('Statistics do not agree. Error in plotTaylor');
            continue;
        end    
        datetxt = string(stats(i).date);
        %user data to construct table
        restxt = sprintf('%s: bias= %.3f; corr= %.3f; ndstd= %.3f with skill S.G= %1.3g',...
                              datetxt,bias(i),corr(i),ndteststd(i),stats(i).global);
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
        hp = polarplot(ax,acoscor(i),ndteststd(i),'+','DisplayName',num2str(i),...           
                 'ButtonDownFcn',@godisplay,'Tag',num2str(i),'UserData',usertxt);
        hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end

    %find centroid of points
    [N,Xedges,Yedges] = histcounts2(asin(corr),ndteststd);
    [~,idi] = max(N, [],'all');
    [idx,idy] = ind2sub(size(N),idi);
    xHist = (Xedges(idx)+Xedges(idx+1))/2;
    yHist = (Yedges(idy)+Yedges(idy+1))/2;
    score = 4*(1+sin(xHist))^skill.n/((yHist+1/yHist)^2*(1+skill.Ro)^skill.n);  

    %plot peak of 2D histogram
    usertxt = sprintf('%d points: bias= %s; corr= %.3f; ndstd= %.3f with skill S.G= %1.3g',...
                                        nrec,'n/a',sin(xHist),yHist,score);
    hp = polarplot(ax,xHist,yHist,'o',...
                'LineWidth',2.0,'MarkerSize',6,...
                'DisplayName','Histogram centroid','ButtonDownFcn',@godisplay,...           
                'Tag','Histogram marker','UserData',usertxt);
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';

    %add point to plot for sample means and text to match
    mn_ndteststd = mean(ndteststd,'All','omitnan');
    mn_corr = mean(corr,'All','omitnan');    
    mn_bias = mean(bias,'All','omitnan');
    mn_gskill = mean([stats(:).global],'All','omitnan');
    mn_acoscor = asin(mn_corr);
    nrec = numel(ndteststd);
    %plot mean value of all points
    usertxt = sprintf('%d points: bias= %.3f; corr= %.3f; ndstd= %.3f with skill S.G= %1.3g',...
                            nrec,mn_bias,mn_corr,mn_ndteststd,mn_gskill);
    hp = polarplot(ax,mn_acoscor,mn_ndteststd,'x',...
                'LineWidth',2.0,'MarkerSize',10,...
                'DisplayName','Sample Means','ButtonDownFcn',@godisplay,...           
                'Tag','Mean marker','UserData',usertxt);
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

    %add legend
    hold(ax,'off')
    legend(ax,'show','Location','northeastoutside');
        % %re-impose suppression of grid lines in the Taylor diagram plot
        % hp = findobj(figax,'Type','Line');
        % hgrd = findobj(hp,'Tag','RMSgrid'); 
        % hp = findobj(hp,'-not','Tag','RMSgrid'); 
        % newhp = vertcat(hp,hgrd(1)); 
end