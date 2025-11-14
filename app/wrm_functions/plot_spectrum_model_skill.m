function ax = plot_spectrum_model_skill(obsobj,modobj,mobj)
    %compute the skill of model v measured spectrum data and produce Taylor
    %plot of timeseries results
    %calls MS_RunParams class and taylor_plot function    
    x = obsobj(1).Spectrum.dir;              %direction
    y = obsobj(1).Spectrum.freq;             %frequency
    skill = getSkillParameters(mobj,x,y);    %get the parameters for skill model

    %get the statistics for the Taylor plot
    stats = get_skill_stats(obsobj,modobj,skill);

    %add the timeseries results to the Taylor plot
    ndteststd = [stats(:).teststd]./[stats(:).refstd]; %normalised std
    rLim = ceil(max(ndteststd));                        %radial limit for the plot
    ax = taylor_plot_figure(rLim);    
    metatxt = {'Measured','Model'};
    ax = taylor_plot_ts(ax,stats,skill,metatxt);
    subtitletxt = sprintf('%s',getModelInputText(modobj(1),0));
    ax.Title.String = obsobj(1).inpData.tsdst.Description;
    ax.Subtitle.String = subtitletxt;  
end

%%
function skill = getSkillParameters(mobj,x,y)
    %extract Skill parameters using MS_RunParams class for input
    robj = MS_RunParams.setInput(mobj);  %default or current values if user cancels

    skill.Ro = robj.maxcorr;
    skill.n  = robj.skillexponent;
    skill.W = robj.skillwindow;
    if skill.W==0
        skill.Inc = false; 
    else
        skill.Inc = true;                   %flag to include skill score
    end
    subdomain = robj.skillsubdomain;
    skill.SD = getSubDomain(x,y,subdomain);
    skill.iter = robj.skilliteration;
end

%%
function sd = getSubDomain(x,y,subdomain)
    %find the subdomain in integer grid indices defined by x,y range
    %subdomain defined as [x0,xN,y0,yN];    
    if isempty(subdomain) || length(subdomain)~=4
        subdomain = [min(x),max(x),min(y),max(y)];
    end
    ix0 = find(x<=subdomain(1),1,'last');
    ixN = find(x>=subdomain(2),1,'first');
    iy0 = find(y<=subdomain(3),1,'last');
    iyN = find(y>=subdomain(4),1,'first');
    sd.x = [ix0,ix0,ixN,ixN];
    sd.y = [iyN,iy0,iy0,iyN];
end