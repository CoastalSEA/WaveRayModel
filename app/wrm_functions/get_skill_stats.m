function stats = get_skill_stats(obsobj,modobj,skill)
    nrec = size(obsobj,1);
    stats = initialiseStruct(nrec);
    hpw = PoolWaitbar(nrec, 'Processing skill statistics');  %and increment(hpw);
    parfor i=1:nrec
        ds1 = obsobj(i).Spectrum.SG;
        ds2 = modobj(i).Spectrum.SG;
        stats(i) = DS_DifferenceStatistics(ds1,ds2);
        stats(i).global = getSkill(skill,stats(i));     
        if skill.Inc      %get local skill score if required
            stats(i).local = getLocalSkill(skill,ds1,ds2);
        else
            stats(i).local = [];
        end
        stats(i).date = obsobj(i).Spectrum.date; 
        increment(hpw);
    end
    delete(hpw)
end
%%
function cfstats = DS_DifferenceStatistics(ds1,ds2)
    %compute the difference between two vectors or arrays of the same size
    if isa(ds1,'dstable')
        ds1 = ds1.DataTable{:,1};
        ds2 = ds2.DataTable{:,1};
    end
    
    if numel(ds1)~=numel(ds2)
        warndlg('Data sets must be of the same size');
        return
    end
    cfstats = initialiseStruct(1);
    cfstats.refstd = std(ds1,0,'All','omitnan'); %requires R2018b or later
    cfstats.teststd = std(ds2,0,'All','omitnan');
    cfstats.refmean = mean(ds1,'All','omitnan');
    cfstats.testmean = mean(ds2,'All','omitnan');
    cfstats.corrcoef = corrcoef(ds1,ds2,'rows','complete');
    cfstats.crmsd = centredRMSD(ds1,ds2,cfstats);
end

%%
function cRMSD = centredRMSD(ds1,ds2,cfstats)
    %compute centred root mean square difference 
    %ds1 is the reference data set and ds2 is the test data set
    cSD = ((ds2-cfstats.testmean)-(ds1-cfstats.refmean)).^2;
    cSD = cSD(~isnan(cSD));
    nrec = length(cSD);
    cMSD = sum(cSD)/nrec;
    cRMSD = sqrt(cMSD);
end

%%

%%
function score = getSkill(skill,cfstats)
    %get the skill score as defined in Taylor, K. E. (2001). 
    %"Summarizing multiple aspects of model performance in a single diagram." 
    %Journal of Geophysical Research - Atmospheres 106(D7): 7183-7192.Eq(4)
    Ro = skill.Ro;
    n = skill.n;
    sigobs = cfstats.teststd/cfstats.refstd;
    R = cfstats.corrcoef(1,2);
    score = 4*(1+R)^n/((sigobs+1/sigobs)^2*(1+Ro)^n);   
end

%%
function score = getLocalSkill(skill,refvar,testvar)
    %iterate over data set based on interval W defined in skill struct
    %iter =  true iterates for i=1:m-2W,
    %iter = false avoids overlaps and iterates over i=1:2W:m--2W
    W = skill.W;
    [m,n] = size(refvar);
    %set up sampling index to either examine a window at every point or a
    %set of windows that do not overlap
    if skill.iter %true iterates over every point, i=1:m-2W
        indx = 1:m-2*W;
        indy = 1:n-2*W;
    else         %false avoids overlaps and iterates over i=1:2W:m-2W
        indx = 1:2*W:m-2*W;
        indy = 1:2*W:n-2*W;
    end
    
    %data is a 2-D array
    ni = 1; nj = 1;
    ss = zeros(length(indx),length(indy));
    for i=indx
        for j=indy
            subref = refvar(i:i+2*W,j:j+2*W);
            subtest = testvar(i:i+2*W,j:j+2*W);
            cfstats = DS_DifferenceStatistics(subref,subtest);
            ss(ni,nj) = getSkill(skill,cfstats);
            nj = nj+1;
        end
        nj = 1;
        ni = ni+1;
    end
    [ms,ns] = size(ss);
    xr = m/ms; yr = n/ns; %ratios to rescale indices based on
    %window used for local skill
    skill.SD.x = round(skill.SD.x/xr);
    skill.SD.y = round(skill.SD.y/yr);

    %apply subdomain mask
    if ~isempty(skill.SD)   
        [X,Y] = meshgrid(1:size(ss,1),1:size(ss,2));
        in = inpolygon(X,Y,round(skill.SD.x),round(skill.SD.y));
        ss(~in') = NaN;
    end

    %subdomain average ('All' only implemented for >2018a)
    score = mean(ss,'All','omitnan');
end

%%
function stats = initialiseStruct(nrec)
    %initialise an empty struct
    stats(nrec).refstd = [];
    stats(nrec).teststd = [];
    stats(nrec).refmean = [];
    stats(nrec).testmean = [];
    stats(nrec).corrcoef = [];
    stats(nrec).crmsd = [];
    stats(nrec).global = [];
    stats(nrec).local = [];
    stats(nrec).date = [];
end

