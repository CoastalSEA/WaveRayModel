classdef SpectralTransfer < muiDataSet                    
%-------class help---------------------------------------------------------
% NAME
%   SpectralTransfer.m
% PURPOSE
%   Class description - Builds the offshore and inshore Transfer Tables from
%   a backward ray tracking data set (class RayTracks) for use in
%   WRM_WaveModel. Also has a method to create plots of the transfer
%   coefficients for a unit wave height.
%
% SEE ALSO
%   muiDataSet, WaveRayModel, RayTracks, WRM_WaveModel.
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        interp                    %struct for interpolation settings
    end
    
    methods (Access={?muiDataSet,?muiStats,?WRM_WaveModel,?ctWaveSpectrum})
        function obj = SpectralTransfer()             
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model to construct spectral transfer table
%--------------------------------------------------------------------------   
        function obj = runModel(mobj,rayobj,sptname)
            %create offshore and inshore spectral transfer tables
            obj = SpectralTransfer;                           
            muicat = mobj.Cases;    
%--------------------------------------------------------------------------
% Model code>
%--------------------------------------------------------------------------
            %select back tracking ray case to use
            if isempty(rayobj)
                promptxt = 'Select Backward Ray Trace Dataset:';
                rayobj = selectCaseObj(muicat,{'backward_model'},[],promptxt); 
                if isempty(rayobj), return; end
            end            
            rayrec = caseRec(muicat,rayobj.CaseIndex);
            %assign the run parameters to the model instance
            setRunParam(obj,mobj,rayrec); %input caserecs passed as varargin 

            [results,indir,inprops] = specTransfer(obj,rayobj);            
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------  
            %inshore celerities for period and water level
            inccg = inprops(4:5);
            dspi = modelDSproperties(obj,true);
            indst = dstable(inccg{:},'DSproperties',dspi);
            indst.Dimensions.Period = rayobj.Data.Dataset.Dimensions.Period; 
            indst.Dimensions.WaterLevel = rayobj.Data.Dataset.Dimensions.WaterLevel; 
            %assign metadata about model
            indst.Source = metaclass(obj).Name;
            indst.MetaData = sprintf('Derived using %s (cid: %d)',...
                          rayobj.Data.Dataset.Description,rayobj.CaseIndex);
            indst.UserData.Location = [inprops{1},inprops{2}];
            indst.UserData.Depths = inprops{3};
            indst.UserData.ShoreAngle = rayobj.RunParam.WRM_RunParams.ShorelineAngle;

            %offshore celerities for direction, period and water level
            dspo = modelDSproperties(obj,false);
            offdst = dstable(results{:},'RowNames',indir,'DSproperties',dspo);
            offdst.Dimensions.Period = rayobj.Data.Dataset.Dimensions.Period; 
            offdst.Dimensions.WaterLevel = rayobj.Data.Dataset.Dimensions.WaterLevel;                        
            %assign metadata about model
            offdst.Source = metaclass(obj).Name;
            offdst.MetaData = sprintf('Derived using %s (cid: %d)',...
                        rayobj.Data.Dataset.Description,rayobj.CaseIndex);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            dst.Inshore = indst;
            dst.Offshore = offdst;
            %save results

            if isempty(sptname)            %single case
                sptname = sprintf('Transfer %s',rayobj.Data.Dataset.Description);
                setDataSetRecord(obj,muicat,dst,'spectral_model',{sptname},false);
                getdialog('Run complete');
            else                           %part of a batch run
                setDataSetRecord(obj,muicat,dst,'spectral_model',{sptname},true);
            end          
        end    

%%
        function [obj,meta] = getSTcase(mobj,isbatch)
            %select a spectral transfer case to use
            if nargin<2, isbatch = false; end
        
            msgtxt = 'Spectral Transfer Table not found';
            meta.source = 'Spectral transfer table';
            if isbatch     %select ST tables for multiple points
                [obj,ok] = selectRecord(mobj.Cases,...
                                            'CaseClass',{'SpectralTransfer'},...
                                            'PromptText', 'Select cases',...
                                            'ListSize', [250,200],...
                                            'SelectionMode', 'multiple' ); 
                if ok<1, getdialog(msgtxt); return; end
                meta.caserecs = obj(1);
            else            %select a single case to use
                promptxt = 'Select a Transfer Table Case to use:';             
                obj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);                
                if isempty(obj), getdialog(msgtxt); return; end
        
                meta.inptxt = sprintf('%s used for spectral transfer',...
                                          obj.Data.Inshore.Description);
                meta.caserecs = caseRec(mobj.Cases,obj.CaseIndex);
            end  
        end

%%
        function runPlotSpectrum(mobj)
            %create a plot of the offshore and inshore 2-D specrum surfaces 
            %for a single wave condition. Uses ctWaveSpectraPlots fucntions
            % but includes transfer to obtain inshore spectrum
            ptype = questdlg('What type of plot','XY Polar','XY','Polar','XY'); 
            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            obj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            if isempty(obj)
                getdialog('Spectral Transfer table not found'); return; 
            end
            
            [tsdst,meta] = waveModels.getInputData(mobj);
            if isempty(tsdst), return; end
            swl = tsdst(1).swl;
            tsdst(1) = removevars(tsdst(1),'swl');          %remove swl from dstable
            
            %select record from dateset get spectrum and plot
            dates = tsdst(1).DataTable.Properties.RowNames;
            ok = 0;
            spselect = [];
            while ok<1 
                irow = listdlg("PromptString",'Select event to plot',...
                         'SelectionMode','single','ListSize',[160,300],...
                         'ListString',dates);
                if isempty(irow), return; end
                spobj = ctWaveSpectraPlots;             %initialise class object
                spobj.spModel.selection = spselect;
                spobj = getSpectrumObject(spobj,meta.inptype,tsdst,irow);
                getPlot(spobj,ptype,'off');
                if ~isempty(spobj.spModel)
                    spselect = spobj.spModel.selection;
                end

                isout = checkWLrange(obj,swl(irow));
                if isout
                    warndlg('Water levels are outside the range of the Transfer Table')
                    return;
                else
                    spobj.inpData.swl = swl(irow);
                end

                inobj = get_inshore_spectrum(obj,spobj(1)); %returns ctWacveSpectrum object
                inobj.Params = wave_spectrum_params(inobj);
                inobj.Plotxt.vtxt = 'Inshore Spectral Energy (m^2s)';
                inobj.Plotxt.ttxt = 'Inshore spectrum';
                spobj(2) = Spectrum2SpectralPlots(inobj);      %convert to ctWaveSpectraPlots
                %call plot function
                getPlot(spobj(2),ptype,'off');
                shoreNorm = obj.Data.Inshore.UserData.ShoreAngle+90;
                addShoreNormal(spobj(2),shoreNorm);

                hf = getMultiPlot(spobj);
                %add button to access wave parameters
                summary = vertcat(spobj(1).Params,spobj(2).Params);
                if height(summary)==2
                    summary.Properties.RowNames = {'Offshore','Inshore'};
                else
                    rnames = spobj(1).Params.Properties.RowNames;
                    summary.Properties.RowNames = [rnames(:)',{'Inshore'}];
                end
                addDataButton(spobj,hf,summary);
            end       
        end

%%
        function runAnimation(mobj)
            %create an animation of the 2-D spectrum surfaces using a
            %timeseries input
            obj = WRM_WaveModel; 
            [tsdst,meta] = obj.getInputData(mobj);
            if isempty(tsdst), return; end   %user cancelled data selection 

            if height(tsdst(1))>5000
                promptxt = sprintf('Times series contains %d records\nThis could take a while to run and genearte large file\nUse time sub-selection to extract shorter time period',...
                                                            height(tsdst(1)));
                answer = questdlg(promptxt,'Time','Continue','Abort','Abort');                                  
                if strcmp(answer,'Abort'), return; end
            end

            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            obj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
             if isempty(obj)
                getdialog('Spectral Transfer table not found'); return; 
             end

            [offobj,inobj] = runWaves(obj,tsdst,meta);

            if strcmp(offobj(1).inpData.input,'Spectrum')
                tsdst(1) = addvars(tsdst(1),tsdst(2).Hs,'NewVariableNames',{'Hs'});
            end   

            wrm_animation(mobj,obj,tsdst(1),offobj,inobj)
        end

    end
%%
    methods
%--------------------------------------------------------------------------
% Model to transfer a wave timeseries or spectral data set
%--------------------------------------------------------------------------         
        function [offobj,inobj] = runWaves(obj,tsdst,meta)
            %run the spectral transfer model for a timeseries of offshore
            %wave conditions and return a table of wave conditions            
            % islog = false;
            % answer = questdlg('Write log of missing dates to file?',...
            %                        'SpectralTransfer','Yes','No','No');
            % if strcmp(answer,'Yes'), islog = true; end
            % filename = sprintf('Sprectra_log_%s.txt',char(datetime,"ddMMMyy_HH-mm"));

            spobj = ctWaveSpectrum;
            spobj.inpData.input = meta.inptype; 
            inptype = meta.inptype; 
            swl = tsdst(1).swl;                          %extract swl as vector
            if strcmp(inptype,'Spectrum')      
                spobj.inpData.output = 'Measured';           %dummy value                
                tsdst = removevars(tsdst(1),'swl');          %remove swl from dstable
                                                             %*needed for getDSTable to work
            else
                spobj = setSpectrumModel(spobj);  %define the model to be used (Jonswap etc)
                if isempty(spobj.spModel), offobj = []; inobj = []; return; end
                spobj.inpData.output = 'Modelled'; %TMA, etc used in get_inshore_spectrum
            end
                
            tsdst(1).DataTable = rmmissing(tsdst(1).DataTable);%remove nans
            nint = height(tsdst(1).DataTable);
            hpw = PoolWaitbar(nint, 'Processing timeseries');
            offobj(nint) = copy(spobj);
            inobj(nint) = copy(spobj);
            parfor i=1:nint                                %parfor loop  
                %for each offshore wave get the inshore results
                offspectrum = copy(spobj);
                tsdstrow = ctWaveSpectrum.getDatasetRow(tsdst,i);
                %set input parameters for selected record
                offspectrum = setInputParams(offspectrum,tsdstrow,inptype);
                offspectrum.inpData.swl = swl(i);
                %get the spectrum data for selected record
                if strcmp(inptype,'Spectrum')
                    offspectrum = getMeasuredSpectrum(offspectrum);
                elseif strcmp(inptype,'Wind')
                    offspectrum = getModelSpectrum(offspectrum);
                else
                    offspectrum = getMultiModalSpectrum(offspectrum);
                end
                offspectrum.Spectrum.date = tsdstrow.RowNames;
                offspectrum.Params = wave_spectrum_params(offspectrum);
                offobj(i) = offspectrum;

                %inshore spectrum and wave parameters
                inspectrum = get_inshore_spectrum(obj,offspectrum);
                inspectrum.Params = wave_spectrum_params(inspectrum);
                inobj(i) = inspectrum;
                coeftable = get_transfer_coefficients(offobj(i),inobj(i));
                inobj(i).Params = [inobj(i).Params,coeftable];
                % if islog
                %     lines = sprintf('%s',tsdst.RowNames(i)); %#ok<PFBNS> 
                %     writelines(lines,filename,WriteMode="append")%v2022a or later
                % end
                increment(hpw);
            end
            delete(hpw)  
        end
%% 
%--------------------------------------------------------------------------
% Plots of model output and utilities
%-------------------------------------------------------------------------- 
        function coefficientsPlot(obj)
            %generate data to plot the coefficents as a function of
            %direction, period and water level
            spobj = ctWaveSpectrum;
            spobj = setSpectrumModel(spobj);      %define the model to be used (Jonswap etc)            
            if isempty(spobj.spModel),return; end
            %force selection of source to be waves
            % if strcmp(spobj.spModel.input,'Wind'), spobj.spModel.input = 'Wave'; end
     
            T = obj.Data.Offshore.Dimensions.Period;
            zwl = obj.Data.Offshore.Dimensions.WaterLevel;   
            Diri = obj.Data.Offshore.RowNames;  %inshore ray directions
            ndir = length(Diri);
            nper = length(T);              
            nwls = length(zwl);
            spobj.inpData = getloopinput(obj,Diri,T,zwl,1,1,1); %dummy inpData
            hpw = PoolWaitbar(ndir, 'Processing transfer tables');
            kw = zeros(ndir,nper,nwls); kt2 = kw; ktp = kw; kd = kw; gm = kw;
            parfor i=1:ndir                     %parfor loop  
                for j=1:nper
                    for k=1:nwls
                        %inputs
                        offobj = copy(spobj);                        
                        offobj.inpData = getloopinput(obj,Diri,T,zwl,i,j,k);
                        %offshore spectrum and wave parameters
                        offobj = getModelSpectrum(offobj);
                        if spobj.spModel.gamma==0
                            gm(i,j,k) = offobj.inpData.gamma;
                        end
                        offobj.Params = wave_spectrum_params(offobj); 
                        %inshore spectrum and wave parameters
                        inobj = get_inshore_spectrum(obj,offobj);                     
                        inobj.Params = wave_spectrum_params(inobj);
                        %resultant transfer coefficients
                        outable = get_transfer_coefficients(offobj,inobj);
                        kw(i,j,k) = outable.kw;
                        kt2(i,j,k) = outable.kt2;
                        ktp(i,j,k) = outable.ktp;
                        kd(i,j,k) = outable.kd;
                        increment(hpw);
                    end
                end
            end
            delete(hpw)
            output = struct('kw',kw,'kt2',kt2,'ktp',ktp,'kd',kd);
            if spobj.spModel.gamma==0
                spobj.spModel.gamma = mean(gm,'all');
            end
            get_coefficientsPlot(obj,Diri,T,zwl,output,spobj);
        end

%%
        function input = getloopinput(~,Diri,T,zwl,i,j,k)
            %define input from arrays for use in parfor loop
            input.Hs = 1.0;       %transfer coefficients for unit wave height
            input.Dir = Diri(i);  %limit examination of mean offshore
            input.Tp = T(j);      %directions to inshore range
            input.swl = zwl(k);
            input.input = 'Wave';
            input.output = 'SpectralTransfer';
            input.date = 'now';
        end    

%%
        function depi = inshoreDepths(obj,spobj,depi,G,swl)   
            %determine the depths to use if the TMA spectrum or depth 
            %saturation is being used (eg in get_inshore_spectrum)
            % obj - SpectralTransfer; spobj - ctWaveSpectrum; 
            % depi - inshore depth; G - ; swl - still water level
            offdst = obj.Data.Offshore;         %offshore properties
            %check limits for valid rays based on depth and shoreline angle  
            %parfor not finding dynamic property, so use table
            idx = offdst.DataTable.depth<=0 | isnan(offdst.DataTable.depth);            
            T = offdst.Dimensions.Period;          %wave periods used in ray model
            fray = 1./T;                           %wave frequencies used in ray model
        
            hmn = offdst.DataTable.mindepth;       %minimum depth along ray      
            hmn(idx) = 0; 


            Tdef = [2,3,4,5,6,7,8,9,10,11.8,13.3,15.4,18.2,22.2,28.6,40];
            %pad the high frequencies with values from the maximum frequency
            idx = find(Tdef<min(T));
            if ~isempty(idx)
                addfray = 1./Tdef(idx);
                nadd = numel(addfray);
                fray = [addfray';fray];  %pad wave ray frequencies to maximum (1/2s)
                hmn = [repmat(hmn(:,1,:),1,nadd),hmn];
            end

            %pad the high frequencies with values from the minimum frequency
            idx = find(Tdef>max(T));
            addfray = 1./Tdef(idx);
            nadd = numel(addfray);
            if ~isempty(idx)
                fray = [fray;addfray'];  %pad wave ray frequencies to minimum (1/40s)
                hmn = [hmn,repmat(hmn(:,1,:),1,nadd)];

            end

            % if min(T)>3
            %     addfray = [1,0.5,0.33];
            %     fray = [addfray';fray];  %pad wave ray frequencies for periods 1-3s
            %     hmn = [repmat(hmn(:,1,:),1,3),hmn];
            % end
        
            %replicate grid vectors to produce grids for interpn using
            %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
            [PFW,XoFoWo] = transferDims(obj,spobj,fray,swl);
            hGmn = interpn(PFW{:},hmn,XoFoWo{:},'linear',0);  
            %for saturation limit use site or dominant rays min depth
            [~,freq] = spectrumDimensions(spobj);
            radint = deg2rad(spobj.Interp.dir);
            hGmn = trapz(radint,abs(trapz(freq,G.*hGmn,2)));%minimum depth
            depi = min([depi,min(hGmn,[],'All','omitnan')]);
        end

%%
        function [PFW,XoFoWo,fw,fowo] = transferDims(obj,spobj,fray,swl)
            %replicate grid vectors to produce grids for interpn using
            %variable number of dimensions (offshore 2 or 3, inshore 1 or 2)
            [dir,freq] = spectrumDimensions(spobj); 
            offdst = obj.Data.Offshore;            %offshore properties
            % T = offdst.Dimensions.Period;          %wave periods used in ray model
            % fray = 1./T;                           %wave frequencies used in ray model
            zwl = offdst.Dimensions.WaterLevel;    %water levels used in ray model
            InDir = offdst.RowNames;  %inshore directions are held in the offshore table 

            if isscalar(zwl)
                [P,F] = ndgrid(InDir,fray);
                PFW = {P,F};  
                fw = {fray};
                %frequency, direction and water level arrays to interpolate to 
                %var(xso,fro)
                [Xso,Fro] = ndgrid(dir,freq);  
                XoFoWo = {Xso,Fro};  
                fowo = {freq};
            else
                [P,F,W]  = ndgrid(InDir,fray,zwl);
                PFW = {P,F,W}; 
                fw = {fray,zwl};
                %frequency, direction and water level arrays to interpolate to 
                %var(xso,fro) for selected water level swl
                [Xso,Fro,Wln] = ndgrid(dir,freq,swl);  
                XoFoWo = {Xso,Fro,Wln}; 
                fowo = {freq,swl};
            end
        end
    
%%
        function tabPlot(obj,src,mobj) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab

            answer = questdlg('Inshore or Offshore results?','Spectrans',...
                                         'Inshore','Offshore','Offshore');
            if strcmp(answer,'Inshore')
                inshore_tt_plot(obj,src,mobj);
                return;
            end

            dst = obj.Data.Offshore;

            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src,mobj);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end

            %extract required data
            options = get_selection(obj);
            T = dst.Dimensions.Period;
            zwl = dst.Dimensions.WaterLevel;
            phi = dst.RowNames;

            %construct Q-Plot      
            if isscalar(T) || isscalar(phi)
                var = dst.(options.var)(:,:,:);
                scalar_tt_plot(obj,ax,T,zwl,phi,var,options); 
            else
                var = dst.(options.var)(:,:,options.ki);
                surface_tt_plot(obj,ax,T,zwl,phi,var,options); 
            end
        end

%%
        function ok = checkWLrange(obj,swl)
            %check that water level input conditions are within the range
            %of the Transfer Table
            zwl = minmax(obj.Data.Offshore.Dimensions.WaterLevel);
            ok = any(swl<zwl(1) | swl>zwl(2));
        end

%%
        function surf_plot(~,T,phi,var,varname,ax)
            %surface plot used in coefficients_plot and wrm_animation
            if nargin<6
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end
            surf(ax,T,phi,var,'Tag','PlotFigSurface');
            view(2);
            shading interp
            axis tight
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = varname{1};
            xlabel(varname{2}); 
            ylabel(varname{3}); 
            if length(varname)>3
                cb.Tag = varname{4};
            else
                cb.Tag = varname{1};
            end
        end

%%
        function polar_plot(~,Period,Phi,var,varname,~)
            %surface polar plot used in wrm_animation
            if nargin<6
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                axes(hf);
            end
            wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
            radrange = [0,25];
            %interpolate var(phi,T) onto plot domain defined by tints,rints
            %NB this should match properties of ctWaveSpectrum
            tints = linspace(0,2*pi,360); %theta intervals
            rints = linspace(1,20,20);    %radual intervals
            [Tq,Rq] = meshgrid(tints,rints); 
            warning('off',wid)
            vq = griddata(deg2rad(Phi),Period,var',Tq,Rq);
            vq(isnan(vq)) = 0;  %fill blank sector so that it plots the period labels
            warning('on',wid)
            color = mcolor('dark blue');
            [X,Y] = polarplot3d(vq,'plottype','surfn','TickSpacing',45,...
                'RadLabels',4,'RadLabelLocation',{20 'top'},...
                'GridColor',color,'TickColor',color,...
                'RadLabelColor',mcolor('dark grey'),...
                'RadialRange',radrange,'polardirection','cw');
            view(2)
            shading interp 
            axis(gca,'off')
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = varname{1};
            text(max(X,[],'all')/2,max(Y,[],'all')*0.95,0,varname{2});

            if length(varname)>3
                cb(1).Tag = varname{4};
            else
                cb(1).Tag = varname{1};
            end
        end

    end
%%    
    methods (Access = private)
        function [spectran,indir,inprops] = specTransfer(~,rayobj)
            %generate the spectral transfer table of offshore directions,
            %depths and celerities
            dst = rayobj.Data.Dataset;
            
            indir = dst.RowNames;  %assign inshore ray directions
            ndir = length(indir);
            nper = length(dst.Dimensions.Period);
            nwls = length(dst.Dimensions.WaterLevel);

            %use inshore point of first ray to define inshore point properties
            %for each period and water level combination
            xi = dst.xr{1,1,1}(1);
            yi = dst.yr{1,1,1}(1);
            ci = zeros(1,nper,nwls); hi = ci; cgi = ci;
            for j=1:nper
                for k=1:nwls     
                    hi(1,j,k) = dst.depth{1,j,k}(1);
                    ci(1,j,k) = dst.celerity{1,j,k}(1);
                    cgi(1,j,k) = dst.cgroup{1,j,k}(1);  
                end
            end         
            hi = squeeze(hi(1,1,:));
            inprops = {xi,yi,hi,ci,cgi};  %inshore properties
            %compile table data for offshore properties
            
            c = zeros(ndir,nper,nwls); cg = c; h = c; offdir = c; 
            hav = c; hmn = c;
            for i=1:ndir
                for j=1:nper
                    for k=1:nwls
                        flag = dst.UserData.flag(i,j,k)>0;
                        if flag>0
                            theta = dst.alpha{i,j,k}(end);
                            offdir(i,j,k) = mod(compass2trig(theta,true),360);
                        else
                            offdir(i,j,k) = NaN;
                        end
                        h(i,j,k) = dst.depth{i,j,k}(end).*flag;
                        hav(i,j,k) = mean(dst.depth{i,j,k}).*flag;
                        hmn(i,j,k) = min(dst.depth{i,j,k}).*flag;
                        c(i,j,k) = dst.celerity{i,j,k}(end).*flag;
                        cg(i,j,k) = dst.cgroup{i,j,k}(end).*flag;
                    end
                end
            end
            spectran = {offdir,h,c,cg,hav,hmn};
        end

%%       
        function options = get_selection(obj)
            %get index of period, water level and variable to use in plots or model
            %   Defined using varargin for the following fields
            %    FigureTitle     - title for the UI figure
            %    PromptText      - text to guide user on selection to make
            %    InputFields     - text prompt for input fields to be displayed
            %    Style           - uicontrols for each input field (same no. as input fields)
            %    ControlButtons  - text for buttons to edit or update selection 
            %    DefaultInputs   - default text or selection lists
            %    UserData        - data assigned to UserData of uicontrol
            %    DataObject      - data object to use for selection
            %    SelectedVar     - index vector to define case,dataset,variable selection  
            %    ActionButtons   - text for buttons to take action based on selection
            %    Position        - poosition and size of figure (normalized units)
            vardesc = obj.Data.Offshore.VariableDescriptions;
            zwl = obj.Data.Offshore.Dimensions.WaterLevel;
            selection = inputgui('FigureTitle','Levels',...
                                 'InputFields',{'Variable','Water level'},...
                                 'Style',{'popupmenu','popupmenu'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{vardesc,string(zwl)},...
                                 'PromptText','Select values to use');
            if isempty(selection)
                options = []; 
            else
                options.var = obj.Data.Offshore.VariableNames{selection{1}};
                options.desc = obj.Data.Offshore.VariableDescriptions{selection{1}};
                options.lab = obj.Data.Offshore.VariableLabels{selection{1}};
                options.ki = selection{2};
            end  
        end
%%
        function inshore_tt_plot(obj,src,mobj)
            %plot the inshore spectral transfer results
            dst = obj.Data.Inshore;

            %set >Figure button and create axes
            if strcmp(src.Tag,'Plot') || strcmp(src.Tag,'FigButton')
                tabcb  = @(src,evdat)tabPlot(obj,src,mobj);            
                ax = tabfigureplot(obj,src,tabcb,false);
                ax.NextPlot = 'add';
            else
                ax = src; %user passing an axis as src rather than a uicontrol
            end

            answer = questdlg('Celerity of Group Celerity results?','Spectrans',...
                                    'Celerity','Group celerity','Celerity');
            if strcmp(answer,'Celerity')
                var = squeeze(dst.celerity);
            else
                var = squeeze(dst.cgroup);
            end
            T = dst.Dimensions.Period;
            zwl = dst.Dimensions.WaterLevel;
            if isscalar(zwl) && isscalar(T)
                if isgraphics(ax.Parent,'figure')
                    delete(ax.Parent)
                end
                msgbox(sprintf('%s %0.2f (m/s)',answer,var));
                return;
            elseif isscalar(zwl) 
                plot(ax,T,var);
                ylabel(sprintf('%s (m/s)',answer)); 
                xlabel('Wave period (s)'); 
            elseif isscalar(T)
                plot(ax,zwl,var);
                ylabel(sprintf('%s (m/s)',answer)); 
                xlabel('Water Level (mOD)'); 
            else
                surf(ax,T,zwl,var');
                view(2);
                shading interp
                %add the colorbar and labels
                cb = colorbar;
                cb.Label.String = sprintf('%s (m/s)',answer);    
                xlabel('Wave period (s)'); 
                ylabel('Water Level (mOD)');                
            end            

            title(sprintf('%s for Case: %s',answer,dst.Description));
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot            
        end

%%
        function scalar_tt_plot(obj,ax,T,zwl,phi,var,options)
            %offshore transfer table data plot when only a single case such
            %that T, zwl or phi are scalar
            dst = obj.Data.Offshore;
            if isscalar(zwl) && isscalar(T) && isscalar(phi)
                if isgraphics(ax.Parent,'figure')
                    delete(ax.Parent)
                end
                msgbox(sprintf('%s %0.2f',options.desc,var));
                return;
            elseif isscalar(zwl) && isscalar(T) 
                plot(ax,phi,var);
                ylabel(options.lab); 
                xlabel('Direction (degTN)'); 
            elseif isscalar(T) && isscalar(phi) 
                plot(ax,zwl,var);
                ylabel(options.lab); 
                xlabel('Water Level (mOD)'); 
            elseif isscalar(zwl) && isscalar(phi)    
                plot(ax,T,var);
                ylabel(options.lab); 
                xlabel('Wave period (s)');
%             elseif isscalar(zwl) 
%                 surface_tt_plot(obj,ax,T,zwl,phi,var,options)             
            elseif isscalar(T) 
                var = squeeze(dst.(options.var)(:,:,:));
                surf(ax,zwl,phi,var); 
                view(2);
                shading interp                
                cb = colorbar;
                cb.Label.String = options.lab;    
                xlabel('Water Level (mOD)'); 
                ylabel('Inshore Direction (degTN)');
            elseif isscalar(phi)      
                var = squeeze(dst.(options.var)(:,:,:));
                surf(ax,T,zwl,var);
                view(2);
                shading interp
                cb = colorbar;
                cb.Label.String = options.lab; 
                xlabel('Wave period (s)'); 
                ylabel('Water Level (mOD)');                
            end

            title(sprintf('%s for Case: %s',options.desc,dst.Description));
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot  
        end
%%
        function surface_tt_plot(obj,ax,T,zwl,phi,var,options)
            %offshore transfer table data plot for selected variable
            %against period and direction
            dst = obj.Data.Offshore;
            % if sum(isnan(var),'all')/numel(var)<0.3
            %     %direction has NaN value when no offshore result
            %     %to infill the holes can use inpaint but maybe misleading
            %     var = inpaint_nans(var,4); %infill NaNs if < 40%
            % end
            surf(ax,T,phi,var);
            view(2);
            shading interp

            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = options.desc;
            xlabel('Wave period (s)'); 
            ylabel('Inshore direction (degTN)'); 
            title(sprintf('%s for water level of %.2g mOD',dst.Description,zwl(options.ki)));
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot

            if strcmp(options.var,'theta')
                %when variable is theta ???????
                mindir = min(var,[],'All');
                mindir = mindir-mod(mindir,30)+30;
                maxdir = max(var,[],'All');
                maxdir = maxdir-mod(maxdir,30);
                nc = mindir:30:maxdir;
                hold(ax,'on')
                [C,h] = contour3(ax,T,phi,var,nc,'-k');
                clabel(C,h,'LabelSpacing',300,'FontSize',8)
                hold(ax,'off')
            end
        end
%%
        function get_coefficientsPlot(obj,Dir,T,zwl,output,spobj)
            %interactive selection to plot the coefficients for range of
            %directions, periods and water levels
            ki = 1;
            if length(zwl)>1
                while ~isempty(ki)
                    ki = inputgui('FigureTitle','Levels',...
                                         'InputFields',{'Water level'},...
                                         'Style',{'popupmenu','popupmenu'},...
                                         'ActionButtons', {'Select','Cancel'},...
                                         'DefaultInputs',{string(zwl)},...
                                         'PromptText','Select values to use');
                    if isempty(ki), return; else, ki = ki{1}; end
                    coefficients_plot(obj,Dir,T,zwl,output,spobj,ki);
                end
            else
                coefficients_plot(obj,Dir,T,zwl,output,spobj,ki);
            end
        end
%%
        function coefficients_plot(obj,Dir,T,zwl,output,spobj,ki)
            %plot the coefficients for range of directions, periods and
            %water levels
            figure('Name','SpecTrans','Tag','PlotFig');
            labelx = 'Wave Period (s)';
            labely = 'Direction (degTN)';
            s1 = subplot(2,2,1);
            var = output.kw(:,:,ki);
            surf_plot(obj,T,Dir,var,{'Transfer coefficient, kw',labelx,labely},s1);
            s2 = subplot(2,2,2);
            var = output.kt2(:,:,ki);
            surf_plot(obj,T,Dir,var,{'Transfer coefficient, kt2',labelx,labely},s2);
            s3 = subplot(2,2,3);
            var = output.ktp(:,:,ki);
            surf_plot(obj,T,Dir,var,{'Transfer coefficient, ktp',labelx,labely},s3);
            s4 = subplot(2,2,4);
            var = output.kd(:,:,ki);
            surf_plot(obj,T,Dir,var,{'Direction shift (deg)',labelx,labely},s4);
            %plot title
            desc = obj.Data.Inshore.Description;
            sg1 = sprintf('Transfer Coefficients(Tp,Dir) for %s at swl=%g mOD',desc,zwl(ki));
            exptxt = 'Coefficients: kw-wave height; ktp-peak period; kt2-mean period';
            spobj = setModelInputText(spobj);
            sgtxt = sprintf('%s\n%s\n%s',sg1,spobj.Plotxt.stxt,exptxt);
            sgtitle(sgtxt,'FontSize',11,'Margin',1);
        end

%%
function dsp = modelDSproperties(~,isin) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            if isin
                dsp.Variables = struct(...  
                    'Name',{'celerity','cgroup'},...
                    'Description',{'Celerity','Group Celerity'},...
                    'Unit',{'m/s','m/s'},...
                    'Label',{'Celerity (m/s)','Group Celerity (m/s)'},...                           
                    'QCflag',repmat({'model'},1,2)); 
                dsp.Row = struct(...
                    'Name',{'-'},...
                    'Description',{'-'},...
                    'Unit',{'-'},...
                    'Label',{'-'},...
                    'Format',{'-'}); 
                dsp.Dimensions = struct(...    
                    'Name',{'Period','WaterLevel'},...
                    'Description',{'Wave Period','Water Level'},...
                    'Unit',{'m','mOD'},...
                    'Label',{'Wave Period','Water Level'},...
                    'Format',{'-','-'});                  
            else
                dsp.Variables = struct(...  
                    'Name',{'theta','depth','celerity','cgroup','avdepth','mindepth'},...
                    'Description',{'Offshore Direction','Offshore depth',...
                                   'Celerity','Group Celerity',...
                                   'Average depth','Minimum depth'},...
                    'Unit',{'degTN','m','m/s','m/s','m','m'},...
                    'Label',{'Direction (degTN)','Water depth (m)',...
                             'Celerity (m/s)','Group Celerity (m/s)',...
                             'Water depth (m)','Water depth (m)'},...                           
                    'QCflag',repmat({'model'},1,6)); 
                dsp.Row = struct(...
                    'Name',{'InDir'},...
                    'Description',{'Inshore Direction'},...
                    'Unit',{'degTN'},...
                    'Label',{'Direction (degTN)'},...
                    'Format',{'-'}); 
                dsp.Dimensions = struct(...    
                    'Name',{'Period','WaterLevel'},...
                    'Description',{'Wave Period','Water Level'},...
                    'Unit',{'m','mOD'},...
                    'Label',{'Wave Period','Water Level'},...
                    'Format',{'-','-'});  
            end
        end
    end    
end