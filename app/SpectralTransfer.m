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
        interp
    end
    
    methods (Access={?muiDataSet,?muiStats,?WRM_WaveModel})
        function obj = SpectralTransfer()             
            %class constructor
            obj.interp.dir = 0.5; %interval used to interpolate directions (deg)
            obj.interp.per = 0.5; %interval used to interpolate periods (s)
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model to construct spectral transfer table
%--------------------------------------------------------------------------   
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = SpectralTransfer;                           
            muicat = mobj.Cases;    
%--------------------------------------------------------------------------
% Model code>
%--------------------------------------------------------------------------
            %select back tracking ray case to use
            promptxt = 'Select Backward Ray Trace Dataset:';
            rayobj = selectCaseObj(muicat,{'backward_model'},[],promptxt);                                                    
            if isempty(rayobj), return; end
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
            setDataSetRecord(obj,muicat,dst,'spectral_model');
            getdialog('Run complete');
        end    
    end
%%
    methods
%--------------------------------------------------------------------------
% Model to transfer a wave timeseries or spectral data set
%--------------------------------------------------------------------------         
        function [Sot,Sit,Dims,output] = runWaves(obj,tsdst,select)
            %run the spectral transfer model for a timeseries of offshore
            %wave conditions and return a table of wave conditions
            Sot = []; Sit = []; Dims = [];
            islog = false;  filename = [];
            if select.issave
                answer = questdlg('Write log of missing dates to file?','Rays','Yes','No','No');
                if strcmp(answer,'Yes'), islog = true; end
                filename = sprintf('Sprectra_log_%s.txt',char(datetime,"ddMMMyy_HH-mm"));
            end

            ndir = 360/obj.interp.dir;         %number of direction intervals
            nper = length(1:obj.interp.per:30);%number of period intervals
            nint = height(tsdst);
            output = table();
            inputable = tsdst.DataTable;       %needed for parallel option
            hpw = PoolWaitbar(nint, 'Processing timeseries');
            blank = zeros(1,ndir,nper);           %fixed intervals assigned in get_inshore_spectrum
            parfor i=1:nint                    %parfor loop  
                %for each offshore wave get the inshore results
                input = inputable(i,:);
                [SGo,SGi,dims] = get_inshore_spectrum(obj,input,select);
                output(i,:) = get_inshore_wave(SGo,SGi,dims,input,select);
                if select.issave
                    if isempty(SGo)
                        Sot(i,:,:) = blank; Sit(i,:,:) = blank;
                        if islog
                            lines = sprintf('%s',tsdst.RowNames(i)); %#ok<PFBNS> 
                            writelines(lines,filename,WriteMode="append")
                        end
                        continue; 
                    end  
                    Sot(i,:,:) = SGo;
                    Sit(i,:,:) = SGi;
                    Fri(i,:) = dims.freq;
                    Dir(i,:) = dims.dir;
                    Depi(i) = dims.depi;    
                end
                increment(hpw);
            end

            if select.issave
                [Sot,Sit,Dims] = packSpectra(obj,Sot,Sit,Fri,Dir,Depi);
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
            select = get_model_selection(obj);  %select spectral form and data type
            if isempty(select), return; end     %user cancelled
            %force selection of source to be waves
            if strcmp(select.source,'Wind'), select.source='Wave'; end
     
            T = obj.Data.Offshore.Dimensions.Period;
            zwl = obj.Data.Offshore.Dimensions.WaterLevel;   
            Diri = obj.Data.Offshore.RowNames;  %inshore ray directions
            ndir = length(Diri);
            nper = length(T);              
            nwls = length(zwl);

            kw = zeros(ndir,nper,nwls); kt2 = kw; ktp = kw; kd = kw;
            parfor i=1:ndir                     %parfor loop  
                for j=1:nper
                    for k=1:nwls
                        input = getloopinput(obj,Diri,T,zwl,i,j,k);
                        [SGo,SGi,Dims] = get_inshore_spectrum(obj,input,select);                                                   
                        outable = get_inshore_wave(SGo,SGi,Dims,input,select);
                        kw(i,j,k) = outable.kw;
                        kt2(i,j,k) = outable.kt2;
                        ktp(i,j,k) = outable.ktp;
                        kd(i,j,k) = outable.kd;
                    end
                end
            end
            output = struct('kw',kw,'kt2',kt2,'ktp',ktp,'kd',kd);
        
            get_coefficientsPlot(obj,Diri,T,zwl,output,select);
        end
%%
        function input = getloopinput(~,Diri,T,zwl,i,j,k)
            %define input from arrays for use in parfor loop
            input.Hs = 1.0;       %transfer coefficients for unit wave height
            input.Dir = Diri(i);  %limit examination of mean offshore
            input.Tp = T(j);      %directions to inshore range
            input.swl = zwl(k);
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
            %check that water leven input conditions are within the range
            %of the Transfer Table
            zwl = minmax(obj.Data.Offshore.Dimensions.WaterLevel);
            ok = any(swl<zwl(1) | swl>zwl(2));
        end
%%
        function  getSpectrumPlot(obj,SGo,SGi,Dims,ins,off,sel)
            %plot an offshore and inshore spectrum for a single wave height
            hfig = figure('Name','O/I Spectrum', 'Units','normalized', ...
                            'Resize','on','HandleVisibility','on','Visible','off',...
                            'Position',[0.38,0.42,0.30,0.42],'Tag','PlotFig');
            figax = axes(hfig);

            answer = questdlg('What type of plot','O?I spectrum','XY','Polar','XY');
            if strcmp(answer,'XY')
                figtype = true;          %Cartesian dir-freq plot
            else
                figtype = false;         %Polar dir-freq plot     
            end
            hfig.Visible = 'on';
            %add offshore and inshore plots of spectra
            [s1,s2] = off_in_plot(obj,1./Dims.freq,Dims.dir,SGo,SGi,figax,figtype);
            
            %offshore titles are source dependent
            if strcmp(sel.source,'Wave')
                sgtxt = sprintf('%s, gamma=%.2g, and %s, n=%d ',sel.form,...
                                        sel.gamma,sel.spread,sel.nspread);
                st1 = title(s1,sprintf('Hso=%.2f m; Tp=%.1f s; Dir=%.3g degTN; swl=%.2f mOD',...
                                off.Hs,off.Tp,off.Dir,off.swl),'Margin',1);             
            elseif strcmp(sel.source,'Wind')
                sgtxt = sprintf('%s, gamma=%.2g, and %s, n=%d ',sel.form,...
                                        sel.gamma,sel.spread,sel.nspread);
                st1 = title(s1,sprintf('Uw=%.2f m/s; Fetch=%0.0f m; swl=%.2f mOD\nHso=%.2f m; Tp=%.1f s; Dir=%.3g degTN',...
                      off.AvSpeed,off.Fetch,off.swl,off.Hs,off.Tp,off.Dir),'Margin',1);
            else
                ofd = off.tsdst;    
                sgtxt = sprintf('Measured spectrum at %s on %s',...
                                   ofd.Description,string(ofd.RowNames));
                p = wave_spectrum_params(SGo,Dims.freq,Dims.dir);
                st1 = title(s1,sprintf('Hso=%.2f m; Tz=%.1f s; Dir=%.3g degTN; swl=%.2f mOD',...
                                ofd.Hs,ofd.Tz,p.Dir0,off.swl),'Margin',1);      
            end
            
            %inshore titles
            sgtitle(sgtxt,'FontSize',12,'Margin',1);
            st2 = title(s2,sprintf('Hsi=%.2f m; Tp=%.1f s; Dir=%.3g degTN; hmin=%.2f m',...
                            ins.Hsi,ins.Tpi,ins.Diri,ins.depi),'Margin',1);   

            %for polar plot adjust figure size and subplot title position
            if ~figtype
                st1.Position(2) = -st1.Position(2)*1.2;
                st2.Position(2) = -st2.Position(2)*1.2;
                hfig.Position(3) = 0.6;
            end
        end
%%
        function [s1,s2] = off_in_plot(obj,T,dir,var0,vari,ax,isXY)
            %plot offshore and inshore spectra
            if nargin<6
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                ax = axes(hf);
            end
            labeli = 'Inshore direction (degTN)';
            label0 = 'Offshore direction (degTN)';
            labelx = 'Wave period (s)';
            shorenorm = obj.Data.Inshore.UserData.ShoreAngle+90;         
            grey = mcolor('light grey');

            if isXY
                s1 = subplot(2,1,1,ax);
                check_plot(obj,T,dir,var0,{'Spectral Energy (m^2s)',labelx,label0,'Off'},s1)
                s2 = subplot(2,1,2);
                check_plot(obj,T,dir,vari,{'Spectral Energy (m^2s)',labelx,labeli,'In'},s2)
                hold on
                    xn = minmax(T); yn = [shorenorm,shorenorm]; 
                    plot(s2,xn,yn,'Color',grey,'LineStyle','--','LineWidth',1);
                hold off
                s2.YLim = s1.YLim;
            else
                s1 = subplot(1,2,1,ax);
                polar_plot(obj,T,dir,var0,{'Spectral Energy (m^2s)',labelx,label0,'Off'},s1)
                s2 = subplot(1,2,2);
                polar_plot(obj,T,dir,vari,{'Spectral Energy (m^2s)',labelx,labeli,'In'},s2)
                hold on
                    ang = compass2trig(shorenorm);
                    xn = [0,max(T)*cos(ang)]; yn = [0,max(T)*sin(ang)];
                    plot(s2,xn,yn,'Color',grey,'LineStyle','--','LineWidth',1);
                hold off                
                s2.YLim = s1.YLim;
            end
        end        
%%
        function spectra = get_model_selection(~,source)
            %get spectrum form, data source, and parameters for wave
            %spectrum and directions spreading functions
            %  Defined using varargin as in above function   
            %  source iswind is used to prioritise selection of wind            
            sp = {'JONSWAP fetch limited','TMA shallow water',...
                   'Pierson-Moskowitz fully developed','Bretschneider open ocean'};
            src = {'Wave','Wind'};
            if nargin>1 && strcmp(source,'Wind')
                src = {'Wind','Wave'};
            end
            
            spr = {'SPM cosine function','Donelan secant function'};
            nsp = string(0:1:10);
            selection = inputgui('FigureTitle','Spectrum',...
                                 'InputFields',{'Wave Spectra','Data type',...
                                    'Spread function','Spread exponent',...
                                    'Jonswap gamma (if req.)'},...
                                 'Style',{'popupmenu','popupmenu',...
                                         'popupmenu','popupmenu','edit'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{sp,src,spr,nsp,'3.3'},...%use nspread=0 if included in wave data
                                 'PromptText','Select values to use');
            if isempty(selection)
                spectra = [];
            else
                spectra.form = sp{selection{1}};
                spectra.source = src{selection{2}};
                spectra.spread = spr{selection{3}};
                spectra.nspread = str2double(nsp{selection{4}});
                spectra.gamma = str2double(selection{5});
                spectra.ismodel = true;
                spectra.issat = false;
                if strcmp(spectra.form,'TMA shallow water')
                    spectra.issat = true;
                end
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
        function [Sot,Sit,Dims] = packSpectra(~,Sot,Sit,Fri,Dir,Depi)
            %pack spectral timeseries to minimise size of arrays        
            idx = find(Fri(:,1)~=0,1,'first'); %find first non-null result
            Dims.freq = squeeze(Fri(idx,:));      %dimensions used for run
            Dims.dir = squeeze(Dir(idx,:));
            Dims.depi = Depi(idx);
        
            %minimise array size - remove directions and frequencies that
            %are not contributing to spectra (S<0.05 m^2/Hz)
            tol = max(Sot,[],'all')/1000;
            [~,idro,idco] = compact3Darray(Sot,1,tol); %pivot dim is time
            [~,idri,idci] = compact3Darray(Sit,1,tol);
            idr = minmax([idro;idri]);
            idc = minmax([idco;idci]);
            Sot = Sot(:,idr(1):idr(2),idc(1):idc(2));
            Sit = Sit(:,idr(1):idr(2),idc(1):idc(2));
            Dims.dir = Dims.dir(idr(1):idr(2));
            Dims.freq = Dims.freq(idc(1):idc(2));    
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
            %offshore transfer table data plot when only a single case such that 
            %T, zwl or phi are scalar
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
            %offhsore transfer table data plot when only a single case such that 
            %T, zwl or phi are scalar
            dst = obj.Data.Offshore;
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
                mindir = min(var,[],'All');
                mindir = mindir-mod(mindir,30)+30;
                maxdir = max(var,[],'All');
                maxdir = maxdir-mod(maxdir,30);
                nc = mindir:30:maxdir;
                hold on
                [C,h] = contour3(ax,T,phi,var,nc,'-k');
                clabel(C,h,'LabelSpacing',300,'FontSize',8)
                hold off
            end
        end
%%
        function get_coefficientsPlot(obj,Dir,T,zwl,output,sel)
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
                    coefficients_plot(obj,Dir,T,zwl,output,sel,ki);
                end
            else
                coefficients_plot(obj,Dir,T,zwl,output,sel,ki);
            end
        end
%%
        function coefficients_plot(obj,Dir,T,zwl,output,sel,ki)
            %plot the coefficients for range of directions, periods and
            %water levels
            figure('Name','SpecTrans','Tag','PlotFig');
            labelx = 'Wave Period (s)';
            labely = 'Direction (degTN)';
            s1 = subplot(2,2,1);
            var = output.kw(:,:,ki);
            check_plot(obj,T,Dir,var,{'Transfer coefficient, kw',labelx,labely},s1);
            s2 = subplot(2,2,2);
            var = output.kt2(:,:,ki);
            check_plot(obj,T,Dir,var,{'Transfer coefficient, kt2',labelx,labely},s2);
            s3 = subplot(2,2,3);
            var = output.ktp(:,:,ki);
            check_plot(obj,T,Dir,var,{'Transfer coefficient, ktp',labelx,labely},s3);
            s4 = subplot(2,2,4);
            var = output.kd(:,:,ki);
            check_plot(obj,T,Dir,var,{'Direction shift (deg)',labelx,labely},s4);
            sg1 = sprintf('Transfer Coefficients(Mean Direction, Period) for swl=%g mOD',zwl(ki));
            if contains({'JONSWAP fetch limited','TMA shallow water'},sel.form)
                sgtxt = sprintf('%s\n%s, gamma=%.2g, and %s, n=%d',sg1,...
                                sel.form,sel.gamma,sel.spread,sel.nspread);
            else
                sgtxt = sprintf('%s\n%s, and %s, n=%d',sg1,...
                                sel.form,sel.spread,sel.nspread);
            end
            sgtitle(sgtxt,'FontSize',12,'Margin',1);
        end
%%
        function check_plot(~,T,phi,var,varname,ax)
            %check plot for data selection
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
            %check plot for data selection
            if nargin<6
                hf = figure('Name','SpecTrans','Tag','PlotFig');
                axes(hf);
            end
            wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
            radrange = [0,25];
            %interpolate var(phi,T) onto plot domain defined by tints,rints
            tints = linspace(0,2*pi,360);  
            rints = linspace(1,25,25);
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