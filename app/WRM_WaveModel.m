classdef WRM_WaveModel < muiDataSet 
%
%-------class help------------------------------------------------------
% NAME
%   WRM_WaveModel.m
% PURPOSE
%    Class for wave refraction using backward ray transfer function.
%    Constructs inshore time series from an offshore timeseries. Also
%    includes methods to plot the offshore and inshore spectra.
% SEE ALSO
%   muiDataSet, WaveRayModel, RayTracks, SpectralTransfer
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2023
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
    end
    
    properties (Hidden)
        ModelType        %model used for the particular instance
    end
    
    methods (Access={?muiDataSet,?muiStats})
        function obj = WRM_WaveModel()                    
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run the wave refraction model.
            obj = WRM_WaveModel;                            
            [dspec,dsprop] = setDSproperties(obj);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            %Note: getInputData calls setRunParam
            [tsdst,sptobj,inputxt,source] = getInputData(obj,mobj);
            if isempty(tsdst), return; end   %user cancelled data selection

            if strcmp(source,'Measured spectra')
                select = WRM_WaveModel.getSprectraConditions();   
                select.freq = tsdst.Dimensions.freq;
                if isempty(select), return; end       %user cancelled
                ModelType = 'transfer spectra';
            else
                select = get_model_selection(sptobj); %select spectral form and data type
                if isempty(select), return; end       %user cancelled
                ModelType = 'transfer timeseries';
            end
            inputxt = sprintf('%s\n%s',inputxt,select.txt); %add spectral selection to meta data

            answer = questdlg('Save the spectra?','Wave model','Yes','No','No');
            if strcmp(answer,'Yes'), select.issave=true; else, select.issave=false; end
            [Sot,Sit,Dims,results] = runWaves(sptobj,tsdst,select);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            time = tsdst.RowNames;
            dst.Properties = dstable(results,'RowNames',time,'DSproperties',dsprop);                      
            %assign metadata about model            
            dst.Properties.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                         ModelType);
            dst.Properties.MetaData = inputxt;   
            %add depths of inshore point for which there are backward rays
            dst.Properties.UserData = sptobj.Data.Inshore.UserData.Depths;
            %check whether wave spectra should also be saved
            sze = 2*getfield(whos('Sot'),'bytes')*9.53674e-7;
            if select.issave
                questxt = sprintf('Save the full wave spectra (arrays are %.1f Mb)?',sze);
                answer = questdlg(questxt,'Wave model','In Full','Minimised','Minimised');
                if strcmp(answer,'Minimised')
                    %reduce to 1 degree spacing and 1 second periods
                    idx = rem(Dims.dir,1)>0;
                    idy = rem(1./Dims.freq,1)>0;
                    Dims.dir(idx) = [];
                    Dims.freq(idy) = [];
                    Sot(:,idx,:) = [];  Sot(:,:,idy) = [];
                    Sit(:,idx,:) = [];  Sit(:,:,idy) = [];
                end
                dst.Spectra = dstable(Sot,Sit,'RowNames',time,'DSproperties',dspec); 
                dst.Spectra.Dimensions.dir = Dims.dir;    %NB order is X,Y and must
                dst.Spectra.Dimensions.freq = Dims.freq;  %match variable dimensions           
                %assign metadata about model
                dst.Spectra.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                    ModelType);
                dst.Spectra.MetaData = inputxt;
            end
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            obj.ModelType = 'Inwave_model';               %added to match ctWaveModel
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');
        end
%--------------------------------------------------------------------------
% Other utilities
%--------------------------------------------------------------------------
%%
        function runSpectrum(mobj)
            %create a plot of the offshore and inshore 2-D specrum surfaces 
            %for a single wave condition
            offdata = WRM_WaveModel.getForcingConditions();   %get the input conditions
            if isempty(offdata); return; end

            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            sptobj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            isout = checkWLrange(sptobj,offdata.swl);
            if isout
                warndlg('Water levels are outside the range of the Transfer Table')
                return;
            end

            if strcmp(offdata.source,'Spectrum')               
                select.form = 'Measured';           %initialise select properties
                select.source = 'Spectrum';         
                select.ismodel = false;
                select.issat = offdata.issat;       %copy to sprectrum selection  
                select.freq = offdata.tsdst.Dimensions.freq;
                intable = offdata.tsdst.DataTable;
                [SGo,SGi,Dims] = get_inshore_spectrum(sptobj,intable,select);
                ins = get_inshore_wave(SGo,SGi,Dims,intable,select);
            else
                select = get_model_selection(sptobj,offdata.source); %select spectral form and data type
                if isempty(select), return; end           %user cancelled

                [SGo,SGi,Dims] = get_inshore_spectrum(sptobj,offdata,select);
                if isempty(SGo), return; end
    
                if strcmp(offdata.source,'Wind')
                    offdata = WRM_WaveModel.addWaveConditions(SGo,Dims,offdata);
                end
                ins = get_inshore_wave(SGo,SGi,Dims,offdata,select);
            end
            

            getSpectrumPlot(sptobj,SGo,SGi,Dims,ins,offdata,select);            
        end
%%
        function runAnimation(mobj)
            %create an animation of the 2-D spectrum surfaces using a
            %timeseries input
            obj = WRM_WaveModel; 
            [tsdst,~,source] = getInputData(obj,mobj);
            if isempty(tsdst), return; end   %user cancelled data selection  
            tsdst.DataTable = rmmissing(tsdst.DataTable);%remove nans

            if height(tsdst)>5000
                promptxt = sprintf('Times series contains %d records\nThis could take a while to run and genearte large file\nUse time sub-selection to extract shorter time period',...
                                                            height(tsdst));
                answer = questdlg(promptxt,'Time','Continue','Abort','Abort');                                  
                if strcmp(answer,'Abort'), return; end
            end

            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            sptobj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            isout = checkWLrange(sptobj,tsdst.swl);
            if isout
                warndlg('Water levels are outside the range of the Transfer Table')
                return;
            end

            if strcmp(source,'Measured spectra')
                select = WRM_WaveModel.getSprectraConditions();   
                select.freq = tsdst.Dimensions.freq;
                if isempty(select), return; end       %user cancelled
            else
                select = get_model_selection(sptobj); %select spectral form and data type
                if isempty(select), return; end       %user cancelled
            end
            
            select.issave = true;
            [SGo,SGi,Dims] = runWaves(sptobj,tsdst,select);
            if isempty(SGo), return; end

            wrm_animation(mobj,sptobj,tsdst,SGo,SGi,Dims)
        end
    end
%%
    methods
        function [tsdst,caserec] = getWaveModelDataset(obj,mobj,type,addnames,caserec)
            %prompt user to select a model wave dataset and add Tp if inshore
            muicat = mobj.Cases;
            if nargin<4
                addnames = {'Tp'};  %default is to add Tp
            end
            %
            if nargin<5   %no caserec to prompt for selection
                [wvobj,wvdst,ok] = selectClassInstance(obj,'ModelType',type);
                if ok<1, tsdst = []; return; end
                caserec = caseRec(muicat,wvobj.CaseIndex);
            else
                wvobj = getCase(muicat,caserec);
                wvdst = wvobj.Data.Dataset;
            end

            %if inshore wave dataset add variables requested otherwise just
            %copy dstable
            tsdst = copy(wvdst);
            if strcmp(type,'Inwave_model')
                inpwavecid = wvobj.RunParam.ctWaveData.caseid;
                inpwaverec = caseRec(muicat,inpwavecid);
                inpdst = getDataset(muicat,inpwaverec,1);
                dstnames = inpdst.VariableNames;
                varnames = tsdst.VariableNames;
                for i=1:length(addnames)                    
                    if any(strcmp(varnames,addnames{i}))
                        continue;
                    elseif any(contains(varnames,addnames{i})) && ...
                                      ~any(strcmp(varnames,addnames{i}))  
                        idx = contains(varnames,addnames{i});
                        tsdst.VariableNames{idx} = addnames{i};
                    elseif any(strcmp(dstnames,addnames{i})) && ...
                                            ~isempty(inpdst.(addnames{i}))
                        tsdst = addvars(tsdst,inpdst.(addnames{i}),...
                                           'NewVariableNames',addnames{i});
                    else
                        warndlg('Variable %s not found so not added to wave dataset',...
                                                           addnames{i});
                    end
                end
            end
        end
%%        
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        function [tsdst,sptobj,inputxt,source] = getInputData(obj,mobj)
            %prompt user to select wave and water level data and return in
            %input dstable of data and metadata for inputs used
            tsdst = []; sptobj = []; inputxt = []; source = [];
            muicat = mobj.Cases;
            promptxt = 'Select input wave data set:';           
            [wv_caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                           'CaseClass',{'ctWaveData','muiUserModel'},'ListSize',[300,100]);                                    
            if ok<1, return; end
            wvdst = getDataset(muicat,wv_caserec,1);  %1 selects first dataset in struct
            wvtime = wvdst.RowNames;
            inputxt = sprintf('%s used for offshore waves',wvdst.Description);

            source = 'Measured waves';
            if isfield(wvdst.Dimensions,'freq')
                source = 'Measured spectra';
            end
            
            varnames = wvdst.VariableNames;
            if strcmp(source,'Measured waves') && ~any(strcmp(varnames,'Hs'))
                wvdst = extract_wave_data(wvdst);
                if isempty(wvdst), return; end
            end

            promptxt = 'Select input water level data set (Cancel to use SWL=0):';           
            [wl_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseClass',{'ctWaterLevelData','ctTidalAnalysis',...
                                              'muiUserModel'},...
                                'ListSize',[300,100]); 

            swl = zeros(size(wvtime));               
            if ok<1 || isempty(wl_crec)
                getdialog('Using SWL=0');
                inputxt = sprintf('%s, 0mOD used for water level',inputxt);
                wl_crec = 0;     %assign a null value if no water level data available
            else
                wldst = getDataset(muicat,wl_crec,1); 
                %check that there is water level data for period of interest
                [idst,idnd] = ts2_endpoints_in_ts1(wvdst,wldst);
                if isempty(idst)
                    getdialog('Data do not overlap. Using SWL=0');
                    inputxt = sprintf('%s, 0mOD used for water level',inputxt);
                    wl_crec = 0; %assign a null value if no water level data available
                else 
                    %select a variable from the water level dataset
                    varnames = wldst.VariableNames;
                    idx = 1;
                    if length(varnames)>1
                        [idx,ok] = listdlg('Name','WL options', ...
                            'PromptString','Select a variable:', ...
                            'SelectionMode','single','ListSize',[200,100],...
                            'ListString',varnames);
                        if ok<1, idx = 1; end
                    end
                    wldata = wldst.(varnames{idx});
                    wltime = wldst.RowNames;
                    swltime = wvtime(idst:idnd);                                        
                    %now interpolate water levels onto wave height times
                    swl(idst:idnd,1) = interp1(wltime,wldata,swltime,'linear','extrap');
                    swl(isnan(swl)) = 0;
                    inputxt = sprintf('%s, %s used for water levels',...
                                                inputxt,wldst.Description);
                end
            end

            [wvdst,timerange] = getSubSet(obj,wvdst);      %allow user to extract a subset 
            swl = swl(timerange);

            tsdst = addvars(wvdst,swl,'NewVariableNames','swl');
            
            if strcmp(source,'Measured spectra')                
                tsprops = getDataset(muicat,wv_caserec,2);  %2 selects Properties dataset from spectra input
                tsprops = removerows(tsprops,find(~timerange));
                tsdst = horzcat(tsdst,tsprops);                
            end
            tsdst = activatedynamicprops(tsdst);

            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:';             
            sptobj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            msgtxt = 'Spectral Transfer Table not found';
            if isempty(sptobj), warndlg(msgtxt); return; end
            inputxt = sprintf('%s, %s used for spectral transfer',inputxt,...
                                      sptobj.Data.Inshore.Description);
                
            %assign the run parameters to the model instance
            spt_caserec = caseRec(mobj.Cases,sptobj.CaseIndex);
            if wl_crec==0
                setRunParam(obj,mobj,wv_caserec,spt_caserec);
            else
                setRunParam(obj,mobj,wv_caserec,wl_crec,spt_caserec); %input caserecs passed as varargin     
            end
        end
%%
function [subdst,timerange] = getSubSet(~,tsdst)
            %
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
            timerange = var2range(tsdst.RowRange);
            selection = inputgui('FigureTitle','Levels',...
                                 'InputFields',{'Time'},...
                                 'Style',{'edit'},...
                                 'ControlButtons',{'Ed'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{timerange},...
                                 'PromptText','Select time range to use');

            if isempty(selection)
                subdst = tsdst;
                timerange = true(size(tsdst.RowNames));
            else
                seltime = range2var(selection{1});
                times = tsdst.RowNames;
                timerange = isbetween(times,seltime{:});
                subdst = getDSTable(tsdst,'RowNames',times(timerange));
            end
            
        end
%%
        function results = callRefraction(~,tsdst,site_params)
            %parse inputs for call to refraction function
            z0 = site_params.OffshoreBedLevel;
            theta0 = site_params.OffshoreAngle;
            Kf = site_params.FrictionCoefficient;
            
            prompt = {'Angle of bed contour at refraction point(s) (degTN)',...
                      'Bed level at refraction point(s) (mOD)'};
            default = {num2str(theta0),num2str(-100)};
            answer = inputdlg(prompt,'Deepwater refraction',1,default);
            if isempty(answer), return; end            
            thetaT = str2double(answer{1});
            zT = str2double(answer{2});
            dep0 = tsdst.swl-z0;   %source water depth from swl to bed level 
            depT = tsdst.swl-zT;   %target water depth from swl to bed level 
            
            %this is a general refraction model so 'isshore' is false
            [Hso,Diro] = refraction(tsdst.Hs,tsdst.Tp,tsdst.Dir,[dep0,depT],...
                                                [theta0,thetaT],Kf,false); 
            results = {Hso,tsdst.Tp,Diro};
        end
%%
        function [dspec,dsprop] = setDSproperties(~)
            %define the variables in the dataset
            %define the metadata properties for the demo data set
            dspec = struct('Variables',[],'Row',[],'Dimensions',[]); dsprop = dspec;    
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dspec.Variables = struct(...
                'Name',{'So','Si'},...
                'Description',{'Offshore spectral density','Inshore spectral density'},...
                'Unit',{'m^2/Hz','m^2/Hz'},...
                'Label',{'Spectral density (m^2/Hz)','Spectral density (m^2/Hz)'},...
                'QCflag',repmat({'raw'},1,2)); 
            dspec.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dspec.Dimensions = struct(...    
                'Name',{'dir','freq'},...
                'Description',{'Direction','Frequency'},...
                'Unit',{'degTN','Hz'},...
                'Label',{'Direction (degTN)','Frequency (Hz)'},...
                'Format',{'',''});    

            dsprop.Variables = struct(...   
                'Name',{'Hsi','T2i','Diri','Tpi','Diripk','kw','kt2','ktp',...
                                                       'kd','swl','depi'},...
                'Description',{'Inshore wave height','Inshore mean period',...
                               'Inshore wave direction','Inshore peak period',...
                               'Inshore peak direction','Wave transfer coefficient',...
                               'Mean period coefficient','Peak period coefficient',...
                               'Mean direction shift','Still water level','Inshore depth'},...                               
                'Unit',{'m','s','degTN','s','degTN','-','-','-','deg','mOD','m'},...
                'Label',{'Wave height (m)','Wave period (s)','Wave direction (degTN)',...
                         'Wave period (s)','Wave direction (degTN)',...
                         'Transfer coefficient, kw','Transfer coefficient, kt2',...
                         'Transfer coefficient, ktp','Direction shift (deg)',...
                         'Water level (mOD)','Water depth (m)'},...
                'QCflag',repmat({'model'},1,11)); 
            dsprop.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsprop.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''}); 
        end   
    end  
%%
    methods (Static, Access=private)
        function inputs = getForcingConditions()
            %get the user input of wave or wind conditions
            inputs = []; 

            answer = questdlg('Wind, Wave or Spectrum input?','Input',...
                                         'Wind','Wave','Spectrum','Wave');            
            if strcmp(answer,'Wave')
                promptxt = {'Wave height (m)','Peak period (s)',...
                                 'Wave direction (degTN)','Still water level'};
                defaults = {'1.1','8.2','185','0.0'};
                inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                if isempty(inpt), return; end  %user cancelled
                inputs.Hs = str2double(inpt{1});
                inputs.Tp = str2double(inpt{2});
                inputs.Dir = str2double(inpt{3});
                inputs.swl = str2double(inpt{4});
            elseif strcmp(answer,'Wind')
                promptxt = {'Wind Speed (m/s)','Wind Direction (degTN)',...
                            'Height above msl (m)','Fetch Length (m)',...
                            'Still water level'};                                 
                defaults = {'20.0','185','10.0','20000','0.0'};
                inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                if isempty(inpt), return; end  %user cancelled
                inputs.AvSpeed = str2double(inpt{1});
                inputs.Dir = str2double(inpt{2});
                inputs.zW = str2double(inpt{3});
                inputs.Fetch = str2double(inpt{4});
                inputs.swl = str2double(inpt{5});
            else
                promptxt = {'Still water level','Include depth saturation (1=Yes,0=No)'};
                defaults = {'0.0','1'};
                inpt = inputdlg(promptxt,'Input conditions',1,defaults);
                if isempty(inpt), return; end  %user cancelled
                inputs.swl = str2double(inpt{1});
                inputs.issat = logical(str2double(inpt{2}));    
                [filename,path,~] = getfiles('MultiSelect','off','PromptText','Select file:');
                if filename==0, inputs = []; return; end  %user cancelled
                varlist = {'',[path,filename]};
                specdst = wave_cco_spectra('getData',varlist{:});
                tsdst = horzcat(specdst.Spectra,specdst.Properties);
                tsdst = activatedynamicprops(tsdst,specdst.Properties.VariableNames);
                tsdst = addvars(tsdst,inputs.swl,'NewVariableNames','swl');
                inputs.tsdst = tsdst;
            end
            inputs.source = answer;
        end
%%
        function inputs = getSprectraConditions()
            %get conditions for getForcingConditions when the input is a
            %measured spectrum and is called for same purpose in runWaves
            inputs = []; 
            promptxt = {'Include depth saturation (1=Yes,0=No)'};
            defaults = {'1'};
            inpt = inputdlg(promptxt,'Input conditions',1,defaults);
            if isempty(inpt), return; end  %user cancelled
            inputs.issat = logical(str2double(inpt{1}));    
            inputs.form = 'Measured';           %initialise select properties
            inputs.source = 'Spectrum';         
            inputs.ismodel = false;
            if inputs.issat
                inputs.txt = 'Measured spectrum and including depth saturation';
            else
                inputs.txt = 'Measured spectrum and excluding depth saturation';
            end
        end
%%
        function off = addWaveConditions(SGo,Dims,off)
            %when using wind input define the offshore wave conditions
            g = 9.81;
            off.Tp = 0.54*g^-0.77*off.AvSpeed.^0.54.*off.Fetch.^0.23;

            dir_int = 0.5;     %interval used to interpolate directions (deg)
            radint = deg2rad(dir_int);
            So = trapz(radint,abs(trapz(Dims.freq,SGo,2))); %integral of offshore spectrum
            off.Hs = 4*sqrt(So);
        end
    end
end    