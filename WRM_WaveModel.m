classdef WRM_WaveModel < muiDataSet 
%
%-------class help------------------------------------------------------
% NAME
%   ctWaveModel.m
% PURPOSE
%    Class for inshore and offshore wave models to be run in CoastalTools 
%    and other ModelUI apps
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        InshoreBedSlope   %bed slope as single value or vector array
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
            %function to run the wave refraction model. isin is true for
            %inshore waves and false for deepwater waves
            obj = WRM_WaveModel;                            
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in CoastalTools
            %NB - water level data not checked because optional
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            %Note: getInputData calls setRunParam
            [tsdst,inputxt] = getInputData(obj,mobj);
            if isempty(tsdst), return; end   %user cancelled data selection
            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            sptobj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            
            results = runWaves(sptobj,mobj,tsdst);

%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            time = tsdst.RowNames;
            dst = dstable(results,'RowNames',time,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model            
            dst.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                         obj.ModelType);
            dst.MetaData = inputxt;            
            %save results
            setDataSetRecord(obj,mobj.Cases,dst,'wave_model');
            getdialog('Run complete');
        end
%%
        function runAnimation(mobj)
            %create an animation of the 2-D spectrum surfaces using a
            %timeseries input
            obj = WRM_WaveModel; 
            [tsdst,~] = getInputData(obj,mobj);
            if isempty(tsdst), return; end   %user cancelled data selection
            %get the refraction transfer table
            promptxt = 'Select a Transfer Table Case to use:'; 
            sptobj = selectCaseObj(mobj.Cases,[],{'SpectralTransfer'},promptxt);
            
            [SGo,SGi,Dims] = runSpectra(sptobj,mobj,tsdst);

            off_in_plot(sptobj,1./Dims.f,Dims.xso,squeeze(SGo(2,:,:)),squeeze(SGi(2,:,:))) 
        end
    end
%%
    methods
        function [tsdst,caserec] = getWaveModelDataset(obj,mobj,type,varnames,caserec)
            %prompt user to select a model wave dataset and add Tp if inshore
            %
            muicat = mobj.Cases;
            if nargin<4
                varnames = {'Tp'};  %default is to add Tp
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
                for i=1:length(varnames)                    
                    if any(strcmp(dstnames,varnames{i})) && ...
                                            ~isempty(inpdst.(varnames{i}))
                        tsdst = addvars(tsdst,inpdst.(varnames{i}),...
                                           'NewVariableNames',varnames{i});
                    else
                        warndlg('Variable %s not found so not added to wave dataset',...
                                                           varnames{i});
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
        function [tsdst,inputxt] = getInputData(obj,mobj)
            %prompt user to select wave and water level data and return in
            %input dstable of data and metadata for inputs used
            tsdst = []; inputxt = [];
            muicat = mobj.Cases;
            promptxt = 'Select input wave data set:';           
            [wv_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                           'CaseClass',{'ctWaveData'},'ListSize',[300,100]);                                    
            if ok<1, return; end
            wvdst = getDataset(muicat,wv_crec,1);
            wvtime = wvdst.RowNames;
            inputxt = sprintf('%s used for offshore waves',wvdst.Description);

            promptxt = 'Select input water level data set (Cancel to use SWL=0):';           
            [wl_crec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseClass',{'ctWaterLevelData','ctTidalAnalysis'},...
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
            tsdst = addvars(wvdst,swl,'NewVariableNames','swl');
            %assign the run parameters to the model instance
            if wl_crec==0
                setRunParam(obj,mobj,wv_crec);
            else
                setRunParam(obj,mobj,wv_crec,wl_crec); %input caserecs passed as varargin     
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
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...   
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
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});  
        end
        
    end    
end    