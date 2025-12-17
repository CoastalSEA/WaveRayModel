classdef WRM_WaveModel < waveModels
%
%-------class help------------------------------------------------------
% NAME
%   WRM_WaveModel.m
% PURPOSE
%    Class for wave refraction using backward ray transfer function.
%    Constructs inshore time series from an offshore timeseries.Inherits 
%    abstract class waveModels which is a subclass of muiDataSet.
% SEE ALSO
%   muiDataSet, waveModels, WaveRayModel, RayTracks, SpectralTransfer
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
    end                  %abstract property required by waveModels
    
    methods (Access={?muiDataSet,?muiStats,?ctWaveData,?ctWaveSpectra})
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
            %function to run the wave refraction spectral transfer model.
            obj = WRM_WaveModel;                            
            [dspec,dsprop] = setDSproperties(obj);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            [tsdst,meta] = getInputData(obj,mobj);            
            if isempty(tsdst), return; end   %user cancelled data selection

            if strcmp(meta.source,'Measured spectra')
                ModelType = 'transfer spectra';
            else
                ModelType = 'transfer timeseries';
            end

            %select a spectral transfer case to use
            [sptobj,sptmeta] = SpectralTransfer.getSTcase(mobj);
            caserecs = [meta.caserecs,sptmeta.caserecs];
            setRunParam(obj,mobj,caserecs{:})     %assign run parameters
            %add spectral transfer selection to meta data
            inputxt = sprintf('%s, %s',meta.inptxt,sptmeta.inptxt);

            [offobj,inobj] = runWaves(sptobj,tsdst,meta);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------    
            [sp,results] = unpackSpectrum(inobj,offobj);
            results = addvars(results,sp.swl,sp.depths,'NewVariableNames',{'swl','depi'});
            dir = offobj(1).Spectrum.dir;
            freq = offobj(1).Spectrum.freq;

            answer = questdlg('Save the spectra?','Wave model','Yes','No','No');
            if strcmp(answer,'Yes'), issave=true; else, issave=false; end
            %check whether wave spectra should also be saved
            Sot = sp.Sot; %#ok<NASGU>
            sze = 2*getfield(whos('Sot'),'bytes')*9.53674e-7;
            if issave
                questxt = sprintf('Save the full wave spectra (arrays are %.1f Mb) or spt format',sze);
                answer = questdlg(questxt,'Wave model','Full','SPT format','SPT format');
                source = sprintf('Class %s, using %s',metaclass(obj).Name,ModelType);                        
                if strcmp(answer,'SPT format')
                    dst.OffshoreSpectra = saveSpectrum(offobj);
                    dst.OffshoreSpectra.Source = source;
                    dst.OffshoreSpectra.MetaData = inputxt; 
                    dst.InshoreSpectra = saveSpectrum(inobj);
                    dst.InshoreSpectra.Description = sprintf('Inshore using %s',dst.OffshoreSpectra.Description);
                    dst.InshoreSpectra.Source = source;
                    dst.InshoreSpectra.MetaData = inputxt;  
                else
                    dst.oiSpectra = dstable(sp.Sot,sp.Sit,'RowNames',sp.time,'DSproperties',dspec); 
                    dst.oiSpectra.Dimensions.dir = dir;    %NB order is X,Y and must
                    dst.oiSpectra.Dimensions.freq = freq;  %match variable dimensions  
                    %assign metadata about model
                    dst.oiSpectra.Source = source;
                    dst.oiSpectra.MetaData = inputxt;                    
                end                    
            end

            dst.Properties = dstable(results,'RowNames',sp.time,'DSproperties',dsprop);                      
            %assign metadata about model            
            dst.Properties.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                         ModelType);
            dst.Properties.MetaData = inputxt;   
            %add depths of inshore point for which there are backward rays
            dst.Properties.UserData = sp.depths;
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            obj.ModelType = 'Inwave_model';           %added to match ctWaveModel
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
            if isempty(sptobj)
                getdialog('Spectral Transfer table not found'); return; 
            end

            isout = checkWLrange(sptobj,offdata.swl);
            if isout
                warndlg('Water levels are outside the range of the Transfer Table')
                return;
            end

            if strcmp(offdata.source,'Spectrum')               
                select.form = 'Measured';             %initialise select properties
                select.source = 'Spectrum';         
                select.ismodel = false;
                select.issat = offdata.issat;         %copy to sprectrum selection  
                select.freq = offdata.tsdst.Dimensions.freq;
                intable = offdata.tsdst.DataTable;
                [SGo,SGi,Dims] = get_inshore_spectrum(sptobj,intable,select);
                ins = get_inshore_wave(SGo,SGi,Dims,intable,select);
            else
                select = get_model_selection(sptobj,offdata.source); %select spectral form and data type
                if isempty(select), return; end       %user cancelled

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
            [tsdst,~,source] = getInputData(obj,mobj,false);
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
             if isempty(sptobj)
                getdialog('Spectral Transfer table not found'); return; 
             end

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

%%
        function runBatchMode(mobj)
            %run spectral transfer for multiple points and load a set of
            %dstables as datasets in a single case
            obj = WRM_WaveModel;                            
            [~,dsprop] = setDSproperties(obj);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            %get the timeseries input data and site parameters
            %Note: getInputData calls setRunParam
            % [tsdst,sptrecs,inputxt,~] = getInputData(obj,mobj,true); %multi-case
            % if isempty(tsdst), return; end   %user cancelled data selection
            muicat = mobj.Cases;
            [tsdst,meta] = getInputData(obj,mobj);            
            if isempty(tsdst), return; end   %user cancelled data selection

            %select a spectral transfer case to use
            [sptrecs,sptmeta] = SpectralTransfer.getSTcase(mobj,true);
            caserecs = [meta.caserecs,sptmeta.caserecs];
            setRunParam(obj,mobj,caserecs{:})     %assign run parameters

            ModelType = 'transfer timeseries';
            % inputxt = meta.inptxt;
            hw = waitbar(0,'Processing point 0');
            npnts = length(sptrecs);
            for i=1:npnts      %NOT parfor because used in runWaves               
                waitbar(i/npnts,hw,sprintf('Processing point %d',i));
                sptobj = getCase(muicat,sptrecs(i));
                [offobj,inobj] = runWaves(sptobj,tsdst,meta);
                %each variable should be an array in the 'results' cell array
                %if model returns single variable as array of doubles, use {results}
                time = tsdst.RowNames;
                [sp,results] = unpackSpectrum(inobj,offobj);
                results = addvars(results,sp.swl,sp.depths,'NewVariableNames',{'swl','depi'});
                adst = dstable(results,'RowNames',time,'DSproperties',dsprop);
                %assign metadata about model            
                adst.Source =  sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                             ModelType);
                adst.MetaData = sprintf('%s, %s used for spectral transfer',meta.inptxt,...
                                              sptobj.Data.Inshore.Description);
                %add depths of inshore point for which there are backward rays
                adst.UserData = sptobj.Data.Inshore.UserData.Depths;
                pname{i} = sprintf('Point%d',i);
                pdst(i) = adst;
            end
            
            for j = 1:npnts
                dst.(pname{j}) = pdst(j);
            end
            delete(hw)
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            obj.ModelType = 'Inwave_model';           %added to match ctWaveModel
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');          
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab            
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 

%%
    methods (Access=private)
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
                'Name',{'Hs','m0','Dir','Sp','Tp','Dp',...
                        'Sfdpk','Tfdpk','Dfdpk','T2',...
                        'kw','kt2','ktp','kd','swl','depi'},...
                'Description',{'Inshore wave height',...
                               'Inshore zero moment',...
                               'Inshore wave direction',...
                               'Inshore peak spectral density'...
                               'Inshore peak period',...
                               'Inshore peak direction',...
                               'Inshore f-d spectral density peak',...
                               'Inshore f-d peak period',...
                               'Inshore f-d peak direction',...
                               'Inshore mean period',...
                               'Wave transfer coefficient',...
                               'Mean period coefficient',...
                               'Peak period coefficient',...
                               'Mean direction shift',...
                               'Still water level',...
                               'Inshore depth'},...                               
                'Unit',{'m','m^2','degTN','m^2/Hz','s','degTN','m^2/Hz',...
                            's','degTN','s','-','-','-','deg','mOD','m'},...
                'Label',{'Wave height (m)','Zero moment (m^2)',...
                         'Wave direction (degTN)','Spectral density (m^2/Hz)',...                
                         'Wave period (s)','Wave direction (degTN)',... 
                         'Spectral density (m^2/Hz)','Wave period (s)',... 
                         'Wave direction (degTN)','Wave period (s)',...
                         'Transfer coefficient, kw','Transfer coefficient, kt2',...
                         'Transfer coefficient, ktp','Direction shift (deg)',...
                         'Water level (mOD)','Water depth (m)'},...
                'QCflag',repmat({'model'},1,16)); 
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