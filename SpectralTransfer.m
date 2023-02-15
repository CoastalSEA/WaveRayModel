classdef SpectralTransfer < muiDataSet                    
%-------class help---------------------------------------------------------
% NAME
%   Model_template.m
% PURPOSE
%   Class description - Class for Model XXXXX to be run as a muitoolbox app
%
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
    end
    
    methods (Access = private)
        function obj = SpectralTransfer()             
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = SpectralTransfer;                           
            dsp = modelDSproperties(obj);
            
%             %now check that the input data has been entered
%             %isValidModel checks the InputHandles defined in XXXMainUI <<%Edit to UI
%             if ~isValidModel(mobj, metaclass(obj).Name)  
%                 warndlg('Use Setup to define model input parameters');
%                 return;
%             end
            muicat = mobj.Cases;    
%--------------------------------------------------------------------------
% Model code>
%--------------------------------------------------------------------------
            %select back tracking ray case to use
            promptxt = 'Select Backward Ray Trace Dataset:';
            [rayobj,rayrec] = selectCaseObj(muicat,{'backward_model'},...
                                                    {'RayTracks'},promptxt);
            %assign the run parameters to the model instance
            setRunParam(obj,mobj,rayrec); %input caserecs passed as varargin 

            [results,indir,intable] = specTransfer(obj,rayobj);            
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(results{:},'RowNames',indir,'DSproperties',dsp);
            dst.Dimensions.Period = rayobj.Data.Dataset.Dimensions.Period; 
            dst.Dimensions.WaterLevel = rayobj.Data.Dataset.Dimensions.WaterLevel;  
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = sprintf('Derived using Case ID %d',rayobj.CaseIndex);
            dst.UserData.Inshore = intable;
            %save results
            setDataSetRecord(obj,muicat,dst,'spectral_model');
            getdialog('Run complete');
        end
    end
%%
    methods
        function tabPlot(obj,src,mobj) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            dst = obj.Data.Dataset;

            %set >Figure button and create axes
            tabcb  = @(src,evdat)tabPlot(obj,src,mobj); 
            ax = tabfigureplot(obj,src,tabcb,false);
            ax.NextPlot = 'add';

            %extract required data
            ki = get_selection(obj);
            T = dst.Dimensions.Period;
            zwl = dst.Dimensions.WaterLevel;
            phi = dst.RowNames;
            var = dst.theta(:,:,ki);

            mindir = min(var,[],'All');
            mindir = mindir-mod(mindir,30)+30;
            maxdir = max(var,[],'All');
            maxdir = maxdir-mod(maxdir,30);
            nc = mindir:30:maxdir;
            %construct Q-Plot
            
            surf(ax,T,phi,var);
            view(2);
            shading interp
            hold on
            [C,h] = contour3(ax,T,phi,var,nc,'-k');
            clabel(C,h,'LabelSpacing',300,'FontSize',8)
            hold off
            %add the colorbar and labels
            cb = colorbar;
            cb.Label.String = 'Offshore direction (degTN)';
            xlabel('Wave period (s)'); 
            ylabel('Inshore direction (degTN)'); 
            title(sprintf('%s for water level of %.2g mOD',dst.Description,zwl(ki)));
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
        end
    end 
%%    
    methods (Access = private)
        function [spectran,indir,intable] = specTransfer(~,rayobj)
            %generate the spectral transfer table of offshore directions,
            %depths and celerities
            dst = rayobj.Data.Dataset;
            
            indir = dst.RowNames;  %assign inshore ray directions
            %use first point of frist ray to define inshore point properties
           
            xi = dst.xr{1,1,1}(1);
            yi = dst.yr{1,1,1}(1);
            hi = dst.depth{1,1,1}(1);
            ci = dst.celerity{1,1,1}(1);
            cgi = dst.cgroup{1,1,1}(1);
            intable = table(xi,yi,hi,ci,cgi);  %inshore properties
            %compile table data for offshore properties
            ndir = length(indir);
            nper = length(dst.Dimensions.Period);
            nwls = length(dst.Dimensions.WaterLevel);
            c = zeros(ndir,nper,nwls); cg = c; h = c; offdir = c; 
            for i=1:ndir
                for j=1:nper
                    for k=1:nwls
                        theta = dst.alpha{i,j,k}(end);
                        offdir(i,j,k) = mod(compass2trig(theta,true),360);
                        h(i,j,k) = dst.depth{i,j,k}(end);
                        c(i,j,k) = dst.celerity{i,j,k}(end);
                        cg(i,j,k) = dst.cgroup{i,j,k}(end);
                    end
                end
            end
            spectran = {offdir,h,c,cg};
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
            
            zwl = obj.Data.Dataset.Dimensions.WaterLevel;
            selection = inputgui('FigureTitle','Levels',...
                                 'InputFields',{'Water level'},...
                                 'Style',{'popupmenu'},...
                                 'ActionButtons', {'Select','Cancel'},...
                                 'DefaultInputs',{string(zwl)},...
                                 'PromptText','Select values to use');
            if isempty(selection)
                options = []; 
            else
                options = selection{1};
            end  
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
                'Name',{'theta','depth','celerity','cgroup'},...
                'Description',{'Offshore Direction',...
                           'Water depth','Celerity','Group Celerity'},...
                'Unit',{'degTN','m','m/s','m/s'},...
                'Label',{'Offshore Direction (degTN)','Water depth (m)',...
                                'Celerity (m/s)','Group Celerity (m/s)'},...                           
                'QCflag',repmat({'model'},1,4)); 
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