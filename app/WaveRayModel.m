classdef WaveRayModel < muiModelUI                       
%
%-------class help---------------------------------------------------------
% NAME
%   WaveRayModel.m
% PURPOSE
%   Main GUI for the ModelSkill interface, which implements the 
%   muiModelUI abstract class to define main menus.
% SEE ALSO
%   Abstract class muiModelUI.m and tools provided in muitoolbox
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
% 
    properties  (Access = protected)
        %implement properties defined as Abstract in muiModelUI
        vNumber = '1.2'
        vDate   = 'Nov 2025'
        modelName = 'WaveRayModel'                       
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = WaveRayModel             
            %constructor function initialises GUI
            isok = check4muitoolbox(obj);
            if ~isok, return; end
            %
            obj = setMUI(obj);             
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus    
            %classes required to run model:
            obj.ModelInputs.RayTracks = {'WRM_RunParams'};
            obj.ModelInputs.WRM_Bathy = {'GD_GridProps'};
            obj.ModelInputs.SpectralTransfer = {''};
            obj.ModelInputs.WRM_WaveModel = {'WRM_RunParams'};
            obj.ModelInputs.WRM_SpectraModel = {'RayTracks'};
            obj.ModelInputs.WRM_Mesh = {'WRM_Mesh'};
            obj.ModelInputs.WRM_SedimentTransport = {'WRM_WaveModel','WRM_SedimentTransport'};
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D','2DT','3DT'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Timeseries','Taylor'};  
            
            modelLogo = 'WRM_logo.jpg';  %default splash figure - edit to alternative
            initialiseUI(obj,modelLogo); %initialise menus and tabs                  
        end    
        
%% ------------------------------------------------------------------------
% Definition of Menu Settings
%--------------------------------------------------------------------------
        function menu = setMenus(obj)
            %define top level menu items and any submenus
            %MenuLabels can any text but should avoid these case-sensitive 
            %reserved words: "default", "remove", and "factory". If label 
            %is not a valid Matlab field name this the struct entry
            %is modified to a valid name (eg removes space if two words).
            %The 'gcbo:' Callback text triggers an additional level in the 
            %menu. Main menu labels are defined in sequential order and 
            %submenus in order following each brach to the lowest level 
            %before defining the next branch.         
                                                              % << Edit menu to suit model 
            MenuLabels = {'File','Tools','Project','Setup','Run',...
                                                        'Analysis','Help'};
            menu = menuStruct(obj,MenuLabels);  %create empty menu struct
            %
            %% File menu --------------------------------------------------
             %list as per muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %list as per muiModelUI.toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Model','Figures','Cases'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,3]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            % submenu for Scenarios
            menu.Project(2).List = {'Edit Description','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,6]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------
            menu.Setup(1).List = {'Input Data','Grid Parameters',...
                                  'Grid Tools','Run Parameters',...
                                  'Data Clean-up','Model Constants'};                                                                      
            menu.Setup(1).Callback = [{'gcbo;'},{@obj.setupMenuOptions},...
                                      {'gcbo;'},{'gcbo;'},{'gcbo;'},...
                                      {@obj.setupMenuOptions}];
            %add separators to menu list (optional - default is off)
            menu.Setup(1).Separator = [repmat({'off'},[1,5]),{'on'}]; %separator preceeds item
            
            % submenu for Import Data (if these are changed need to edit
            % loadMenuOptions to be match)
            menu.Setup(2).List = {'Bathymetry','Waves','Water levels','Winds'};
            menu.Setup(2).Callback = repmat({'gcbo;'},[1,4]);
            % submenu for Gridded and Timeseries Data 
            nitems = 4;
            for j=1:nitems  %add standard submenu to all import menu items
                menu.Setup(j+2).List = {'Load','Add','Delete','Quality Control'};                                   
                menu.Setup(j+2).Callback = repmat({@obj.loadMenuOptions},[1,4]);
            end
            % submenu for Grid Tools
            menu.Setup(7).List = {'Translate Grid','Rotate Grid',...
                                  'Re-Grid','Sub-Grid',...
                                  'Combine Grids','Add Surface','Infill Surface',...
                                  'To curvilinear','From curvilinear',... 
                                  'Display Dimensions','Difference Plot',...
                                  'Plot Sections','Grid Image','Digitise Line',...
                                  'Export xyz Grid','User Function'};                                                                         
            menu.Setup(7).Callback = repmat({@obj.gridMenuOptions},[1,16]);
            menu.Setup(7).Separator = [repmat({'off'},[1,7]),...
                       {'on','off','on','off','off','off','on','on','on'}]; %separator preceeds item  
            % submenu for Run Parameters
            menu.Setup(8).List = {'Run Conditions','Forward Tracking',...
                                          'Backward Tracking','Batch Start Points'};
            menu.Setup(8).Callback = repmat({@obj.setupMenuOptions},[1,4]);
            menu.Setup(8).Separator = [repmat({'off'},[1,3]),{'on'}];
            % submenu for Data clean-up
            menu.Setup(9).List = {'Concatenate two timeseries',...
                            'Resample timeseries','Patch timeseries',...                
                            'Trim timeseries'};
            menu.Setup(9).Callback = repmat({@obj.datacleanup},[1,4]);
            menu.Setup(9).Separator = {'off','off','off','off'};
            
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Check Start Points','Forward Rays',...
                                'Check Start Depth','Backward Rays',...
                                'Batch Run Rays','Transfer Table',...
                                'Run Wave Timeseries','Batch Run Waves',...
                                'Batch Sediment Transport',...
                                'Create Mesh','Test Grid','Derive Output'};
            menucall = repmat({@obj.runMenuOptions},[1,12]);            
            menu.Run(1).Callback = menucall;
            menu.Run(1).Separator = [repmat({'off'},[1,2]),{'on','off','off',...
                                   'on','off','off','on','on','off','on'}];
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics','Plot Mesh',...
                                     'Ray Plots','Spectral Plots',...
                                     'Multi-point Plots'};
            menu.Analysis(1).Callback = [repmat({@obj.analysisMenuOptions},[1,4]),...
                                        {'gcbo;'},{@obj.analysisMenuOptions}];
            menu.Analysis(1).Separator = {'off','off','on','off','off','on'};
            
            %submenu for Spectral Plots
            menu.Analysis(2).List = {'Transfer Table','Transfer Coefficients',...
                                     'Spectrum Plots','O/I Spectrum','O/I Animation'};
            menu.Analysis(2).Callback = repmat({@obj.analysisMenuOptions},[1,5]);
            menu.Analysis(2).Separator = {'off','off','on','off','off'};
            
            %% Help menu --------------------------------------------------
            menu.Help.List = {'Documentation','Manual'};
            menu.Help.Callback = repmat({@obj.Help},[1,2]);
            
        end
        
%% ------------------------------------------------------------------------
% Definition of Tab Settings
%--------------------------------------------------------------------------
        function [tabs,subtabs] = setTabs(obj)
            %define main tabs and any subtabs required. struct field is 
            %used to set the uitab Tag (prefixed with sub for subtabs). 
            %Order of assignment to struct determines order of tabs in figure.
            %format for tabs: 
            %    tabs.<tagname> = {<tab label>,<callback function>};
            %format for subtabs: 
            %    subtabs.<tagname>(i,:) = {<subtab label>,<callback function>};
            %where <tagname> is the struct fieldname for the top level tab.
            tabs.Data  = {'   Data  ',@obj.refresh}; 
            tabs.Rays  = {'   Rays  ',@obj.refresh}; 
            tabs.Transfer  = {' Transfer ',@obj.refresh}; 
            tabs.Inputs = {'  Inputs  ',@obj.InputTabSummary};
            tabs.Plot   = {'  Q-Plot  ',@obj.getTabData};
            tabs.Stats = {'   Stats   ',@obj.setTabAction};
            subtabs = [];
        end
       
%%
        function props = setTabProperties(~)
            %define the tab and position to display class data tables
            %props format: {class name, tab tag name, position, ...
            %               column width, table title}
            % position and column widths vary with number of parameters
            % (rows) and width of input text and values. Inidcative
            % positions:  top left [0.95,0.48];    top right [0.95,0.97]
            %         bottom left [0.45, 0.48]; bottom rigth [0.45,0.97]
            props = {...
                'WRM_RunParams','Inputs',[0.94,0.97],{160,70},'Wave model parameters:'; ...
                'WRM_FT_Params','Inputs',[0.94,0.50],{160,90},'Forward tracking parameters:'; ...
                'WRM_BT_Params','Inputs',[0.62,0.50],{160,90},'Backward tracking parameters:'; ...
                'GD_GridProps','Inputs',[0.35,0.50],{160,90}, 'Grid parameters:';...
                'WRM_Bathy','Inputs',[0.35,0.97],{160,70}, 'Bathymetry parameters:';...
                'WRM_Mesh','Inputs',[0.50,0.97],{160,70}, 'Mesh parameters:';...
                'WRM_SedimentTransport','None',[1,1],{0,0},'None'}; %not displayed on a tab
        end    
 %%
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            switch src.Tag                                  
                case 'Plot' 
                     if isa(cobj,'RayTracks') || isa(cobj,'SpectralTransfer')
                        tabPlot(cobj,src,obj);
                     else
                        tabPlot(cobj,src);
                     end
                case 'Stats'                    
                    lobj = getClassObj(obj,'mUI','Stats');
                    if isempty(lobj), return; end
                    tabStats(lobj,src);     
            end
        end      
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu -----------------------------------------------------
        %use default menu functions defined in muiModelUI
                
        %% Project menu ---------------------------------------------------
        %use default menu functions defined in muiModelUI           

        %% Setup menu -----------------------------------------------------
        function setupMenuOptions(obj,src,~)
            %callback functions for data input
            switch src.Text
                case 'Grid Parameters'
                    GD_GridProps.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);
                case 'Run Conditions'                         
                    WRM_RunParams.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);
                case 'Forward Tracking'                         
                    WRM_FT_Params.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);
                case 'Backward Tracking'                         
                    WRM_BT_Params.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);  
                case 'Batch Start Points'
                    WRM_BT_Params.batchPoints(obj); 
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
            end
        end  
%%
        function gridMenuOptions(obj,src,~)
            %callback functions for grid tools options
            gridclasses = {'GD_ImportData','WRM_Bathy'}; %add other classes if needed
            GD_ImportData.gridMenuOptions(obj,src,gridclasses);
            DrawMap(obj);
        end
%%
        function datacleanup(obj,src,~)
            %all cleanup options call the same function
            ct_data_cleanup(obj.Cases,src);
            DrawMap(obj);
        end
%%
        function loadMenuOptions(obj,src,~)
            %callback functions to import timeseries data
            switch src.Parent.Text
                case 'Bathymetry'
                    classname = 'GD_ImportData';
                case 'Waves'
                    classname = 'ctWaveData';
                case 'Water levels'
                    classname = 'ctWaterLevelData';
                case 'Winds'
                    classname = 'ctWindData';
            end
            %
            switch src.Text
                case 'Load'
                    fname = sprintf('%s.loadData',classname);
                    callStaticFunction(obj,classname,fname); 
                case 'Add'
                    useCase(obj.Cases,'single',{classname},'addData');
                case 'Delete'
                    useCase(obj.Cases,'single',{classname},'deleteGrid');
                case 'Quality Control'
                    useCase(obj.Cases,'single',{classname},'qcData');
            end
            DrawMap(obj);
        end

        %% Run menu -------------------------------------------------------
        function runMenuOptions(obj,src,~)
            %callback functions to run model            
            switch src.Text     
                case 'Check Start Points'
                    RayTracks.checkStart(obj);
                case 'Forward Rays'
                    RayTracks.runModel(obj,src,[],[]); 
                case 'Check Start Depth'
                    RayTracks.startDepth(obj);
                case 'Backward Rays'                    
                    RayTracks.runModel(obj,src,[],[]);
                case 'Batch Run Rays'
                    RayTracks.batchRun(obj);
                case 'Transfer Table'
                    SpectralTransfer.runModel(obj,[],[]);
                case 'Run Wave Timeseries'
                    WRM_WaveModel.runModel(obj);
                case 'Batch Run Waves'
                    WRM_WaveModel.runBatchMode(obj);
                case 'Batch Sediment Transport'
                    WRM_SedimentTransport.runModel(obj);
                case 'Create Mesh'
                    WRM_Mesh.runModel(obj);    
                case 'Test Grid'
                    WRM_Bathy.runModel(obj);
                case 'Derive Output'
                    obj.mUI.ManipUI = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            %callback functions for analysis menu
            promptxt = 'Select a Case to use:'; 
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
                case 'Plot Mesh'
                    [cobj,~] = selectCaseObj(obj.Cases,[],{'WRM_Mesh'},promptxt);
                    if isempty(cobj), getdialog('No mesh available'); return; end
                    hf = figure('Name','Mesh','Tag','PlotFig');
                    ax = axes(hf);
                    tabPlot(cobj,ax);
                case 'Ray Plots'                    
                    [cobj,~] = selectCaseObj(obj.Cases,[],{'RayTracks'},promptxt);
                    if isempty(cobj), getdialog('No rays available');return; end
                    hf = figure('Name','RayTracks','Tag','PlotFig');
                    ax = axes(hf);
                    tabPlot(cobj,ax,obj);
                case 'Transfer Table'
                    [cobj,~] = selectCaseObj(obj.Cases,[],{'SpectralTransfer'},promptxt);
                    if isempty(cobj), getdialog('No Transfer tables available');return; end
                    hf = figure('Name','SpecTrans','Tag','PlotFig');
                    ax = axes(hf);
                    tabPlot(cobj,ax,obj);
                case 'Transfer Coefficients'
                    [cobj,~] = selectCaseObj(obj.Cases,[],{'SpectralTransfer'},promptxt);
                    if isempty(cobj), getdialog('No Transfer tables available'); return; end
                    coefficientsPlot(cobj);
                case 'Spectrum Plots'
                    ctWaveSpectra.getPlotOption(obj);
                case 'O/I Spectrum'
                    WRM_WaveModel.runSpectrum(obj);
                case 'O/I Animation'
                    WRM_WaveModel.runAnimation(obj);
                case 'Multi-point Plots'
                    WRM_SedimentTransport.transportPlots(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,src,~)
            %menu to access online documentation and manual pdf file
            switch src.Text
                case 'Documentation'
                    doc waveraymodel   %must be name of html help file  
                case 'Manual'
                    wrm_open_manual;
            end
        end 
        %% Check that toolboxes are installed------------------------------
        function isok = check4muitoolbox(~)
            %check that dstoolbox and muitoolbox have been installed
            fname = 'dstable.m';
            dstbx = which(fname);
        
            fname = 'muiModelUI.m';
            muitbx = which(fname);
        
            if isempty(dstbx) && ~isempty(muitbx)
                warndlg('dstoolbox has not been installed')
                isok = false;
            elseif ~isempty(dstbx) && isempty(muitbx)
                warndlg('muitoolbox has not been installed')
                isok = false;
            elseif isempty(dstbx) && isempty(muitbx)
                warndlg('dstoolbox and muitoolbox have not been installed')
                isok = false;
            else
                isok = true;
            end
        end   
%% ------------------------------------------------------------------------
% Overload muiModelUI.MapTable to customise Tab display of records (if required)
%--------------------------------------------------------------------------     
        function MapTable(obj,ht)
            %create tables for Data and Model tabs - called by DrawMap
            % load case descriptions
            muicat = obj.Cases;
            caserec = find(tabSubset(obj,ht.Tag));
            caseid = muicat.Catalogue.CaseID(caserec);
            casedesc = muicat.Catalogue.CaseDescription(caserec);
            caseclass = muicat.Catalogue.CaseClass(caserec);

            %intiailise blank row of data
            if strcmp(ht.Tag,'Rays')
                headers = {'ID','Case Description','Type','Dir','Nray','Nper','Nwl'};
                cwidth = {25 230 70 60 60 60 60};
                cdata = {'0','Description of individual cases','Type','#','#','#','#'};
            elseif strcmp(ht.Tag,'Transfer')
                headers = {'ID','Case Description','Dir','Nray','Nper','Nwl'};
                cwidth = {25 250 80 70 70 70};
                cdata = {'0','Description of individual cases','#','#','#','#'};
            else
                headers = {'ID','Data Class','Data Description','Nrec', 'Start','End'};
                cwidth = {25 100 230 50 80 80};
                cdata = {'0','Type','Description of individual cases','#','#','#'};
            end
            % cdata = {'0','Type','Description of individual cases','#','#','#'};
            %ecdata = {'0','Type','Description of individual cases','#','#','xxx'};
            irec = 1;
            for i=1:length(caseid)
                case_id = num2str(caseid(i));
                if ~isfield(muicat.DataSets,caseclass{i}) || ...
                                  isempty(muicat.DataSets.(caseclass{i}))
                    type = 'New';
                else
                    type = caseclass{i};
                end
                %
                if strcmp(ht.Tag,'Rays') 
                    cobj = getCase(muicat,caserec(i));
                    dst = cobj.Data.Dataset;
                    raytype = dst.RowDescription;
                    nray = height(dst);
                    if strcmp(raytype,'Forward ray')
                        raytype = 'Forward';
                        dir = num2str(cobj.RunParam.WRM_FT_Params.dir0TN);
                    else
                        raytype = 'Backward';
                        dir = num2str(cobj.RunParam.WRM_BT_Params.DirectionRange);
                    end 
                    nper = cobj.RunParam.WRM_RunParams.nPeriod;
                    nwl  = cobj.RunParam.WRM_RunParams.nWaterLevel;
                    cdata(irec,:) = {case_id,char(casedesc{i}),raytype,dir,nray,nper,nwl};
                    
                elseif strcmp(ht.Tag,'Transfer')
                    cobj = getCase(muicat,caserec(i));
                    dst = cobj.Data.Offshore;
                    nray = height(dst);
                    dir = num2str([dst.RowRange{:}]);
                    nper = length(dst.Dimensions.Period);
                    nwl = length(dst.Dimensions.WaterLevel);
                    cdata(irec,:) = {case_id,char(casedesc{i}),dir,nray,nper,nwl};

                else
                    dst = getDataset(muicat,caserec(i),1);
                    if isempty(dst), continue; end
                    reclen = num2str(height(dst.DataTable));
                    range = dst.RowRange;
                    if ~isempty(range) && isdatetime(range{1})                 
                        stdate = char(range{1},'dd-MMM-yyyy'); %control ouput format
                        endate = char(range{2},'dd-MMM-yyyy');  
                    else
                        stdate = '';
                        endate = '';
                    end
                    cdata(irec,:) = {case_id,type,char(casedesc{i}),reclen,stdate,endate};
                    
                end
                irec = irec+1;
            end
            
            % draw table of case descriptions
            tc=uitable('Parent',ht,'Units','normalized',...
                'CellSelectionCallback',@obj.caseCallback,...
                'Tag','cstab');
            tc.ColumnName = headers;
            tc.RowName = {};
            tc.Data = cdata;
            tc.ColumnWidth = cwidth;
            tc.RowStriping = 'on';
            tc.Position(3:4)=[0.935 0.8];    %may need to use tc.Extent?
            tc.Position(2)=0.9-tc.Position(4);
        end   
%%
        function subset = tabSubset(obj,srctxt)  
            %get the cases of a given CaseType and return as logical array
            %in CoastalTools seperate data from everything else
            % srctxt - Tag for selected tab (eg src.Tag)
            % Called by MapTable. Separate function so that it can be 
            % overloaded from muiModelUI version.
            caseclass = obj.Cases.Catalogue.CaseClass;
            switch srctxt
                case 'Rays'
                    subset = contains(caseclass,'RayTracks');  
                case 'Transfer'
                    subset = contains(caseclass,'SpectralTransfer');  
                otherwise
                    subset = ~contains(caseclass,{'RayTracks','SpectralTransfer'});  
            end
        end
    end
end    
    
    
    
    
    
    
    
    
    
    