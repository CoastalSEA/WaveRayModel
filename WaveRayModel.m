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
        vNumber = '1.0'
        vDate   = 'Jan 2023'
        modelName = 'WaveRayModel'                       
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = WaveRayModel             
            %constructor function initialises GUI
            obj = setMUI(obj);             
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus    
            %classes required to run model, format:
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            %                                        % << Edit to model and input parameters classnames 
            obj.ModelInputs.RayTracks = {'WRM_RunParams'};
            obj.ModelInputs.WRM_Bathy = {'GD_GridProps'};
            obj.ModelInputs.SpectralTransfer = {''};
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
            menu.Setup(1).List = {'Input Data','Grid Parameters','Grid Tools',...
                                  'Run Parameters','Model Constants'};                                    
            menu.Setup(1).Callback = [{'gcbo;'},{@obj.setupMenuOptions},...
                              {'gcbo;'},{'gcbo;'},{@obj.setupMenuOptions}];
            %add separators to menu list (optional - default is off)
            menu.Setup(1).Separator = [repmat({'off'},[1,4]),{'on'}]; %separator preceeds item
            
            % submenu for Import Data (if these are changed need to edit
            % loadMenuOptions to be match)
            menu.Setup(2).List = {'Gridded Data','Timeseries'};
            menu.Setup(2).Callback = {'gcbo;','gcbo;'};
            % submenu for Gridded and Timeseries Data 
            nitems = 2;
            for j=1:nitems  %add standard submenu to all import menu items
                menu.Setup(j+2).List = {'Load','Add','Delete','Quality Control'};                                   
                menu.Setup(j+2).Callback = repmat({@obj.loadMenuOptions},[1,4]);
            end
            % submenu for Grid Tools
            menu.Setup(5).List = {'Translate Grid','Rotate Grid',...
                                  'Re-Grid','Sub-Grid',...
                                  'Combine Grids','Add Surface',...
                                  'To curvilinear','From curvilinear',... 
                                  'Display Dimensions','Difference Plot',...
                                  'Plot Sections','Digitise Line',...
                                  'Export xyz Grid','User Function'};                                                                         
            menu.Setup(5).Callback = repmat({@obj.gridMenuOptions},[1,14]);
            menu.Setup(5).Separator = [repmat({'off'},[1,6]),...
                             {'on','off','on','off','off','on','on','on'}]; %separator preceeds item  
            % submenu for Run Parameters
            menu.Setup(6).List = {'Run Conditions','Forward Tracking','Backward Tracking'};
            menu.Setup(6).Callback = repmat({@obj.setupMenuOptions},[1,3]);
            menu.Setup(6).Separator = repmat({'off'},[1,3]);

            
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Check Start Points','Forward Rays',...
                                'Backward Rays','Spectral Transfer',...
                                              'Test Grid','Derive Output'};
            menucall = repmat({@obj.runMenuOptions},[1,6]);
            menucall{4} = 'gcbo;';
            menu.Run(1).Callback = menucall;
            menu.Run(1).Separator = [repmat({'off'},[1,4]),{'on','on'}];

            % submenu for Spectral Transfer
            menu.Run(2).List = {'Transfer Table','Run Timeseries'};
            menu.Run(2).Callback = repmat({@obj.runMenuOptions},[1,2]);
            menu.Run(2).Separator = repmat({'off'},[1,2]);
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics','Ray Plots','Spectral Plots'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,4]);
            menu.Run(2).Separator = {'off','off','on','off'};
            
            %% Help menu --------------------------------------------------
            menu.Help(1).Callback = {@obj.Help}; %make model specific?
            
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
            tabs.Cases  = {'   Cases  ',@obj.refresh};        % << Edit tabs to suit model 
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
                'WRM_Bathy','Inputs',[0.35,0.97],{160,70}, 'Bathymetry parameters:'};
        end    
 %%
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            switch src.Tag                                  
                case 'Plot' 
                     tabPlot(cobj,src,obj);
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
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
            end
        end  
        %%
        function gridMenuOptions(obj,src,~)
            %callback functions for grid tools options
            gridclasses = {'GD_ImportData'}; %add other classes if needed
            GD_ImportData.gridMenuOptions(obj,src,gridclasses);
            DrawMap(obj);
        end
%%
        function loadMenuOptions(obj,src,~)
            %callback functions to import timeseries data
            switch src.Parent.Text
                case 'Gridded Data'
                    classname = 'GD_ImportData';
                case 'Timeseries'
                    classname = 'muiUserData';
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
            muicat = obj.Cases;   %handle to muiCatalogue
            promptxt = 'Select a Case to use:'; 
            gridclasses = {'GD_ImportData'};
            switch src.Text     
                case 'Check Start Points'
                    RayTracks.checkStart(obj);
                case 'Forward Rays'
                    RayTracks.runModel(obj,src); 
                case 'Backward Rays'                    
                    RayTracks.runModel(obj,src);
                case 'Transfer Table'
%                     [cobj,classrec] = selectCaseObj(muicat,[],gridclasses,promptxt);
%                     if isempty(cobj), return; end
%                     getUserTools(cobj,classrec,muicat);
                    
                    SpectralTransfer.runModel(obj);
                case 'Run Timeseries'
                    SpectralTransfer.runWaves(obj);
                case 'Test Grid'
                    WRM_Bathy.runModel(obj);
                case 'Derive Output'
                    obj.mUI.Manip = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,~,~)
            doc modelskill                             
        end
%% ------------------------------------------------------------------------
% Overload muiModelUI.MapTable to customise Tab display of records (if required)
%--------------------------------------------------------------------------     
%         function MapTable(obj,ht)
%             %create tables for Record display tabs - called by DrawMap
%             % ht - tab handle
%         end
    end
end    
    
    
    
    
    
    
    
    
    
    