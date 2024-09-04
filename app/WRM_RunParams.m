classdef WRM_RunParams < muiPropertyUI        
%
%-------class help---------------------------------------------------------
% NAME
%   WRM_RunParams.m
% PURPOSE
%   Class for run parameters to the WaveRayModel
% NOTE
%   
% USAGE
%   obj = WRMrunparams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Range of Wave Periods (s)',...
                          'Number of Wave Periods (-)',...
                          'Range of Water Levels (mOD)',...
                          'Number of Waver Levels (-)',...
                          'Ray cutoff depth (m)',...
                          'Shoreline angle (degTN)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        %Minimum of one direction, one period and one water level
        %required for model to run  
        PeriodRange = [4 18]        %Range of Wave Periods (s)
        nPeriod = 8                 %Number of Wave Periods
        WaterLevelRange = [-1 1]    %Range of Water Levels (mOD)
        nWaterLevel = 3             %Number of Waver Levels
        hCutOff = 0.1               %wave ray cut-off water depth (m)
        ShorelineAngle              %angle of contours from north (degrees TN)
    end    

%%   
    methods (Access=protected)
        function obj = WRM_RunParams(mobj)             
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'WRM_RunParams';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = WRM_RunParams(mobj);       
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end   
end