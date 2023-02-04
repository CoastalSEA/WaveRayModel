classdef WRM_BT_Params < muiPropertyUI        
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
        PropertyLabels = {'BT Range of Directions (degTN)',...
                          'BT Number of Directions (-)',...
                          'BT start point [x,y]'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        %FT used in Forward Tracking model and BT used in Backward Tracking
        %model. Minimum of one direction, one period and one water level
        %required for model to run
        DirectionRange = [90 270]   %BT Range of Direction (degTN)
        nDirections = 18            %BT Number of Directions    
        StartPoint = [1 1]          %BT start point [x,y]
    end    

%%   
    methods (Access=protected)
        function obj = WRM_BT_Params(mobj)             
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
            classname = 'WRM_BT_Params';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = WRM_BT_Params(mobj);       
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