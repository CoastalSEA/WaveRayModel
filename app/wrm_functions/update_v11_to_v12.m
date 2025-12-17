function update_v11_to_v12(obj)
%
%-------header-------------------------------------------------------------
% NAME
%   update_v11_to_v12.m 
% PURPOSE
%   update saved models from v1.1 to v1.2.
% USAGE
%   update_v11_to_v12(obj)
% INPUTS
%   obj - instance of model
% RESULTS
%   saved model updated from v1.1 to v1.2
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
% To use from command line, open WaveRayModel using:
% >>wm = WaveRayModel;     open project file and then call
% >>update_v11_to_v12(wm)
%
% Author: Ian Townend
% CoastalSEA (c) Dec 2025 
%--------------------------------------------------------------------------
%

    %update label for WRM_RunParams properties
     if isfield(obj.Inputs,'WRM_RunParams')
         obj.Inputs.WRM_RunParams.PropertyLabels = {'Range of Wave Periods (s)',...
                          'Number of Wave Periods (-)',...
                          'Range of Water Levels (mOD)',...
                          'Number of Waver Levels (-)',...
                          'Ray cut-off depth (m)',...
                          'Distance tolerance',...
                          'Shoreline angle (degTN)'};
     end
     getdialog('Project updated from v1.1 to v1.2')
end