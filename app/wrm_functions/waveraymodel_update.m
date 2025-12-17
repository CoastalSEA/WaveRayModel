function waveraymodel_update(obj,oldV,newV)
%
%-------header-------------------------------------------------------------
% NAME
%   waveraymodel_update.m 
% PURPOSE
%   update saved models to newer versions of WaveRayModel
% USAGE
%   waveraymodel_update(oldV,newV) 
% INPUTS
%   obj - instance of model
%   oldV - old version number as a character string
%   newV - new version number as a character string
% RESULTS
%   saved model updated to new version. If this is called from WaveRayModel
%   this will be the version that is being run at the time.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
%
% Author: Ian Townend
% CoastalSEA (c) Dec 2025 
%--------------------------------------------------------------------------
%
    if str2double(oldV)<1.2
        update_v11_to_v12(obj);
    else
        warndlg(sprintf('No update for version %s to version %s', oldV,newV))
    end
end