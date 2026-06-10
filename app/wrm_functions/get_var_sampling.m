function selection = get_var_sampling(inpsel,isstats)    
%
%-------function help------------------------------------------------------
% NAME
%   get_var_sampling.m
% PURPOSE
%   UI to prompt user to select Seasons, Direction Sampling interval and
%   recurrence period for binning data
% USAGE
%   selection = get_var_sampling(isstats);
% INPUTS
%   inpsel - input selection: default is {4,2,1,1} = month, year, all, all
%   isstats - logical include statistical methods selection if true
%             (default is false)
% OUTPUT
%   selection - vector array of user selection indices for
%               1 - Season; 2 - Direction; 3 - bin interval; 4 - bin period
% SEE ALSO
%   binned_variable.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2026
%--------------------------------------------------------------------------
%   
    if nargin<1, inpsel = {4,2,1,1}; isstats = false; 
    elseif nargin<2, isstats = false; end

    if ~iscell(inpsel), inpsel = num2cell(inpsel); end

    lst1 = {'All','Year', 'Quarter', 'Month', 'Week', 'Day', 'Hour'};
    lst2 = {'All','Winters','Summers'};
    lst3 = {'All','+ve only','-ve only'};
    deflists = {lst1,lst1,lst2,lst3};
    fields = {'Interval (bin size):','Recurrence period:','Seasons:','Direction:'};
    if isstats
        statlist = {'All','sum','mean','std','peclet','mode','median','95pct','5pct'};
        deflists = [deflists,statlist];
    end

    ok = 0;
    while ok<1
        varargin =  {'FigureTitle','Select variables',...
                     'PromptText','Select any subsampling of data required',...
                     'InputFields',fields,...
                     'InputOrder',{'',''},...
                     'Style',repmat({'popupmenu'},1,numel(fields)),...
                     'ControlButtons',{},...                            
                     'DefaultInputs',deflists,...
                     'UserData',inpsel,...%for popupmenu cell array used to set initial values
                     'DataObject',[],...
                     'SelectedVar',{},...
                     'ActionButtons',{'Select','Close'},...
                     'Position',[]};
        selection = inputUI.getUI(varargin{:});
        if isempty(selection) || selection{1}>=selection{2}
            ok = 1; 
        else
            getdialog('Bin interval must be shorter than or equal to Recurruence period')
        end
    end  
end
