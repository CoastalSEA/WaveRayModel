function wrm_runmovie(obj,src,~)
%
%-------function help------------------------------------------------------
% NAME
%   wrm_runmovie.m
% PURPOSE
%   callback function for animation figure buttons and slider
%   modified from muiPlots to handle two subplots
% USAGE
%   wrm_runmovie(obj,src,~)
% INPUTS
%   obj - muiPlots class object
%   src - callback handle
%   ~   - unused event handle
% OUTPUT
%   animation figure
%
% Author: Ian Townend
% CoastalSEA (c) March 2023
%--------------------------------------------------------------------------
%
    hfig = src.Parent;
    idm = hfig.Number==[obj.ModelMovie{:,1}];
    if strcmp(src.Tag,'runMovie')       %user pressed run button
        if license('test', 'Image_Processing_Toolbox')   %tests whether product is licensed (returns 1 if it is)
            implay(obj.ModelMovie{idm,2});
        else
            hmf = figure('Name','Animation', 'Units','normalized', ...
            'Resize','on','HandleVisibility','on','Visible','on',...
            'Position',[0.38,0.42,0.30,0.42],'Tag','PlotFig');
            if ~obj.MetaData
                hmf.Position(3) = 0.6;   %wider figure for polar plot
            end
            movie(hmf,obj.ModelMovie{idm,2});
        end
    elseif strcmp(src.Tag,'saveMovie')  %user pressed save button 
        saveanimation2file(obj.ModelMovie{idm,2});
    else                                %user moved slider
        val = ceil(src.Value);          %slider value 
        %get figure axis, extract variable and refresh plot                
        s1 = findobj(hfig,'Tag','PlotFigAxes1'); 
        s2 = findobj(hfig,'Tag','PlotFigAxes2');                                                  
        hp1 = findobj(s1.Children,'Type','surface');
        hp2 = findobj(s2.Children,'Type','surface');

        var = s1.UserData.Z;       
        var1 = squeeze(var.off(val,:,:)); %#ok<NASGU> 
        refreshdata(hp1,'caller')
        var2 = squeeze(var.in(val,:,:)); %#ok<NASGU> 
        refreshdata(hp2,'caller') 

        if obj.MetaData             %XY plot - add markers
            hp3 = findobj(s1.Children,'Tag','DirPk');
            hp4 = findobj(s2.Children,'Tag','DirPk');

            xyz = s1.UserData.xyz;
            var3 = xyz.off(val,1); %#ok<NASGU>
            var4 = xyz.off(val,2); %#ok<NASGU>
            refreshdata(hp3,'caller')
            var5 = xyz.in(val,1); %#ok<NASGU>
            var6 = xyz.in(val,2); %#ok<NASGU>
            refreshdata(hp4,'caller')
        end

        %update title
        time = s1.UserData.T(val);   %time slice selected
        w = obj.Data.Waves;
        sg = findobj(s1.Parent.Children,'Tag','PlotFigTitle');
        sg.String = sprintf('%s \nTime = %s, Hs=%.3g; Tp=%.3g; Dir=%.3g\n',...
            obj.Title,string(time),w(val,1),w(val,2),w(val,3));
        drawnow;
        %update slider selection text
        stxt = findobj(hfig,'Tag','FrameTime');
        stxt.String = string(time);
    end
end