classdef WRM_SedimentTransport < muiPropertyUI & muiDataSet & matlab.mixin.Copyable 
%
%-------class help---------------------------------------------------------
% NAME

% USAGE
%   obj = SedimentTransport.setInput(mobj); %mobj is a handle to Main UI
%   obj = SedimentTransport.runModel(mobj);
% SEE ALSO
%   inherits muiPropertyUI and muiDataSet
%
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2023
%--------------------------------------------------------------------------
% 
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Upper beach slope (1:s - enter value for s)',...
                          'Bed level 1km out from SWL (mOD) [or y,z]',...
                          'Sediment grain size (m)',...
                          'Drift coefficient, Kc',...
                          'Shoreline angle (degTN)',...
                          'Distance between points (m)'};        
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed         
    end

    properties
        UpperBeachSlope = 10      %Upper beach slope (1:s - enter value for s)   
        BedLevelat1km = -8        %Bed level 1km out from SWL (mOD) [or y,z]
        GrainSize = 0.001         %Sediment grain size, D50 (m)
        DriftCoefficient = 1      %Drift coefficient, Kc
        ShorelineAngle = NaN      %angle of contours from north (degrees TN)
        PointDistances = NaN      %distance between backtracking points (m)
    end

%%   
    methods (Access=protected)
        function obj = WRM_SedimentTransport(mobj)             
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
            obj = WRM_SedimentTransport(mobj);    
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs','WRM_SedimentTransport',obj);
        end 

%%
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            promptxt = 'Select nearshore wave timeseries to use';
            wavobj = selectCaseObj(mobj.Cases,[],{'WRM_WaveModel'},promptxt); 
            if isempty(wavobj), return; end
            dst = wavobj.Data;
            pntnames = fieldnames(dst);
            npnts = length(pntnames);
%--------------------------------------------------------------------------
% Model code 
%-------------------------------------------------------------------------- 
            obj = WRM_SedimentTransport.setInput(mobj); %set class properties
            if isnan(obj.ShorelineAngle)
                obj = getShoreDir(obj);
                if isnan(obj.ShorelineAngle), return; end
            else
                obj.ShorelineAngle = repmat(obj.ShorelineAngle(1),1,npnts);
                obj.PointDistances = (0:1:npnts)*obj.PointDistances(1);
            end         
            %
            driftmodel = WRM_SedimentTransport.getDriftMethod();
            if isempty(driftmodel), return; end
            setRunParam(obj,mobj);  

            %run model and add results to a struct array of dstables
            dst = getDrift(obj,mobj,dst,driftmodel);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------  
            %save results
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete'); 
        end

%%
        function transportPlots(mobj)
            %derived plots from the sediment transport data at multiple alongshore
            %points
            promptxt = 'Sediment Tranport Plots - select Case to use:';
            listxt = {'Annual Mean Drift','Monthly Mean Drift',...
                      'Summary Point Drift','Summary Shore Drift',...
                      'Monthly Peclet Ratio','Cluster Peclet Ratio',...
                      'Wave-Drift Tables'};
                ok = 1;
                while ok>0
                    [obj,~] = selectCaseObj(mobj.Cases,[],{'WRM_SedimentTransport'},promptxt);
                    if isempty(obj), ok = 0; continue; end
                    selection = listdlg("ListString",listxt,"PromptString",...
                                        'Select option:','SelectionMode','single',...
                                        'ListSize',[150,200],'Name','EDBtools');
                    if isempty(selection), ok = 0; continue; end
                    wrm_transport_plots(obj,mobj,listxt{selection});
                end
        end
    end

%%
    methods
        function dst = getDrift(obj,mobj,dst,driftmodel)
            %calculate the littoral drift, drift gradient and x-shore 
            %transport for all points along the shore 
            pntnames = fieldnames(dst);
            npnts = length(pntnames);
            %properties for bed slope within surf zone (half depth of inshore wave point)  
            ubs = obj.UpperBeachSlope;
            z1km = obj.BedLevelat1km; 
            d50 = obj.GrainSize;
            Kc = obj.DriftCoefficient;

            %call drift model and add longhshore drift to wave time series
            g = mobj.Constants.Gravity;
            rhw = mobj.Constants.WaterDensity;
            rhs = mobj.Constants.SedimentDensity;
            vsc = mobj.Constants.KinematicViscosity;

            %initialise metatdata
            wv = dst.(pntnames{1});
            zi = mean((wv.swl-wv.depi),'omitnan');
            theta = obj.ShorelineAngle;
            mtime = wv.RowNames;
            nrec = length(mtime);
            dsp = WRM_SedimentTransport.setDSproperties();

            hw = waitbar(0,'Processing point 0');
            Qs = zeros(nrec,npnts); Qx = Qs; dQdx = Qs;
            for i=1:npnts                
                wv = dst.(pntnames{i});
                bs = profileslope(wv.depi/2,wv.swl,z1km,ubs); %first argument is depth
                Qall = littoraldrift(wv.Hsi,wv.Tpi,wv.Diri,wv.depi,...
                                            theta(i),bs,d50,0.0006,g,rhs,rhw,vsc);   
                Qs(:,i) = Qall(:,driftmodel.id)*Kc;
                % dQdx(:,1) = diff(Qs(:,1))/obj.PointDistance(i);
                % dQdx = [dQdx;dQdx(end,1)];    %pad to make same length as Qs
                %Note the current formulation dose NOT use Tp, Dir and theta
                Qx(:,i) = xshore_bailard(wv.Hsi,wv.Tpi,wv.Diri,wv.depi,...
                                            theta(i),bs,d50,g,rhw,rhs,vsc);
            end
            
            for j=1:nrec
                if any(Qs(j,:)~=0) %only call if Qs~=0
                    dQdx(j,:) = central_differences(Qs(j,:),obj.PointDistances);
                end
            end
            
            %save results to dstables - one for each point
            for i=1:npnts   
                waitbar(i/npnts,hw,sprintf('Processing point %d',i));
                adst = dstable(Qs(:,i),dQdx(:,i),Qx(:,i),'RowNames',mtime,'DSproperties',dsp);
                srctxt = sprintf('Class %s, at %s',metaclass(obj).Name,pntnames{i});
                adst.Source = srctxt;   
                bs = mean(bs,'omitnan');
                mtxt1 = sprintf('Drift using %s; Theta=%g; d50=%g; Kc=%g; Beach slope=1:%.1f; Zi=%g',...
                                    driftmodel.name,theta(i),d50,Kc,bs,zi);    
                mtxt2 = sprintf('Using %s case for wave input',wv.Description);
                inptxt = sprintf('%s\n%s',mtxt1,mtxt2);
                adst.MetaData = inptxt;
                dst.(pntnames{i}) = adst;
            end
            
            delete(hw)  
        end
%%
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end

%%
        function obj = getShoreDir(obj)
            %load the backtracking point coordinates and calculate
            %shoreline angle
            [fname,path,nfiles] = getfiles('MultiSelect','off','FileType','*.txt;',...
                                           'PromptText','Select backtracking start points file:');
            if nfiles<1, return; end
            xypnts = readmatrix([path,fname]); 
            plines = gd_lines2points([xypnts;NaN,NaN]);
            [plinedir,~,cumlen] = gd_curvelineprops({plines},1);
            plinedir = [plinedir{1}(1),plinedir{1}(1:end-1),plinedir{1}(end-1)];
            for j=2:length(plinedir)-1
                theta(j-1) = sum(plinedir(j-1:j+1))/3*180/pi(); %#ok<AGROW>
            end
            obj.ShorelineAngle = theta;
            % dist = diff(cumlen{1}(1:end-1)); 
            % obj.PointDistance = [dist,dist(end)];   %pad to make same length  
            obj.PointDistances = cumlen{1}(1:end-1);  %remove traling NaN
        end

    end


    methods (Static, Access=private)


%%
        function driftmodel = getDriftMethod()
            %prompt user for transport formul to use
            tlist = {'CERC formula (SPM, 1994)',...
                     'Dynamics of Marine Sands, Soulsby',...
                     'Kamphuis formula',...
                     'Damgaard & Soulsby (shingle)'};
            [driftmodel.id,ok] = listdlg('Name','Plot profile', ...
                                 'PromptString','Select formula', ...
                                 'ListSize',[200,80], ...
                                 'SelectionMode','single', ...
                                 'ListString',tlist);
            if ok==0, driftmodel = []; return; end 
            driftmodel.name = tlist{driftmodel.id};
        end

%%
        function dsprop = setDSproperties(~)
            %define the variables in the dataset
            %define the metadata properties for the demo data set
            dsprop= struct('Variables',[],'Row',[],'Dimensions',[]); 
            
            %struct entries are cell arrays and can be column or row vectors
            dsprop.Variables = struct(...
                'Name',{'Qs','dQdx','Qx'},...
                'Description',{'Alongshore drift rate potential',...
                               'Alongshore drift gradient'...
                               'Cross-shore transport rate'},...
                'Unit',{'m^3/s','m^3/s/m','m^3/s'},...
                'Label',{'Transport rate (m^3/s)','Flux gradient (m^3/s/m)',...
                                               'Transport Rate (m^3/s)'},...
                'QCflag',repmat({'model'},1,3));
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
end