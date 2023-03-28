function wrm_open_manual()
%find the location of the asmita app and open the manual
appname = 'WaveRayModel';
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'WaveRayModel'));
path = appinfo(idx(1)).location;
if isfolder([path,filesep,appname])
    %Matlab installs the App as a subfolder of the App folder if there
    %are folders included that are on the same level (ie not subfolders)
    path = [path,filesep,appname];
end

fpath = [path,[filesep,'doc',filesep,'WaveRayModel manual.pdf']];
open(fpath)
