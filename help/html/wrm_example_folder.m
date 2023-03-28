function wrm_example_folder()
%find the location of the example folder and open it
appname = 'WaveRayModel';
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},appname));
path = appinfo(idx(1)).location;
if isfolder([path,filesep,appname])
    %Matlab installs the App as a subfolder of the App folder if there
    %are folders included that are on the same level (ie not subfolders)
    path = [path,filesep,appname];
end

fpath = [path,[filesep,'example']];
try
    winopen(fpath)
catch
    msg = sprintf('The examples can be found here:\n%s',fpath);
    msgbox(msg)
end