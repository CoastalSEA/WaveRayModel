function wrm_example_folder()
%find the location of the example folder and open it
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'WaveRayModel'));
fpath = [appinfo(idx(1)).location,[filesep,'example']];
try
    winopen(fpath)
catch
    msg = sprintf('The examples can be found here:\n%s',fpath);
    msgbox(msg)
end