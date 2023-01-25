function open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'AppName'));
fpath = [appinfo(idx(1)).location,[filesep,'doc',filesep,'AppName manual.pdf']];
open(fpath)
