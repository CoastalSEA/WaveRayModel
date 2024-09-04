function wrm_open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'WaveRayModel'));
fpath = [appinfo(idx(1)).location,[filesep,'WaveRayModel',filesep,'app',...
                        filesep,'doc',filesep,'WaveRayModel manual.pdf']];
open(fpath)
