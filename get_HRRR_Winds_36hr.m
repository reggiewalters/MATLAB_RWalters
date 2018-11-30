function [wind_speed, wind_gust, time] = get_HRRR_Winds_36hr(sLat, sLon)
%
% function pulls wind speed and gust data from the High-Resolution Rapid
% Refresh (HRRR) model. Resolution: 3km, hourly
% HRRR integrates to 36 hours for 00/06/12/18z UTC cycles and for 18 hours 
% for all other cycles.
% this function will force a download of the most recent 36-hour cycle and
% return time, wind speed, and maximum gust arrays
% r. walters, hhwp, november 2018
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% INPUT:
%       sLat: site latitude in decimal degrees, singular or array
%             minimum: 21.1405 deg, maximum: 52.6133 deg)
%
%       sLon: site longitude in decimal degrees (negative values W of 
%             Greenwich), same dimension as sLat
%             minimum: -134.0955 deg, maximum: -60.9365 deg
%
%       * lat/lon pairs must be within given ranges
%
% % %
% OUTPUT:
%       wind_speed in miles per hour
%       wind_gust in miles per hour
%       time: hourly timestamps, corrected to local time zone
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% % %

if length(sLat) ~= length(sLon)
    disp('input latitude and longitude arrays must be of same size!');
    return
end
    
vars = {'wind10m', 'gustsfc'};
conv = 3600/.3048/5280;                     % conversion factor (meters/sec to miles/hr)

timeZone = 'America/Los_Angeles';           % hard-wired time zone locale
T = datetime('today','TimeZone',timeZone);
[dT,dST] = tzoffset(T); dT = hours(dT); dST = hours(dST);
if dST == 1
    dst_str = 'PDT';
else
    dst_str = 'PST';
end

date = datestr(now,'yyyymmdd');
hrC = datevec(now); hrC = hrC(4);
runC = hrC-dT;
if hrC > (23+dT)
    date = datestr(now+1,'yyyymmdd');
end

stop_flag = 0;
runV = [18 12 6 0];
runD = runC - runV;
runI = find(runD >= 0, 1);
runC  = runV(runI);
disp('---');
disp('--- checking for most current 36-hr HRRR forecast ---');
disp('--- please wait ... ');

while stop_flag == 0
    try
        if runC < 10
            R = ['0',num2str(runC)];
        else
            R = num2str(runC);
        end
        fullpath = ['http://nomads.ncep.noaa.gov/dods/hrrr/hrrr',date,'/hrrr_sfc.t',R,'z'];
        time = ncread(fullpath,'time',[1],[Inf])+365;
        stop_flag = 1;
    catch
        runC = runC-6;
    end
    if runC < 0
        stop_flag = 1;
        warning('no fx file found, please check path syntax');
    end
end

valid_str = ['forecast found! valid ' R 'Z ' datestr(datenum(date,'yyyymmdd'))];
disp(valid_str)
disp('Pulling Lat/Lon/Time Info...');
lat  = ncread(fullpath,'lat',[1],[Inf]);
lon  = ncread(fullpath,'lon',[1],[Inf]);
disp('Pulling Lat/Lon/Time Info... DONE');
disp('---');

% get the min/max lat/lon values to get an output matrix spanning the range
latStart = floor(min(sLat)*10)/10;
lonStart = floor(min(sLon)*10)/10;
latEnd   = ceil(max(sLat)*10)/10;
lonEnd   = ceil(max(sLon)*10)/10;
dfLat = diff(lat);   dfLat = dfLat(end);
dfLon = diff(lon);   dfLon = dfLon(end);
latN     = ceil((latEnd-latStart)/dfLat);
lonN     = ceil((lonEnd-lonStart)/dfLon);
latDiff = abs(lat - latStart);  iLat = find(latDiff == min(latDiff));
lonDiff = abs(lon - lonStart);  iLon = find(lonDiff == min(lonDiff));

% store longitude/latitude vectors corresponding to the above indices
Lat = lat(iLat:iLat+latN-1);
Lon = lon(iLon:iLon+lonN-1);

%=========== read the data ================================================
stop_flag = 0;  wait_Time = 4;   max_attempts = 5;   n=0;
disp('--- Pulling wind forecast data ... ---');
while stop_flag == 0
    disp(['attempt #: ' num2str(n+1)]);
    try      
        pause(0.5);
        wnd = ncread(fullpath,vars{1},[iLon iLat 1],[lonN latN Inf]);   
        pause(0.5);
        gst = ncread(fullpath,vars{2},[iLon iLat 1],[lonN latN Inf]); 
        disp('--- forecast retrieved! ---');
        stop_flag = 1;
    catch
        pause(wait_Time);
    end
    if n> max_attempts
        disp('** unable to retrieve HRRR forecast! **');
        return
    end
end

wnd = permute(wnd,[2,1,3]);         % permute such that matrix is rxc
gst = permute(gst,[2,1,3]);         % permute such that matrix is rxc

% % % get each point's time series
for i = 1:length(sLat)
    dLat = abs(Lat - sLat(i));      dLon = abs(Lon - sLon(i));
    Ri = find(dLat==min(dLat));     Ci = find(dLon==min(dLon));
    g(:,i)  = gst(Ri,Ci,:);         w(:,i) = wnd(Ri,Ci,:);
end

wind_speed = w * conv;
wind_gust  = g * conv;
time       = time + (dT/24);


function o = dPad(i, d)
%DPAD Appends '0's to an input string in order to reach length d
o = i;
    if numel(o) < d
        while numel(o) < d
            o = ['0' o];
        end
    end
end

end
