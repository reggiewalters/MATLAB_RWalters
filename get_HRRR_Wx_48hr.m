function [Precip, Temp, Snow, RelH, time] = get_HRRR_Wx_48hr(sLat, sLon)
%
% function pulls precipitation/temperature/snow data from the High-Resolution Rapid
% Refresh (HRRR) model. Resolution: 3km, hourly
% HRRR integrates to 48 hours for 00/06/12/18z UTC cycles and for 18 hours
% for all other cycles.
% this function will force a download of the most recent 48-hour cycle and
% return time, precip, 2m air temp and snowfall
% r. walters, hhwp, november 2018
% updated september 2019 to add relative humidity
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
%       Precip: surface total precipitation [inches]
%       Temp: 2 m above ground temperature [F]
%       Snow: surface total snowfall [inches]
%       RH:   2 m relative humidity [%]
%       time: hourly timestamps, corrected to local time zone
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% % %

if length(sLat) ~= length(sLon)
    disp('input latitude and longitude arrays must be of same size!');
    return
end

vars = {'apcpsfc', 'asnowsfc', 'tmp2m', 'rh2m'};    % {precip [kg/m2], snowfall [m], temp [k], RH [%]}

timeZone = 'America/Los_Angeles';           % hard-wired time zone locale
T = datetime('today','TimeZone',timeZone);
[dT,dST] = tzoffset(T); dT = hours(dT); dST = hours(dST);
if dST == 1
    dst_str = 'PDT';
else
    dst_str = 'PST';
end

date = datestr(now,'yyyymmdd');
hrC = hour(datetime(datestr(now)));
runC = hrC-dT;
dFlag = 0;
if hrC > (23+dT)
    date = datestr(now+1,'yyyymmdd');
    dFlag = 1;
end

stop_flag = 0;
runV = [18 12 6 0];
runD = runC - runV;
runI = find(runD >= 0, 1);
runC  = runV(runI);
disp('---');
disp('--- searching for most current 48-hr HRRR forecast ---');
disp('--- please wait ... ');

while stop_flag == 0
    try
        if runC < 10
            R = ['0',num2str(runC)];
        else
            R = num2str(runC);
        end
        fullpath = ['https://nomads.ncep.noaa.gov/dods/hrrr/hrrr',date,'/hrrr_sfc.t',R,'z'];
        time = ncread(fullpath,'time',[1],[Inf]);
        time = time + 365;
%         dV = datevec(time);
%         dV(:,1) = dV(:,1) + 1;
%         t = datenum(dV);
%         delta_t = mode(diff(t));
%         t = t(1) : delta_t : t(end);
%         time = t(1:length(time));
        stop_flag = 1;
    catch
        runC = runC-6;
    end
    if runC < 0 && dFlag == 0
        stop_flag = 1;
        warning('no fx file found, please check path syntax');
        return
    elseif runC < 0 && dFlag == 1
        try
            R = '18';
            date = datestr(now,'yyyymmdd');
            fullpath = ['https://nomads.ncep.noaa.gov/dods/hrrr/hrrr',date,'/hrrr_sfc.t',R,'z'];
            time = ncread(fullpath,'time',[1],[Inf]);
            time = time + 365;
%             dV = datevec(time);
%             dV(:,1) = dV(:,1) + 1;
%             t = datenum(dV);
%             delta_t = mode(diff(t));
%             t = t(1) : delta_t : t(end);
%             time = t(1:length(time));
            stop_flag = 1;
        catch
            stop_flag = 1;
            warning('no fx file found, please check path syntax');
        end
    end
    
end

valid_str = ['forecast found! valid ' R 'Z ' datestr(datenum(date,'yyyymmdd'))];
disp(valid_str)
disp('Pulling Lat/Lon Info...');
n = 0;
while n<5
    try
        lat  = ncread(fullpath,'lat',[1],[Inf]);
        lon  = ncread(fullpath,'lon',[1],[Inf]);
        n=5;
    catch
        pause(2);
        n = n+1;
    end
end
disp('Pulling Lat/Lon Info... DONE');
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
stop_flag = 0;  wait_Time = 3;   max_attempts = 7;   n=0;
disp('--- Pulling precip/temperature forecast data ... ---');
while stop_flag == 0
    n = n+1;
    disp(['attempt #: ' num2str(n)]);
    try
        
        pause(1.0);
        
        pcp = ncread(fullpath,vars{1},[iLon iLat 1],[lonN latN Inf]);
        PPt = pcp./25.4;                    % convert km/m2 to in.
        PPt = permute(PPt,[2,1,3]);         % permute to correct matrix size (lat=rows,lon=cols)
        pause(0.5);
        
        sno = ncread(fullpath,vars{2},[iLon iLat 1],[lonN latN Inf]);
        sno = sno.*39.3701;                 % convert m to in.
        sno = permute(sno,[2,1,3]);         % same permutation as above
        pause(0.5);
        
        Ta = ncread(fullpath,vars{3},[iLon iLat 1],[lonN latN Inf]);
        Ta = (Ta-273.15).*9./5+32;          % convert K to F
        Ta = permute(Ta,[2,1,3]);           % same permutation as for PPt
        
        rh = ncread(fullpath,vars{4},[iLon iLat 1],[lonN latN Inf]);
        rh = permute(rh,[2,1,3]);           % same permutation
        
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


% % % get each point's time series
for i = 1:length(sLat)
    dLat = abs(Lat - sLat(i));      dLon = abs(Lon - sLon(i));
    Ri = find(dLat==min(dLat));     Ci = find(dLon==min(dLon));
    Precip(:,i)  = PPt(Ri,Ci,:);
    Temp(:,i)    = Ta(Ri,Ci,:);
    Snow(:,i)    = sno(Ri,Ci,:);
    RelH(:,i)    = rh(Ri,Ci,:);
end

time       = time + (dT/24);

