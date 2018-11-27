% This function will retrieve Global Forecast System (GFS) 0.25 degree
% temperature and precipitation data from the current day, looking out 
% 9 days into the future, including the current date. 
% The most current 06Z UTC forecast is retrieved in order to encompass a
% full day for the current day's forecast (assuming that the UTC offset is
% at least 6 hours - an ok assumption for Pacific/Mountain Time Zones). 
% The function receives a lat/lon pair as input and retrieves the
% forecasted 0.25 grid cell point for which it is the nearest neighbor.
% R. Walters, HHWP, Sept 2018
% Updated September 4, 2018 to include variable argument input
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% *NOTE: User time zone is hard-wired and must be edited when using outside
%  of Pacific Time Zone (line 51)
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%% USAGE:
% >> get_GFS_025(Lat, Lon, T_offset, varargin)
%
%%% INPUTS:
%   'Lat':      scalar or nx1 array of decimal latitude values
%   'Lon':      scalar or nx1 array of decimal longitude values
%
%%% OUTPUTS:
% **First 4 Outputs are Computed Daily Variables:
%   'Tmax':     maximum daily air temperature   [degF]
%   'Tmin':     minimum daily air temperature   [degF]
%   'Ptot':     total daily accumulated precip  [inches]
%   'DT':       daily datetime value (12:00)
% **Next 3 Outputs are at the model temporal resolution (3 hours)
%   'T_air':    3-hourly air temperature        [degF]
%   'Precip':   3-hourly precipitation          [inches]
%   'time:'     3-hourly datetime value
%
%   The function allows for an additional set of lat/lon pairs to be
%   entered as variable argument inputs, to be used for precip outputs that
%   are discretely different than those of the temperature points. This
%   functionality could easily be achieved by separate calls to get_GFS_025
%   with differing site lat/lon. However, this comes at the expense of
%   additional calls to ncread.m which can add computing time.
%
%   EXAMPLE:
% >> [Tmax, ~, Ptot, DT] = get_GFS_025(37.975, -119.916);
% >> plot(DT,Tmax, 'o-', 'lineWidth',2);    datetick('x','ddd');
% >> grid on;   set(gca,'fontSize',14); ylabel('T_a [\circF]');
% retrieves and plots forecasted daily maximum air temperature and 
% accumulated precip for the GFS 1/4-deg grid cell containing the Cherry 
% Valley Met (CVM) weather station 
%
%
function [Tmax,Tmin,Ptot,DT, T_air, Precip, time] = get_GFS_025(T_Lat, T_Lon, varargin)

% set time offset, dependent upon timezone/DST
timeZone = 'America/Los_Angeles';
T = datetime('today','TimeZone',timeZone);

if iscolumn(T_Lat)
    input_LatT = T_Lat;
    input_LonT = T_Lon;
elseif ~iscolumn(T_Lat)
    input_LatT = T_Lat';                % force column vector
    input_LonT = T_Lon';
end

if nargin>2
    input_LatP = varargin{1};
    input_LonP = varargin{2};
else
    input_LatP = T_Lat;
    input_LonP = T_Lon;
end

[dT, dST] = tzoffset(T);    dT = hours(dT); dST = hours(dST);
if isdst(T)
    T_offset = dT;
elseif ~isdst(T)
    T_offset = dT - dST;
end

T_offset = T_offset/24*-1;

date = datestr(now,'yyyymmdd');

% initialize
stop_flag = 0;          R = '06';   % 06Z forecast        
                             
% set path to GFS 1/4deg model
fullpath = ['http://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs',date,'/gfs_0p25_',R,'z'];

% % %
wait_Time = 5;              % time to wait btwn ncep fetch attempts 
max_attempts = 10;          % max # of attempts
n=0;
while stop_flag == 0
    n=n+1;
    try
        disp(['searching for forecast for ' date ...
              ', attempt number: ' num2str(n)]);     
        pause(1);
        time = ncread(fullpath,'time',[1],[Inf])+365 - T_offset;
        disp('forecast found!');
        stop_flag = 1;
    catch
        pause(wait_Time);
    end
    if n > max_attempts
        disp('unable to retrieve GFS forecast!');
        return                
    end
end
T_utc = datevec(time + T_offset);    % save utc time vector
hrs   = T_utc(:,4);

% build lat/lon arrays (always identical, no need to download)
lat = -90: 0.25 : 90;   lon = 0: 0.25 : 360-0.25;
lonInds = find(lon>180);    
lon(lonInds) = lon(lonInds)-360;

% get the indices within lat/lon from station location information arrays
latStart = floor(min([input_LatT; input_LatP])*2)/2;
lonStart = floor(min([input_LonT; input_LonP])*2)/2;
latEnd   = ceil(max([input_LatT; input_LatP])*2)/2;
lonEnd   = ceil(max([input_LonT; input_LonP])*2)/2;
latN     = (latEnd-latStart)/0.25 + 1;
lonN     = (lonEnd-lonStart)/0.25 + 1;
latDiff = abs(lat-latStart);    iLat = find(latDiff == min(latDiff));
lonDiff = abs(lon-lonStart);    iLon = find(lonDiff == min(lonDiff));

% store longitude/latitude vectors corresponding to the above indices
Lat = lat(iLat:iLat+latN-1);
Lon = lon(iLon:iLon+lonN-1);

% get T_air and Precip matrices
disp('getting Wx data matrices ...');
Ta = ncread(fullpath, 'tmp2m', [iLon iLat 1], [lonN latN Inf]);
Ta = (Ta-273.15).*9./5+32;  % convert from K to F
Ta = permute(Ta,[2,1,3]);   % permute to correct matrix size (lat=rows,lon=cols)

PPt = ncread(fullpath, 'apcpsfc', [iLon iLat 1], [lonN latN Inf]);
PPt = PPt./25.4;            % convert kgm^-2 to in.
PPt = permute(PPt,[2,1,3]); % same permutation as for Ta
disp('done!');

% get daily lapse indices
[~, ~, fx_days, fx_hrs] = datevec(time);                
da_inds = [1];
for i = 1:length(fx_days)-1
    if fx_days(i+1) ~= fx_days(i)
        % da_inds becomes starting index of each new day in forecast
        da_inds = [da_inds i+1];
    end
end

% % % get each point's time series and collect min/max Ta and total PPt
for i = 1:length(input_LatT)
    dLat = abs(Lat - input_LatT(i));    dLon = abs(Lon - input_LonT(i));
    Ri = find(dLat==min(dLat));         Ci = find(dLon==min(dLon));
    t = Ta(Ri,Ci,:);    
    T_air(:,i) = t(:);
end

for i = 1:length(input_LatP)
    dLat = abs(Lat - input_LatP(i));    dLon = abs(Lon - input_LonT(i));
    Ri = find(dLat==min(dLat));         Ci = find(dLon==min(dLon));
    p =PPt(Ri,Ci,:);
    p = p(:);
    % every hr divisible by 6 is accum. precip for t-6, whereas hrs
    % divisible by 3 (but not 6) are for t-3. difference such that all
    % points are scaled to t-3 accumulation periods:
    for q = 3:2:length(hrs)
        p(q) = p(q) - p(q-1);
    end
    Precip(:,i) = p;
end

dfX = diff(fx_hrs);
fxTempRes = dfX(find(dfX>0,1,'first'));
% get the index that starts the first full day in the forecast
fullD1 = find(diff(da_inds) == 24/fxTempRes, 1, 'first');
Tmax = [];  Tmin = [];  Ptot = [];  DT = [];
for i = fullD1:length(da_inds)-1
    dWin    = da_inds(i):da_inds(i+1)-1;        % daily inds
    timeWin = time(dWin);   
    [y,m,d] = datevec(timeWin(1));
    tmpWin  = T_air(dWin,:);
    pcpWin  = Precip(dWin,:);
    Tmax    = [Tmax; max(tmpWin)];              % daily maxima
    Tmin    = [Tmin; min(tmpWin)];              % daily minima
    Ptot    = [Ptot; sum(pcpWin)];              % daily totals
    DT      = [DT; datenum([y m d 12 0 0])];    % forecast day
end
end
