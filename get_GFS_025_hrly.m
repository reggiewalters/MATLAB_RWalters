% This function will retrieve Global Forecast System (GFS) 0.25 degree
% HOURLY
% temperature and precipitation data from the current day
%
% The function receives a lat/lon pair as input and retrieves the
% forecasted 0.25 grid cell point for which it is the nearest neighbor.
% R. Walters, HHWP, Jan 2019
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% *NOTE: User time zone is hard-wired and must be edited when using outside
%  of Pacific Time Zone (line 51)
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%% USAGE:
% >> get_GFS_025_hrly(Lat, Lon, varargin)
%
%%% INPUTS:
%   'Lat':      scalar or nx1 array of decimal latitude values
%   'Lon':      scalar or nx1 array of decimal longitude values
%
%%% OUTPUTS:
%   'T_air':    hourly air temperature        [degF]
%   'Precip':   hourly precipitation          [inches]
%   'time:'     hourly datetime value
%
%   The function allows for an additional set of lat/lon pairs to be
%   entered as variable argument inputs, to be used for precip outputs that
%   are discretely different than those of the temperature points. This
%   functionality could easily be achieved by separate calls to get_GFS_025
%   with differing site lat/lon. However, this comes at the expense of
%   additional calls to ncread.m which can add computing time.
%
%   EXAMPLE:
% >> [Tmp, Ppt, dt] = get_GFS_025_hrly(37.975, -119.916);
% >> plot(dt,Tmp, 'o-', 'lineWidth',2);    datetick('x','ddd');
% >> grid on;   set(gca,'fontSize',14); ylabel('T_a [\circF]');
% >> yyaxis right;  bar(dt,Ppt);    ylabel('Precip [in]');
% retrieves and plots forecasted hourly air temperature and
% for the GFS 1/4-deg grid cell containing the Cherry
% Valley Met (CVM) weather station
%
%
function [T_air, Precip, time] = get_GFS_025_hrly(T_Lat, T_Lon, varargin)

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
stop_flag = 0;

hrC = datevec(now); hrC = hrC(4);
runC = hrC+7;
if hrC > 17
    date = datestr(now+1,'yyyymmdd');
end

stop_flag = 0;
runV = [18 12 6 0];
runD = runC - runV;
runI = find(runD >= 0, 1);
runC  = runV(runI);
disp('--- checking for most current GFS hourly forecast ---');
disp('--- please wait ...');

while stop_flag == 0
    try
        R = dPad(num2str(runC),2);
        fullpath = ['http://nomads.ncep.noaa.gov/dods/gfs_0p25_1hr/gfs',date,'/gfs_0p25_1hr_',R,'z'];
        time = ncread(fullpath,'time',1,inf) + 365 - T_offset;
        stop_flag = 1;
        valid_str = ['forecast found! valid ' R 'Z ' datestr(datenum(date,'yyyymmdd'))];
        disp(valid_str);
    catch
        runC = runC-6;
    end
    if runC < 0
        stop_flag = 1;
        warming('no fx file found ...');
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
max_attempts = 5;
disp('getting Wx data matrices ...');

n = 0;  stop_flag = 0;
while stop_flag == 0
    n = n+1;
    disp(['--- air temperature attempt #' num2str(n) ' ---']);
    try
        Ta = ncread(fullpath, 'tmp2m', [iLon iLat 1], [lonN latN Inf]);
        Ta = (Ta-273.15).*9./5+32;  % convert from K to F
        Ta = permute(Ta,[2,1,3]);   % permute to correct matrix size (lat=rows,lon=cols)
        disp('--- air temperature retrieved ---');
        stop_flag = 1;
    catch
        disp('---');
    end
    if n>=max_attempts
        stop_flag = 1;
        disp('unable to get air temp!');
    end
end

n = 0;  stop_flag = 0;
while stop_flag == 0
    n = n+1;
    disp(['--- precip attempt #' num2str(n) ' ---']);
    try
        PPt = ncread(fullpath, 'apcpsfc', [iLon iLat 1], [lonN latN Inf]);
        PPt = PPt./25.4;            % convert kgm^-2 to in.
        PPt = permute(PPt,[2,1,3]); % same permutation as for Ta
        disp('--- precip retrieved ---');
        stop_flag = 1;
    catch
        disp('---');
    end
    if n>=max_attempts
        stop_flag = 1;
        disp('unable to get precip!');
    end
end


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

hrInds = 1:6:length(hrs);
for i = 1:length(input_LatP)
    dLat = abs(Lat - input_LatP(i));    dLon = abs(Lon - input_LonT(i));
    Ri = find(dLat==min(dLat));         Ci = find(dLon==min(dLon));
    p =PPt(Ri,Ci,:);
    p = p(:);
    
    %     % every hr divisible by 6 is accum. precip for t-6
    P = p;  P(2:end) = zeros(1,length(P)-1);
    
    for q = 1:length(hrInds)-1
        iWin = hrInds(q)+1:hrInds(q+1);
        pWin = p(iWin);
        pcp(1) = pWin(1);
        for k = 2:length(pWin)
            pcp(k) = pWin(k) - pWin(k-1);
        end
        P(iWin) = pcp;
    end
    
    Precip(:,i) = P;
    
end
disp('--- done! ---');

end
