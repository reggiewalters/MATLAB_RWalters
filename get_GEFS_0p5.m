% This function will retrieve Global Ensemble Forecast System (GEFS)
% version 12.0 0.5-deg precipitation & temperature forecasts from the current day,
% cycling out 15 days, including the current date.
% The most current 00Z UTC forecast is retrieved in order to encompass a
% full day for the current day's forecast (assuming that the UTC offset is
% at least 6 hours - an ok assumption for Pacific/Mountain Time Zones).
% The function receives a lat/lon pair as input and retrieves the
% forecasted 0.5-deg grid cell point for which it is the nearest neighbor.
% User can alternatively enter a UTC date-time
% R. Walters, HHWP, Sept 2020
%
% DEPENDENCY:
% This function requires nctoolbox
% (https://github.com/nctoolbox/nctoolbox/releases)
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% *NOTE: User time zone is hard-wired and must be edited when using outside
%  of Pacific Time Zone (line xxx)
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%% USAGE:
% >> get_GEFS_0p5(Lat, Lon, varargin)
%
%%% INPUTS:
%   'Lat':      scalar or nx1 array of decimal latitude values
%   'Lon':      scalar or nx1 array of decimal longitude values
%   'varargin': variable argument input allows for date-time of chosen
%               forecast (e.g., '10-01-2019 12:00') [Time is in UTC]
%               - note: hour must be on [0, 6, 12, 18]
%
%%% OUTPUTS:
% **First 4 Outputs are Computed Daily Variables:
% **Size:   [M x N x P], where M = # of days, N = ensembles and P = forecast
%           point (one for each lat/long input pair)
%
%   'Tmax':     maximum daily air temperature   [degF]
%   'Tmin':     minimum daily air temperature   [degF]
%   'Ptot':     total daily accumulated precip  [inches]
%   'DT':       daily datetime value (12:00)
%
% **Next 3 Outputs are at the model temporal resolution (6 hours)
% **Size:   [M x N x P], where M = # time-steps (6-hour temporal
%           resolution), N = ensembles and P = forecast point
%   'T_air':    6-hourly air temperature        [degF]
%   'Precip':   6-hourly precipitation          [inches]
%   'time:'     6-hourly datetime value
%
%%% EXAMPLE:
% >> [Tmax, Tmin, PTOT, DT, pcp, ta, dt] = get_GEFS_0p5(37.945, -119.783);
%    gets daily max/min air temp, daily precip, along with 6-hourly 
%    precip and air temperature matrices for Hetch Hetchy
%    Outputs are matrices with each column representing a GEFS ensemble
%    member (column #1 is considered the GEFS control member)
%    If more than one unique forecast point (spaced by more than 1.0 deg)
%    is entered, the outputs will be datacubes of dimension-3 nLat
%
function [Tmax, Tmin, Ptot, DT, Precip, T_air, time] = get_GEFS_0p5(T_Lat, T_Lon, varargin)

Run = '00';

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

[dT, dST] = tzoffset(T);    dT = hours(dT); dST = hours(dST);
if isdst(T)
    T_offset = dT;
elseif ~isdst(T)
    T_offset = dT - dST;
end

T_offset = T_offset/24*-1;

date = datestr(now,'yyyymmdd');


if isempty(varargin)
    
    % initialize
    stop_flag = 0;
    
    % set path to GEFS 0.5-deg model
    fullpath = ['http://nomads.ncep.noaa.gov/dods/gefs/gefs' date '/gefs_pgrb2ap5_all_' Run 'z'];
    
    % % %
    wait_Time = 2.5;                % time to wait btwn ncep fetch attempts
    max_attempts = 10;              % max # of attempts
    n=0;
    while stop_flag == 0
        n=n+1;
        try
            disp(['searching for GEFS forecast for ' date ...
                ', attempt number: ' num2str(n)]);
            pause(0.25);
            time = ncread(fullpath, 'time', 1,inf) - T_offset;
            %             lat  = ncread(fullpath, 'lat', 1, inf);
            %             lon  = ncread(fullpath, 'lon', 1, inf);
            time = time + 365;
%             dV = datevec(time);
%             dV(:,1) = dV(:,1) + 1;
%             t = datenum(dV);
%             delta_t = mode(diff(t));
%             t = t(1) : delta_t : t(end);
%             time = t(1:length(time));
            disp(['forecast found!  ::::  ' datestr(time(1)+T_offset, 'mm/dd/yyyy HH:MM') ' UTC']);
            stop_flag = 1;
        catch
            pause(wait_Time);
        end
        if n > max_attempts
            disp('unable to retrieve GEFS forecast!');
            return
        end
    end
    
else
    
    stop_flag = 0;
    nAttempts = 10;  q=1;
    theDate = varargin{1};
    Run    = datestr(theDate,'HH');
    date = datestr(theDate,'yyyymmdd');
    while stop_flag == 0
        disp(['searching for GEFS forecast for ' theDate]);
        disp(['attempt #: ' num2str(q)]);
        try
            fullpath = ['http://nomads.ncep.noaa.gov/dods/gefs/gefs' date '/gefs_pgrb2ap5_all_' Run 'z'];
            time = ncread(fullpath,'time',1,inf) - T_offset;
            %             lat  = ncread(fullpath, 'lat', 1, inf);
            %             lon  = ncread(fullpath, 'lon', 1, inf);
            time = time + 365;
%             dV = datevec(time);
%             dV(:,1) = dV(:,1) + 1;
%             time = datenum(dV);
            stop_flag = 1;
            valid_str = ['forecast found! valid ' R 'Z ' datestr(datenum(date,'yyyymmdd'))];
            disp(valid_str);
        catch
            q = q+1;
        end
        if q > nAttempts-2
            disp(['unable to find GEFS forecast for ' R 'Z!']);
            return
        end
    end
    
end

% build lat/lon arrays (always identical, no need to download)
hRes = 0.5;                         % horizontal resolution (degrees)
lat = -90: hRes : 90;   lon = 0: hRes : 360-hRes;
lonInds = find(lon>180);
lon(lonInds) = lon(lonInds)-360;

T_utc = datevec(time + T_offset);    % save utc time vector
hrs   = T_utc(:,4);

% get the indices within lat/lon from input location arrays
latStart = floor(min(T_Lat));   latEnd = ceil(max(T_Lat));
lonStart = floor(min(T_Lon));   lonEnd = ceil(max(T_Lon));
latDiff = abs(lat - latStart);  iLat = find(latDiff == min(latDiff));
lonDiff = abs(lon - lonStart);  iLon = find(lonDiff == min(lonDiff));
latN    = max(1, round((latEnd - latStart)/mode(diff(lat)))) + 1;
lonN    = max(1, round((lonEnd - lonStart)/mode(diff(lon)))) + 1;

% store lat/lon vectors
Lat = lat(iLat:iLat+latN-1);
Lon = lon(iLon:iLon+lonN-1);

% get the ensemble precip data cube
% dimensions: [Lon x Lat x Time x Ens]
disp('getting precipitation & temperature ensemble array ... ');

n = 0;
while n < 10
    nEns = 31;              % GEFS v12.0 -> 31 members 
    nTs  = 65;              % 65 timesteps (6h/ea)
    try
%         Ppt = ncread(fullpath, 'apcpsfc', [iLon iLat 1 1], [latN lonN nTs nEns]);
%         Ppt = Ppt./25.4;                % convert kgm^-2 to in.
%         Ppt = permute(Ppt, [2,1,3,4]);  % permute to row(lat) x column(lon)
%         
%         disp('... got precip ...');
%         pause(0.2);
%         
%         Ta = ncread(fullpath, 'tmp2m', [iLon iLat 1 1], [latN lonN nTs nEns]);
%         Ta = (Ta-273.15).*9./5+32;      % convert from K to F
%         Ta = permute(Ta, [2,1,3,4]);    % same permutation as above
        ds  = ncgeodataset(fullpath);
        
        Ppt = ds.data('apcpsfc', [1 1 iLat iLon], [nEns nTs (iLat + latN-1) (iLon + lonN-1)]);
        Ppt = Ppt./25.4;                    % convert kgm^-2 to in.
        Ppt = permute(Ppt, [4, 3, 2, 1]);   % permute to row(lat) x column(lon) x DT x Ens
        
        disp('... got precip ...');
        pause(0.2);
        
        Ta  = ds.data('tmp2m',   [1 1 iLat iLon], [nEns nTs (iLat + latN-1) (iLon + lonN-1)]);
        Ta = (Ta-273.15).*9./5+32;          % convert from K to F
        Ta  = permute(Ta,  [4, 3, 2, 1]);   % same permutation as for precipitation
                
        n = 10;
    catch
        n = n+1;
        if n==10
            return
        end
        pause(2.5);
    end
end
disp('getting precipitation & temperature ensemble array ... DONE');

Precip = zeros(size(Ppt,3), size(Ppt,4), length(T_Lat));
T_air  = Precip;
for i = 1:length(T_Lat)
    dLat = abs(Lat - T_Lat(i));     dLon = abs(Lon - T_Lon(i));
    Ri = find(dLat==min(dLat));     Ci = find(dLon==min(dLon));
    Precip(:,:,i) = Ppt(Ri,Ci,:,:);
    T_air(:,:,i)  = Ta(Ri,Ci,:,:);
end

% get daily lapse indices
[~, ~, fx_days, fx_hrs] = datevec(time);
[~,da_inds] = unique(fx_days);
da_inds = sort(da_inds);

% % % get each point's time series and collect min/max Ta and total PPt
dfX = diff(fx_hrs);
fxTempRes = dfX(find(dfX>0,1,'first'));
% get the index that starts the first full day in the forecast
fullD1 = find(diff(da_inds) == 24/fxTempRes, 1, 'first');
Tmax = [];  Tmin = [];  Ptot = [];  DT = [];

for i = fullD1:length(da_inds)-1
    dWin    = da_inds(i):da_inds(i+1)-1;        % daily inds
    timeWin = time(dWin);
    [y,m,d] = datevec(timeWin(1));
    tmpWin  = T_air(dWin,:,:);
    pcpWin  = Precip(dWin,:,:);
    Tmax    = [Tmax; max(tmpWin)];              % daily maxima
    Tmin    = [Tmin; min(tmpWin)];              % daily minima
    Ptot    = [Ptot; nansum(pcpWin)];              % daily totals
    DT      = [DT; datenum([y m d 12 0 0])];    % forecast day
end

%     function [P, T] = gefs_member_get(date, Run, iLon, iLat, latN, lonN)
%         
%         nEns = 31;              % GEFS v12.0 -> 31 members 
%         nTs  = 65;              % 65 timesteps (6h/ea)
%         % pre-allocation for P and T arrays
%         P    = nan(latN, lonN, 65, nEns);
%         T    = P;
%         
%         for j = 1:nEns
%             dStr = ['Getting GEFS v12.0 Ensemble Member #: ' num2str(j)];
%             disp(dStr);
%             if j == 1           % control Run
%                 iPath = ['http://nomads.ncep.noaa.gov/dods/gefs/gefs' ...   
%                     date '/gec00_' Run 'z_pgrb2a'];  
%             else  
%                 iPath = ['http://nomads.ncep.noaa.gov/dods/gefs/gefs' ...
%                     date '/gep' dPad(num2str(j-1),2) '_' Run 'z_pgrb2a'];     
%             end
%             
%             Tmp = ncread(iPath, 'tmp2m', [iLon iLat 1 1], [latN lonN nTs 1]);
%             Tmp = (Tmp-273.15).*9./5+32;
%             Tmp = permute(Tmp, [2,1,3]);        % permute to row(lat) x column(lon)
%             
%             Pp  = ncread(iPath, 'apcpsfc', [iLon iLat 1 1], [latN lonN nTs 1]);
%             Pp  = Pp./25.4;                     % convert kgm^-2 to in.
%             Pp  = permute(Pp, [2,1,3]);         % same permutation as above
%             
%             T(:,:,:,j) = Tmp;                   % add as 4th (ens) dimension to hypercube
%             P(:,:,:,j) = Pp;                    % same, for precip 
%             
%         end
%     end

end
%%% End of Function
