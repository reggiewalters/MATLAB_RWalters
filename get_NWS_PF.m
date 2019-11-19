% This function will retrieve National Weather Service Hourly Point
% forecast information. Requires coordinate pair (lat/long) inputs.
% Typical forecast length is 7 days
% R. Walters, HHWP, Oct 2019
% Updated Nov 18, 2019 to adapt to SFPUC firewall changes and associated
%       compatibility issues with Matlab certificate
%
%
%%% USAGE:
%   >> get_NWS_PF(Lat, Lon)
%
%%% INPUTS:
%   latitude & longitude: numeric decimal lat/lon pair(s), can be entered
%   as Nx2 array of coordinate points (compute time increases linearly)
%
%%% OUTPUTS:
%   'Ta':       Hourly air temperature (deg F)
%   'precip':   Hourly quantitative precipitation forecast (inches)
%   'RH':       Relative Humidity (%)
%   'WS':       Hourly Wind Speed (mph)
%   'GS':       Hourly Wind Gust (mph)
%   'dt':       Matlab serial date array coinciding with each hourly data point
%
%   'Tmax':     Maximum daily air temperature
%   'Tmin':     Minimum daily air temperature
%   'Ptot':     Daily qpf 
%   'DT':       Daily serial date coinciding with each daily data point (12:00 time stamped)
%   
%%% EXAMPLE:
% >> [T_air, PPT, ~, Wind_spd, ~, dt] = get_NWS_PF(38.200, -119.983)
% gets hourly air temp, precip, wind speed and date array for a coordinate
% point corresponding to Pincrest Met Station

function [Ta, precip, RH, WS, GS, dt, Tmax, Tmin, Ptot, DT] = get_NWS_PF(Lat, Lon)

fx_Len = 168;   % standard forecast look-ahead window

% pre-allocate arrays
precip = nan.*ones(fx_Len,length(Lat));
Ta= precip;  WS = precip; GS = precip; RH = precip; dt = precip(:,1);

for j = 1:length(Lat)
    
    j/length(Lat)
    stopFlag = 0;
    n = 0;
    while stopFlag == 0
        try
            fURL = ['https://forecast.weather.gov/MapClick.php?lat=' num2str(Lat(j)) '&lon=' num2str(Lon(j))  '&FcstType=digitalDWML'];
            s = webread(fURL, weboptions('CertificateFilename', ''));
            stopFlag = 1;
        catch
            n = n+1;
            if n>10
                stopFlag = 1;
                disp(['Unable to fetch data for location ' num2str(j)]);
            end
        end
    end
    
    if j == 1
        % get time array
        dt = extractBetween(s,'<start-valid-time>','</end-valid-time>');
        dt = extractBetween(dt,'','<'); dt = dt(:,1);
        dt = datenum(dt, 'yyyy-mm-ddTHH:MM');
    end
    
    % qpf
    pp = extractBetween(s,'<hourly-qpf','</hourly-qpf>');
    pp = cellfun(@str2double,extractBetween(pp,'<value>','</value>'));
    precip(1:length(pp),j) = pp;
    
    % get air temp
    tt = extractBetween(s,'<temperature type="hou','</temperature>');
    tt = cellfun(@str2double,extractBetween(tt,'<value>','</value>'));
    Ta(:,j) = tt;
    
    % get wind speed and gust
    ws = extractBetween(s,'<wind-speed type="sust','</wind-speed>');
    ws = cellfun(@str2double,extractBetween(ws,'<value>','</value>'));
    WS(:,j) = ws;
    
    wg = extractBetween(s,'<wind-speed type="gust','</wind-speed>');
    wg = extractBetween(wg,'<value','/');
    WG = nan.*ones(size(ws));
    for i = 1:length(wg)
        if strfind(wg{i},'nil')
        else
            WG(i) = cellfun(@str2double,regexp(wg{i},'\d*','Match'));
        end
    end
    GS(:,j) = WG;
    
    % get rh
    rh = extractBetween(s,'<humidity type="relative"','</humidity>');
    rh = cellfun(@str2double,extractBetween(rh,'<value>','</value>'));
    RH(1:length(rh),j) = rh;
    
end

% % % compute daily stats/totals
% get daily lapse indices
[~, ~, fx_days, fx_hrs] = datevec(dt);
[~,da_inds] = unique(fx_days);
da_inds = sort(da_inds);

% % % get each point's time series and collect min/max Ta and total PPt
dfX = diff(fx_hrs);
fxTempRes = dfX(find(dfX>0,1,'first'));
% get the index that starts the first full day in the forecast
fullD1 = find(diff(da_inds) == 24/fxTempRes, 1, 'first');
Tmax = [];  Tmin = [];  Ptot = [];  DT = [];

for i = 1:length(da_inds)-1
    dWin    = da_inds(i):da_inds(i+1)-1;        % daily inds
    timeWin = dt(dWin);
    [y,m,d] = datevec(timeWin(1));
    tmpWin  = Ta(dWin,:,:);
    pcpWin  = precip(dWin,:,:);
    Tmax    = [Tmax; max(tmpWin)];              % daily maxima
    Tmin    = [Tmin; min(tmpWin)];              % daily minima
    Ptot    = [Ptot; nansum(pcpWin)];           % daily totals
    DT      = [DT; datenum([y m d 12 0 0])];    % forecast day
end

end

