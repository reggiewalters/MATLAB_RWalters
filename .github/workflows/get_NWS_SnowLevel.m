function [SL, DT, sL, dt] = get_NWS_SnowLevel(Lat, Lon)
% function retrieves snow level forecast information from NOAA WRH tabular
% web table, given a site latitude and longitude
% r. walters, apr 2020, hhwp
% re-written  mar 2021 to retrieve data via NOAA API (deprecated data
%                                   discovery via WRH tabular arrays)
%
%%%
% INPUTS:
%   Lat:    site latitude [decimal degrees]
%   Lon:    site longitude [decimal degrees]
%
% OUTPUTS:
%   SL:     daily snow level [feet]
%   DT:     daily time array (matlab numeric datetime)
%   sL:     hourly snow level [feet], interpolated from 3-hourly after approx. 48 hours
%   dt:     hourly time array (matlab numeric datetime)
%%%
% EXAMPLE:
%   >> [SnowLine, day, snowline, time] = get_NWS_SnowLevel(37.975, -119.916);
%   >> figure(11);  clf;    hold on;
%   >> plot(day, SnowLine, '-o', 'lineWidth',2);
%   >> plot(time, snowline, 'lineWidth',1.5);
%   >> datetick;   grid on;
% The above plots daily average and hourly snow level forecast for Cherry Valley
%%%
%
opts = weboptions('CertificateFilename', '', 'ContentType', 'auto');
api  = 'https://api.weather.gov/points/';
nLoc = length(Lat);
snowLvar = 'snowLevel';
dtFmt    = 'yyyy-mm-ddTHH:MM:SS';
m2ft     = 3.28084;                     % meters to feet conversion
nAtt     = 5;                           % # attempts at each URL before breaking t/c
ki = 0;
if length(Lat)>1 && length(Lat)<1e2
    ki = 1;
    nL = length(Lat);
end

% first iteration, get the starting time stamp and initial snow level array
i = 1;

n = 0;
while n < nAtt
    
    try
        url = [api num2str(Lat(i)) ',' num2str(Lon(i))];
        S   = webread(url, opts);
        % extract the gridId, gridX, and gridY data for the lat/lon pair
        gid = S.properties.gridId;  % forecast office
        gX  = S.properties.gridX;
        gY  = S.properties.gridY;
        
        % reset api/url definitions for x/y grid identifiers
        API = 'https://api.weather.gov/gridpoints/';
        URL = [API gid '/' num2str(gX), ',', num2str(gY)];
        D       = webread(URL, opts);
        Props   = D.properties;
        
        evalStr = ['struct2cell(Props.' snowLvar '.values)'];
        Dat     = eval(evalStr);
        
        tt = datenum(Dat(1,:), dtFmt);
        sl = (cell2mat(Dat(2,:)) .* m2ft)';
        
        tmp_dt = tt(1) : 1/24 : tt(end);         % template time array
        
        slev    = interp1(tt, sl, tmp_dt);
        sL(:,i) = slev';
        
        break
        
    catch
        n = n+1;
    end
    
end

if n >= nAtt
    disp('snow level retrieval failure, please check matlab routine or nws api');
    return
end

sL = [sL nan(size(sL,1), nLoc-1)];

if nLoc>1
    
    for i = 2:nLoc
        
        n = 0;
        while n < nAtt
            
            try
                url = [api num2str(Lat(i)) ',' num2str(Lon(i))];
                S   = webread(url, opts);
                % extract the gridId, gridX, and gridY data for the lat/lon pair
                gid = S.properties.gridId;  % forecast office
                gX  = S.properties.gridX;
                gY  = S.properties.gridY;
                
                % reset api/url definitions for x/y grid identifiers
                API = 'https://api.weather.gov/gridpoints/';
                URL = [API gid '/' num2str(gX), ',', num2str(gY)];
                D       = webread(URL, opts);
                Props   = D.properties;
                
                evalStr = ['struct2cell(Props.' snowLvar '.values)'];
                Dat     = eval(evalStr);
                
                tt = datenum(Dat(1,:), dtFmt);
                sl = (cell2mat(Dat(2,:)) .* m2ft)';
                
                slev    = interp1(tt, sl, tmp_dt);
                sL(:,i) = slev';
                
                break
                
            catch
                n = n+1;
            end
            
        end
        
        if ki
            disp([num2str(round(i/nL*100)) ' % complete']);
        end
        
    end
    
end

dt = tmp_dt;


% sample daily arrays and time stamps
for i = 1:nLoc
    
    [x, t] = downsample_ts(sL(:,i), tmp_dt, 'day');
    
    % truncate last (partial) day, append to SL array
    SL(:,i) = x(1:end-1);
end

DT = floor(t(1:end-1)) + 0.5;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---  previous code archive  ---
%
% %%
% DT = nan(7, length(Lat));
% SL = DT;
% sl = nan(7*4, length(Lat));
% dt = sl;
% sInd = 1e3;
% ki = 0;
% if length(Lat)>1 && length(Lat)<1e2
%     ki = 1;
%     nL = length(Lat);
% end
% for k = 1:length(Lat)
%     if ki
%         disp([num2str(round(k/nL*100)) ' % complete']);
%     end
%     fURL = ['https://www.wrh.noaa.gov/forecast/wxtables/index.php?' ...
%         'lat=' num2str(Lat(k)) '&lon=' num2str(Lon(k)) ...
%         '&table=custom&duration=10&interval=12'];
%     opts = weboptions;
%     opts.CertificateFilename = '';
%     opts.Timeout = 15;
%     sF = 0; n = 0;
%     while sF == 0
%         try
%             s  = webread(fURL, opts);
%             if contains(s, 'Snow Level')
%                 sF = 1;
%             else
%                 n = n+1;
%                 % some times a negligibly small lat/lon change will bring up a
%                 % proper URL string (shouldn't change the HRAP pixel loc)
%                 fURL = ['https://www.wrh.noaa.gov/forecast/wxtables/index.php?' ...
%                         'lat=' num2str(Lat(k) + randn*.001) '&lon=' num2str(Lon(k) + randn*.001) ...
%                         '&table=custom&duration=10&interval=12'];
%                     if n>10
%                         disp('Unable to retrieve snow level!');
%                         return
%                     end
%             end
%         catch
%             n = n+1;
%             % some times a negligibly small lat/lon change will bring up a
%             % proper URL string (shouldn't change the HRAP pixel loc)
%             fURL = ['https://www.wrh.noaa.gov/forecast/wxtables/index.php?' ...
%                         'lat=' num2str(Lat(k) + randn*.001) '&lon=' num2str(Lon(k) + randn*.001) ...
%                         '&table=custom&duration=10&interval=12'];
%             if n>10
%                 disp('Unable to retrieve snow level!');
%                 return
%             end
%         end
%     end
%
%     % retrieve snow level
%     dta = strfind(s, 'Snow Level');
%     sdt = s(dta:dta+3200);
%     b1 = '<td bgcolor="#d8e2f2" class="tdbody">';
%     b2 = '</td>';
%     sLev = cellfun(@str2num,extractBetween(sdt, b1, b2));
%
%     % get time array
%     nDays = round(length(sLev)/28)*28 / 4;
%     dtStart = datevec(floor(now));
%     dta = strfind(s, '6-Hour');
%     sdt = s(dta:dta+3200);
%     b1 = '<td align="center" class="hhead">';
%     t = extractBetween(sdt, b1, b2);
%     noon_inds = find(contains(t, 'Noon'));
%     mdnt_inds = find(contains(t, 'Mdnt'));
%     t(noon_inds) = {'12pm'};
%     t(mdnt_inds) = {'12am'};
%
%     for i = 1:length(t)
%         T(i) = sscanf(t{i}(2:end), '%d');
%     end
%     deltaT = abs(mode(diff(T)));
%     AorP = t{1}(end-1);
%     if strcmpi(AorP,'p')
%         hOffset = 12;
%     else
%         hOffset = 0;
%     end
%     dtStart(4) = T(1)+hOffset;
%     dtt = datenum(dtStart) : deltaT/24 : datenum(dtStart) + ceil(length(sLev)/4)+1;
%
%     exp_sL_len = nDays * 4;
%     sL_diff    = exp_sL_len - length(sLev);
%     dt = dtt(sL_diff+1 : length(sLev)+sL_diff);
%     dt = dt(:);
%
%     % construct daily arrays
%     [S, D] = downsample_ts(sLev, dt, 'day');
%     DT(1:length(D),k) = floor(D) + 0.5;
%     SL(1:length(S),k) = S;
%
%     sL(1:length(sLev),k) = sLev;
%     sind = find(isnan(sL(:,k)),1);
%     if sind < sInd
%         sInd = sind;
%     end
% end
%
% DT = DT(:,1);
% try
%     sL = sL(1:sInd-1, :);
% catch
%     sL = sL(1:end, :);
% end

