% This function will retrieve data from the California Data Exchange Center
% (CDEC) into Matlab arrays. Requires that user knows sensor numbers and
% associated duration codes.
% R. Walters, Hetch Hetchy Water & Power, June 2018
% Updated Aug 2018 with additional error traps for sensor number alignment
%
%%% USAGE:
%   >> get_CDEC(station_ID, dur_code, sensor_Num, StartDate, EndDate)
%
%%% INPUTS:
%   'station_ID': three-letter station identification (from CDEC)
%   'dur_code':   duration code [e.g., 'd'=daily data, 'e'=event (15-minute)]
%   'sensor_Num': one- or two-digit sensor number (from CDEC)
%   'StartDate':  beginning date in the following format: mm/dd/yyyy
%   'EndDate':    ending date in same format as StartDate
%                 or enter 'now' for today's date
%
%%% OUTPUTS:
% 'Data':         data output array as Nx1 column vector
% 'date':         Matlab serial date array conciding with each 'Data' entry
%
%%% EXAMPLES:
% >> [t_pp, day] = get_CDEC('TUM', 'd', '45', '10/01/2010', '9/30/2011');
% gets daily incremental precipitation for the Tuolumne Meadows Met Station
% for the 2011 water year
%
% >> [Ta_moc, dt] = get_CDEC('mhh', 'e', '4', '10/01/2017', 'now');
% gets event (15-minute) air temperature for the Moccasin Met Station from
% the beginning of WY 2018 through the most current available entry

function [Data, date] = get_CDEC(station_ID, dur_code, sensor_Num, StartDate, EndDate)

EndDate = lower(EndDate);
if ( strncmpi('now', EndDate, 3) ) == 1
    furl = ['http://cdec.water.ca.gov/cgi-progs/queryCSV?station_id=', ...
        station_ID,'&sensor_num=',sensor_Num,'&dur_code=',dur_code, ...
        '&start_date=',datestr(StartDate,'mm/dd/yyyy'), ...
        '&end_date=Now'];
else
    furl = ['http://cdec.water.ca.gov/cgi-progs/queryCSV?station_id=', ...
        station_ID,'&dur_code=',dur_code,'&sensor_num=',sensor_Num, ...
        '&start_date=',datestr(StartDate,'mm/dd/yyyy'), ...
        '&end_date=',datestr(EndDate,'mm/dd/yyyy')];
end
fname = 'temp.csv';
try
    websave(fname,furl);  
catch
    catch_str = ['**cannot find cdec vars with specified parameters** \n', ...
        '**please check syntax or try again later** \n'];
    fprintf(catch_str);
    return
end

fid = fopen(fname);
fSpec = '%s %s %f';
C = textscan(fid,fSpec,'HeaderLines',2);
fclose(fid);    delete(fname);
C = C{1};

l1 = strsplit(C{1},',');    l2 = strsplit(C{end},',');
try
    st = datenum([l1{1} l1{2}],'yyyymmddHHMM');
    en = datenum([l2{1} l2{2}],'yyyymmddHHMM');
catch
    fail_str = ['**no timing information found - check to be sure ',upper(station_ID),' includes the specified sensor number** \n'];
    fprintf(fail_str);
    return
end

% expected dates entries (to scan for missing data)
if dur_code == 'e'
    refVec = st:(1/24/4):en;            % event (15-minute)
elseif dur_code == 'h'
    refVec = st:(1/24):en;              % hourly
elseif dur_code == 'd'
    refVec = st:(24/24):en;             % daily
end

for i = 1:length(C)
    li = strsplit(C{i},',');
    date(i) = datenum([li{1} li{2}],'yyyymmddHHMM');
    Data(i) = str2double(li{3});
end

refInds = ismembertol(refVec,date,1e-10);
newDateVec = nan(numel(refVec),1);  newDataVec = newDateVec;
newDateVec(refInds) = date;

newDataVec(refInds) = Data;

Data = newDataVec;
date = newDateVec;

nanx = isnan(date); t = 1:numel(date);
date(nanx) = interp1(t(~nanx), date(~nanx), t(nanx));