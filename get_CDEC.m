% This function will retrieve data from the California Data Exchange Center
% (CDEC) into Matlab arrays. Requires that user knows sensor numbers and
% associated duration codes.
% R. Walters, HHWP, June 2018
% with snippets borrowed from R. Picklum
% Updated Jul 2018 with additional error traps for sensor number alignment
% Updated Aug 27 2018 to account for CDEC url changes (dynamicapp)
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
% grabs cdec data and fills missing data with NaN values
% r. walters, hetch hetchy water and power, july 2018
% all inputs in single quotes. use 'now' for EndDate to run thru current


if ( strncmpi('now', EndDate, 3) ) == 1
    furl = ['http://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=', ...
        station_ID,'&SensorNums=',sensor_Num,'&dur_code=',dur_code, ...
        '&Start=',datestr(StartDate,'yyyy-mm-dd'), ...
        '&end_date=Now'];
else
    furl = ['http://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=', ...
        station_ID,'&SensorNums=',sensor_Num,'&dur_code=',dur_code, ...
        '&Start=',datestr(StartDate,'yyyy-mm-dd'), ...
        '&end_date=',datestr(EndDate,'yyyy-mm-dd')];
end

catch_str = ['**cannot find cdec vars with specified parameters** \n', ...
        '**please check syntax or try again later** \n'];
try
    s = urlread(furl);
catch
    fprintf(catch_str)
    return
end

if length(s) < 100
    fprintf(catch_str)
    return
end

express = ['\w{3},' upper(dur_code) ',[^,]*,[^,]*,([^,]*),[^,]*,([^,]*)'];
tok = regexp(s, express, 'tokens');

date = zeros(1, numel(tok));    Data = date;

for i = 1:length(tok)
    date(i) = datenum(tok{i}{1});
    Data(i) = str2double(tok{i}{2});
end

