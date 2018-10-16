% Function get_USGS.mat retrieves station data from USGS Nat'l Water
% Information System (NWIS) and stores into Matlab arrays. Requires that
% user knows station ID number for station of interest. As of current
% version, the function only retrieves elevation and discharge/storage
% data.
% R. Walters, HHWP, Aug 2018
% with snippet(s) from R. Picklum
% Updated 10/16/2018 to fix time array issue
%
%%% USAGE:
%   >> get_USGS(siteNum, StartDate, EndDate)
%
%%% INPUTS:
%   'siteNum':      8-digit USGS station number
%   'StartDate':    beginning date
%   'EndDate':      ending date
%
%%% OUTPUTS:
% 'elev':           river gage height or reservoir elevation [ft]
% 'QV':             river discharge [cfs] or reservoir storage [acre-ft]
% 'dates':          Matlab serial date array
% 'qual':           Data-value qualification code (Approved/Provisional)
%
%%% EXAMPLE:
%   >> [Stage, Flow, dT, qual] = get_USGS(11276500, '10/01/2016', '9/30/2017');
%   gets gage height and discharge for USGS Tuolumne R nr Hetch Hetchy and
%   stores a corresponing dates vector along with a quality code array.
%
function [elev, QV, dates, qual] = get_USGS(siteNum, StartDate, EndDate)

if ~ischar(siteNum)     siteNum = num2str(siteNum);     end
if strcmpi(EndDate, 'now') == 1     EndDate = now;      end
    
dFormat  = 'yyyy-mm-dd';
startStr = datestr(StartDate, dFormat);
endStr   = datestr(EndDate, dFormat);


% River Case:
cb1 = '00060';    % discharge [cfs]
cb2 = '00065';    % stage [ft]
furl = ['https://nwis.waterdata.usgs.gov/ca/nwis/uv/?cb_' cb1 '=on&&cb_' ...
    cb2 '=on&format=rdb&period&begin_date=' startStr '&end_date=' endStr ...
    '&site_no=' siteNum];
a = urlread(furl);
if numel(a) > 250
    
    hLines = numel(strfind(a, '#')) + 2;    % number of headers
    s = textscan(a, '%s %s %s %s %s %f %s %f %s', 'headerLines', hLines);
    
% Reservoir Case:
else
    cb1 = '00054';    % storage [acre-ft]
    cb2 = '62614';    % elevation [ft]
    furl = ['https://nwis.waterdata.usgs.gov/ca/nwis/uv/?cb_' cb1 '=on&&cb_' ...
    cb2 '=on&format=rdb&period&begin_date=' startStr '&end_date=' endStr ...
    '&site_no=' siteNum];
    a = urlread(furl);
    hLines = numel(strfind(a, '#')) + 2;    % number of headers
    s = textscan(a, '%s %d %s %s %s %f %s %f %s', 'headerLines', hLines);
end

% timing vector array
allDates = s{3};    allTimes = s{4};
DM = [allDates allTimes];
dates = zeros(1,size(DM,1));
for i = 1:size(DM,1)
    d = strjoin(DM(i,:));
    dates(i) = datenum(d);
end
dates = dates';

% store output arrays
QV =    s{6};
elev =  s{8};
qual =  s{end};
