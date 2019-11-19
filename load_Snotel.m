% This function will load in a SNOTEL .csv file and parse it into a matrix
% whose columns represent water years of the period of record (to date) and 
% whose rows represent the day of year, beginning on October 1.
% Leap year data are omitted from the time series.
% R. Walters, Idaho Power, Feb 2016.
% Updated Oct. 2016 to truncate potential zeros from end of data matrix.
% 
%%% USAGE:
%   >> load_Snotel(state_name, site_num, prec_or_swe)
%
%%% INPUTS: 
%   'state_name': two-letter state abbreviation in single quotes
%
%  'site_num': integer site number for the SNOTEL location of interest
%
%  'prec_or_swe': binary variable (0 or 1)
%                   0 = Accumulated Precip
%                   1 = Snow Water Equivalent
%
%%% OUTPUTS:
% 'DM': data matrix with columns representing water years and rows 
%       representing days of the respective year.
%       Only full water years are included with the exception of the
%       current, presumably incomplete water year.
%
% 'st_yr':  first complete water year of period of record
% 'end_yr': last year of period of record (current water year)
% 'n_yrs':  number of water years contained in period of record, including
%           the current, presumably incomplete water year.
%
%%% EXAMPLE:
% >> [DM,st_yr] = load_Snotel('ID',324,1);
% gives Period of Record output in matrix form for Bear Saddle, ID SWE
% and the starting year of the period of record.
%

function [DM,st_yr,end_yr,n_yrs] = load_Snotel(state_name, site_num, prec_or_swe)
opts = weboptions;
opts.CertificateFilename = '';
% switch for precip or SWE - go to the appropriate NRCS URL and download.
% prec_or_swe == 0 for precip and == 1 for snow water equivalent
if prec_or_swe == 0
    furl = ['https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/',num2str(site_num),':',state_name,':SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/PREC::value'];
    s = webread(furl,opts);
elseif prec_or_swe == 1
    furl = ['https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/',num2str(site_num),':',state_name,':SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/WTEQ::value'];
    s = webread(furl,opts);
end

% read s
M = textscan(s,'%s %s','Delimiter',',','EmptyValue',NaN);

% parse into separate columns (date and data)
c1 = M{1};  c2 = M{2};

% find first instance of Oct 1 in date column
s = strfind(c1,'-10-01');
s1 = find(~cellfun(@isempty,s),1);

% get the starting year for the first full water year
% (add one since water year is denoted by latter calendar year)
st_yr = c1(s1);     st_yr = cell2mat(st_yr);    st_yr = str2double(st_yr(1:4)) + 1;

% generate vector of all time series starting at first Oct 1 index s1
c1 = c1(s1:end);        c2 = c2(s1:end);

% find leap year indices
s = strfind(c1,'-02-29');
li = find(~cellfun(@isempty,s));            % leap year indices

% omit leap year instances from time series
c1(li) = [];
c2(li) = [];

% covert c2 time series from string to double-precision values
for i = 1:length(c2)
    dd(i) = str2double(c2{i});
end
dd = dd';

% get ending year and number of years
e_yr = cell2mat(c1(end));       end_yr = str2num(e_yr(1:4));
e_month = cell2mat(c1(end));    end_month = str2num(e_month(6:7)); 
if end_month>=10
    end_yr = end_yr + 1;        % ending yr might be partial water/fiscal year (e.g. after Oct 1 but before Dec 30)
end

n_yrs = end_yr - st_yr + 1;

% pre-allocate container for final data matrix 'DM'
DM = nan.*ones(365,n_yrs);

% parse data into a matrix whose columns represent the water year and whose
% rows represent the day of year, beginning on Oct-01
for i = 1:n_yrs-1
    di = dd(1+365*(i-1):365*i);
    DM(:,i) = di;
end
d_end = dd(365*i+1:end);
DM(1:length(d_end),end) = d_end;

% truncate DM to 365 rows if it happens to be larger
DM = DM(1:365,:);

