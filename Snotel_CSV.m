% Snotel_CSV.m
% This function will load in daily SNOTEL precip and SWE data, parsing into 
% matrices whose columns represent water years of the period of record and
% whose rows represent the day of year, beginning on October 1.
% Leap year data are omitted from the time series.
%
% The function will produce a .csv file for use in MS Excel with two
% separate sheets for precipitation and SWE, respectively. The .csv file
% is given a unique name based on the site name and download date. The 
% output file will be placed in the same working directory from which the 
% function is called.
%
% Reggie Walters, Idaho Power, April 2016.
% Updated Oct. 2016 to truncate potential zeros from end of data matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USAGE:
%   >> Data = load_Snotel(state_name, site_num)
%
%%% INPUTS:
%  'state_name': two-letter state abbreviation in single quotes
%
%  'site_num':   integer site number for the SNOTEL location of interest
%
%%% OUTPUTS:
% 'Site_Name_MonthDay_Year.csv': 
%       .csv file for use in MS Excel containing Precipitation and SWE in
%       two separate sheets. This file will be written to the current
%       working directory.
%
% 'Data': Matlab structure containing data matrices with columns 
%       representing water years and rows
%       representing days of the respective year.
%       Only full water years are included with the exception of the
%       current, potentially incomplete water year.
%
%%% EXAMPLE:
% >> [Data] = load_Snotel('ID',324);
%    gives Period of Record output in matrix form for Bear Saddle, ID SWE
%    and writes a .csv file to the current working directory.
%

function Data = Snotel_CSV(state_name, site_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  PRECIPITATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go to the appropriate NRCS URL and download precip:

furl = ['https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/',num2str(site_num),':',state_name,':SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/PREC::value'];
file_name = 'temp.csv';
urlwrite(furl,file_name);

% read the temp.csv file
fid = fopen(file_name);
M = textscan(fid,'%s %s','Delimiter',',','EmptyValue',NaN);
fclose(fid);
delete(file_name);

% print status
fprintf('Formatting Precipitation Data\n');

% parse data into separate columns (date and data)
c1 = M{1};  c2 = M{2};

%%%%%%%%%%%%%%% Header Information routine %%%%%%%%%%%%%
sF = 0; i=0;
while sF==0
        i = i+1;
    cN = str2num(c1{i});            % try to convert to numeric
    if ~isempty(cN)                 % test if numeric
        sF=1;
    end
end

% get the site info from the headers
hStop = i-1;                        % get ending row of headers data
Headers = M{1}; Headers = Headers(1:hStop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find first instance of Oct 1 in date column
s = strfind(c1,'-10-01');

% get first index where the Oct 1 date occurs
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

% get ending year and produce a year vector
e_yr = cell2mat(c1(end));       end_yr = str2num(e_yr(1:4));
e_month = cell2mat(c1(end));    end_month = str2num(e_month(6:7)); 
if end_month>=10
    end_yr = end_yr + 1;        % ending yr might be partial water/fiscal year (e.g. after Oct 1 but before Dec 30)
end

n_yrs = end_yr - st_yr + 1;     yr_vec = st_yr:end_yr;

% pre-allocate container for final data matrix 'DM'
DM = nan.*ones(365,n_yrs);

% parse data into a matrix whose columns represent the water year and whose
% rows represent the day of year, beginning on Oct-01
for i = 1:n_yrs-1
    di = dd(1+365*(i-1):365*i)';
    DM(:,i) = di;
end
d_end = dd(365*i+1:end);
DM(1:length(d_end),end) = d_end;

% Error trap to delete columns composed completely of NaN values:
stop_flag = 0;                              % flag for while loop
i = 1;                                      % initialize while iterator
while stop_flag == 0
    dCol = DM(:,i);                         % get column i
    if any(dCol)==0                         % check for all NaN's
        i = i+1;                            % if so, increment i
    else                                    % if not, then exit
        stop_flag = 1;          
    end
end
DM = DM(:,i:end);                           % truncate if necessary
yr_vec = yr_vec(i:end);                     % truncate year vector, if necessary
st_yr = st_yr + i-1;                        % increment start year, if necessary                        

% Insert header row for year
DM = [yr_vec;DM];

DM = DM(1:366,:);

% store the precip data into data structure
Data.Precip = DM;

clearvars -except Data site_num state_name Headers hStop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  SNOW WATER EQUIVALENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go to the appropriate NRCS URL and download SWE:
furl = ['https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/',num2str(site_num),':',state_name,':SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/WTEQ::value'];
file_name = 'temp.csv';
urlwrite(furl,file_name);

% read the temp.csv file
fid = fopen(file_name);
M = textscan(fid,'%s %s','Delimiter',',','EmptyValue',NaN);
fclose(fid);
delete(file_name);

% print status
fprintf('Formatting Snow Water Equivalent Data\n');

% parse into separate columns (date and data)
c1 = M{1};  c2 = M{2};

% find first instance of Oct 1 in date column
s = strfind(c1,'-10-01');

% get first index where the Oct 1 date occurs
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

% convert c2 time series from string to double-precision values
for i = 1:length(c2)
    dd(i) = str2double(c2{i});
end
dd = dd';

% get ending year and produce a year vector
e_yr = cell2mat(c1(end));       end_yr = str2num(e_yr(1:4));
e_month = cell2mat(c1(end));    end_month = str2num(e_month(6:7)); 
if end_month>=10
    end_yr = end_yr + 1;        % ending yr might be partial water/fiscal year (e.g. after Oct 1 but before Dec 30)
end

n_yrs = end_yr - st_yr + 1;     yr_vec = st_yr:end_yr;

% pre-allocate container for final data matrix 'DM'
DM = nan.*ones(365,n_yrs);

% parse data into a matrix whose columns represent the water year and whose
% rows represent the day of year, beginning on Oct-01
for i = 1:n_yrs-1
    di = dd(1+365*(i-1):365*i)';
    DM(:,i) = di;
end
d_end = dd(365*i+1:end);
DM(1:length(d_end),end) = d_end;

% Error trap to delete columns composed completely of NaN values:
stop_flag = 0;                              % flag for while loop
i = 1;                                      % initialize while iterator
while stop_flag == 0
    dCol = DM(:,i);                         % get column i
    if any(dCol)==0                         % check for all NaN's
        i = i+1;                            % if so, increment i
    else                                    % if not, then exit
        stop_flag = 1;          
    end
end
DM = DM(:,i:end);                           % truncate if necessary
yr_vec = yr_vec(i:end);                     % truncate year vector, if necessary
st_yr = st_yr + i-1;                        % increment start year, if necessary    

%DM = DM(1:365,:);                           % truncate if necessary

% Insert header row for year
DM = [yr_vec;DM];

DM = DM(1:366,:);

% store the SWE data into data structure
Data.SWE = DM;

clearvars -except Data end_yr site_num state_name Headers hStop;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DATA WRITING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write data to appropriate sheets in .csv (using xlswrite)

%build a string for the file name
h1 = Headers{hStop-6};    hE = find(h1 == '(');
sN = h1 (3:hE-2);   hB = find(sN == ' ');
hD = Headers{hStop-2};    hE = find(hD == ':',1);
hD = hD(hE+2:end);
if any(hB) sN(hB)='_'; end
file_str = [sN,'_',hD(1:3),hD(5:end),'_',num2str(end_yr)];
fName = [file_str,'.csv'];

% build a string for dates column
d_st = datenum([00 10 01 0 0 0]);   d_en = datenum([01 09 30 0 0 0]);
dS = d_st:d_en;
for i = 1:length(dS)    dd = datestr(dS(i));    dC{i} = dd(1:6);    end
dC = dC';   dC = [{''};dC];

% generate blank cells to make headers lines identical dimensions to DM
blks_Precip = repmat({''},length(Headers(end-6:end)),size(Data.Precip,2));
blks_SWE    = repmat({''},length(Headers(end-6:end)),size(Data.SWE,2));

% error trap to make sure file is not currently existing and open
check_Open = ddeinit('Excel',fName);
if check_Open ~=0
    disp(['WARNING: Failed to write file -- "',fName,'" is open. Please close for writing/overwriting']);
    % Truncate the Matlab Structure to remove the year vectors
    Data.Precip = Data.Precip(2:end,:);
    Data.SWE    = Data.SWE(2:end,:);
    return
end

% print status
disp(['Writing data to "',pwd,'\',fName,'"'])

% write the header lines, blanks, and data matrices to the .csv (xls) file
warning('off','MATLAB:xlswrite:AddSheet');
xlswrite(fName,[Headers(end-6:end) blks_Precip; dC num2cell(Data.Precip)] ,'PRECIP');
xlswrite(fName,[Headers(end-6:end) blks_SWE; dC num2cell(Data.SWE)],'SWE');

% set working path and build the Excel File objects
fPath = pwd;
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(fullfile(fPath, fName));
% Delete the default sheets 1-3:
try     % Throw an error if the sheets do not exist.
      objExcel.ActiveWorkbook.Worksheets.Item('Sheet1').Delete;
      objExcel.ActiveWorkbook.Worksheets.Item('Sheet2').Delete;
      objExcel.ActiveWorkbook.Worksheets.Item('Sheet3').Delete;
catch
       % Do nothing otherwise
end

% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;

% Truncate the Matlab Structure to remove the year vectors
Data.Precip = Data.Precip(2:end,:);
Data.SWE    = Data.SWE(2:end,:);

fprintf('Complete!\n');

end