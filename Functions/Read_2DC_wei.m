% function [SD_tot_max, SD_tot_area, bin_endpoints, dlogDp] = Read_2DC_wei(Aircraft_data_folder,Aircraft_data_filename)
function [time_2DC_wei, SD_2DC_wei_max, SD_2DC_wei_area, bin_endpoints_2DC_wei, bin_midpoints_2DC_wei, dlogDp] = Read_2DC_wei(SD_wei_path,Aircraft_data_filename)
% Reads NCAR netCDF-file

filename = [SD_wei_path,Aircraft_data_filename];
finfo = ncinfo(filename);
attributes = finfo.Attributes;

%%

% dN/dD
c2DS_all_max = ncread(filename,'conc_all')'; % 2DS Concentration (per cell), cm-4 
c2DS_all_area = ncread(filename,'conc_AreaR')';

% c2DS_drop = ncread(filename,'conc_AreaR_droplet')';
% c2DS_ice = ncread(filename,'conc_AreaR_ice')'; 
% c2DS_drop_2 = ncread(filename,'conc_liquid')';
% c2DS_ice_2 = ncread(filename,'conc_ice')'; 
% F2DC = squeeze(F2DC(:,1,:))';
% F2DC(:,1) = []; % remove first bin

% cm4 = cm3 * cm = 10-3 L * 10+4 um = 10 L um

bin_midpoints_2DC_wei = double(ncread(filename, 'bin_mid'))'*1000; % mm = 1000 um
bin_start = double(ncread(filename, 'bin_min'))'*1000;
bin_end = double(ncread(filename, 'bin_max'))'*1000;

bin_endpoints_2DC_wei = [bin_start, bin_end(end)];

dD = double(ncread(filename,'bin_dD'))'*1000;
% dDp = diff(bin_endpoints); % consistency check: this is the same

% Time
hhmmss = double(ncread(filename,'time')); % HHMMSS

%% Time to datenum

% get date (clunky)
date = strsplit(Aircraft_data_filename,'_');
date = strsplit(date{2},'.');
date = str2num ( date{2} ); % yyyymmdd

year = floor(date/10000); % first 4 digits
date = date - year*10000;
month = floor(date/100); % next 2
day = date - month * 100; % last 2

% get time
hour = floor(hhmmss/10000);
mmss = hhmmss - hour * 10000;
minute = floor(mmss/100);
second = mmss - minute*100;

% combine
time_2DC_wei = datenum(year, month, day, hour, minute, second);

%% Make SD
dlogDp = diff(log10(bin_endpoints_2DC_wei));
SD_2DC_wei_max = c2DS_all_max.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)
SD_2DC_wei_area = c2DS_all_area.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)

% SD_tot_max = [];
% SD_tot_max(1,:) = [0 0 bin_midpoints];
% SD_tot_max(2:length(time_2DC_wei)+1,1) = time_2DC_wei;
% SD_tot_max(2:length(time_2DC_wei)+1,2) = nansum(c2DS_all_max.*dD,2)/10; % 1/cm4 * um = 10-4 cm / cm4 = 10-4 /cm3 = 10-1 /L
% SD_tot_max(2:end,3:end) = c2DS_all_max.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)
% 
% SD_tot_area = [];
% SD_tot_area(1,:) = [0 0 bin_midpoints];
% SD_tot_area(2:length(time_2DC_wei)+1,1) = time_2DC_wei;
% SD_tot_area(2:length(time_2DC_wei)+1,2) = nansum(c2DS_all_area.*dD,2)/10; % 1/cm4 * um = 10-4 cm / cm4 = 10-4 /cm3 = 10-1 /L
% SD_tot_area(2:end,3:end) = c2DS_all_area.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)

% SD_ice = [];
% SD_ice(1,:) = [0 0 bin_midpoints];
% SD_ice(2:length(time)+1,1) = time;
% SD_ice(2:length(time)+1,2) = nansum(c2DS_ice.*dD,2)/10; % 1/cm4 * um = 10-4 cm / cm4 = 10-4 /cm3 = 10-1 /L
% SD_ice(2:end,3:end) = c2DS_ice.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)
% 
% SD_drop = [];
% SD_drop(1,:) = [0 0 bin_midpoints];
% SD_drop(2:length(time)+1,1) = time;
% SD_drop(2:length(time)+1,2) = nansum(c2DS_drop.*dD,2)/10; % 1/cm4 * um = 10-4 cm / cm4 = 10-4 /cm3 = 10-1 /L
% SD_drop(2:end,3:end) = c2DS_drop.*dD./dlogDp /10;   % dN/dlogD = D * dN/dD (since dlogD/dD = 1/D)

%%




end

