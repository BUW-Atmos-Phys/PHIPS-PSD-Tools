function [SD_2DS_RAF, time_2DS_RAF, bin_endpoints_2DS_RAF] = Read_2DS_RAF(Aircraft_data_folder,Aircraft_data_filename)
% Reads NCAR netCDF-file

filename = [Aircraft_data_folder,Aircraft_data_filename];
finfo = ncinfo(filename);
attributes = finfo.Attributes;

%%
% date
index =structfind(attributes,'Name','FlightDate');
C = {attributes.Value};
date = datenum(C{index},'mm/dd/yyyy');

% dN
SD_2DS_RAF = ncread(filename,'C2DSA_2H'); % 2DS Concentration (per cell), #/L (_R for only round
SD_2DS_RAF = squeeze(SD_2DS_RAF(:,1,:))';
SD_2DS_RAF(:,1) = []; % remove first bin

bin_endpoints_2DS_RAF = ncreadatt(filename,'C2DSA_2H','CellSizes');

% Time
time_2DS_RAF = double(ncread(filename,'Time'));

%% Time to datenum
time_2DS_RAF = date + datenum(0,0,0,0,0,time_2DS_RAF);


%% Size bins
bin_midpoints = bin_endpoints_2DS_RAF(1:end-1)+diff(bin_endpoints_2DS_RAF)./2;

% Calculate the bin widths
dlogDp = [];
for i = 1:length(bin_midpoints)
    dlogDp(i) = log10(bin_endpoints_2DS_RAF(i+1))-log10(bin_endpoints_2DS_RAF(i));
end

SD_2DS_RAF = SD_2DS_RAF./repmat(dlogDp,size(SD_2DS_RAF,1),1);

% %% Make SD
% SD = [];
% SD(1,:) = [0 0 bin_midpoints];
% SD(2:length(time)+1,1) = time;
% SD(2:length(time)+1,2) = nansum(SD_F2DC,2);
% SD(2:end,3:end) = SD_F2DC./repmat(dlogDp,size(SD_F2DC,1),1);

end

