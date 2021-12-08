function [SD_2DC_RAF, time_2DC_RAF, bin_endpoints_2DC_RAF] = Read_F2DC_RAF(Aircraft_data_folder,Aircraft_data_filename)
% Reads NCAR netCDF-file

filename = [Aircraft_data_folder,Aircraft_data_filename];
finfo = ncinfo(filename);
attributes = finfo.Attributes;


%% 'CONC2DCA_RWOI' - all - total conc
% 'CONC2DCR_RWOI' - round
% 'PLWC2DCA_RWOI' and R_ - liq water content
% 'C2DCA_RWOI' - and _A - concentration per channel

% date
index =structfind(attributes,'Name','FlightDate');
C = {attributes.Value};
date = datenum(C{index},'mm/dd/yyyy');

% dN
%C2DCA_RWOI
SD_2DC_RAF = ncread(filename,'C2DCA_RWOI'); % F2DC Concentration (per cell), #/L (_R for only round
SD_2DC_RAF = squeeze(SD_2DC_RAF(:,1,:))';
SD_2DC_RAF(:,1) = []; % remove first bin

bin_endpoints_2DC_RAF = ncreadatt(filename,'C2DCA_RWOI','CellSizes'); 

% Time
time_2DC_RAF = double(ncread(filename,'Time'));

%% Time to datenum
time_2DC_RAF = date + datenum(0,0,0,0,0,time_2DC_RAF);

%% Size bins
bin_midpoints = bin_endpoints_2DC_RAF(1:end-1)+diff(bin_endpoints_2DC_RAF)./2;

% Calculate the bin widths
dlogDp = [];
for i = 1:length(bin_midpoints)
    dlogDp(i) = log10(bin_endpoints_2DC_RAF(i+1))-log10(bin_endpoints_2DC_RAF(i));
end

SD_2DC_RAF = SD_2DC_RAF./repmat(dlogDp,size(SD_2DC_RAF,1),1);

%% Make SD
% SD = [];
% SD(1,:) = [0 0 bin_midpoints];
% SD(2:length(time_2DC_RAF)+1,1) = time_2DC_RAF;
% SD(2:length(time_2DC_RAF)+1,2) = nansum(SD_2DC_RAF,2);
% SD(2:end,3:end) = SD_2DC_RAF./repmat(dlogDp,size(SD_2DC_RAF,1),1);

end

