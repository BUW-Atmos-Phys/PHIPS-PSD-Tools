function [time_PHIPS, ShatteringFlag, SD_PHIPS_ice, SD_PHIPS_drop, N_ice, N_drop, N_ice_uncertainty, N_drop_uncertainty, counts_ice, counts_drop, ...
    bin_endpoints_PHIPS, bin_midpoints_PHIPS] = Read_PHIPS_SD(savepath, tstep);

% folder = [particleopticspath,filesep, 'PHIPS Results', filesep, campaign, filesep, flight, filesep, 'SD'];
% C:\Users\wa9929\Desktop\SD_Emma\PHIPS Results\Campaigns\SOCRATES\RF08\SD

% folder = 'C:\Users\Fritz\Desktop\RF02\SD\'

% tstep = 1;

%% Load PHIPS SD
cd(savepath)
listings = dir(['*', num2str(tstep), 's_droplet.sum']); 
if isempty(listings)
    disp('No SD found!')
    return 
end
filename = listings(end).name;
filename = [savepath,'\',filename];
SD_drop_raw = dlmread(filename);

listings = dir(['*', num2str(tstep), 's_ice.sum']); 
filename = listings(end).name;
filename = [savepath,'\',filename];
SD_ice_raw = dlmread(filename);

%%
% consistency check: are SD files the same size?
if size(SD_drop_raw) ~= size(SD_ice_raw)
    disp('ERROR! FILE SIZE MISMATCH')
    return
end

time_PHIPS = SD_ice_raw(2:end,1);

ShatteringFlag = SD_ice_raw(2:end,2);
SD_PHIPS_drop = SD_drop_raw(2:end,5:end);
N_drop = SD_drop_raw(2:end,3);
N_drop_uncertainty = SD_drop_raw(2:end,4);
SD_PHIPS_ice = SD_ice_raw(2:end,5:end);
N_ice = SD_ice_raw(2:end,3);
N_ice_uncertainty = SD_ice_raw(2:end,4);
%SD_tot = SD_drop + SD_ice;

%% consistency check: diameter
bin_midpoints_PHIPS = SD_drop_raw(1, 5:end);

% bin_endpoints_PHIPS = [20,45,70,100,150,200,250,300,350,400,500,600,700];
% bin_endpoints_PHIPS = [20,40,60,80,100,125,150,200,250,300,350,400,500,600,700];
% bin_endpoints_PHIPS = [35,50,60,80,100,125,150,200,250,300,350,400,500,600,700];
bin_endpoints_PHIPS = [30,60,100,150,200,250,300,350,400,500,600,700];

% problem: you cant calculate bin_endpoints from midpoints
% idea:  at least check the other way, if it fits
bin_midpoints_consitency_check = [];
for i = 1:length(bin_endpoints_PHIPS)-1
    bin_midpoints_consitency_check(i) = (bin_endpoints_PHIPS(i) + bin_endpoints_PHIPS(i+1))/2;
    if bin_midpoints_consitency_check(i) ~= bin_midpoints_PHIPS(i)
       disp('ERROR! BIN DIAMETERS ARE WRONG!') 
       return
    end
end

%dlogDp_PHIPS = diff(log10(bin_endpoints_PHIPS));

%% Load PHIPS SD COUNTS
cd(savepath)
listings = dir(['*', num2str(tstep), 's_droplet_counts.sum']); 
if isempty(listings)
    disp('No COUNTS SD found!')
    return 
end
filename = listings(end).name;
filename = [savepath,'\',filename];
counts_drop_raw = dlmread(filename);

listings = dir(['*', num2str(tstep), 's_ice_counts.sum']); 
filename = listings(end).name;
filename = [savepath,'\',filename];
counts_ice_raw = dlmread(filename);

counts_drop = counts_drop_raw(2:end,4:end);
counts_ice = counts_ice_raw(2:end,4:end);

end







