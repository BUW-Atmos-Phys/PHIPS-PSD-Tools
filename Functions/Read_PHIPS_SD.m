function [time_PHIPS, ShatteringFlag, SD_PHIPS_ice, SD_PHIPS_ice_only_images, SD_PHIPS_drop, ...
    N_ice, N_ice_only_images, N_drop, ...
    N_ice_uncertainty, N_ice_uncertainty_only_images, N_drop_uncertainty, ...
    counts_ice, counts_ice_only_images, counts_drop, ...
    bin_endpoints_PHIPS, bin_midpoints_PHIPS,bin_midpoints_area_PHIPS, ...
    SD_PHIPS_ice_area, SD_PHIPS_ice_area_only_images,A_ice] = Read_PHIPS_SD(savepath, tstep)


%% Load PHIPS PSD - Number
cd(savepath)
listings = dir(['*', num2str(tstep), 's_droplet_v1.sum']); % PSDs based on images + scattering
if isempty(listings)
    disp('No SD found!')
    return 
end
filename = listings(end).name;
filename = [savepath,filesep,filename];
SD_drop_raw = dlmread(filename);

listings = dir(['*', num2str(tstep), 's_ice_v1.sum']); % PSDs based on images + scattering
filename = listings(end).name;
filename = [savepath,filesep,filename];
SD_ice_raw = dlmread(filename);

listings = dir(['*', num2str(tstep), 's_ice_only_images.sum']); % PSDs based on ONLY images (no shattering)
filename = listings(end).name;
filename = [savepath,filesep,filename];
SD_ice_raw_only_images = dlmread(filename);

%% Load PHIPS PSD - Area
listings = dir(['*', num2str(tstep), 's_ice_area.sum']); % PSDs based on images
if isempty(listings)
    disp('No SD found!')
    return 
end
filename = listings(end).name;
filename = [savepath,filesep,filename];
SD_ice_area_raw = dlmread(filename);

listings = dir(['*', num2str(tstep), 's_ice_area_only_images.sum']); % PSDs based on ONLY images (no shattering)
if isempty(listings)
    disp('No SD found!')
    return 
end
filename = listings(end).name;
filename = [savepath,filesep,filename];
SD_ice_area_raw_only_images = dlmread(filename);


%%
% consistency check: are SD files the same size?
if size(SD_drop_raw) ~= size(SD_ice_raw)
    disp('ERROR! FILE SIZE MISMATCH')
    return
end

time_PHIPS = SD_ice_raw(2:end,1);
ShatteringFlag = SD_ice_raw(2:end,2);

% PSD counts
SD_PHIPS_drop = SD_drop_raw(2:end,5:end);
N_drop = SD_drop_raw(2:end,3);
N_drop_uncertainty = SD_drop_raw(2:end,4);
SD_PHIPS_ice = SD_ice_raw(2:end,5:end);
N_ice = SD_ice_raw(2:end,3);
N_ice_uncertainty = SD_ice_raw(2:end,4);
SD_PHIPS_ice_only_images = SD_ice_raw_only_images(2:end,5:end);
N_ice_only_images = SD_ice_raw_only_images(2:end,3);
N_ice_uncertainty_only_images = SD_ice_raw_only_images(2:end,4);

% PSD area 
SD_PHIPS_ice_area = SD_ice_area_raw(2:end,5:end);
A_ice = SD_ice_area_raw(2:end,3); % total area [um2 L-1]
SD_PHIPS_ice_area_only_images = SD_ice_area_raw_only_images(2:end,5:end);
A_ice_only_images = SD_ice_area_raw_only_images(2:end,3); % total area [um2 L-1]


%% consistency check: diameter
bin_midpoints_PHIPS = SD_drop_raw(1, 5:end);
bin_midpoints_area_PHIPS = pi .* (bin_midpoints_PHIPS./2).^2;

bin_endpoints_PHIPS = [15,30,60,100,150,200,250,300,350,400,500,600,700];

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
filename = [savepath,filesep,filename];
counts_drop_raw = readmatrix(filename,'FileType','text');

listings = dir(['*', num2str(tstep), 's_ice_counts.sum']); 
filename = listings(end).name;
filename = [savepath,filesep,filename];
counts_ice_raw = readmatrix(filename,'FileType','text');

listings = dir(['*', num2str(tstep), 's_ice_only_images_counts.sum']); 
filename = listings(end).name;
filename = [savepath,filesep,filename];
counts_ice_raw_only_images = readmatrix(filename,'FileType','text');

counts_drop = counts_drop_raw(2:end,4:end);
counts_ice = counts_ice_raw(2:end,4:end);
counts_ice_only_images = counts_ice_raw_only_images(2:end,4:end);

end







