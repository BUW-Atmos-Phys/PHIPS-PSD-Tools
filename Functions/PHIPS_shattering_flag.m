function ShatteringFlag = PHIPS_shattering_flag(SD_2DS_path,SD_2DC_path, flight, taxis)


disp('Calculate shattering Flag...')
%% load 2DS SD from Wei (UO)
% Aircraft_data_folder = [SD_2DS_path, filesep, flight, filesep];
SD_wei_path = SD_2DS_path;

cd(SD_2DS_path)
listings = dir(['*', flight, '.2DS.H.nc']); 
if ~isempty(listings)
    disp('Load 2DS data')
    Aircraft_data_filename = listings(end).name;
    %[time_2DS_wei, SD_2DS_wei_max, SD_2DS_wei_area, bin_endpoints_2DS_wei, dlogDp_2DS] = Read_2DS_wei(SD_2DS_path,Aircraft_data_filename);
    [time_2DS_wei, SD_2DS_wei_max, SD_2DS_wei_area_liq, SD_2DS_wei_area_ice, SD_2DS_wei_area, ...
    bin_endpoints_2DS_wei, bin_midpoints_2DS_wei, dlogDp, TWC_2DS_wei] = Read_2DS_wei(SD_wei_path,Aircraft_data_filename);
%     [time_2DS_wei, SD_2DS_wei_max, SD_2DS_wei_area, bin_endpoints_2DS_wei, dlogDp] = Read_2DS_wei(SD_wei_path,Aircraft_data_filename);
    exist_2DS_Wei = true; 
else
    disp('No 2DS data available')
    exist_2DS_Wei = false; % this marker is needed because we cant calculate the SF when there's no 2DS data
end

SD_wei_path = SD_2DC_path;
cd(SD_2DC_path)
listings = dir(['*', flight, '.F2DC.nc']); 
if ~isempty(listings)
    disp('Load 2DC data')
    exist_2DC_Wei = true;
    Aircraft_data_filename = listings(end).name;
    % [time_2DC_wei, SD_2DC_wei_max, SD_2DC_wei_area, bin_endpoints_2DC_wei, dlogDp_2DC] = Read_2DC_wei(SD_2DC_path,Aircraft_data_filename);
    [time_2DC_wei, SD_2DC_wei_max, SD_2DC_wei_area, bin_endpoints_2DC_wei, bin_midpoints_2DC_wei, dlogDp] = Read_2DC_wei(SD_wei_path,Aircraft_data_filename);

else 
    disp('No 2DC data available')
    exist_2DC_Wei = false;
end
    
%% Calculate Shattering Flag
% SF for droplet is always 1

ratio_size_threshold_2DS = 800; % particles >X are "large" particles that are prone to shattering
ratio_SF_threshold_2DS = 10; %segments with ratio > X% are flagged as shattering
ratio_size_threshold_2DC = 800;
ratio_SF_threshold_2DC = 10; 



if exist_2DS_Wei == 1
    disp('Calculate 2DS SF')
    ShatteringFlag_2DS_UO = calculate_SF_2DX(time_2DS_wei, SD_2DS_wei_max, bin_endpoints_2DS_wei, ratio_size_threshold_2DS, ratio_SF_threshold_2DS, taxis);
else 
    disp('No 2DS data from Wei, cannot use this for SF')
end
if exist_2DC_Wei == 1
    disp('Calculate 2DC SF')
    ShatteringFlag_2DC_UO = calculate_SF_2DX(time_2DC_wei, SD_2DC_wei_max, bin_endpoints_2DC_wei, ratio_size_threshold_2DC, ratio_SF_threshold_2DC, taxis);
else
    disp('No 2DC data from Wei, cannot use this for SF')
end

if exist_2DS_Wei == 1;
    ShatteringFlag = ShatteringFlag_2DS_UO;
        if exist_2DC_Wei == 1;
            disp('Combine both, if 2DS == NaN, then use 2DC')
            idx = find(isnan(ShatteringFlag));
            ShatteringFlag(idx) = ShatteringFlag_2DC_UO(idx); % when 2DS = NaN, use 2DC
        end
else
    ShatteringFlag = ShatteringFlag_2DC_UO; % if there is no 2DS AT ALL, then use 2DC
end

ShatteringFlag(isnan(ShatteringFlag)) = 1; % quick fix of the QF: SF = NaN means SF = 1

end



%% function that calculates the indiviual SD of 2DS and 2DC (so we dont have to do it twice)
function ShatteringFlag_2DX = calculate_SF_2DX(time_ref,SD_tot_max,bin_endpoints, ratio_size_threshold, ratio_SF_threshold, taxis)
% calculate SF for 2DS or 2DC

ShatteringFlag_2DX = [];

for i=1:length(taxis)-1 % loop over all timesteps
    t1 = taxis(i)-datenum(0,0,0,0,0,1/1000); %1ms tolerance, because sometimes the timestamps are off by 1/100 ms
    t2 = taxis(i+1)+datenum(0,0,0,0,0,1/1000);

    % calculate 
    ratio = calculate_ratio_SF(time_ref,SD_tot_max,bin_endpoints,ratio_size_threshold,t1,t2);

    if isnan(ratio) 
        ShatteringFlag_2DX(i,1) = NaN;
    elseif ratio > ratio_SF_threshold
        ShatteringFlag_2DX(i,1) = 0;
    else
        ShatteringFlag_2DX(i,1) = 1;
    end

end

end


%% function to calculate ratio
function ratio = calculate_ratio_SF(time_ref, SD_ref, bin_endpoints, ratio_size_threshold, t1, t2)
% size_threshold = 1000
% bin_endpoints_2DC = [25,50,100,150,200,250,300,350,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000]; 


    lower_threshold = 200;

    bin_midpoints = bin_endpoints(1:end-1)+diff(bin_endpoints)./2;
    % Calculate the bin widths
    dlogDp = [];
    for j = 1:length(bin_midpoints)
        dlogDp(j) = log10(bin_endpoints(j+1))-log10(bin_endpoints(j));
    end

    % find 2DS data in this region
    idx = find(time_ref >= t1 & time_ref < t2);
    if isempty(idx) % no data in this interval, QF = NaN
        ratio = NaN;
    else
        SD_roi = nansum(SD_ref(idx,:),1); % sum over all seconds in this interval (downwards summation)
        N_tot = SD_roi .* dlogDp;


        % calculate ratio
        idx = find(bin_endpoints >= ratio_size_threshold);
        idx(end) = [];
        SD_sum_large = nansum(N_tot(:,idx),2);

        idx = find(bin_endpoints >= lower_threshold);
        idx(end) = [];
        SD_sum_all = nansum(N_tot(:,idx),2);
        ratio = SD_sum_large/SD_sum_all * 100;
    end
end








