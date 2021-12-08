%% Revision history
% 17.05.2021, FW, adjust for updated PO structure and added ASF from lvl 0

function [SDice,SDdroplet] = PHIPS_SD(particleopticspath,SD_2DS_path,SD_2DC_path,campaign,flight,tstep,save_status,start_time,end_time)
% PHIPS size distribution from scattering data
% Input arguments: POpath -> path to ParticleOptics
%                  SD_2DS_path -> path to Wei's 2DS and 2DC data
%                  flight -> flight name
%                  campaign -> campaign name, e.g. 'SOCRATES'
%                  tstep -> time step in seconds
%                  t_start/end -> optional

% ADJUST BINS ALSO IN Read_PHIPS_SD.m for translate_SD_sum_to_nc.m
% bin_endpoints = [20,40,60,80,100,125,150,200,250,300,350,400,500,600,700];
bin_endpoints = [30,60,100,150,200,250,300,350,400,500,600,700];

%savepath = [particleopticspath,filesep,'PHIPS Results', filesep, 'Campaigns', filesep,campaign,filesep,flight, filesep, 'SD'];
% folder = 'C:\Users\Fritz\Desktop\PHIPS SD test sens area\';
% savepath = [folder, filesep, num2str(bin_endpoints(1)), '_new2', filesep, flight, filesep];
if ~isdir(savepath)
    mkdir(savepath)
end

%% Calibration coefficients
if strcmp(campaign,'ACLOUD')
    % calibration coefficient
    % Droplets
    a_droplet = 0.86098;
    b_droplet = 0.5;
    % Ice
    a_ice = 1.4605;
    b_ice = 0.5;
    
    BG_offset = 236.2868;
    
elseif strcmp(campaign,'SOCRATES')
    % Droplets
    a_droplet = 1.4441; % 1.4118
    b_droplet = 0.5;
    % Ice
    a_ice = 1.4167; % 1.5137
    b_ice = 0.5;
    
    BG_offset = 238.188;
    
elseif strcmp(campaign,'CIRRUS-HL')
    % Droplets
    a_droplet = 0.61138;
    b_droplet = 0.5;
    % Ice
    a_ice = 0.85673;
    b_ice = 0.5;
    
    BG_offset = 254.6728;
    disp('Preliminary calibration coefficients!')
    
else
    disp('Calculate Calibration coefficients first');
    return
end

% Folder to PHIPS data
folder = [particleopticspath,filesep,'Phips',filesep,campaign,filesep,flight];
folder = [particleopticspath,filesep,'PHIPS Results', filesep, 'Campaigns', filesep,campaign,filesep,flight, filesep, 'Level Files', filesep];

%% Define bins

bin_midpoints = bin_endpoints(1:end-1)+diff(bin_endpoints)./2;

dlogDp = zeros(1,length(bin_endpoints)-1);
for i=1:length(bin_endpoints)-1
    dlogDp(i) = log10(bin_endpoints(i+1))- log10(bin_endpoints(i));  % =/= log dDp !
end

%% PHIPS sensitive area
% Seansitive area function, a*x^b+c
% new, FRED simulation

A_sens_parameters_ice.a = 0.091077;
A_sens_parameters_ice.b = 0.026101;
A_sens_parameters_ice.c = -0.09853;

A_sens_parameters_drop.a = -0.037236;
A_sens_parameters_drop.b = -0.28694;
A_sens_parameters_drop.c = 0.013069;


%%
ignore_above_saturation = 1;
ignore_below_background = 0;

%%  ------- Initialization ended ------------


%% Make SD

if exist('a_droplet')
    %% Load PhipsData Lvl 0 for ASF
    cd(folder)
    % Level 0 data
    listings = dir('*level_0.csv');
    csv_name = listings.name;
    filename = [folder,filesep,csv_name];
    PhipsData_level_0 = import_phips(filename); %use the level_0 file later for deadtime calculation
    
    %% get Droplet Flag from lvl 2
    listings = dir('*level_2.csv');
    if isempty(listings)
        disp('Level File not found!')
        SDice = NaN;
        SDdroplet = NaN;
    else
        csv_name = listings(end).name;
        filename = [folder,filesep,csv_name];
        PhipsData = import_phips(filename);
        
        %% get ASF from lvl_0
        PhipsData = replace_ASF_with_lvl0(folder, PhipsData);
        
        %% Remove saturates
        % if more than 4 channels are saturated/background, remove the particle (DropletFlag = NaN)
        
        PhipsArray = table2array(PhipsData);
        
        sat_tresh = 4; %saturated if >= X chanels are saturated
        bg_tresh = 4;
        
        num_sat = 0; num_BG = 0;
        for i = 1:size(PhipsArray,1)
            if sum(PhipsArray(i, end-16:end) >= 2047) >= sat_tresh && ignore_above_saturation == 1
                PhipsData.DropletFlag(i) = NaN; %erase those particles from Drop and Ice, but they still appear in Total
                num_sat = num_sat + 1;
            end
            
            if sum(PhipsArray(i, end-16:end) <= 5) >= bg_tresh && ignore_below_background == 1
                PhipsData.DropletFlag(i) = NaN; %erase those particles from Drop and Ice, but they still appear in Total
                num_BG = num_BG + 1;
            end
        end
        
        disp([num2str(num_sat), ' particles above saturation removed'])
        %  disp([num2str(num_BG), ' particles below BG removed'])
        
        %% Load Aircraft Data
        if strcmp(campaign,'ACLOUD')
            Aircraft_data_folder = [particleopticspath, '/SID3/ACLOUD/Flight data/'];
            flight_date = strsplit(flight, ' ');
            flight_date = flight_date{2};
            Aircraft_data_filename = ['AircraftData_', flight_date, '.dat'];
        elseif strcmp(campaign,'SOCRATES')
            Aircraft_data_folder = [particleopticspath,'/Phips/SOCRATES/Flight data/'];
            Aircraft_data_filename = ['SOCRATES',lower(flight),'.nc'];
        elseif strcmp(campaign,'CIRRUS-HL')
            Aircraft_data_folder = [particleopticspath,'/Phips/CIRRUS-HL/Flight data/'];
            %  'E:\Phips\CIRRUS-HL\Flight data\QL-CIRRUS-HL_F01_20210624a_ADLR_BAHAMAS_v2.nc';
            flightnum = strsplit(flight,'RF');
            flightnum = flightnum{2};
            listings = dir([Aircraft_data_folder, filesep, 'QL-CIRRUS-HL_F', flightnum, '_*_BAHAMAS_v1.nc']);
            Aircraft_data_filename = listings(end).name;
        end
        
        [time,airspeed,height,lat,lon,hFig1,hFig2,T] = Aircraft_data(campaign,Aircraft_data_folder,Aircraft_data_filename);
        close all
        
        
        %% Integrate Scattering Intensity
        % Get index of starting angle (42 degrees)
        A = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle42'));
        % Integrate
        ang_temp = [42:8:170];
        SPF_temp = table2array(PhipsData(:,A:end));
        integrated_intensity = trapz(ang_temp, SPF_temp,2);
        
        PhipsData = addvars(PhipsData,integrated_intensity,'NewVariableNames','integrated_intensity');
        
        %% B.2. Generate histograms
        if nargin == 9 %careful when debugging this step by step: nargin command does only work if used as function
            start_time = start_time + datenum(0,0,0,0,0,tstep/2); %to compensate for the +/-tstep/2 in the next step
            end_time = end_time - datenum(0,0,0,0,0,tstep/2);
            taxis = [start_time-datenum(0,0,0,0,0,tstep/2):datenum(0,0,0,0,0,tstep):end_time+datenum(0,0,0,0,0,tstep/2)];
        else %if no start/end time is given in the input of the function
            %% Start from next full minute minus half of timestep
            temp = datevec(PhipsData_level_0.RealTimeStamp(1));
            taxis = [datenum(temp(1),temp(2),temp(3),temp(4),temp(5),0+tstep/2):datenum(0,0,0,0,0,tstep):PhipsData_level_0.RealTimeStamp(end)];
        end
        
        %% Shattering correction based on interarrivaltime
        thresh_interarrivaltime = 0.5/1000; % = 0.5 ms
        [PhipsData] = PHIPS_remove_shattered(PhipsData,thresh_interarrivaltime);
        
        %% Calculate Shattering Flag
        if strcmp(campaign, 'SOCRATES')
            ShatteringFlag = PHIPS_shattering_flag(SD_2DS_path,SD_2DC_path, flight, taxis);
        else
            ShatteringFlag = zeros(length(taxis)-1, 1)+1;
        end
        
        %% calculate number of particles per size bin
        
        dp_ice = abs(a_ice .*(PhipsData.integrated_intensity(PhipsData.DropletFlag==0)-BG_offset).^(b_ice));
        dp_drop = abs(a_droplet .*(PhipsData.integrated_intensity(PhipsData.DropletFlag==1)-BG_offset).^(b_droplet));
        
        PhipsData.dp (PhipsData.DropletFlag==0) = dp_ice;
        PhipsData.dp (PhipsData.DropletFlag==1) = dp_drop;
        PhipsData.dp (isnan(PhipsData.DropletFlag)) = NaN;
        
        disp('Calculate SD (this may take a while...)')
        disp('... for ice...')
        [particles_per_bin_ice, conc_ice, V_ice, V_err_ice] = generate_histogram(PhipsData(PhipsData.DropletFlag==0,:),...
            PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_ice,airspeed,time);
        disp('... for droplets...')
        [particles_per_bin_drop, conc_drop, V_drop, V_err_drop] = generate_histogram(PhipsData(PhipsData.DropletFlag==1,:),...
            PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_drop,airspeed,time);
        conc_drop(:,1) = NaN; % delete the first bin for droplets as the lower limit is 50
        
        %% Produce PHIPS_counts SD
        SDice_counts = [0 0 0 bin_midpoints;...
            taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag sum(particles_per_bin_ice,2) particles_per_bin_ice];
        SDdroplet_counts = [0 0 0 bin_midpoints;...
            taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 sum(particles_per_bin_drop,2) particles_per_bin_drop]; % SF for droplet is always 1
        
        %% Make SD files
        
        % Total concentration
        Ntot_ice = 1000.*nansum(conc_ice,2); % [L^-1]
        Ntot_droplet = 1000.*nansum(conc_drop,2);
        Ntot_err_droplet = Ntot_ice.* NaN; % Ntot_droplet.*nanmean(nanmean(V_err_droplet./V_droplet));
        Ntot_err_ice = Ntot_ice.* NaN; % Ntot_ice.*nanmean(nanmean(V_err_ice./V_ice)); % THIS NEEDS TO BE CHANGED
        
        % dNdlogDp
        dNdlogDp_ice = zeros(size(conc_ice));
        for i=1:size(conc_ice(:,1),1)
            dNdlogDp_ice(i,:) = 1000.*(conc_ice(i,:)./dlogDp);
        end
        
        dNdlogDp_droplet = zeros(size(conc_drop));
        for i=1:size(conc_drop(:,1),1)
            dNdlogDp_droplet(i,:) = 1000.*(conc_drop(i,:)./dlogDp);
        end
        
        SDice = [0 0 0 0 bin_midpoints;...
            taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice Ntot_err_ice dNdlogDp_ice];
        SDdroplet = [0 0 0 0 bin_midpoints;...
            taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 Ntot_droplet Ntot_err_droplet dNdlogDp_droplet]; % SF for droplet is always 1
        
%         %% fill NaN for time where PHIPS was not operating yet / shut down already
%         % not needed anymore if we produce the SD only for time limits
%         %  based on the lvl0 file instead of from aircraft data
%         t1 = PhipsData_level_0.RealTimeStamp(1);
%         t2 = PhipsData_level_0.RealTimeStamp(end);
%         
%         idx1 = find(taxis< t1)
        
        
        %% check how many segments are SF = 0
        idx = find(SDice(:, 3) > 0); % Ntot > 0
        SF0 = length(find(SDice(idx,2) == 0));
        SF1 = length(find(SDice(idx,2) == 1));
        
        disp([num2str(SF0), ' segments are SF=0, ', num2str(SF1), ' are SF = 1, that means that ', num2str(SF0/(SF0+SF1)*100), '% are rejected)'])
        
        %% Save SD
        % for debugging:
        % input_name = 'Ice';
        % input = SDice;
        % a = a_ice;
        % b = b_ice;
        disp('save...')
        save_SD(SDice, SDice_counts, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bin_endpoints);
        save_SD(SDdroplet, SDdroplet_counts, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bin_endpoints);
        
        translate_SD_sum_to_nc(savepath, campaign, flight, tstep);
        
    end % end of "if exists level file"
end % end of "if exists a_droplet"

disp('..all done.')
end


%% FUNCTIONS 

%%
function [particles_per_bin, conc, V, V_err] = generate_histogram(PhipsData, PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters ,airspeed,time) 
% Function for Generating histograms for PHIPS SD

% A_sens_parameters = A_sens_parameters_ice;

%% Make threshist and calculate volume
V_bin = [];
for i=1:length(taxis)-1
    t1 = taxis(i)-datenum(0,0,0,0,0,1/1000); %1ms tolerance, because sometimes the timestamps are off by 1/100 ms
    t2 = taxis(i+1)+datenum(0,0,0,0,0,1/1000);
    
    %% calculate SD
    % Select time
    % deadtime calculation (including ForceTriggers)
    idx= PhipsData_level_0.RealTimeStamp>=t1 & PhipsData_level_0.RealTimeStamp<t2;
    number_of_triggers(i)=sum(idx); %histcounts(data,binlimits) sorts data into the number of bins determined by binlimits (intensity)
    qdead = (1-sum(number_of_triggers(i),2)./tstep*12*10^-6);% Fraction of missed volumne due to el. dead time of 12*10^-6s
    
    % size of particles in this time
    idx = PhipsData.RealTimeStamp>=t1 & PhipsData.RealTimeStamp<t2;
    dp = PhipsData.dp(idx); % size of each particle in this time
    
    % corresponding A_sens of EACH PARTICLE in this time frame
    sensitive_area = A_sens_parameters.a.*dp.^ A_sens_parameters.b + A_sens_parameters.c ; % cm2 % error 27.5%
    % Do not allow negative sensitive areas
  %  sensitive_area(sensitive_area<0) = min(sensitive_area(sensitive_area>0));
    
    % Calculate probed volume (corresponding to each particle)
    idx = time>=t1 & time<t2;  %aircraft_time
    V_per_particle = sensitive_area .* nanmedian(airspeed(idx)).*100.*tstep*qdead; % The sampling Volume [cm-3]
    
    % % OLD
    sensitive_area_bin = A_sens_parameters.a.*bin_endpoints.^ A_sens_parameters.b + A_sens_parameters.c ; % cm2 % error 27.5%
    % Do not allow negative sensitive areas
    %  sensitive_area(sensitive_area<0) = min(sensitive_area(sensitive_area>0));
    
    % Calculate probed volume (corresponding to each particle)
    idx = time>=t1 & time<t2;  %aircraft_time
    V_bin(i,:) = sensitive_area_bin .* nanmedian(airspeed(idx)).*100.*tstep*qdead;
    
%% sort in bins

    % loop over all bins
    for j = 1:length(bin_endpoints)-1
        
        idx = find(dp > bin_endpoints(j) & dp <= bin_endpoints(j+1));
        
        particles_per_bin(i,j) = length(idx); % number of particles per bin
        % old: conc = [sum ( N ) ] / V , where V was const for one bin
        % new: conc = sum ( N / V ) = 1 / sum(V), where N = 1 (since we
        % count every particle individually and V is the Volume of each particle)
        conc(i,j) = sum(1./V_per_particle(idx)); % recipreocely sumated volume per particle
        
        % V not needed as output...
        V(i,j) = 1 / (conc(i,j) / particles_per_bin(i,j)); %  effective sensitive volume of each bin (weighted by size of each particle) = conc / num of particles
        V_err(i,j) = 28.5/100 * V(i,j); %an error of 27.5% for the Area and 1% for the airspeed.  % this def. needs improvement...
        
        % % OLD 
%        conc(i,j) = sum(particles_per_bin(i,j))/V_bin(i,j+1);

        % if no particles, then no concentration
        if particles_per_bin(i,j) == 0
            conc(i,j) = 0;
        end
    end % end loop size bins
    
end % end loop time

end % end function


