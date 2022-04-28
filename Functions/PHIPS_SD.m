%% Revision history
% 17.05.2021, FW, adjust for updated PO structure and added ASF from lvl 0
% 08.12.2021, FW, sens volume is calculated for each particle (not just once per bin)
% 07.02.2022, FW, adjusted paths and improved documentation

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
% dN_ice = -999; dN_drop = -999;

savepath = [particleopticspath,filesep,'PHIPS Results', filesep, 'Campaigns', filesep,campaign,filesep,flight, filesep, 'SD', filesep];
% savepath = 'C:\Users\Fritz\Desktop\PHIPS SD test flight speed correction\new\';

% savepath = [folder, filesep, num2str(bin_endpoints(1)), '_new2_1s', filesep, flight, filesep];
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

% % Define bins

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

% %
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
    
    % % get Droplet Flag from lvl 2
    listings = dir('*level_2.csv');
    %%
    if isempty(listings)
        disp('Level File not found!')
        SDice = NaN;
        SDdroplet = NaN;
    else
        csv_name = listings(end).name;
        filename = [folder,filesep,csv_name];
        PhipsData = import_phips(filename);
        
        % % get ASF from lvl_0
        PhipsData = replace_ASF_with_lvl0(folder, PhipsData);
        
        % % Remove saturates
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
        
        % % Integrate Scattering Intensity
        % Get index of starting angle (42 degrees)
        A = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle42'));
        % Integrate
        ang_temp = [42:8:170];
        SPF_temp = table2array(PhipsData(:,A:end));
        integrated_intensity = trapz(ang_temp, SPF_temp,2);
        
        PhipsData = addvars(PhipsData,integrated_intensity,'NewVariableNames','integrated_intensity');
        
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
        
        % % Calculate Shattering Flag
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
        
%         % % estimate uncertainty 
        % ice
        dN_minus =  [ -24.9479  -24.0072  -18.7372  -17.7478  -14.5033  -11.3761   -7.8475   -8.1237   -5.9329   -5.2310   -3.8175];
        dN_plus =   [45.5047   44.6152   51.4139   48.3331   32.1553   23.5189   23.0540   17.5292   12.1717    6.3060    7.0891];
        dN_ice = (dN_plus - dN_minus)./2; % [%]
        % droplets
        dN_minus = [ -39.8952   -9.7177  -10.0118  -10.9244   -9.2650  -13.1825  -11.8750   -3.4833   -8.8068       NaN       NaN];
        dN_plus =   [42.5975   12.1592   12.7038   17.2285   14.0938   22.1102   35.3615   15.7782    9.8641       NaN       NaN];
        dN_drop = (dN_plus - dN_minus)./2;
        dN_drop(end-1:end) = dN_drop(end-2); % set the last 2 bins, which dont have enough data, on the last valid point
        
        % % Make SD files
        % Total concentration
        Ntot_ice = 1000.*nansum(conc_ice,2); % [L^-1]
        Ntot_droplet =   1000.*nansum(conc_drop,2);
        
        Ntot_err_ice = 1000.*nansum(conc_ice.*dN_ice./100,2); % [L^-1]   
        Ntot_err_droplet = 1000.*nansum(conc_drop.*dN_drop./100,2);
        
        
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
        
        
        % % check how many segments are SF = 0
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
        
        % save SD as .nc file
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
    %%
    
%     i = 1293
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
    TAS = nanmedian(airspeed(idx));
    
    % correction flight speed        
    % https://amt.copernicus.org/articles/9/5135/2016/
    % -25 % for D<70, -5% for D>70
    
    xi = 0.99 + 2.55e-4 * TAS - 3.3e-6 * TAS .^2;
%     mu = [];
%     if length(sensitive_area)>0
%         for j = 1:length(sensitive_area)
%             if dp(j) < 70
%                 mu(j,1) = 0.8;
%             else
%                 mu(j,1) = 1;
%             end
%         end
%     else
%         mu = 1;
%     end
    mu = 1;
    corr_flight_speed = xi * mu; 
    % corr_flight_speed = 1;
    
    V_per_particle = sensitive_area * TAS .*100.*tstep*qdead ./corr_flight_speed; % The sampling Volume [cm-3]
    
    % % OLD
    sensitive_area_bin = A_sens_parameters.a.*bin_endpoints.^ A_sens_parameters.b + A_sens_parameters.c ; % cm2 % error 27.5%
    % Do not allow negative sensitive areas
    %  sensitive_area(sensitive_area<0) = min(sensitive_area(sensitive_area>0));
    
    % Calculate probed volume (corresponding to each particle)
    idx = time>=t1 & time<t2;  %aircraft_time
    V_bin(i,:) = sensitive_area_bin * TAS .*100.*tstep*qdead;
    
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

%% get the raw ASF from the lvl 0 file
function PhipsData = replace_ASF_with_lvl0(folder, PhipsData);
% 
% produce level 1 new but keep image assignment so we dont have to do it again
% folder = path, e.g. 
% folder = [particleopticspath, '\PHIPS Results\Campaigns\', campaign, filesep, flight, filesep, 'Level Files', filesep];
% level = 1 or 2 or ...

% initialize
% particleopticspath = 'P:\';
% campaign = 'SOCRATES';
% flight = 'RF11';
% folder = [particleopticspath, '\PHIPS Results\Campaigns\', campaign, filesep, flight, filesep, 'Level Files', filesep];
cd(folder)

%% load lvl 0
listings = dir('*_level_0.csv');
filename_0 = listings(end).name;
Level_0 = readtable(filename_0);
Level_0 = Level_0(Level_0.ForcedParticleFlag == 0,:); % remove forcetriggers
Level_0.RealTimeStamp = datenum(Level_0.RealTimeStamp);

%% load lvl 1
% listings = dir(['*_level_', num2str(level), '.csv']);
% filename_1 = listings(end).name;
% Level_1_old = readtable(filename_1);
% Level_1_old.RealTimeStamp = datenum(Level_1_old.RealTimeStamp);

Level_1_old = PhipsData;

%% clear saturated channels

Level_0_varnames = Level_0.Properties.VariableNames;

Level_0a = table2array(Level_0);
scatt_channels = Level_0a(:,end-30:end);
scatt_channels(scatt_channels >= 2047) = NaN;
Level_0a(:,end-30:end) = scatt_channels;

Level_0_cleared = array2table(Level_0a, 'VariableNames', Level_0_varnames);

%% Do channle to angle order

csv_filename = dir('*level_0.csv'); 
csv_filename = csv_filename(end).name;
% Make a filename
C = strsplit(csv_filename,'_');
save_filename = [C{1},'_',C{2},'_level_1.csv'];

% Date
PHIPS_date = floor(datenum(C{2},'yyyymmdd-HHMM'));

if PHIPS_date>=datenum(2019,6,4)
    v = 5;
elseif PHIPS_date>=datenum(2017,9,22) & PHIPS_date<datenum(2019,6,4)
    v = 4;
elseif PHIPS_date>=datenum(2017,2,23) & PHIPS_date<datenum(2017,9,22)
    v = 3;
    Level_0_cleared.Properties.VariableNames{25} = 'Channel13'; % there was a mixup for the naming of ARISTO
    Level_0_cleared.Properties.VariableNames{27} = 'Channel15';
    Level_0.Properties.VariableNames{25} = 'Channel13'; % there was a mixup for the naming of ARISTO
    Level_0.Properties.VariableNames{27} = 'Channel15';
elseif PHIPS_date>=datenum(2017,2,14) & PHIPS_date<datenum(2017,2,23)
    v = 2;
else
    v = 1;
end
Level_0_ordered = PHIPS_PMTchannel_to_angle(Level_0_cleared,v);
Level_0_ordered = PHIPS_PMTchannel_to_angle(Level_0,v);

%% put the ASF based on the cleared lvl0 back into the lvl1
Level_1_new = Level_1_old;
Level_1_new(:, end-19:end) = Level_0_ordered(:, end-19:end);

PhipsData = Level_1_new;



end

%% remove shattered
function [PhipsData,shattering_percentage] = PHIPS_remove_shattered(PhipsData,threshold)
% Removes particles that are around interarrival times smaller than
% threshold (threshold can be calculated by shattering_threshold.m)
% else, 0.5 ms is a good guess


%angle = [1	2	3	4	5	6	7	8	9	10	18	26	34	42	50	58	66	74	82	90	98	106	114	122	130	138	146	154	162	170];


%% Length of raw_data
% length_old_data = height(PhipsData(PhipsData.ForcedParticleFlag==0,:)); %
length_old_data = height(PhipsData); %because in the new lvl2+ files there is no FPF

%% Interarrival times
% Find force triggers
time_in_s = PhipsData.ParticleTimeStamp./48e6;
interarrivaltime = diff(time_in_s); % in seconds
ind = find(interarrivaltime<threshold);   

%% Histogram of interarrival times
figure,
histogram(interarrivaltime.*1000,logspace(-3,5,100)); hold on
set(gca,'FontSize',14,'linewidth',2,'xscale','log','yscale','log')
xlabel('Interarrival time [ms]','fontsize',16)
title('Histogram of Interarrval Times','fontsize',16)
yl =ylim; plot([threshold.*1000 threshold.*1000],[1 yl(2)],'r','linewidth',2)
grid on


%% Mark entries before and after short interarrival times as NaN
PhipsData.ParticleTimeStamp(ind,1) = NaN;
PhipsData.ParticleTimeStamp(ind+1,1) = NaN;
shattered_raw_data = PhipsData(isnan(PhipsData.ParticleTimeStamp),:);
PhipsData(isnan(PhipsData.ParticleTimeStamp),:) = [];
% length_new_data = height(PhipsData(PhipsData.ForcedParticleFlag==0,:));
length_new_data = height(PhipsData);
disp([num2str(100.*(1-round(length_new_data./length_old_data,3))),' % of data removed due to shattering'])

% disp([num2str(num_prev-num_corr), ' out of ', num2str(num_corr), ' = ', num2str((num_prev-num_corr)/num_corr*100), '% of particles removed due to shattering'])


%% How many PARTICLES (not data points) were removed? (e.g. 3 successive shattering events are assigned as 1 particle
%about 30% of shattering events shatter in more than 2 fragments

% num = 0;
% for i = 1:length(ind)-1
%     if ind(i+1)~=ind(i)+1
%         num=num+1; %number of removed particles
%     end
% end
% 
% shattering_percentage = num/length_old_data*100; % x percent of actual particles were removed due to shattering --> this is the correction factor
% disp([num2str(round(shattering_percentage,2)),' % of actual particles removed due to shattering (= correction factor)'])


%% Figure
% figure,
% semilogy(angle(end-19:end),nanmean(table2array(shattered_raw_data(:,end-19:end)),1),'sr')
% set(gca,'fontsize',14,'linewidth',2)
% hold on
% semilogy(angle(end-19:end),nanmean(table2array(PhipsData(:,end-19:end)),1),'sk')
% ylim([1 2e3]); legend('Shattered SPF','Accepted SPF'); title([num2str(100.*(1-round(length_new_data./length_old_data,3))),' % of data removed due to shattering'],'fontsize',18)

    

end



%% PHIPS_PMTchannel_to_angle function
function [TABLE] = PHIPS_PMTchannel_to_angle(PhipsData,v)
% Converts PMT channels to PHIPS angle
% v is version number
% v1 -> -13.2.2017
% v2 -> 14.2.2017-22.2.2017
% v3 -> 23.2.2017-21.9.2017
% v4 -> 22.9.2017-03.06.2019
% v5 -> 04.06.2019-08.11.2021
% v6 -> 09.11.2021-

%% Make an empty table
TABLE  = cell2table(cell(size(PhipsData,1),20+13), 'VariableNames', {'DataSet', 'RealTimeStamp', 'ParticleTimeStamp', 'ElapsedParticleTime', 'ParticleTimeOfFlight', 'CameraTrigger', 'ForcedParticleFlag', 'MultiplexerValue', 'ImageCamera1','TimeStampCamera1', 'ImageCamera2','TimeStampCamera2','TriggerIntensity','ScatteringAngle18','ScatteringAngle26','ScatteringAngle34','ScatteringAngle42','ScatteringAngle50','ScatteringAngle58','ScatteringAngle66','ScatteringAngle74','ScatteringAngle82','ScatteringAngle90','ScatteringAngle98','ScatteringAngle106','ScatteringAngle114','ScatteringAngle122','ScatteringAngle130','ScatteringAngle138','ScatteringAngle146','ScatteringAngle154','ScatteringAngle162','ScatteringAngle170'});
TABLE.DataSet = PhipsData.DataSet;
TABLE.RealTimeStamp = PhipsData.RealTimeStamp;
TABLE.ParticleTimeStamp = PhipsData.ParticleTimeStamp;
TABLE.ElapsedParticleTime = PhipsData.ElapsedParticleTime;
TABLE.ParticleTimeOfFlight = PhipsData.ParticleTimeOfFlight;
TABLE.CameraTrigger = PhipsData.CameraTrigger;
TABLE.ForcedParticleFlag = PhipsData.ForcedParticleFlag;
TABLE.MultiplexerValue = PhipsData.MultiplexerValue;
TABLE.ImageCamera1 = PhipsData.ImageCamera1;
TABLE.ImageCamera2 = PhipsData.ImageCamera2;
TABLE.TimeStampCamera1 = PhipsData.TimeStampCamera1;
TABLE.TimeStampCamera2 = PhipsData.TimeStampCamera2;
TABLE.TriggerIntensity = PhipsData.TriggerIntensity;

if v == 1
TABLE.ScatteringAngle18 = PhipsData.Channel17;
TABLE.ScatteringAngle26 = PhipsData.Channel19;
TABLE.ScatteringAngle34 = PhipsData.Channel2;
TABLE.ScatteringAngle42 = PhipsData.Channel3;
TABLE.ScatteringAngle50 = PhipsData.Channel4;
TABLE.ScatteringAngle58 = PhipsData.Channel5;
TABLE.ScatteringAngle66 = PhipsData.Channel6;
TABLE.ScatteringAngle74 = PhipsData.Channel7;
TABLE.ScatteringAngle82 = PhipsData.Channel8;
TABLE.ScatteringAngle90 = PhipsData.Channel9;
TABLE.ScatteringAngle98 = PhipsData.Channel10;
TABLE.ScatteringAngle106 = PhipsData.Channel11;
TABLE.ScatteringAngle114 = PhipsData.Channel12;
TABLE.ScatteringAngle122 = PhipsData.Channel13;
TABLE.ScatteringAngle130 = PhipsData.Channel14;
TABLE.ScatteringAngle138 = PhipsData.Channel15;
TABLE.ScatteringAngle146 = PhipsData.Channel23;
TABLE.ScatteringAngle154 = PhipsData.Channel20;
TABLE.ScatteringAngle162 = PhipsData.Channel21;
TABLE.ScatteringAngle170 = PhipsData.Channel22;

elseif v == 2
TABLE.ScatteringAngle18 = PhipsData.Channel25+PhipsData.Channel26;
TABLE.ScatteringAngle26 = PhipsData.Channel24;
TABLE.ScatteringAngle34 = PhipsData.Channel27+PhipsData.Channel28;
TABLE.ScatteringAngle42 = PhipsData.Channel29+PhipsData.Channel30;
TABLE.ScatteringAngle50 = PhipsData.Channel31+PhipsData.Channel32;
TABLE.ScatteringAngle58 = PhipsData.Channel11;
TABLE.ScatteringAngle66 = PhipsData.Channel8;
TABLE.ScatteringAngle74 = PhipsData.Channel6;
TABLE.ScatteringAngle82 = PhipsData.Channel4;
TABLE.ScatteringAngle90 = PhipsData.Channel2;
TABLE.ScatteringAngle98 = PhipsData.Channel3;
TABLE.ScatteringAngle106 = PhipsData.Channel5;
TABLE.ScatteringAngle114 = PhipsData.Channel7;
TABLE.ScatteringAngle122 = PhipsData.Channel9;
TABLE.ScatteringAngle130 = PhipsData.Channel10;
TABLE.ScatteringAngle138 = PhipsData.Channel12;
TABLE.ScatteringAngle146 = PhipsData.Channel13;
TABLE.ScatteringAngle154 = PhipsData.Channel14;
TABLE.ScatteringAngle162 = PhipsData.Channel15;
TABLE.ScatteringAngle170 = PhipsData.Channel23;

elseif v == 3
TABLE.ScatteringAngle18 = PhipsData.Channel25+PhipsData.Channel26;
TABLE.ScatteringAngle26 = PhipsData.Channel20;
TABLE.ScatteringAngle34 = PhipsData.Channel27+PhipsData.Channel28;
TABLE.ScatteringAngle42 = PhipsData.Channel29+PhipsData.Channel30;
TABLE.ScatteringAngle50 = PhipsData.Channel31+PhipsData.Channel32;
TABLE.ScatteringAngle58 = PhipsData.Channel11;
TABLE.ScatteringAngle66 = PhipsData.Channel8;
TABLE.ScatteringAngle74 = PhipsData.Channel6;
TABLE.ScatteringAngle82 = PhipsData.Channel4;
TABLE.ScatteringAngle90 = PhipsData.Channel2;
TABLE.ScatteringAngle98 = PhipsData.Channel3;
TABLE.ScatteringAngle106 = PhipsData.Channel5;
TABLE.ScatteringAngle114 = PhipsData.Channel7;
TABLE.ScatteringAngle122 = PhipsData.Channel9;
TABLE.ScatteringAngle130 = PhipsData.Channel10;
TABLE.ScatteringAngle138 = PhipsData.Channel12;
TABLE.ScatteringAngle146 = PhipsData.Channel13;
TABLE.ScatteringAngle154 = PhipsData.Channel14;
TABLE.ScatteringAngle162 = PhipsData.Channel15;
TABLE.ScatteringAngle170 = PhipsData.Channel23;

elseif v == 4
TABLE.ScatteringAngle18 = PhipsData.Channel25+PhipsData.Channel26;
TABLE.ScatteringAngle26 = PhipsData.Channel20;
TABLE.ScatteringAngle34 = PhipsData.Channel17;
TABLE.ScatteringAngle42 = PhipsData.Channel29+PhipsData.Channel30;
TABLE.ScatteringAngle50 = PhipsData.Channel31+PhipsData.Channel32;
TABLE.ScatteringAngle58 = PhipsData.Channel11;
TABLE.ScatteringAngle66 = PhipsData.Channel8;
TABLE.ScatteringAngle74 = PhipsData.Channel6;
TABLE.ScatteringAngle82 = PhipsData.Channel4;
TABLE.ScatteringAngle90 = PhipsData.Channel2;
TABLE.ScatteringAngle98 = PhipsData.Channel3;
TABLE.ScatteringAngle106 = PhipsData.Channel5;
TABLE.ScatteringAngle114 = PhipsData.Channel7;
TABLE.ScatteringAngle122 = PhipsData.Channel9;
TABLE.ScatteringAngle130 = PhipsData.Channel10;
TABLE.ScatteringAngle138 = PhipsData.Channel12;
TABLE.ScatteringAngle146 = PhipsData.Channel13;
TABLE.ScatteringAngle154 = PhipsData.Channel14;
TABLE.ScatteringAngle162 = PhipsData.Channel15;
TABLE.ScatteringAngle170 = PhipsData.Channel23;

elseif v == 5
    % empty: 14, 22, 24(defect), 26-31(free for pol.) 32(defect)
TABLE.ScatteringAngle18 = PhipsData.Channel25;
TABLE.ScatteringAngle26 = PhipsData.Channel23;
TABLE.ScatteringAngle34 = PhipsData.Channel21;
TABLE.ScatteringAngle42 = PhipsData.Channel20;
TABLE.ScatteringAngle50 = PhipsData.Channel19;
TABLE.ScatteringAngle58 = PhipsData.Channel18;
TABLE.ScatteringAngle66 = PhipsData.Channel17;
TABLE.ScatteringAngle74 = PhipsData.Channel16;
TABLE.ScatteringAngle82 = PhipsData.Channel15;

TABLE.ScatteringAngle90 = PhipsData.Channel3;
TABLE.ScatteringAngle98 = PhipsData.Channel4;
TABLE.ScatteringAngle106 = PhipsData.Channel5;
TABLE.ScatteringAngle114 = PhipsData.Channel6;
TABLE.ScatteringAngle122 = PhipsData.Channel7;
TABLE.ScatteringAngle130 = PhipsData.Channel8;
TABLE.ScatteringAngle138 = PhipsData.Channel9;
TABLE.ScatteringAngle146 = PhipsData.Channel10;
TABLE.ScatteringAngle154 = PhipsData.Channel11;
TABLE.ScatteringAngle162 = PhipsData.Channel12;
TABLE.ScatteringAngle170 = PhipsData.Channel13;

elseif v == 6
    % 2, 3, 14, 20, 32 (faulty); 28, 29 free but could be used for cross
    % pol. with 13, 15; 26, 27, 30, 31 (unused) 
TABLE.ScatteringAngle18  =  PhipsData.Channel25;
TABLE.ScatteringAngle26  =  PhipsData.Channel24;
TABLE.ScatteringAngle34  =  PhipsData.Channel23;
TABLE.ScatteringAngle42  =  PhipsData.Channel22;
TABLE.ScatteringAngle50  =  PhipsData.Channel21;
TABLE.ScatteringAngle58  =  PhipsData.Channel19;
TABLE.ScatteringAngle66  =  PhipsData.Channel18;
TABLE.ScatteringAngle74  =  PhipsData.Channel17;
TABLE.ScatteringAngle82  =  PhipsData.Channel16;

TABLE.ScatteringAngle90  =  PhipsData.Channel4;
TABLE.ScatteringAngle98  =  PhipsData.Channel5;
TABLE.ScatteringAngle106 =  PhipsData.Channel6;
TABLE.ScatteringAngle114 =  PhipsData.Channel7;
TABLE.ScatteringAngle122 =  PhipsData.Channel8;
TABLE.ScatteringAngle130 =  PhipsData.Channel9;
TABLE.ScatteringAngle138 =  PhipsData.Channel10;
TABLE.ScatteringAngle146 =  PhipsData.Channel11;
TABLE.ScatteringAngle154 =  PhipsData.Channel12;
TABLE.ScatteringAngle162 =  PhipsData.Channel13;
TABLE.ScatteringAngle170 =  PhipsData.Channel15;

elseif v == 7
    % 2, 3, 14, 15, 20, 32 (faulty); 28, 29 free but could be used for cross
    % pol. with 13, 15; 26, 30, 31 (unused) 
TABLE.ScatteringAngle18  =  PhipsData.Channel25;
TABLE.ScatteringAngle26  =  PhipsData.Channel24;
TABLE.ScatteringAngle34  =  PhipsData.Channel23;
TABLE.ScatteringAngle42  =  PhipsData.Channel22;
TABLE.ScatteringAngle50  =  PhipsData.Channel21;
TABLE.ScatteringAngle58  =  PhipsData.Channel19;
TABLE.ScatteringAngle66  =  PhipsData.Channel18;
TABLE.ScatteringAngle74  =  PhipsData.Channel17;
TABLE.ScatteringAngle82  =  PhipsData.Channel16;

TABLE.ScatteringAngle90  =  PhipsData.Channel4;
TABLE.ScatteringAngle98  =  PhipsData.Channel5;
TABLE.ScatteringAngle106 =  PhipsData.Channel6;
TABLE.ScatteringAngle114 =  PhipsData.Channel7;
TABLE.ScatteringAngle122 =  PhipsData.Channel8;
TABLE.ScatteringAngle130 =  PhipsData.Channel9;
TABLE.ScatteringAngle138 =  PhipsData.Channel10;
TABLE.ScatteringAngle146 =  PhipsData.Channel11;
TABLE.ScatteringAngle154 =  PhipsData.Channel12;
TABLE.ScatteringAngle162 =  PhipsData.Channel13;
TABLE.ScatteringAngle170 =  PhipsData.Channel27;

end

end
