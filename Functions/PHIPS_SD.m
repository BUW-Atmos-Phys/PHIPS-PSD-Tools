%% Revision history
% 12.03.2025, EJ, added mass-equivalent spherical radius 
% 04.03.2025, EJ, added PSD for only imaged particles
% 05.04.2024, EJ, added one bin for 15 to 30 micron.
% 01.06.2023, EJ, added extinction coefficient, IWC and effective radius.
% 31.05.2023, EJ, added PSD based on projected area.
% 30.05.2023, EJ, added final CIRRUS-HL calib. coeffs.
% 24.05.2023, EJ, removed the replacement of Level0 scattering data
% 17.05.2021, FW, adjust for updated PO structure and added ASF from lvl 0
% 08.12.2021, FW, sens volume is calculated for each particle (not just once per bin)
% 07.02.2022, FW, adjusted paths and improved documentation

function [SDice,SDdroplet,SDice_v1,SDdroplet_v1,SDice_area,Bext,IWC] = ...
    PHIPS_SD(Aircraft_data_folder, folder, savepath, SD_2DS_path,SD_2DC_path,...
    campaign,flight,tstep,save_status,start_time,end_time)
% PHIPS size distribution from scattering data
% Input arguments: Aircraft_data_folder -> path to aircraft data
%                  folder -> folder to PHIPS level files
%                  savepath -> path where PSDs are saved
%                  SD_2DS_path -> path to Wei's 2DS and 2DC data
%                  flight -> flight name
%                  campaign -> campaign name, e.g. 'SOCRATES'
%                  tstep -> time step in seconds
%                  t_start/end -> optional

% ADJUST BINS ALSO IN Read_PHIPS_SD.m for translate_SD_sum_to_nc.m
bin_endpoints = [15,30,60,100,150,200,250,300,350,400,500,600,700];
% dN_ice = -999; dN_drop = -999;

%savepath = [particleopticspath,filesep,'PHIPS Results', filesep, 'Campaigns', filesep,campaign,filesep,flight, filesep, 'SD', filesep];
%savepath = [particleopticspath,filesep,flight, filesep, 'SD', filesep]; % POISTA
if ~isfolder(savepath)
    mkdir(savepath)
end

% Shattering threshold
if strcmp(campaign,'CIRRUS-HL')
    thresh_interarrivaltime = 0.1/1000; % = 0.1 ms
else
    thresh_interarrivaltime = 0.5/1000; % = 0.5 ms
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
    a_droplet = 1.212;
    b_droplet = 0.5;
    % Ice
    a_ice = 0.21589;
    b_ice = 0.56378;
    a_ice_area = 0.020561;
    b_ice_area = 1.1746;
    
    BG_offset = 250.6858;
    
else
    disp('Calibration coefficients unknown. PSDs could not be produced.');
    return
end

%% Folder to PHIPS data
%folder = [particleopticspath,filesep,'Phips',filesep,campaign,filesep,flight];
%folder = [particleopticspath,filesep,'PHIPS Results',filesep,flight]; % POISTA


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

%% 0) Define bins
bin_midpoints = bin_endpoints(1:end-1)+diff(bin_endpoints)./2;

% % dlogDp
dlogDp = zeros(1,length(bin_endpoints)-1);
for i=1:length(bin_endpoints)-1
    dlogDp(i) = log10(bin_endpoints(i+1))- log10(bin_endpoints(i));  % =/= log dDp !
end


%% A.1) Load Data
% PhipsData Lvl 0, Lvl 3 and Lvl 5
% Level 0 is needed for dead time calculation
cd(folder)

% Level 5
listings = dir('*level_5.csv');   
if isempty(listings)
    disp('Level 5 File not found!')
    % Check if level 3 exists
    listings = dir('*level_3.csv');   
    if isempty(listings)
        disp('Level 3 File not found! PSDs cannot be produced.')
        SDice = NaN;
        SDdroplet = NaN;
        return
    end
    SDice_only_images = NaN;
    Reff_only_images = NaN;
end
csv_name = listings(end).name; filename = [folder,filesep,csv_name];
PhipsData = readtable(filename); PhipsData.RealTimeStamp = datenum(PhipsData.RealTimeStamp);

% Level 0
listings = dir('*level_0.csv');
csv_name = listings(end).name; filename = [folder,filesep,csv_name];
PhipsData_level_0 = readtable(filename); PhipsData_level_0.RealTimeStamp = datenum(PhipsData_level_0.RealTimeStamp);

% Replace NaNs in scattering data with 2048
A = find(strcmp(PhipsData.Properties.VariableNames, 'ScatteringAngle18'));
for i = A: width(PhipsData)
    PhipsData.(i)(isnan(PhipsData.(i))) = 2048;
end  

% Load Aircraft Data
[time,airspeed] = get_airspeed(campaign,flight,Aircraft_data_folder);
       
% % 1.b) Integrate Scattering Intensity
% Get index of starting angle (42 degrees)
A = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle42'));
B = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle170'));
% Integrate
ang_temp = [42:8:170];
SPF_temp = table2array(PhipsData(:,A:B));
integrated_intensity = trapz(ang_temp, SPF_temp,2);
PhipsData = addvars(PhipsData,integrated_intensity,'NewVariableNames','integrated_intensity');

% % 1.c) Get diameter and projecter area from images
[PhipsData] = get_image_size(PhipsData);
        

%% A.2) Remove saturates
% Duplicate DropletFlag
PhipsData = addvars(PhipsData, PhipsData.DropletFlag, 'Before', 'DropletFlag', 'NewVariableNames', 'DropletFlag_image');

% if more than 4 channels are saturated/background, remove the particle (DropletFlag = NaN)
A = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle42'));
B = find(strcmpi(PhipsData.Properties.VariableNames,'ScatteringAngle170'));
PhipsArray = table2array(PhipsData(:,A:B));

sat_tresh = 4; % saturated if >= X chanels are saturated
bg_tresh = 4;

num_sat = 0; num_BG = 0;
for i = 1:size(PhipsArray,1)
    if sum(PhipsArray(i, :) >= 2047) >= sat_tresh && ignore_above_saturation == 1
        PhipsData.DropletFlag(i) = NaN; %erase those particles from Drop and Ice, but they still appear in Total
        num_sat = num_sat + 1;
    end
    
    if sum(PhipsArray(i, :) <= 5) >= bg_tresh && ignore_below_background == 1
        PhipsData.DropletFlag(i) = NaN; %erase those particles from Drop and Ice, but they still appear in Total
        num_BG = num_BG + 1;
    end
end
disp([num2str(num_sat), ' particles above saturation were removed.'])
%  disp([num2str(num_BG), ' particles below BG removed'])

% Mark DropletFlag_image as NaN where DropletFlag = NaN and dp_image = NaN
PhipsData.DropletFlag_image(isnan(PhipsData.DropletFlag) & isnan(PhipsData.dp_image)) = NaN;        
num_sat = length(find(isnan(PhipsData.DropletFlag_image)));     
disp([num2str(num_sat), ' particles above saturation without image size were removed.'])


%% B.1. Generate histograms
% % B.1.a generate time axis
if nargin == 11 %careful when debugging this step by step: nargin command does only work if used as function
    start_time = start_time + datenum(0,0,0,0,0,tstep/2); %to compensate for the +/-tstep/2 in the next step
    end_time = end_time - datenum(0,0,0,0,0,tstep/2);
    taxis = [start_time - datenum(0,0,0,0,0,tstep/2):datenum(0,0,0,0,0,tstep):end_time+datenum(0,0,0,0,0,tstep/2)];
else %if no start/end time is given in the input of the function
    % Start from next full minute minus half of timestep
    temp = datevec(PhipsData.RealTimeStamp(1));
    taxis = [datenum(temp(1),temp(2),temp(3),temp(4),temp(5),0+tstep/2):datenum(0,0,0,0,0,tstep):PhipsData.RealTimeStamp(end)];
end
        
% % B.1.b Shattering correction based on interarrivaltime
[PhipsData] = PHIPS_remove_shattered(PhipsData,thresh_interarrivaltime);

% % Calculate Shattering Flag
if strcmp(campaign, 'SOCRATES')
    ShatteringFlag = PHIPS_shattering_flag(SD_2DS_path,SD_2DC_path, flight, taxis);
else
    ShatteringFlag = zeros(length(taxis)-1, 1)+1;
end
        
% % B.1.c add scattering size and area
dp_ice = abs(a_ice .*(PhipsData.integrated_intensity(PhipsData.DropletFlag==0)-BG_offset).^(b_ice));
dp_drop = abs(a_droplet .*(PhipsData.integrated_intensity(PhipsData.DropletFlag==1)-BG_offset).^(b_droplet));
A_ice = abs(a_ice_area .*(PhipsData.integrated_intensity(PhipsData.DropletFlag==0)-BG_offset).^(b_ice_area));
PhipsData.dp (PhipsData.DropletFlag==0) = dp_ice;
PhipsData.dp (PhipsData.DropletFlag==1) = dp_drop;
PhipsData.dp (isnan(PhipsData.DropletFlag)) = NaN;
PhipsData.A (PhipsData.DropletFlag==0) = A_ice;
PhipsData.A (PhipsData.DropletFlag==1) = NaN;
PhipsData.A (isnan(PhipsData.DropletFlag)) = NaN;

% % B.1.d make histograms
disp('Calculate SD (this may take a while...)')
disp('... for ice...')
% Ice histogram based on scattering size
[particles_per_bin_ice, conc_ice, ~, ~] = generate_histogram(PhipsData(PhipsData.DropletFlag==0,:),...
    PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_ice,airspeed,time);
% Droplet histogram based on scattering size
[particles_per_bin_drop, conc_drop, ~, ~] = generate_histogram(PhipsData(PhipsData.DropletFlag==1,:),...
    PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_drop,airspeed,time);
conc_drop(:,1) = NaN; % delete the first bin for droplets as the lower limit is 50
% Ice histogram based on imaged size
num1 = length(find(~isnan(PhipsData.dp))); num2 = length(find(~isnan(PhipsData.dp_image)));
PhipsData.dp(~isnan(PhipsData.dp_image)) = PhipsData.dp_image(~isnan(PhipsData.dp_image)); % replace scattering dp with image dp
PhipsData.A(~isnan(PhipsData.A_image)) = PhipsData.A_image(~isnan(PhipsData.A_image)); % replace scattering projected area with image projected area
disp([num2str(round(num2/num1*100,2)),' % particles had image size'])
[particles_per_bin_ice_v1, conc_ice_v1, dA, dm] = generate_histogram(PhipsData(PhipsData.DropletFlag_image==0,:),...
    PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_ice,airspeed,time);
disp('... for droplets...')
% Droplet histogram based on imaged size
[particles_per_bin_drop_v1, conc_drop_v1, ~, ~] = generate_histogram(PhipsData(PhipsData.DropletFlag_image==1,:),...
    PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_ice,airspeed,time);
% Ice histogram ONLY FOR IMAGED PARTICLES
if ~exist("SDice_only_images")
    PhipsData = PhipsData(PhipsData.Shattering==0 & PhipsData.Multiple==0,:);
    [particles_per_bin_ice_only_images, conc_ice_only_images, dA_only_images, dm_only_images] = generate_histogram(PhipsData(PhipsData.Droplet==0,:),...
        PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters_ice,airspeed,time);
end


        
%% C.1 Produce PHIPS_counts SD
% v1 is based on image size
SDice_counts = [0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag sum(particles_per_bin_ice,2) particles_per_bin_ice];
SDice_counts_v1 = [0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag sum(particles_per_bin_ice_v1,2) particles_per_bin_ice_v1];
SDdroplet_counts = [0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 sum(particles_per_bin_drop,2) particles_per_bin_drop]; % SF for droplet is always 1
SDdroplet_counts_v1 = [0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 sum(particles_per_bin_drop_v1,2) particles_per_bin_drop_v1]; % SF for droplet is always 1
if ~exist("SDice_only_images")
SDice_counts_only_images = [0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag sum(particles_per_bin_ice_only_images,2) particles_per_bin_ice_only_images];
end        

% % estimate uncertainty (FOR THE FIRST BIN - THIS NEEDS TO BE DONE!)
% ice
dN_minus =  [NaN -24.9479  -24.0072  -18.7372  -17.7478  -14.5033  -11.3761   -7.8475   -8.1237   -5.9329   -5.2310   -3.8175];
dN_plus =   [NaN 45.5047   44.6152   51.4139   48.3331   32.1553   23.5189   23.0540   17.5292   12.1717    6.3060    7.0891];
dN_ice = (dN_plus - dN_minus)./2; % [%]
% droplets
dN_minus = [NaN -39.8952   -9.7177  -10.0118  -10.9244   -9.2650  -13.1825  -11.8750   -3.4833   -8.8068       NaN       NaN];
dN_plus =   [NaN 42.5975   12.1592   12.7038   17.2285   14.0938   22.1102   35.3615   15.7782    9.8641       NaN       NaN];
dN_drop = (dN_plus - dN_minus)./2;
dN_drop(end-1:end) = dN_drop(end-2); % set the last 2 bins, which don't have enough data, on the last valid point
        
% Total concentration
Ntot_ice = 1000.*nansum(conc_ice,2); % [L^-1]
Ntot_ice_v1 = 1000.*nansum(conc_ice_v1,2); % [L^-1]
if ~exist("SDice_only_images")
    Ntot_ice_only_images = 1000.*nansum(conc_ice_only_images,2); % [L^-1]
    Ntot_ice_area_only_images = 1000.*nansum(dA_only_images,2); % [um2 L^-1]
end
Ntot_ice_area = 1000.*nansum(dA,2); % [um2 L^-1]
Ntot_ice_mass = 1000.*nansum(dm,2); % [mg L^-1]
Ntot_droplet =   1000.*nansum(conc_drop,2);
Ntot_droplet_v1 =   1000.*nansum(conc_drop_v1,2);

Ntot_err_ice = 1000.*nansum(conc_ice.*dN_ice./100,2); % [L^-1]   
Ntot_err_ice_v1 = 1000.*nansum(conc_ice_v1.*dN_ice./100,2); % [L^-1]  
Ntot_err_ice_only_images = 1000.*nansum(conc_ice_only_images.*dN_ice./100,2); % [L^-1]  
Ntot_err_droplet = 1000.*nansum(conc_drop.*dN_drop./100,2);
Ntot_err_droplet_v1 = 1000.*nansum(conc_drop_v1.*dN_drop./100,2);

% dNdlogDp
dNdlogDp_ice = zeros(size(conc_ice));
for i=1:size(conc_ice(:,1),1)
    dNdlogDp_ice(i,:) = 1000.*(conc_ice(i,:)./dlogDp);
end
dNdlogDp_ice_v1 = zeros(size(conc_ice_v1));
for i=1:size(conc_ice_v1(:,1),1)
    dNdlogDp_ice_v1(i,:) = 1000.*(conc_ice_v1(i,:)./dlogDp);
end
if ~exist("SDice_only_images")
    dNdlogDp_ice_only_images = zeros(size(conc_ice_only_images));
    for i=1:size(conc_ice_only_images(:,1),1)
        dNdlogDp_ice_only_images(i,:) = 1000.*(conc_ice_only_images(i,:)./dlogDp);
    end
end
dAdlogDp_ice = zeros(size(dA));
for i=1:size(dA(:,1),1)
    dAdlogDp_ice(i,:) = 1000.*(dA(i,:)./dlogDp);
end
if ~exist("SDice_only_images")
    dAdlogDp_ice_only_images = zeros(size(dA_only_images));
    for i=1:size(dA_only_images(:,1),1)
        dAdlogDp_ice_only_images(i,:) = 1000.*(dA_only_images(i,:)./dlogDp);
    end
end
dmdlogDp_ice = zeros(size(dm));
for i=1:size(dm(:,1),1)
    dmdlogDp_ice(i,:) = 1000.*(dm(i,:)./dlogDp);
end

dNdlogDp_droplet = zeros(size(conc_drop));
for i=1:size(conc_drop(:,1),1)
    dNdlogDp_droplet(i,:) = 1000.*(conc_drop(i,:)./dlogDp);
end
dNdlogDp_droplet_v1 = zeros(size(conc_drop_v1));
for i=1:size(conc_drop_v1(:,1),1)
    dNdlogDp_droplet_v1(i,:) = 1000.*(conc_drop_v1(i,:)./dlogDp);
end

       
% % Make SD files
only_images = 0;
SDice = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice Ntot_err_ice dNdlogDp_ice];
SDdroplet = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 Ntot_droplet Ntot_err_droplet dNdlogDp_droplet]; % SF for droplet is always 1
SDice_v1 = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice_v1 Ntot_err_ice_v1 dNdlogDp_ice_v1];
SDdroplet_v1 = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag.*0+1 Ntot_droplet_v1 Ntot_err_droplet_v1 dNdlogDp_droplet_v1]; % SF for droplet is always 1
SDice_area = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice_area Ntot_err_ice_v1.*NaN dAdlogDp_ice]; % Error estimation not yet performed
SDice_mass = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice_mass Ntot_err_ice_v1.*NaN dmdlogDp_ice]; % Error estimation not yet performed
if ~exist("SDice_only_images")
    SDice_only_images = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice_only_images Ntot_err_ice_only_images dNdlogDp_ice_only_images];
    SDice_area_only_images = [0 0 0 0 bin_midpoints;...
    taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) ShatteringFlag Ntot_ice_area_only_images Ntot_err_ice_only_images.*NaN dAdlogDp_ice_only_images]; % Error estimation not yet performed
    only_images = 1;
end
            
% % check how many segments are SF = 0
idx = find(SDice(:, 3) > 0); % Ntot > 0
SF0 = length(find(SDice(idx,2) == 0));
SF1 = length(find(SDice(idx,2) == 1));

disp([num2str(SF0), ' segments are SF=0, ', num2str(SF1), ' are SF = 1, that means that ', num2str(SF0/(SF0+SF1)*100), '% are rejected)'])
        

%% D) Calculate extinction coefficient, IWC, Reff, Rm
b_ext = 2.*nansum(dA,2) .* 1e-18 ./ 1e-15; % 1/km
if only_images == 1
    b_ext_only_images = 2.*nansum(dA_only_images,2) .* 1e-18 ./ 1e-15; % 1/km
end

% Calculate IWC
iwc = nansum(dm,2); % mg / cm3
iwc = 1e-3./1e-6.*iwc; % g / m3
if only_images == 1
    iwc_only_images = nansum(dm_only_images,2); % mg / cm3
    iwc_only_images = 1e-3./1e-6.*iwc_only_images; % g / m3
end

% Calculate effective radius
rho = 0.917; % density of ice in g cm-3
rho = rho.*1e6; % density of ice in g m-3
b_ext_temp = b_ext .* 1e-3; % m-1
reff = 3/2 .* iwc ./ rho ./ b_ext_temp; % [g m-3 / (g m-3) / m-1 -> m]
reff = reff .* 1e6; % um
% For only images
if only_images == 1
    b_ext_temp = b_ext_only_images .* 1e-3; % m-1
    reff_only_images = 3/2 .* iwc_only_images ./ rho ./ b_ext_temp; % [g m-3 / (g m-3) / m-1 -> m]
    reff_only_images = reff_only_images .* 1e6; % um
end

% Calculate mass-equivalent spherical radius
rho_ice = 917;               % Density of ice in kg/m³
iwc_kg = iwc    .* 1e-3; % Convert g/m³ to kg/m³
N = Ntot_ice_v1 .* 1e3;  % Convert L⁻¹ to m⁻³
% Initialize r_m with NaN
r_m = NaN(size(iwc_kg));
% Find valid indices where IWC > 0 and N > 0 to avoid division errors
valid_indices = (iwc_kg > 0) & (N > 0);
% Compute mean ice crystal radius
r_m(valid_indices) = (3 * iwc_kg(valid_indices) ./ (4 * pi * N(valid_indices) * rho_ice)).^(1/3);

% Calculate mass-equivalent spherical radius, for only imaged particles
if only_images == 1
    iwc_kg = iwc_only_images    .* 1e-3; % Convert g/m³ to kg/m³
    N = Ntot_ice_only_images    .* 1e3;  % Convert L⁻¹ to m⁻³
    % Initialize r_m with NaN
    r_m_only_images = NaN(size(iwc_kg));
    % Find valid indices where IWC > 0 and N > 0 to avoid division errors
    valid_indices = (iwc_kg > 0) & (N > 0);
    % Compute mean ice crystal radius
    r_m_only_images(valid_indices) = (3 * iwc_kg(valid_indices) ./ (4 * pi * N(valid_indices) * rho_ice)).^(1/3);
end

Bext = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) b_ext];
IWC = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) iwc];
Reff = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) reff];
Rm = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) r_m.*1e6]; % in micron
if only_images == 1
    Reff_only_images = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) reff_only_images];
    Rm_only_images = [taxis(1:end-1)'+datenum(0,0,0,0,0,tstep/2) r_m_only_images.*1e6]; % in micron
end


%% E) Save SD
if save_status == 1
disp('save...')
    save_SD(SDice, SDice_counts, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bin_endpoints);
    save_SD(SDice_v1, SDice_counts_v1, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bin_endpoints);
    save_SD(SDdroplet, SDdroplet_counts, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bin_endpoints);
    save_SD(SDdroplet_v1, SDdroplet_counts_v1, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bin_endpoints);
    save_SD(SDice_area, SDice_counts_v1, save_status, savepath, campaign, flight, tstep, a_ice_area, b_ice_area, bin_endpoints);
    
    save_SD(SDice_only_images, SDice_counts_only_images, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bin_endpoints);
    save_SD(SDice_area_only_images, SDice_counts_only_images, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bin_endpoints);

    save([savepath,'Bext_',num2str(tstep),'s.txt'],'Bext','-ascii','-double')
    save([savepath,'IWC_',num2str(tstep),'s.txt'],'IWC','-ascii','-double')
    save([savepath,'Reff_',num2str(tstep),'s.txt'],'Reff','-ascii','-double')
    save([savepath,'Reff_only_images_',num2str(tstep),'s.txt'],'Reff_only_images','-ascii','-double')
    save([savepath,'Rm_',num2str(tstep),'s.txt'],'Reff','-ascii','-double')
    save([savepath,'Rm_only_images_',num2str(tstep),'s.txt'],'Rm_only_images','-ascii','-double')

    % save SD as .nc file
    % translate_SD_sum_to_nc(savepath, campaign, flight, tstep);
end
disp('..all done.')

end


%% FUNCTIONS 

%% Get Airspeed
function [time,airspeed] = get_airspeed(campaign,flight,Aircraft_data_folder)
    if strcmp(campaign,'ACLOUD')
        %Aircraft_data_folder = [particleopticspath, '/SID3/ACLOUD/Flight data/'];
        flight_date = strsplit(flight, ' ');
        flight_date = flight_date{2};
        Aircraft_data_filename = ['AircraftData_', flight_date, '.dat'];
        [time,airspeed,~,~,~,~,~,~] = Aircraft_data(campaign,Aircraft_data_folder,Aircraft_data_filename);
 
    elseif strcmp(campaign,'SOCRATES')
        %Aircraft_data_folder = [particleopticspath,'/Phips/SOCRATES/Flight data/'];
        Aircraft_data_filename = ['SOCRATES',lower(flight),'.nc'];
        [time,airspeed,~,~,~,~,~,~] = Aircraft_data(campaign,Aircraft_data_folder,Aircraft_data_filename);
 
    elseif strcmp(campaign,'CIRRUS-HL')
        %Aircraft_data_folder = [particleopticspath,filesep,'Phips',filesep,'CIRRUS-HL',filesep,'Flight data',filesep];
        flightnum = strsplit(flight,'RF');
        flightnum = flightnum{2};
        listings = dir([Aircraft_data_folder, '*CIRRUSHL_F', flightnum, '_*_BAHAMAS_v1.nc']);
        Aircraft_data_filename = listings(end).name;

        % Read TIME and TAS variable from the netCDF file
        time_seconds = ncread([Aircraft_data_folder,Aircraft_data_filename], 'TIME');
        time = datenum(datetime(time_seconds, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC')); % Convert seconds since 1970-01-01 00:00:00 UTC to datetime
        airspeed = ncread([Aircraft_data_folder,Aircraft_data_filename], 'TAS');
    else
        disp('Configure flight data path or data name first!')
        return
    end
    
end

%% Get image size
function [PhipsData] = get_image_size(PhipsData)
    vbs = PhipsData.Properties.VariableNames;
    Index_size = find(contains(vbs,'diameter'));
    Index_area = find(contains(vbs,'proj_area'));
    psize_c1c2 = table2array(PhipsData(:,Index_size));
    parea_c1c2 = table2array(PhipsData(:,Index_area));

    hFig = figure;
    scatter(parea_c1c2(:,1),parea_c1c2(:,2),50,'k'); hold on;
    set(gca,'fontsize',14,'linewidth',1.5); xlabel('Image diameter C1'); ylabel('Image diameter C2')

    % Delete data with too large difference between C1 and C2 area
    % For example dust particles
    parea_c1c2_new = parea_c1c2;
    erotus =  abs(parea_c1c2_new(:,1) - parea_c1c2_new(:,2)) ./ max(parea_c1c2_new, [], 2, "omitnan"); 
    idx = erotus>0.8; % Tolerance
    idx_C1 = parea_c1c2_new(:,1)<parea_c1c2_new(:,2); % find where C1 is small
    idx_C2 = parea_c1c2_new(:,1)>parea_c1c2_new(:,2); % find where C2 is small
    parea_c1c2_new(idx & idx_C1,1) = NaN;
    parea_c1c2_new(idx & idx_C2,2) = NaN;
    psize_c1c2(idx & idx_C1,1) = NaN;
    psize_c1c2(idx & idx_C2,2) = NaN;
    scatter(parea_c1c2_new(:,1),parea_c1c2_new(:,2),50,'r')

    % Make averages
    A_image = mean(parea_c1c2_new,2,"omitnan");
    dp_image = mean(psize_c1c2,2,"omitnan");

    % Fill to PhipsData
    PhipsData = addvars(PhipsData,dp_image,'NewVariableNames','dp_image');
    PhipsData = addvars(PhipsData,A_image,'NewVariableNames','A_image');
        
end

%% Generate histogram
function [particles_per_bin, dN, dA, dm] = generate_histogram(PhipsData, PhipsData_level_0, taxis,tstep,bin_endpoints,A_sens_parameters ,airspeed,time) 
% Function for Generating histograms for PHIPS SD
% A_sens_parameters = A_sens_parameters_ice;
% Mass-dimensional relationship from Baker & Lawson, 2006 is used
% m = 0.115*A^1.218

V_bin = [];
iwc = 0.*ones(length(taxis)-1,1);
for i=1:length(taxis)-1
    t1 = taxis(i)-datenum(0,0,0,0,0,1/1000); %1ms tolerance, because sometimes the timestamps are off by 1/100 ms
    t2 = taxis(i+1)+datenum(0,0,0,0,0,1/1000);
    
    %% calculate SD
    % Select time
    % deadtime calculation (including ForceTriggers)
    idx = PhipsData_level_0.RealTimeStamp>=t1 & PhipsData_level_0.RealTimeStamp<t2;
    number_of_triggers(i) = sum(idx); %histcounts(data,binlimits) sorts data into the number of bins determined by binlimits (intensity)
    qdead = (1-sum(number_of_triggers(i),2)./tstep*12*10^-6);% Fraction of missed volumne due to el. dead time of 12*10^-6s
    
    % size of particles in this time
    idx = PhipsData.RealTimeStamp>=t1 & PhipsData.RealTimeStamp<t2;
    dp = PhipsData.dp(idx); % size of each particle in this time
    A = PhipsData.A(idx); % area of each particle in this time, um2
    m = 0.115.*(A.*1e-6).^1.218; % mass of each particle in this time, mg
    
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
    
    V_per_particle = sensitive_area .* TAS .* 100 .* tstep .* qdead ./ corr_flight_speed; % The sampling Volume [cm-3]
    
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
        area_per_bin(i,j) = nansum(A(idx)); % area of particles per bin, um2
        mass_per_bin(i,j) = nansum(m(idx)); % mass of particles per bin, mg
        % old: conc = [sum ( N ) ] / V , where V was const for one bin
        % new: conc = sum ( N / V ) = 1 / sum(V), where N = 1 (since we
        % count every particle individually and V is the Volume of each particle)
        conc(i,j) = sum(1./V_per_particle(idx)); % recipreocely sumated volume per particle

        % Calculate area (mean of individual areas)
        V = mean(sensitive_area(idx),'omitnan') .* TAS .* 100 .* tstep .* qdead ./ corr_flight_speed; % The sampling Volume [cm-3]

        % Calculate number concentration
        dN(i,j) = particles_per_bin(i,j) ./ V; % cm-3

        % Calculate area concentration
        dA(i,j) = area_per_bin(i,j) ./ V; % um2 / cm3

        % Calculate mass concentration
        dm(i,j) = mass_per_bin(i,j) ./ V; % mg / cm3
        
        % if no particles, then no concentration
        if particles_per_bin(i,j) == 0
            conc(i,j) = 0;
            dN(i,j) = 0;
            dA(i,j) = 0;
            dm(i,j) = 0;
        end
    end % end loop size bins
    
end % end loop time
    
end % end function


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
