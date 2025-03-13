%% translate .sum SD file into .nc
clearvars
addpath('/Users/emma/Documents/GitHub/PHIPS-Scattering-Data-Analysis-Tools/Aircraft Data Functions')
addpath('/Users/emma/Documents/GitHub/PHIPS-Scattering-Data-Analysis-Tools/PHIPS Main Functions')
addpath('/Users/emma/Documents/GitHub/PHIPS-PSD-Tools/Functions/')

particleopticspath = '/Users/emma/ParticleOptics/';
campaign = 'CIRRUS-HL';
tstep = 1; % in seconds

flightno = [2:23]; % CIRRUS-HL
%flightno = [527,602,604,605,608,614,616,617,618]; % ACLOUD




%% loop over all flights
for loop =  1:length(flightno)
    % Make flight name
    flight = make_flight_name(campaign,flightno(loop));

    % Folder where SD files are
    sourcepath = ['/Users/emma/Documents/Field campaigns/CIRRUS-HL/PHIPS Results',filesep, flight, filesep,'SD',filesep];
    savepath = ['/Users/emma/Documents/Field Campaigns/',campaign,'/PHIPS/PHIPS nc files/',num2str(tstep),'s',filesep];
    
    % Folder where aircraft data is
    Aircraft_data_folder = '/Users/emma/Documents/Field campaigns/CIRRUS-HL/BAHAMAS/';

    translate_SD_sum_to_nc_function(sourcepath, campaign, flight, tstep, savepath, Aircraft_data_folder);
end






%% Function

function flight = make_flight_name(campaign,flight_num)
    if strcmp(campaign, 'SOCRATES') || strcmp(campaign, 'CIRRUS-HL')
        if flight_num < 10
            flight = ['RF0', num2str(flight_num)];
        else
            flight = ['RF', num2str(flight_num)];
        end
    elseif strcmp(campaign, 'ACLOUD')
        flight = ['Flight 170', num2str(flight_num)];
    elseif strcmp(campaign, 'IMPACTS')
        flight = ['SF0', num2str(flight_num)];
    end
    disp(flight)
end

function translate_SD_sum_to_nc_function(sourcepath, campaign, flight, tstep, savepath, Aircraft_data_folder)
% #######
% careful, this does not work if the save directory already has an nc file
% with the same name
% #####
  
   if ~isfolder(savepath)
       mkdir(savepath)
   end
   
   %% Read PHIPS data
   % PSDs based on particle images (when not available to scattering data)
    [time_PHIPS, QualityFlag, SD_PHIPS_ice, SD_PHIPS_ice_only_images, SD_PHIPS_drop, ...
    N_ice, N_ice_only_images, N_drop, ...
    N_ice_uncertainty, N_ice_uncertainty_only_images, N_drop_uncertainty, ...
    counts_ice, counts_ice_only_images, counts_drop, ...
    bin_endpoints_PHIPS, bin_midpoints_PHIPS,bin_midpoints_area_PHIPS, ...
    SD_PHIPS_ice_area, SD_PHIPS_ice_area_only_images, A_ice] = ...
    Read_PHIPS_SD(sourcepath, tstep);
    
    % IWC, Reff and Bext
    IWC = load([sourcepath,'IWC_',num2str(tstep),'s.txt']); IWC = IWC(:,2);
    Bext = load([sourcepath,'Bext_',num2str(tstep),'s.txt']); Bext = Bext(:,2);
    Reff = load([sourcepath,'Reff_',num2str(tstep),'s.txt']); Reff = Reff(:,2);
    Reff_only_images = load([sourcepath,'Reff_only_images_',num2str(tstep),'s.txt']); Reff_only_images = Reff_only_images(:,2);
    Rm = load([sourcepath,'Reff_',num2str(tstep),'s.txt']); Rm = Rm(:,2);
    Rm_only_images = load([sourcepath,'Rm_only_images_',num2str(tstep),'s.txt']); Rm_only_images = Rm_only_images(:,2);
    
    % convert datenum to seconds after midnight
    t_vec = datevec(time_PHIPS(:));
    day = t_vec(1,3);
    % seconds + minutes + hours + delta(day)
    time_PHIPS_sec_af_midn = t_vec(:,6) + t_vec(:,5) *60 + t_vec(:,4)*60*60 + (t_vec(:,3)-day)*60*60*24; 


    %% Read aircraft data
    if strcmp(campaign,'ACLOUD')
        %Aircraft_data_folder = [particleopticspath, '/SID3/ACLOUD/Flight data/'];
        flight_date = strsplit(flight, ' ');
        flight_date = flight_date{2};
        Aircraft_data_filename = ['AircraftData_', flight_date, '.dat'];
    elseif strcmp(campaign,'SOCRATES')
        %Aircraft_data_folder = [particleopticspath,'/Phips/SOCRATES/Flight data/'];
        Aircraft_data_filename = ['SOCRATES',lower(flight),'.nc'];
    elseif strcmp(campaign,'CIRRUS-HL')
        %Aircraft_data_folder = [particleopticspath,filesep,'Phips',filesep,'CIRRUS-HL',filesep,'Flight data',filesep];
        flightnum = strsplit(flight,'RF');
        flightnum = flightnum{2};
        listings = dir([Aircraft_data_folder, '*CIRRUSHL_F', flightnum, '_*_BAHAMAS_v1.nc']);
        Aircraft_data_filename = listings(end).name;
    end
        
    [time,airspeed,height,lat,lon,~,~,T,S_ice,w,RH_liq,P] = Aircraft_data(campaign,Aircraft_data_folder,Aircraft_data_filename); close all;
    height = interp1(time, height, time_PHIPS); % interpolate to PHIPS time
    T = interp1(time, T, time_PHIPS); % interpolate to PHIPS time
    P = interp1(time, P, time_PHIPS);
    lat = interp1(time, lat, time_PHIPS);
    lon = interp1(time, lon, time_PHIPS);


    %% Manually remove periods with shattering
    % Manually correct for a timeshift in ACLOUD flight 170616
    if strcmp(flight,'Flight 170616')
        time_PHIPS_sec_af_midn = time_PHIPS_sec_af_midn - 3365;
        time_PHIPS = time_PHIPS - datenum(0,0,0,0,0,3365);
    end
    
    [QualityFlag,quality_string] = manually_remove_shattering(time_PHIPS,QualityFlag,campaign,flight);

    
    %% 1. create file with variable
    filename = [savepath, 'PHIPS_SD_', campaign,'_', flight, '_', num2str(tstep), 's', '.nc'];
    
    % Check if file exists and delete it
    if exist(filename, 'file') == 2
        delete(filename);
    end

    %nccreate(filename,'time_PHIPS', 'Dimensions', {'time',size(time_PHIPS,1),'row',size(time_PHIPS,2)}, 'FillValue','disable');
    nccreate(filename,'time_sec_af_midn', 'Dimensions', {'time',size(time_PHIPS_sec_af_midn,1),'row',size(time_PHIPS_sec_af_midn,2)}, 'FillValue','disable');
    nccreate(filename,'T', 'Dimensions', {'time',size(T,1),'y1',size(T,2)}, 'FillValue','disable');
    nccreate(filename,'P', 'Dimensions', {'time',size(T,1),'y1',size(T,2)}, 'FillValue','disable');
    nccreate(filename,'altitude', 'Dimensions', {'time',size(height,1),'y1',size(height,2)}, 'FillValue','disable');
    nccreate(filename,'lat', 'Dimensions', {'time',size(T,1),'y1',size(T,2)}, 'FillValue','disable');
    nccreate(filename,'lon', 'Dimensions', {'time',size(T,1),'y1',size(T,2)}, 'FillValue','disable');
    nccreate(filename,'QualityFlag', 'Dimensions', {'time',size(QualityFlag,1),'row',size(QualityFlag,2)}, 'FillValue','disable');
    
    nccreate(filename,'PSD_ice', 'Dimensions', {'time',size(SD_PHIPS_ice,1),'size_bins',size(SD_PHIPS_ice,2)}, 'FillValue','disable');
    nccreate(filename,'PSD_ice_image', 'Dimensions', {'time',size(SD_PHIPS_ice_only_images,1),'size_bins',size(SD_PHIPS_ice_only_images,2)}, 'FillValue','disable');
    nccreate(filename,'PSD_drop', 'Dimensions', {'time',size(SD_PHIPS_drop,1),'size_bins',size(SD_PHIPS_drop,2)}, 'FillValue','disable');
    nccreate(filename,'PAD_ice', 'Dimensions', {'time',size(SD_PHIPS_ice,1),'size_bins',size(SD_PHIPS_ice,2)}, 'FillValue','disable');
    nccreate(filename,'PAD_ice_image', 'Dimensions', {'time',size(SD_PHIPS_ice_only_images,1),'size_bins',size(SD_PHIPS_ice_only_images,2)}, 'FillValue','disable');
    nccreate(filename,'PAD_drop', 'Dimensions', {'time',size(SD_PHIPS_drop,1),'size_bins',size(SD_PHIPS_drop,2)}, 'FillValue','disable');
    
    nccreate(filename,'N_ice', 'Dimensions', {'time',size(N_ice,1),'row',size(N_ice,2)}, 'FillValue','disable');
    nccreate(filename,'N_drop', 'Dimensions', {'time',size(N_drop,1),'row',size(N_drop,2)}, 'FillValue','disable');
    nccreate(filename,'A_ice', 'Dimensions', {'time',size(N_ice,1),'row',size(N_ice,2)}, 'FillValue','disable');
    nccreate(filename,'A_drop', 'Dimensions', {'time',size(N_drop,1),'row',size(N_drop,2)}, 'FillValue','disable');
    nccreate(filename,'counts_ice', 'Dimensions', {'time',size(counts_ice,1),'size_bins',size(counts_ice,2)}, 'FillValue','disable');
    nccreate(filename,'counts_drop', 'Dimensions', {'time',size(counts_drop,1),'size_bins',size(counts_drop,2)}, 'FillValue','disable');
    
    nccreate(filename,'IWC', 'Dimensions', {'time',size(IWC,1),'y1',size(IWC,2)}, 'FillValue','disable');
    nccreate(filename,'Reff', 'Dimensions', {'time',size(Reff,1),'y1',size(Reff,2)}, 'FillValue','disable');
    nccreate(filename,'Reff_image', 'Dimensions', {'time',size(Reff,1),'y1',size(Reff_only_images,2)}, 'FillValue','disable');
    nccreate(filename,'b_ext', 'Dimensions', {'time',size(Bext,1),'y1',size(Bext,2)}, 'FillValue','disable');
    nccreate(filename,'Rm', 'Dimensions', {'time',size(Rm,1),'y1',size(Rm,2)}, 'FillValue','disable');
    nccreate(filename,'Rm_image', 'Dimensions', {'time',size(Rm,1),'y1',size(Rm,2)}, 'FillValue','disable');
    
    nccreate(filename,'bin_endpoints', 'Dimensions', {'row',size(bin_endpoints_PHIPS,1),'size_bin_limits',size(bin_endpoints_PHIPS,2)}, 'FillValue','disable');
    nccreate(filename,'bin_midpoints', 'Dimensions', {'row',size(bin_midpoints_PHIPS,1),'size_bins',size(bin_midpoints_PHIPS,2)}, 'FillValue','disable');
    

    %% 2. write in file
    
    % % 1. write Aircraft data
    ncwrite(filename,'time_sec_af_midn',time_PHIPS_sec_af_midn)
        ncwriteatt(filename,'time_sec_af_midn','unit','second after midnight');
        ncwriteatt(filename,'time_sec_af_midn','long name','Second After Midnight');
    ncwrite(filename,'T',T)
        ncwriteatt(filename,'T','unit','degree Celcius');
        ncwriteatt(filename,'T','long name','Temperature');
    ncwrite(filename,'P',P)
        ncwriteatt(filename,'P','unit','hPa');
        ncwriteatt(filename,'P','long name','Static Pressure');
    ncwrite(filename,'altitude',height)
        ncwriteatt(filename,'altitude','unit','meter');
        ncwriteatt(filename,'altitude','long name','Flight Altitude');
    ncwrite(filename,'lat',lat)
        ncwriteatt(filename,'lat','unit','degree');
        ncwriteatt(filename,'lat','long name','Latitude');
    ncwrite(filename,'lon',lon)
        ncwriteatt(filename,'lat','unit','degree');
        ncwriteatt(filename,'lat','long name','Longitude');

    % % 1. write PSD data
    ncwrite(filename,'QualityFlag',QualityFlag)
        ncwriteatt(filename,'QualityFlag','long name','Indicates shattering influenced periods (QualityFlag = 0).');
    ncwrite(filename,'PSD_ice',SD_PHIPS_ice)
        ncwriteatt(filename,'PSD_ice','unit','dN/dlogDp in #/L');
        ncwriteatt(filename,'PSD_ice','long name','Particle size distribution – Ice');
    ncwrite(filename,'PSD_ice_image',SD_PHIPS_ice_only_images)
        ncwriteatt(filename,'PSD_ice_image','unit','dN/dlogDp in #/L');
        ncwriteatt(filename,'PSD_ice_image','long name','Particle size distribution – Ice. PSD based on imaged particles only that are manually classified. Shattered particles and multiple particles in one frame are excluded.');
    ncwrite(filename,'PSD_drop',SD_PHIPS_drop)
        ncwriteatt(filename,'PSD_drop','unit','dN/dlogDp in #/L');
        ncwriteatt(filename,'PSD_drop','long name','Particle size distribution – Droplet');
    ncwrite(filename,'PAD_ice',SD_PHIPS_ice_area)
        ncwriteatt(filename,'PAD_ice','unit','dA/dlogDp in um^2/L');
        ncwriteatt(filename,'PAD_ice','long name','Particle area distribution – Ice');
    ncwrite(filename,'PAD_ice_image',SD_PHIPS_ice_area_only_images)
        ncwriteatt(filename,'PAD_ice','unit','dA/dlogDp in um^2/L');
        ncwriteatt(filename,'PAD_ice','long name','Particle area distribution – Ice. PAD based on imaged particles only that are manually classified. Shattered particles and multiple particles in one frame are excluded.');
    ncwrite(filename,'N_ice',N_ice)
        ncwriteatt(filename,'N_ice','unit','particles per liter');
        ncwriteatt(filename,'N_ice','long name','Total concentration – Ice');
    ncwrite(filename,'N_drop',N_drop)
        ncwriteatt(filename,'N_drop','unit','particles per liter');
        ncwriteatt(filename,'N_drop','long name','Total concentration – Droplet');
    ncwrite(filename,'A_ice',A_ice)
        ncwriteatt(filename,'A_ice','unit','area per liter');
        ncwriteatt(filename,'A_ice','long name','Total area – Ice');
    ncwrite(filename,'counts_ice',counts_ice)
        ncwriteatt(filename,'counts_ice','unit','integer number');
        ncwriteatt(filename,'counts_ice','long name','Particle counts per bin - Ice');
    ncwrite(filename,'counts_drop',counts_drop)
        ncwriteatt(filename,'counts_drop','unit','integer number');
        ncwriteatt(filename,'counts_drop','long name','Particle counts per bin - Droplet');
    ncwrite(filename,'IWC',IWC)
        ncwriteatt(filename,'IWC','unit','g m-3');
        ncwriteatt(filename,'IWC','long name','Ice water content using mass-dimensional relationship from Baker & Lawson, 2006. m = 0.115*A^1.218 where A is projected area.');
    ncwrite(filename,'Reff',Reff)
        ncwriteatt(filename,'Reff','unit','micron');
        ncwriteatt(filename,'Reff','long name','Effective radius');
    ncwrite(filename,'Reff_image',Reff_only_images)
        ncwriteatt(filename,'Reff_image','unit','micron');
        ncwriteatt(filename,'Reff_image','long name','Effective radius based on imaged particles only that are manually classified. Shattered particles and multiple particles in one frame are excluded.');
    ncwrite(filename,'Rm',Rm)
        ncwriteatt(filename,'Rm','unit','micron');
        ncwriteatt(filename,'Rm','long name','Mass-equivalent spherical radius. For definition see Baran et al., 2005.');
    ncwrite(filename,'Rm_image',Rm_only_images)
        ncwriteatt(filename,'Rm_image','unit','micron');
        ncwriteatt(filename,'Rm_image','long name','Mass-equivalent spherical radius based on imaged particles only that are manually classified. Shattered particles and multiple particles in one frame are excluded.');
    ncwrite(filename,'b_ext',Bext)
        ncwriteatt(filename,'counts_drop','unit','km-1');
        ncwriteatt(filename,'counts_drop','long name','Extinction coefficient');
    ncwrite(filename,'bin_endpoints',bin_endpoints_PHIPS)
        ncwriteatt(filename,'bin_endpoints','unit','microns');
        ncwriteatt(filename,'bin_endpoints','long name','Bin End-point Diameters');
    ncwrite(filename,'bin_midpoints',bin_midpoints_PHIPS)
        ncwriteatt(filename,'bin_midpoints','unit','microns');
        ncwriteatt(filename,'bin_midpoints','long name','Bin Mid-point Diameters');

    % % 3. write attributes
    fileattrib(filename,"+w");
    ncwriteatt(filename,"/","note",'Please contact E. Jaervinen before publication, jaervinen@uni-wuppertal.de');
    ncwriteatt(filename,"/","creation_date",datestr(now));
    ncwriteatt(filename,"/","time_resolution",[num2str(tstep),'s']);
    ncwriteatt(filename,"/","reference",'https://doi.org/10.5194/amt-14-3049-2021');
    ncwriteatt(filename,"/","quality_control",quality_string);
    
    ncdisp(filename)
end 