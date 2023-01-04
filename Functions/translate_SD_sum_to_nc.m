%% translate .sum SD file into .nc
clearvars
campaign = 'ACLOUD';
flight_num = 617;
tstep = 10; % in seconds

% Make flight name
flight = make_flight_name(campaign,flight_num)

% Folder where SD files are
sourcepath = ['/Users/emma/ParticleOptics/PHIPS Results/Campaigns/',campaign,filesep', flight, filesep,'SD',filesep];
savepath = ['/Users/emma/Documents/Field Campaigns/ACLOUD/PHIPS/PHIPS nc files/',num2str(tstep),'s',filesep];

translate_SD_sum_to_nc_function(sourcepath, campaign, flight, tstep, savepath);






%% Function

function flight = make_flight_name(campaign,flight_num)
    if strcmp(campaign, 'SOCRATES')
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

function translate_SD_sum_to_nc_function(sourcepath, campaign, flight, tstep, savepath)
% #######
% careful, this does not work if the save directory already has an nc file
% with the same name
% #####
  
   if ~isfolder(savepath)
       mkdir(savepath)
   end
   
   %% Read PHIPS data
    [time_PHIPS, QualityFlag, SD_PHIPS_ice, SD_PHIPS_drop, N_ice, N_drop, ... 
        N_ice_uncertainty, N_drop_uncertainty, counts_ice, counts_drop, ...
        bin_endpoints_PHIPS, bin_midpoints_PHIPS] = Read_PHIPS_SD(sourcepath, tstep);
    
    % convert datenum to seconds after midnight
    t_vec = datevec(time_PHIPS(:));
    day = t_vec(1,3);
    % seconds + minutes + hours + delta(day)
    time_PHIPS_sec_af_midn = t_vec(:,6) + t_vec(:,5) *60 + t_vec(:,4)*60*60 + (t_vec(:,3)-day)*60*60*24; 

    %% Manually remove periods with shattering
    % Manually correct for a timeshift in ACLOUD flight 170616
    if strcmp(flight,'Flight 170616')
        time_PHIPS_sec_af_midn = time_PHIPS_sec_af_midn - 3365;
        time_PHIPS = time_PHIPS - datenum(0,0,0,0,0,3365);
    end
    
    [QualityFlag,quality_string] = manually_remove_shattering(time_PHIPS,QualityFlag,campaign,flight);

    
    %% 1. create file with variable

    % careful! you cannot create a nc file with the same name that already exists
        
    filename = [savepath, 'PHIPS_SD_', flight, '_', num2str(tstep), 's', '.nc'];
    nccreate(filename,'time_PHIPS', 'Dimensions', {'time',size(time_PHIPS,1),'row',size(time_PHIPS,2)}, 'FillValue','disable');
    nccreate(filename,'time_PHIPS_sec_af_midn', 'Dimensions', {'time',size(time_PHIPS_sec_af_midn,1),'row',size(time_PHIPS_sec_af_midn,2)}, 'FillValue','disable');
    nccreate(filename,'QualityFlag', 'Dimensions', {'time',size(QualityFlag,1),'row',size(QualityFlag,2)}, 'FillValue','disable');
    nccreate(filename,'SD_PHIPS_ice', 'Dimensions', {'time',size(SD_PHIPS_ice,1),'size_bins',size(SD_PHIPS_ice,2)}, 'FillValue','disable');
    nccreate(filename,'SD_PHIPS_drop', 'Dimensions', {'time',size(SD_PHIPS_drop,1),'size_bins',size(SD_PHIPS_drop,2)}, 'FillValue','disable');
    
    nccreate(filename,'N_ice', 'Dimensions', {'time',size(N_ice,1),'row',size(N_ice,2)}, 'FillValue','disable');
    nccreate(filename,'N_drop', 'Dimensions', {'time',size(N_drop,1),'row',size(N_drop,2)}, 'FillValue','disable');
    nccreate(filename,'counts_ice', 'Dimensions', {'time',size(counts_ice,1),'size_bins',size(counts_ice,2)}, 'FillValue','disable');
    nccreate(filename,'counts_drop', 'Dimensions', {'time',size(counts_drop,1),'size_bins',size(counts_drop,2)}, 'FillValue','disable');
    
    nccreate(filename,'bin_endpoints_PHIPS', 'Dimensions', {'row',size(bin_endpoints_PHIPS,1),'size_bin_limits',size(bin_endpoints_PHIPS,2)}, 'FillValue','disable');
    nccreate(filename,'bin_midpoints_PHIPS', 'Dimensions', {'row',size(bin_midpoints_PHIPS,1),'size_bins',size(bin_midpoints_PHIPS,2)}, 'FillValue','disable');
    
    % % 2. write in file
    ncwrite(filename,'time_PHIPS',time_PHIPS)
    ncwrite(filename,'time_PHIPS_sec_af_midn',time_PHIPS_sec_af_midn)
    ncwrite(filename,'QualityFlag',QualityFlag)
    ncwrite(filename,'SD_PHIPS_ice',SD_PHIPS_ice)
    ncwrite(filename,'SD_PHIPS_drop',SD_PHIPS_drop)
    ncwrite(filename,'N_ice',N_ice)
    ncwrite(filename,'N_drop',N_drop)
    ncwrite(filename,'counts_ice',counts_ice)
    ncwrite(filename,'counts_drop',counts_drop)
    ncwrite(filename,'bin_endpoints_PHIPS',bin_endpoints_PHIPS)
    ncwrite(filename,'bin_midpoints_PHIPS',bin_midpoints_PHIPS)

    % % 3. write attributes
    fileattrib(filename,"+w");
    ncwriteatt(filename,"/","note",'Please contact E. Jaervinen before publication, emma.jaervinen@kit.edu');
    ncwriteatt(filename,"/","creation_date",datestr(now));
    ncwriteatt(filename,"/","time_resolution",[num2str(tstep),'s']);
    ncwriteatt(filename,"/","reference",'https://doi.org/10.5194/amt-14-3049-2021');
    ncwriteatt(filename,"/","quality_control",quality_string);
    
    ncdisp(filename)
end 