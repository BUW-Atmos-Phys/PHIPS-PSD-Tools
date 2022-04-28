


% translate .sum SD file into .nc
% read as e.g. % a = ncread(filename, 'time_PHIPS');
% campaign = 'SOCRATES';
% flight_num = 1:15;
% FN = 2; % for debugging
% 
% mainfolder = 'C:\Users\Fritz\Desktop\SD_Emma\PHIPS Results\Campaigns\SOCRATES\';
% %%
% for FN = 1:size(flight_num,2)
%     
%     %%
%     
%     if strcmp(campaign, 'SOCRATES')
%         if flight_num(FN) < 10
%             flight = ['RF0', num2str(flight_num(FN))];
%         else
%             flight = ['RF', num2str(flight_num(FN))];
%         end
%     elseif strcmp(campaign, 'ACLOUD')
%         flight = ['Flight 170', num2str(flight_num(FN))];
%     elseif strcmp(campaign, 'IMPACTS')
%         flight = ['SF0', num2str(flight_num(FN))];
%     end
%     
%     disp(flight)
%     
%     folder = [mainfolder, flight, '\SD\'];

function translate_SD_sum_to_nc(sourcepath, campaign, flight, tstep);
% #######
% careful, this does not work if the save directory already has an nc file
% with the same name
% #####

  %  savepath = folder; % then it saves to the same directory
  savepath = sourcepath;
  savepath = 'C:\Users\Fritz\Desktop\PHIPS SD test flight speed correction\new\'; % saves to a specific directory  
  
   if ~isdir(savepath)
       mkdir(savepath)
   end
   
        % tstep = 1;
        time_resolution = tstep;
        [time_PHIPS, ShatteringFlag, SD_PHIPS_ice, SD_PHIPS_drop, N_ice, N_drop, N_ice_uncertainty, N_drop_uncertainty, counts_ice, counts_drop, ...
    bin_endpoints_PHIPS, bin_midpoints_PHIPS] = Read_PHIPS_SD(sourcepath, tstep);
        
        %% convert datenum to seconds after midnight

        t_vec = datevec(time_PHIPS(:));
        day = t_vec(1,3);
        % seconds + minutes + hours + delta(day)
        time_PHIPS_sec_af_midn = t_vec(:,6) + t_vec(:,5) *60 + t_vec(:,4)*60*60 + (t_vec(:,3)-day)*60*60*24; 
        
        %% 1. create file with variable

        % careful! you cannot create a nc file with the same name that already exists
        
        % the xn and y15 etc are "dimesion names". they are necessairy, though I
        % have no clue what they are actually used for. dimension names must NOT be
        % different size for different vars
        
        filename = [savepath, 'PHIPS_SD_', flight, '_', num2str(tstep), 's', '.nc'];
        nccreate(filename,'time_resolution', 'Dimensions', {'x1',size(time_resolution,1),'y1',size(time_resolution,2)}, 'FillValue','disable');
        nccreate(filename,'time_PHIPS', 'Dimensions', {'xn',size(time_PHIPS,1),'y1',size(time_PHIPS,2)}, 'FillValue','disable');
        nccreate(filename,'time_PHIPS_sec_af_midn', 'Dimensions', {'xn',size(time_PHIPS_sec_af_midn,1),'y1',size(time_PHIPS_sec_af_midn,2)}, 'FillValue','disable');
        nccreate(filename,'ShatteringFlag', 'Dimensions', {'xn',size(ShatteringFlag,1),'y1',size(ShatteringFlag,2)}, 'FillValue','disable');
        nccreate(filename,'SD_PHIPS_ice', 'Dimensions', {'xn',size(SD_PHIPS_ice,1),'y14',size(SD_PHIPS_ice,2)}, 'FillValue','disable');
        nccreate(filename,'SD_PHIPS_drop', 'Dimensions', {'xn',size(SD_PHIPS_drop,1),'y14',size(SD_PHIPS_drop,2)}, 'FillValue','disable');
        
        nccreate(filename,'N_ice', 'Dimensions', {'xn',size(N_ice,1),'y1',size(N_ice,2)}, 'FillValue','disable');
        nccreate(filename,'N_drop', 'Dimensions', {'xn',size(N_drop,1),'y1',size(N_drop,2)}, 'FillValue','disable');
        nccreate(filename,'counts_ice', 'Dimensions', {'xn',size(counts_ice,1),'y14',size(counts_ice,2)}, 'FillValue','disable');
        nccreate(filename,'counts_drop', 'Dimensions', {'xn',size(counts_drop,1),'y14',size(counts_drop,2)}, 'FillValue','disable');
        
        nccreate(filename,'bin_endpoints_PHIPS', 'Dimensions', {'x1',size(bin_endpoints_PHIPS,1),'y15',size(bin_endpoints_PHIPS,2)}, 'FillValue','disable');
        nccreate(filename,'bin_midpoints_PHIPS', 'Dimensions', {'x1',size(bin_midpoints_PHIPS,1),'y14',size(bin_midpoints_PHIPS,2)}, 'FillValue','disable');
        
        % % 2. write in file
        ncwrite(filename,'time_resolution',time_resolution)
        ncwrite(filename,'time_PHIPS',time_PHIPS)
        ncwrite(filename,'time_PHIPS_sec_af_midn',time_PHIPS_sec_af_midn)
        ncwrite(filename,'ShatteringFlag',ShatteringFlag)
        ncwrite(filename,'SD_PHIPS_ice',SD_PHIPS_ice)
        ncwrite(filename,'SD_PHIPS_drop',SD_PHIPS_drop)
        ncwrite(filename,'N_ice',N_ice)
        ncwrite(filename,'N_drop',N_drop)
        ncwrite(filename,'counts_ice',counts_ice)
        ncwrite(filename,'counts_drop',counts_drop)
        ncwrite(filename,'bin_endpoints_PHIPS',bin_endpoints_PHIPS)
        ncwrite(filename,'bin_midpoints_PHIPS',bin_midpoints_PHIPS)
        % ncwrite(filename,'time_PHIPS',time_PHIPS)
        
        nini = ncinfo(filename);
 %   end % loop tstep
end % loop flighgts
