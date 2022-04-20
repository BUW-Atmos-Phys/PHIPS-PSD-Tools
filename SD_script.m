

%%

clear all
close all
%%
if ismac
    particleopticspath = '/Users/emma/ParticleOptics';
%elseif isunix
    % Code to run on Linux platform
elseif ispc
    particleopticspath = 'P:';
    particleopticspath = '\\os.lsdf.kit.edu\kit\imk-aaf\projects\ParticleOptics\';
%     SD_2DS_path = ['C:\Users\wa9929\Desktop', filesep, '2DS SD', filesep];
    SD_2DS_path = [particleopticspath, '\Other Probes\SOCRATES\2DS\'];
    SD_2DC_path = [particleopticspath, '\Other Probes\SOCRATES\2DC\'];
else
    disp('Platform not supported')
end

campaign = 'SOCRATES';
flightno = 1:15;
flight = 'RF11'; % for debugging only
%
%campaign = 'ACLOUD';
%flightno =  [527, 530, 602, 604, 605, 608, 613, 614, 616, 617, 618, 620, 623, 6261, 6262];


% flight = 'Flight 170602';
% campaign = 'CIRRUS-HL';
% flightno = 2:17;

save_status = 1;

% Time Resolution
tstep = 1;

loop = 2; %for debugging

%%
addpath('C:\Users\Fritz\Documents\GitHub\PHIPS-PSD-Tools\Functions')

% addpath([particleopticspath,'/Software/PHIPS analysis/Size Distribution/functions'])
% addpath([particleopticspath,'/Software/PHIPS analysis/PHIPS inversion 2\Functions'])
% addpath([particleopticspath,'\Software\SID3_MATLAB\SID3 analysis\functions']) % for read aircraft data
% addpath([particleopticspath,'\Other Probes\Software Read other probes\']) % for read 2DS

addpath('C:\Users\Fritz\Documents\GitHub\PHIPS-Scattering-Data-Analysis-Tools\Aircraft Data Functions\') % aircraft data
addpath('C:\Users\Fritz\Documents\GitHub\PHIPS-Scattering-Data-Analysis-Tools\PHIPS Main Functions') % import phips

%% loop over all flights

for loop = 1:length(flightno)
    if strcmp(campaign, 'SOCRATES') || strcmp(campaign, 'CIRRUS-HL') %name the correct flight
        if flightno(loop) < 10
            flight = ['RF0', num2str(flightno(loop))];
        else
            flight = ['RF', num2str(flightno(loop))];
        end
    elseif strcmp(campaign, 'ACLOUD')
        flight = ['Flight 170', num2str(flightno(loop))];
        if flightno(loop) == 6261
            flight = 'Flight 170626a'; %because we have a problem with the 'a'
        elseif flightno(loop) == 6262
            flight = 'Flight 170626b';
        end
    end

    %%
    


% if the plot should go only over a particular time range
specific_interval_only = 0;
disp(['Calculate SD for flight ', flight])

if specific_interval_only == 1
    start_time  = datenum('18-Feb-2018 03:28:00');
    end_time = datenum('18-Feb-2018 03:38:00');
    tstep = etime(datevec(end_time),datevec(start_time));
    [SDice,SDdroplet] = PHIPS_SD(particleopticspath,SD_2DS_path,SD_2DC_path,campaign,flight,tstep,save_status,start_time,end_time);
else
    [SDice,SDdroplet] = PHIPS_SD(particleopticspath,SD_2DS_path,SD_2DC_path,campaign,flight,tstep,save_status);
end


disp(['Finished with SD for flight ', flight])




%%
end

disp(['Finished with SD for all flights.'])

%% Save Files
% save_SD(SDice, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bins)
% save_SD(SDdroplet, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bins)


%% Loop over all Flights

% campaign = 'ACLOUD';
% flightno = [527 530 602 604 605 608 618];
% % flightno = [527 530 602 604 605 608 613 614 616 617 618 620 623 625 626];

% flight_num = [4:1:15];
% 
% for FN = 1:size(flight_num,2)
% 
%     if strcmp(campaign, 'SOCRATES')
%         if flight_num(FN) < 10
%             flight = ['RF0', num2str(flight_num(FN))];
%         else
%             flight = ['RF', num2str(flight_num(FN))];
%         end
%     elseif strcmp(campaign, 'ACLOUD')
%         flight = ['Flight 170', num2str(flight_num(FN))];
%         if flightno(loop) == 626
%             flight = 'Flight 170626a'; %because we have a problem with the 'a'
%         end
%     end
% 
%     savepath = [particleopticspath, filesep, 'PHIPS Results',filesep, campaign, filesep, flight, filesep, 'SD']; 
% %    savepath_quicklooks = [particleopticspath, filesep, 'PHIPS Results',filesep, campaign, filesep, flight, filesep, 'Quicklooks']; 
% 
%     tstep = 1;
%     [SDice,SDdroplet] = PHIPS_SD(POpath,campaign,flight,tstep) %,start_time,end_time)
%     save_SD(SDice, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bins)
%     save_SD(SDdroplet, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bins)
% 
%     tstep = 10;
%     [SDice,SDdroplet] = PHIPS_SD(POpath,campaign,flight,tstep) %,start_time,end_time)
%     save_SD(SDice, save_status, savepath, campaign, flight, tstep, a_ice, b_ice, bins)
%     save_SD(SDdroplet, save_status, savepath, campaign, flight, tstep, a_droplet, b_droplet, bins)
%     
% end

%%





