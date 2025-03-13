% Master Script for generating PSDs for an entire campaign
clear all; close all; clc;
addpath('Functions/')
addpath('/Users/emma/Documents/GitHub/PHIPS-Scattering-Data-Analysis-Tools/Aircraft Data Functions')
addpath('/Users/emma/Documents/GitHub/PHIPS-Scattering-Data-Analysis-Tools/PHIPS Main Functions')

% campaign = 'SOCRATES';
% flightno = 1:15;

% campaign = 'ACLOUD';
% flightno =  [616];

campaign = 'CIRRUS-HL';
flightno = 2:23;

% Time Resolution
tstep = 1; % in seconds
save_status = 1;


%% Generate paths
Aircraft_data_folder = '/Users/emma/Documents/Field campaigns/CIRRUS-HL/BAHAMAS/';
particleopticspath = [];

% Path to 2DS / CIP
SD_2DS_path = [particleopticspath, filesep,'Other Probes',filesep,campaign,'ACLOUD_CDP_CIP_PIP_data',filesep];
SD_2DC_path = [particleopticspath, filesep,'Other Probes',filesep,campaign,'ACLOUD_CDP_CIP_PIP_data',filesep];


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
    disp(['Calculate SD for flight ', flight])
    % Folder for PHIPS level files
    folder = ['/Users/emma/Documents/Field campaigns/',campaign,'/PHIPS Results/',flight,'/Level Files/'];
    % Folder for saving data
    savepath = ['/Users/emma/Documents/Field campaigns/',campaign,'/PHIPS Results/',flight,'/SD/'];

    [SDice,SDdroplet,SDice_v1,SDdroplet_v1,SDice_area,Bext,IWC] = ...
    PHIPS_SD(Aircraft_data_folder, folder, savepath, SD_2DS_path,SD_2DC_path,...
    campaign,flight,tstep,save_status);
    disp(['Finished with SD for flight ', flight])


end
disp('Finished with SD for all flights.')






