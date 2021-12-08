function BG_offset = BG_offset_calib(campaign)
% load BG file for SD calibration offset
% mean over all flights of integrated BG intensity
% runtime ~.2sec
%%

if ismac
    particleopticspath = '/Users/emma/ParticleOptics';
%elseif isunix
    % Code to run on Linux platform
elseif ispc
    particleopticspath = 'P:';
%    particleopticspath = 'C:\Users\wa9929\Desktop';
else
    disp('Platform not supported')
end

if strcmp(campaign, 'SOCRATES')
    flight_num = [1:15];
elseif strcmp(campaign, 'ACLOUD')
    flight_num = [527 530 602 604 605 608 613 614 616 617 618 620 623 625 626];
    
elseif strcmp(campaign, 'ARISTO2017')
    
elseif strcmp(campaign, 'CIRRUS-HL')
    flight_num = 2:17;
else
    disp('UNKNOWN CAMPAIGN')
    
end


BG_mean = []; BG_3sigma = [];

angle = 18:8:170;

for FN = 1:size(flight_num,2)

    if strcmp(campaign, 'SOCRATES') || strcmp(campaign, 'CIRRUS-HL')
        if flight_num(FN) < 10
            flight = ['RF0', num2str(flight_num(FN))];
        else
            flight = ['RF', num2str(flight_num(FN))];
        end
    elseif strcmp(campaign, 'ACLOUD')
        flight = ['Flight 170', num2str(flight_num(FN))];
        if flight_num(FN) == 626
            flight = 'Flight 170626a'; %because we have a problem with the 'a'
        end
    end
    
 %   folder = [particleopticspath,filesep,'Phips ',filesep,campaign,filesep,flight];
    folder = [particleopticspath, '\Phips Results\Campaigns\', campaign, '\', flight, filesep, 'Level Files', filesep];


    cd(folder)
    listings = dir('*BG.txt'); 
    if isempty(listings)
        disp('Level File for this flight does not exist!')
        BG_mean(FN) = NaN;
    else
    csv_name = listings.name;
    filename = [folder,filesep,csv_name];
    BG = dlmread(filename,',');

    BG_mean(FN) =  trapz(angle, BG(1,:));
%    BG_3sigma(FN) = BG(2,:);
    end

end

BG_offset = nanmean(BG_mean);


