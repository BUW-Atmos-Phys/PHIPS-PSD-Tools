%% Size Calc from SPF vs Size image
% partially based on SD_from_SPF_main.m

clear all
% close all
%% Initialize
save_status = 0;
doplot = 1;

campaign = 'SOCRATES';
flight = 'RF08';
%campaign = 'ACLOUD';
%flight = 'Flight 170602';

particleopticspath = 'P:';
%particleopticspath = '/Users/emma/ParticleOptics';
folder = [particleopticspath, filesep, 'Phips', filesep, campaign, filesep, flight];
savepath = [particleopticspath, filesep, 'PHIPS Results',filesep, campaign, filesep, flight, filesep, 'SD']; 
savepath_quicklooks = [particleopticspath, filesep, 'PHIPS Results',filesep, campaign, filesep, flight, filesep, 'Quicklooks']; 
savepath = 'B:\';

numbins = 20; %number of size bins
tstep_sec = 10; %how many seconds = 1 timestep?
tstep = tstep_sec *2.5 / 3600/60; %conversion seconds to datenum

ignore_above_saturation = 1;
ignore_below_background = 0;

angle = [42:8:170];

addpath([particleopticspath,filesep,'Software',filesep,'SID3_MATLAB',filesep,'SID3 analysis',filesep,'functions']) %here is the function for aircraft_data analysis
addpath([particleopticspath,filesep,'Software',filesep,'PHIPS analysis',filesep,'PHIPS inversion 2',filesep,'Functions'])
addpath([particleopticspath,filesep,'Software',filesep,'PHIPS analysis',filesep,'Size Distribution'])

%%

flight_num = [2 4 6 7 8 9 10 11 12 13]; 
dw_ice = [];
dw_ice_image = [];
dw_drop = [];
dw_drop_image = [];
int_ice = [];
int_drop = [];

trig_drop = [];
trig_ice = [];

for fn = 1:size(flight_num,2) %loop over all flights
   if flight_num(fn)  >= 1 && flight_num(fn)< 10
      flight = strcat('RF0', num2str(flight_num(fn)));
   elseif flight_num(fn)  >= 10
       flight = ['RF', num2str(flight_num(fn))];
   end
   
   message = ['Flight ', num2str(flight)];
   disp (message)

  
%% Load PhipsData
disp('Loading PhipsData...')

folder = [particleopticspath, filesep, 'Phips', filesep, campaign, filesep, flight];

% load lvl 1 file
cd(folder)
listings = dir('*level_3.csv'); 
csv_name = listings.name;
filename = [folder,filesep,csv_name];

% Date 
date = strsplit(csv_name,'_');
date = strsplit(date{2},'-');

PhipsData = import_phips(filename); 
disp('...done.')

start_time = PhipsData.RealTimeStamp(1, 1);
end_time = PhipsData.RealTimeStamp(end-1, 1);
%start_time = datenum('17-Feb-2018 03:00:00'); %RF11
%end_time = datenum('17-Feb-2018 05:00:00');
%start_time = datenum('24-Feb-2018 04:00:00'); %RF15
%end_time = datenum('24-Feb-2018 07:00:00');

if flight == 'RF06'
    for i=1:size(PhipsData,1)
       if PhipsData.ImageCamera1(i) < 2228 && PhipsData.ImageCamera1(i) >= 1 %to leave in the NaN
          PhipsData.dw_C1(i) = NaN;
       end
    end
end

%% Remove shattered
thresh_interarrivaltime = 0.5/1000; % = 0.5 ms
[PhipsData] = PHIPS_remove_shattered_new(PhipsData,thresh_interarrivaltime);

%% which particles have channels saturated or in background?

PhipsArray = table2array(PhipsData);
%    PhipsArray = PhipsArray((isnan(PhipsData.DropletFlag_Algorithm)==0),:);
num_p = size(PhipsArray,1); %number of particles
num_c = 20-4; %starting at channel end-X, angle 4 = 42°
num_sat = zeros(num_p,3); %number of satureated chanels of each particle (i.e. 5/20 chanels are sat.)
num_bg = zeros(num_p,3); %particles in BG
num_sat(:,3) = PhipsArray(:,2) ; %= Time
num_bg(:,3) = PhipsArray(:,2) ; 

for i = 1:num_p
    for j = 1:num_c
        if PhipsArray(i,end-num_c+j) > 2030
            num_sat(i,1) = num_sat(i,1) + 1;
        end
        if isnan(PhipsArray(i,end-num_c+j))
            num_bg(i,1) = num_bg(i,1) + 1;
        elseif PhipsArray(i,end-num_c+j) < 5
            num_bg(i,1) = num_bg(i,1) + 1;
        end
    end
    if num_sat(i,1) >= 4 && ignore_above_saturation == 1
        num_sat(i,2) = 1; %marker: more than 3 chanels are saturated
        PhipsData.DropletFlag(i) = NaN; %erase those particles from Drop and Ice, but they still appear in Total
    end
    if num_bg(i,1) >= 4  && ignore_below_background == 1%particles below BG
        num_bg(i,2) = 1;
        PhipsData.DropletFlag(i) = NaN;
    end
end    


%% Calculate SD
for loop=1:3 %loop = for drop and ice
    if loop==1
        type = 'Total';
    elseif loop==2
        type = 'Droplet';
    elseif loop == 3
        type = 'Ice';
    end

%% Integrating
[int_SPF, dw_C1, dw_C2] = integrate_SPF(PhipsData, angle, type);

% calibration: diam = a * I^0.5, 
% ACLOUD: a = 0.82897
% SOCTRATES:  a = 1.3709
if strcmp(campaign, 'ACLOUD')
    a = 0.82897;
elseif strcmp(campaign, 'SOCRATES')
    a = 1.3709;
end



if strcmp(type, 'Ice')
    a = 0.3281;
    b = 0.6340;
    dw_calc = a .* int_SPF.^(b);
    
    dw_ice = [dw_ice; dw_calc];
    dw_ice_image = [dw_ice_image; dw_C1];
    int_ice = [int_ice; int_SPF];
    trig = PhipsData.TriggerIntensity(PhipsData.DropletFlag==0);
    trig_ice = [trig_ice; trig]; %Trigger intensity
%    RTS_Ice = RTS;
%    SD_Ice = SD;
elseif strcmp(type, 'Droplet')
    a = 1.3709;
    b = 0.5;
    dw_calc = a .* int_SPF.^(b);
    
    dw_drop = [dw_drop; dw_calc];
    dw_drop_image = [dw_drop_image; dw_C1];
    int_drop = [int_drop; int_SPF];
    trig = PhipsData.TriggerIntensity(PhipsData.DropletFlag==1);
    trig_drop = [trig_drop; trig]; 
%    RTS_Drop = RTS;
%    SD_Drop = SD;
% elseif strcmp(type, 'Total')
%     dw_Total = dw_calc;
%     dw_Total_image = dw_mess;
%     RTS_Total = RTS;
%     SD_Total = SD;
%     SAT_Total = SAT;
end

end %END OF LOOP over 3 types

    clear drop_intensity
    clear DropletData
    clear PhipsData

end %END of loop over all flights

A_ice = pi /4 .* dw_ice_image;
A_drop = pi /4 .* dw_drop_image;

%% ICE 
type = 'ice';

num_bins_ice = 75;
bins_ice = logspace(1,3,num_bins_ice+1); 

bins_SPF = logspace(2,5.3, num_bins_ice+1);

%Sort Data into bins
array = zeros(num_bins_ice,num_bins_ice);
tic


for num = 1:size(dw_ice,1)
    for i = 1:num_bins_ice %rows = calc
        for j = 1:num_bins_ice %coloumns = images
            %if dw_ice(num) >= bins_ice(i) && dw_ice(num) < bins_ice(i+1)&& dw_ice_image(num) >= bins_ice(j)&& dw_ice_image(num) < bins_ice(j+1)
            if int_ice(num) >= bins_SPF(i) && int_ice(num) < bins_SPF(i+1)&& dw_ice_image(num) >= bins_ice(j)&& dw_ice_image(num) < bins_ice(j+1)
                array(i,j) = array(i,j) + 1;
            end
        end
    end
end
toc

%% Plot heatmap Ice
clear mid_bin
clear mid_SPF
for i=1:num_bins_ice
%    mid_bin(i) = bins_ice(i+1);
    mid_bin(i) = mean(bins_ice(i:i+1));
    mid_SPF(i) = mean(bins_SPF(i:i+1));
end

num_particles_ice = size(dw_ice_image,1) - sum(isnan(dw_ice_image));

% % find max
pos = find(array == max(max(array)));

pos_row = mod(pos,num_bins_ice);
pos_col = (pos-pos_row)/num_bins_ice +1;
pos2 = [pos_row, pos_col];

% %lin fit through max and farfield
% x1 = 500; %should fit X=Y for big sizes
% y1 = 500;
% x2 = mid_bin(pos_col);
% y2 = mid_bin(pos_row);
% m = (y1-y2)/(x1-x2);
% c = y1 - m *x1;
% 
% 
% % alternativ:
% xx = [x1,x2];
% yy = [y1,y2];
% %xx = [22, 26, 34, 400, 150]
% %yy = [32, 45, 59, 400, 160]
% p = polyfit(xx,yy,1)

x = 10:1000; %sizerange

figure(74)
pcolor(mid_bin, mid_SPF,array)
shading interp
hand = gca;
colormap(jet(250));
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel ('Surface equiv. Diameter from C1 Images - [log(D)]')
ylabel ('Calculated Diameter from SPF - [log(D)]')
xticks ([20,50 , 100, 200, 500])
yticks (xticks)
title(['SOCRATES, All Flights, Ice, ', num2str(num_particles_ice), ' particles' ]) 
%colorbar

hold on
plot(x, x, 'color', 'white', 'linewidth', 2)

y = 0.99* x + 10;
%plot(x, y, 'color', 'red', 'linewidth', 2)

plot(x,polyval(p,x),'color', 'white', 'linewidth', 2)

txt = 'Linear Fit through max \rightarrow';
text(25,130,txt, 'color', 'white', 'FontSize',16')
txt = ['y = ' num2str(p(1)), ' * x + ', num2str(p(2))];
text(25,110,txt, 'color', 'white', 'FontSize',12')
txt = '\leftarrow X = Y Line';
text(80,70,txt, 'color', 'white', 'FontSize',16)

 X = [0.15 0.15];
 Y = [0.15   0.3];
annotation('doublearrow',X,Y, 'linewidth', 2, 'color', 'white');
txt = '20\mum';
text(12,20,txt, 'color', 'white', 'FontSize',16)

hold off
grid on

if save_status ==1
    save_filename = [filesep,'D_image_vs_D_SPF_SOCRATES_all_flights_', type, '.png']; 
    print(strcat (savepath, save_filename) ,'-dpng')    
end


%% corrected

type = 'ice_corrected';

num_bins_ice = 75;
bins_ice = logspace(1,3,num_bins_ice+1); 

array_corrected = zeros(num_bins_ice,num_bins_ice);

dw_corr = (dw_ice-c)/m;

tic
for num = 1:size(dw_corr,1)
    for i = 1:num_bins_ice %rows = calc
        for j = 1:num_bins_ice %coloumns = images
            if dw_corr(num) >= bins_ice(i) && dw_corr(num) < bins_ice(i+1)&& dw_ice_image(num) >= bins_ice(j)&& dw_ice_image(num) < bins_ice(j+1)
                array_corrected(i,j) = array_corrected(i,j) + 1;
            end
        end
    end
end
toc

% % Plot heatmap Ice
clear mid_bin
for i=1:num_bins_ice
%    mid_bin(i) = bins_ice(i+1);
    mid_bin(i) = mean(bins_ice(i:i+1));
end

num_particles_ice = size(dw_ice_image,1) - sum(isnan(dw_ice_image));
x = 10:1000; %sizerange

figure(86)
pcolor(mid_bin, mid_bin,array_corrected)
shading interp
hand = gca;
colormap(jet(250));
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel ('Surface equiv. Diameter from C1 Images - [log(D)]')
ylabel ('Corrected Calculated Diameter from SPF (corrected) - [log(D)]')
xticks ([20,50 , 100, 200, 500])
yticks (xticks)
title(['SOCRATES, All Flights, Ice, ', num2str(num_particles_ice), ' particles, Corrected by Lin. Fit' ]) 
%colorbar

hold on
plot(x, x, 'color', 'white', 'linewidth', 2)

hold off
grid on

if save_status ==1
    save_filename = [filesep,'D_image_vs_D_SPF_SOCRATES_all_flights_', type, '.png']; 
    print(strcat (savepath, save_filename) ,'-dpng')    
end

%% same for droplets

type = 'droplet';

num_bins = 100;
bins = logspace(log10(40), log10(600),num_bins+1); % from 10^1 to 10^3 = 10:1000

array_drop = zeros(num_bins,num_bins);

tic
for num = 1:size(dw_drop,1)
    for i = 1:num_bins
        for j = 1:num_bins
            if dw_drop(num) >= bins(i) && dw_drop(num) < bins(i+1)&& dw_drop_image(num) >= bins(j)&& dw_drop_image(num) < bins(j+1)
                array_drop(i,j) = array_drop(i,j) + 1;
            end
        end
    end
end
toc

% Plot heatmap
clear mid_bin
for i=1:num_bins
    mid_bin(i) = bins(i);
    %mid_bin(i) = mean(bins(i:i+1));
end

num_particles_drop = size(dw_drop_image,1) - sum(isnan(dw_drop_image));

figure(85)
pcolor(mid_bin, mid_bin,array_drop)
shading interp
hand = gca;
colormap(jet(250));
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel ('Surface equiv. Diameter from C1 Images  - [log(D)]')
ylabel ('Calculated Diameter from SPF - [log(D)]')
xticks ([20,50 , 100, 200, 500])
yticks (xticks)
title(['SOCRATES, All Flights, Droplets, ', num2str(num_particles_drop), ' particles' ]) 
%colorbar

hold on
plot(10:1000, 10:1000, 'color', 'white', 'linewidth', 2)
hold off
grid on

if save_status ==1
    save_filename = [filesep,'D_image_vs_D_SPF_SOCRATES_all_flights_', type, '.png']; 
    print(strcat (savepath, save_filename) ,'-dpng')    
end

%% Normal plot (not heatmap)

% % Plot Image Diam vs SPF Diam
figure(7)
plot(dw_ice_image,dw_ice, '+')
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
xlim([10,700])
ylim(xlim)
xlabel D-SPF
ylabel D-Image
hold on
idx = find(~isnan(dw_ice_image));
p = polyfit(dw_ice_image(idx),dw_ice(idx),1);
x1 = linspace(0,700);
y1 = polyval(p,x1);
plot(x1,y1, 'linewidth', 2, 'color', 'red')

plot (1:700, 1:700, 'linewidth', 2)
plot (1:700, 0.95* [1:700]+18.9, 'linewidth', 2)

%legend ('RF08 Ice Data', 'X=Y')

hold off
grid on




