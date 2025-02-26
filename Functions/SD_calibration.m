%%
% Load PhipsData from all flights and use dw_C1 size from image data to
% calibrate for size distribution with dw = a * I^0.5 fit
% using only droplets

% ACLOUD: a = 0.82897
% SOCTRATES:  a = 1.3709 ( = 1.2819 if you dont exclude the lower mode of RF06)

clear all; close all; clc;
%%

% campaign = 'ACLOUD';
% flightno = [527 530 602 604 605 608 613 614 616 617 618 620 623 6261 6262];
% campaign = 'SOCRATES';
% flightno = 1:15;
campaign = 'CIRRUS-HL';
flightno = 2;%2:23;
%campaign = 'IMPACTS';
%flightno = [118, 125, 201, 207, 213, 218, 220, 225];

if ismac
    particleopticspath = '/Users/emma/ParticleOptics';
%elseif isunix
    % Code to run on Linux platform
elseif ispc
    particleopticspath = 'P:';
else
    disp('Platform not supported')
end
save_path = '/Users/emma/Documents/Field Campaigns/CIRRUS-HL/PHIPS/PSD Calibration';
if ~isdir(save_path)
    mkdir(save_path)
end

doplot = 1;
save_status = 0;

addpath([particleopticspath, '/Software/PHIPS analysis/PHIPS inversion 2/Functions/'])
addpath([particleopticspath, '/Software/SID3_MATLAB/SID3 analysis/functions'])
addpath([particleopticspath, '/Software/PHIPS analysis/Quicklook Ice Concentration/Functions/']);
addpath([particleopticspath, '/Software/PHIPS analysis/Size Distribution/Functions/'])

angle = [42:8:170];
ang = size(angle,2);

%sizelim = 100;
%tresh = 10;

PhipsData_all = [];

% load all data
for loop=1:length(flightno)
    
    %% Load PhipsData
    
    if strcmp(campaign, 'SOCRATES') || strcmp(campaign, 'CIRRUS-HL') %name the correc flight
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
    elseif strcmp(campaign, 'IMPACTS')
        flight = ['Flight_200', num2str(flightno(loop))];
    end
    
    disp(['Processing ', flight])
    folder = [particleopticspath, filesep,'PHIPS Results',filesep,'Campaigns',filesep, campaign, filesep, flight, filesep, 'Level Files', filesep];
    
    cd(folder)
    listings = dir('*level_3.csv');
    if isempty(listings)
        disp('Level File for this flight does not exist!')
    else
        csv_name = listings(end).name;
        filename = [folder,filesep,csv_name];
        
        PhipsData = import_phips(filename);
        
        if strcmp(flight, 'RF06') && strcmp(campaign, 'SOCRATES') %delete the lower mode of RF06
            for i=1:size(PhipsData,1)
                if PhipsData.ImageCamera1(i) < 2228 && PhipsData.ImageCamera1(i) > 1
                    PhipsData.DropletFlag(i) = NaN;
                end
            end
        end
        
        
        %% Remove shattered
        thresh_interarrivaltime = 0.5/1000; % = 0.5 ms
        [PhipsData] = PHIPS_remove_shattered(PhipsData,thresh_interarrivaltime);
        
        %% Remove saturates
        % if more than 4 channels are saturated/background, remove the particle (DropletFlag = NaN)
        ignore_above_saturation = 1;
        ignore_below_background = 0;
        
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
        disp([num2str(num_BG), ' particles below BG removed'])
        
        
        PhipsData_all = [PhipsData_all; PhipsData];
    end % end of "if listings exists"
    
end %end of loop over all flights
disp('... done importing all flights.')
close all

PhipsData_all_uncorrected = PhipsData_all;
% return

%%

A = PhipsData_all;
idx = find(A.DataSet >0);

%% C2 correction
% correction changes nothing for the actual calibration if the data is binned
% for non-binned, the change is 0.2%
% input_data = PhipsData_all_uncorrected(PhipsData_all_uncorrected.DropletFlag == 1,:);
% particle_type = 'droplet';
% PhipsData_all = SD_calibration_correction_C2(input_data, campaign, particle_type, save_status, save_path);
% 
% input_data = PhipsData_all_uncorrected(PhipsData_all_uncorrected.DropletFlag == 0,:);
% particle_type = 'ice';
% PhipsData_all = SD_calibration_correction_C2(input_data, campaign, particle_type, save_status, save_path);

input_data = PhipsData_all_uncorrected;
particle_type = 'all';
PhipsData_all = SD_calibration_correction_C2(input_data, campaign, particle_type, save_status, save_path);


%% Get image diameter and projected area
PhipsData_all = get_image_size(PhipsData_all);


%% SD Calibration ICE
% % sort in size bins
bins = [30,60,100,150,200,250,300,350,400,500,600,700];
bins_area = pi.*bins.^2./4;
binned_status = 0; % Fit through bins (1) or throug all data (0)
particle_type = 'Ice';
input_data = PhipsData_all(PhipsData_all.DropletFlag == 0,:);

fitresult = do_plot_calibration(input_data, bins, bins_area, campaign, binned_status, particle_type, save_status, save_path);
a_ice = fitresult.a;
ci = confint(fitresult); a_ice_err = fitresult.a - ci(1);


%% SD Calibration Droplet
% % sort in size bins
bins = [30,60,100,150,200,250,300,350,400,500,600,700];
bins_area = pi.*bins.^2./4;
binned_status = 1;
particle_type = 'Droplets';
input_data = PhipsData_all(PhipsData_all.DropletFlag == 1,:);

fitresult = do_plot_calibration(input_data, bins, bins_area, campaign, binned_status, particle_type, save_status, save_path);
a_drop = fitresult.a;
ci = confint(fitresult); a_drop_err = fitresult.a - ci(1);


%% FUNCTIONS


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

%%
function fitresult = do_plot_calibration(input_data, bins, bins_area, campaign, binned_status, particle_type, save_status, save_path);

dp_midpoint = bins(1:end-1)+diff(bins)./2;
dp_midpoint_area = bins_area(1:end-1)+diff(bins_area)./2;

angle = [42:8:170];
ang = size(angle,2);

SPF = table2array(input_data(:,end-ang+1:end));
num = size(SPF,1); % = anzahl_zeilen = Anzahl der Partikel

image_diameter = input_data.dp_image;
image_area = input_data.A_image;

integrated_intensity = zeros (num, 1);
for i = 1:num
    ang_temp = angle;
    SPF_temp = SPF(i,:);
    nans = find(isnan(SPF_temp)); %find NaNs
    ang_temp(nans) = []; % remove those NaNs
    SPF_temp(nans) = [];
    integrated_intensity(i, 1) = trapz(ang_temp, SPF_temp);
end
disp('Done with Integrating')

% median over each size bin
med_SPF = []; num_bin = [];
for i = 1:length(bins)-1
    idx = find(image_diameter > bins(i) & image_diameter <= bins(i+1));
    num_bin(i) = length(idx);
    med_SPF(i) = nanmedian(integrated_intensity(idx));
end

if strcmp(campaign, 'IMPACTS')
    BG_offset = 238.120; disp('BG OFFSET IS NOT CALCULATED PROPERLY !')
else
    [BG_offset,BG_offset_std] = BG_offset_calib(campaign);
end
%


lastbin = length(dp_midpoint);
idx = 1:lastbin-1; %which datapoints are taken into account?
num_bin_thresh = 50; %exclude bins with less than X particles
idx = find(num_bin >= num_bin_thresh);
if binned_status == 1
    [xData, yData] = prepareCurveData(abs(med_SPF(idx)-BG_offset),dp_midpoint(idx)); % binned
    [xData_area, yData_area] = prepareCurveData(abs(med_SPF(idx)-BG_offset),dp_midpoint_area(idx)); % binned
    filename = 'binned_corrected_fixed_c';
elseif binned_status == 0
    [xData, yData] = prepareCurveData(abs(integrated_intensity-BG_offset), image_diameter); % unbinned, through all data points
    [xData_area, yData_area] = prepareCurveData(abs(integrated_intensity-BG_offset), image_area); % unbinned, through all data points
    filename = 'non_binned_corrected_fixed_c';
end
% [xData, yData] = prepareCurveData(med_SPF(idx),dp_midpoint(idx)); % without subtracted BG
% [xData, yData] = prepareCurveData(log10(abs(med_SPF(idx)-BG_offset)),log10(dp_midpoint(idx)));
% [xData, yData] = prepareCurveData(log10(abs(integrated_intensity_ice-BG_offset)), log10(image_diameter_ice));

if strcmp(campaign, 'IMPACTS') || strcmp(campaign, 'CIRRUS-HL')
    fixed_b_status = 0;
    disp('b is set free as fit parameter')
else
    fixed_b_status = 1;
    disp('b is fixed to b = 0.5 (diameter) and to b = 1 (area)')
end
    
    
if fixed_b_status == 1
    % Set up fittype and options.
    ft = fittype( 'a * (x).^0.5',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a'});
    opts = fitoptions( 'Method', 'NonLinearLeastSquares');
    opts.Robust = 'LAR';
    
    % Diameter
    fitresult = fit( xData, yData, ft, opts);
    a_ice = fitresult.a;
    ci = confint(fitresult);
    a_ice_err = fitresult.a - ci(1);
    b_ice = 0.5; %fitresult.b;

    % Set up fittype and options.
    ft = fittype( 'a * (x).^1',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a'});
    opts = fitoptions( 'Method', 'NonLinearLeastSquares');
    opts.Robust = 'LAR';

    % Area
    fitresult = fit( xData_area, yData_area, ft, opts);
    a_ice_area = fitresult.a; ci = confint(fitresult);
    a_ice_err_area = fitresult.a - ci(1);
    b_ice_area = 0.5; %fitresult.b;
else
    % Set up fittype and options.
    ft = fittype( 'a * (x).^b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a', 'b'});
    opts = fitoptions( 'Method', 'NonLinearLeastSquares', 'StartPoint', [0.5, 0.5]);
    opts.Robust = 'LAR';
    
    % Diamater
    fitresult = fit( xData, yData, ft, opts);
    a_ice = fitresult.a; ci = confint(fitresult);
    a_ice_err = fitresult.a - ci(1);
    b_ice = fitresult.b;

    % Area
    opts = fitoptions( 'Method', 'NonLinearLeastSquares', 'StartPoint', [0.1, 1]);
    opts.Robust = 'LAR';
    fitresult = fit( xData_area, yData_area, ft, opts);
    a_ice_area = fitresult.a; ci = confint(fitresult);
    a_ice_err_area = fitresult.a - ci(1);
    b_ice_area = fitresult.b;
end
disp (['Diameter = a * [Intensity - BG]^b, a = ', num2str(a_ice), ' \pm ', num2str(a_ice_err), ', b = ', num2str(b_ice),...
    ', BG = ', num2str(BG_offset), ' \pm ', num2str(BG_offset_std)])
disp (['Area = a * [Intensity - BG]^b, a = ', num2str(a_ice_area), ' \pm ', num2str(a_ice_err_area), ', b = ', num2str(b_ice_area),...
    ', BG = ', num2str(BG_offset), ' \pm ', num2str(BG_offset_std)])

%% plot: diameter
fontsize = 16;
color_grey = [0.5 0.5 0.5];
% position and size of the plots
% pos = [x1 y1 dx dy];
x1 = 0.1;
y1 = 0.15;
h = 1 - y1 - 0.03; %0.85;
dx = 0.5 - x1 - 0.015;

pos1 = [x1 y1 dx h];
pos2 = [0.5+x1 y1 dx h];

f = figure(23);
subplot('Position', pos1)
% subplot(1,2,1)
f1 = plot (integrated_intensity, image_diameter, '.', 'color', 'k', 'linewidth', .2, 'markersize', 5, 'DisplayName', 'Single-particle measurements');
hold on
if binned_status == 1
    idx = find(num_bin < num_bin_thresh);
    if size(idx,1) >=1
  %  idx = [1:idx(1)-1,idx(end)+1:lastbin]; %first and last is grey
%   idx = [1:idx(1)-1,idx(end)+1:lastbin]; 
  f33 = plot(med_SPF(idx), dp_midpoint(idx), '+', 'color', [0.8, 0.8, 0.8], 'linewidth', 4, 'markersize', 15, 'DisplayName', 'Bin medians');
    end
    % idx = [1:lastbin];
    idx = find(num_bin > num_bin_thresh);
    f2 = plot(med_SPF(idx), dp_midpoint(idx), '+', 'color', 'r', 'linewidth', 4, 'markersize', 15, 'DisplayName', 'Bin medians');
end
xx = 10^2:5*10^6;
if binned_status == 1
    f3 = plot(xx, a_ice * (xx - BG_offset).^b_ice , 'linewidth', 3, 'color', 'r', 'DisplayName', 'Fit through bin medians');
else
    f3 = plot(xx, a_ice * (xx - BG_offset).^b_ice , 'linewidth', 3, 'color', 'r', 'DisplayName', 'Fit through all points');
end

% plot bin edges
XL = [xx(1), xx(end)];
for j = 1:length(bins)
    plot(XL, [bins(j), bins(j)], '--', 'color', color_grey)
end
hold off

if binned_status == 1
    %    legend([f1, f2, f3], 'Location','southeast')
    legend([f1, f2, f3], 'Location','southeast')
else
    legend([f1, f3], 'Location','southeast')
end

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
% xlabel ( 'Partial scattering cross section' )
% ylabel ('Image area equivalent diameter [\mum]')
xlabel ('\sigma^{partial}_{scatt} [\mum]')
ylabel ('D^{geom}_p [\mum]')

xlim(XL)
ylim ([10, 1000])
yticks (([10, 20, 50, 100, 200, 500, 1000]))
xticks ([100,1000,10000,100000,1e6])
%xlim ([10, 20000])

% (['Diameter = a * [Intensity - BG]^b, a = ', num2str(a_ice), ...
%    ', b = ', num2str(b_ice), ', BG = ', num2str(BG_offset)], 'FontSize', fontsize+4)
%sgtitle ([campaign, ', all flights, ', particle_type], 'FontSize', fontsize+4)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.FontSize = fontsize-2;
ax.FontWeight = 'bold';
ax.YLabel.FontSize = fontsize+4;
ax.XLabel.FontSize = fontsize+4;

% % pcolor heatmap
x = a_ice * (integrated_intensity - BG_offset).^b_ice;
y = image_diameter;

bins_colormap = logspace(1,3,100);
dp_midpoint_colormap = bins_colormap(1:end-1)+diff(bins_colormap)./2;

[N,Xedges,Yedges] = histcounts2(abs(x),y,bins_colormap,bins_colormap);

subplot('Position', pos2)
% subplot(1,2,2)
pcolor(dp_midpoint_colormap,dp_midpoint_colormap,N')
hold on
plot(bins_colormap,bins_colormap,'Color', 'white', 'linewidth', 2)
hold off
colormap(hot) %jet
shading interp
colorbar
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
% ylabel ('Image area equivalent diameter [\mum]')
% xlabel ('Scattering equivalent diameter [\mum]', 'Position', [140, 7.8])

xlabel ('D^{scatt}_p [\mum]', 'Position', [140, 7.8])
ylabel ('D^{geom}_p [\mum]')

xticks ([10, 20, 50, 100, 200, 500, 1000])
yticks ([10, 20, 50, 100, 200, 500, 1000])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.FontSize = fontsize-2;
ax.FontWeight = 'bold';
ax.YLabel.FontSize = fontsize+4;
ax.XLabel.FontSize = fontsize+4;

set(f, 'Units', 'normalized', 'Position', [0.3, 0.1, 0.5, 0.6]); %size


save_name = ['heatmap_', particle_type, '_', filename, '_', campaign];
if save_status == true
    print([save_path, filesep,save_name],'-dpng')
end

%% plot: area
fontsize = 16;
color_grey = [0.5 0.5 0.5];
% position and size of the plots
% pos = [x1 y1 dx dy];
x1 = 0.1;
y1 = 0.15;
h = 1 - y1 - 0.03; %0.85;
dx = 0.5 - x1 - 0.015;

pos1 = [x1 y1 dx h];
pos2 = [0.5+x1 y1 dx h];

f2 = figure(24);
subplot('Position', pos1)
% subplot(1,2,1)
f1 = plot (integrated_intensity, image_area, '.', 'color', 'k', 'linewidth', .2, 'markersize', 5, 'DisplayName', 'Single-particle measurements');
hold on
if binned_status == 1
    idx = find(num_bin < num_bin_thresh);
    if size(idx,1) >=1
  %  idx = [1:idx(1)-1,idx(end)+1:lastbin]; %first and last is grey
%   idx = [1:idx(1)-1,idx(end)+1:lastbin]; 
  f33 = plot(med_SPF(idx), dp_midpoint_area(idx), '+', 'color', [0.8, 0.8, 0.8], 'linewidth', 4, 'markersize', 15, 'DisplayName', 'Bin medians');
    end
    % idx = [1:lastbin];
    idx = find(num_bin > num_bin_thresh);
    f2 = plot(med_SPF(idx), dp_midpoint_area(idx), '+', 'color', 'r', 'linewidth', 4, 'markersize', 15, 'DisplayName', 'Bin medians');
end
xx = 5*10^2:10^7;
if binned_status == 1
    f3 = plot(xx, a_ice_area * (xx - BG_offset).^b_ice_area , 'linewidth', 3, 'color', 'r', 'DisplayName', 'Fit through bin medians');
else
    f3 = plot(xx, a_ice_area * (xx - BG_offset).^b_ice_area , 'linewidth', 3, 'color', 'r', 'DisplayName', 'Fit through all points');
end

% plot bin edges
XL = [xx(1), xx(end)];
for j = 1:length(bins_area)
    plot(XL, [bins_area(j), bins_area(j)], '--', 'color', color_grey)
end
hold off

if binned_status == 1
    legend([f1, f2, f3], 'Location','southeast')
else
    legend([f1, f3], 'Location','southeast')
end

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on
xlabel ('\sigma^{partial}_{scatt} [\mum]')
ylabel ('projected area [\mum^2]')

xlim(XL)
ylim ([20, 1e6])
yticks (([20, 50, 100, 200, 500, 1000,1e4,1e5,1e6]))
xticks ([100,1000,10000,100000,1e6,1e7])
%xlim ([10, 20000])

% (['Diameter = a * [Intensity - BG]^b, a = ', num2str(a_ice), ...
%    ', b = ', num2str(b_ice), ', BG = ', num2str(BG_offset)], 'FontSize', fontsize+4)
%sgtitle ([campaign, ', all flights, ', particle_type], 'FontSize', fontsize+4)

ax = gca;
ax.YAxis(1).Color = 'k';
ax.FontSize = fontsize-2;
ax.FontWeight = 'bold';
ax.YLabel.FontSize = fontsize+4;
ax.XLabel.FontSize = fontsize+4;

% % pcolor heatmap
x = a_ice_area * (integrated_intensity - BG_offset).^b_ice_area;
y = image_area;

bins_colormap = logspace(1,7,100);
dp_midpoint_colormap = bins_colormap(1:end-1)+diff(bins_colormap)./2;

[N,Xedges,Yedges] = histcounts2(abs(x),y,bins_colormap,bins_colormap);

subplot('Position', pos2)
% subplot(1,2,2)
pcolor(dp_midpoint_colormap,dp_midpoint_colormap,N')
hold on
plot(bins_colormap,bins_colormap,'Color', 'white', 'linewidth', 2)
hold off
colormap(hot) %jet
shading interp
colorbar
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
% ylabel ('Image area equivalent diameter [\mum]')
% xlabel ('Scattering equivalent diameter [\mum]', 'Position', [140, 7.8])

xlabel ('A^{scatt}_p [\mum^2]', 'Position', [140, 7.8])
ylabel ('A^{geom}_p [\mum^2]')

xticks ([10, 20, 50, 100, 200, 500, 1000, 1e4, 1e5, 1e6, 1e7])
yticks ([10, 20, 50, 100, 200, 500, 1000, 1e4, 1e5, 1e6, 1e7])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.FontSize = fontsize-2;
ax.FontWeight = 'bold';
ax.YLabel.FontSize = fontsize+4;
ax.XLabel.FontSize = fontsize+4;

set(f2, 'Units', 'normalized', 'Position', [0.3, 0.1, 0.5, 0.6]); %size


save_name = ['heatmap_area_', particle_type, '_', filename, '_', campaign];
if save_status == true
    print([save_path, filesep,save_name],'-dpng')
end

end % end of function
%%


