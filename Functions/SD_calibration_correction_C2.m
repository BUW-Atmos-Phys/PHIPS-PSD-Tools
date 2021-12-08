function input_data = SD_calibration_correction_C2(input_data, campaign, particle_type, save_status, save_path);
% correction changes nothing for the actual calibration if the data is binned
% for non-binned, the change is 0.2%


% save_status = false;
% input_data = IceData;

%% C2 correction

for ii = 1:2 %(uncorrected and corrected)
    if ii == 1
        x = input_data.diameter_C2;
        y = (input_data.diameter_C2-input_data.diameter_C1)./input_data.diameter_C2*100;
    elseif ii==2
        x = C2;
        y = (C2-C1)./C2 * 100;
    end

x(x>200) = NaN; x(x<20) = NaN;
y(y<-100) = NaN; y(y>100) = NaN;

% % pcolor heatmap
bins_x = linspace(10,200,100);
midpoint_x = bins_x(1:end-1)+diff(bins_x)./2;
bins_y = -50:2:110;
% bins_y = linspace(-50,100,200);
midpoint_y = bins_y(1:end-1)+diff(bins_y)./2;

[N,Xedges,Yedges] = histcounts2(abs(x),y,bins_x,bins_y);

bins_C2 = [10:2:50,55:5:200]; 
midpoint_C2 = bins_C2(1:end-1)+diff(bins_C2)./2;
% median over each size bin
med_C2 = []; % num_bin = [];
for i = 1:length(bins_C2)-1
    idx = find(x > bins_C2(i) & x <= bins_C2(i+1));
%    num_bin(i) = length(idx);
    med_C2(i) = nanmedian(y(idx));
end

med_C2(midpoint_C2<30) = NaN;

[xData, yData] = prepareCurveData(midpoint_C2, med_C2);

% Set up fittype and options.
ft = fittype('1 + a * exp(-(x)./b)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b'});
opts = fitoptions( 'Method', 'NonLinearLeastSquares');
opts.StartPoint = [10, 10];%,30,0]; % wenn x = C1
opts.Upper = [1000, 100];%, 50, 0.1];
opts.Lower = [0,0];%,1,-0.1];
%opts.StartPoint = [100, 2,30,0];

fitresult_C1C2 = fit( xData, yData, ft, opts); 
a = fitresult_C1C2.a;
b = fitresult_C1C2.b;
xx = 10:200;
yy = a * exp(-(xx)./b);

color_light_grey = [0.7 0.7 0.7];
f = figure(77);
subplot(1,2,ii)
pcolor(midpoint_x,midpoint_y,N')
    hold on
plot(midpoint_C2,med_C2,'+', 'Color',color_light_grey, 'LineWidth',2)
plot(midpoint_x, midpoint_x*0,'Color', 'white','LineWidth',2)
if ii==1
    plot(xx, yy,'Color', 'r','LineWidth',2)
    plot(1, 1,'Color', 'white','LineWidth',2) % just for the legend
    plot(1, 1,'Color', 'white','LineWidth',2)
    plot(1, 1,'Color', 'white','LineWidth',2)
end

    hold off
colormap(jet)
shading interp
colorbar
xlabel ('Diameter C2')
ylabel ('Relative Difference (C2-C1)/C2 [%]')

    if ii==1
        title ([campaign, ', ', particle_type, ', Uncorrected'])
        legend('Ice', 'Bin Medians', 'Fit through medians', 'y = a * exp( - x / b)', ['a = ', num2str(a)],['b = ', num2str(b)])
    elseif ii==2
        title ([campaign, ', ', particle_type, ', Corrected'])
        legend('Ice', 'Bin Medians')
    end

if ii==1
    % correction C2 for small sizes
    % a = 60.15 ;   b = 23.94 ;
    C1 = input_data.diameter_C1;
    C2 = input_data.diameter_C2;
    % C2-C1/C2 * 100 = a * exp(-C2/b) 
    % 1-C1/C2 = a/100 * exp
    % C1 = -C2 * ( a/100 * exp - 1)
    C2 = -C2 .* (a./100 .* exp(-C2/b) - 1); % correction for small sizes
end

end

set(f, 'Units', 'normalized', 'Position', [0.3, 0.1, 0.5, 0.7]); %size 

save_name = ['\ratio_heatmap_C1_C2_',particle_type ,'_', campaign];
if save_status == true
   print([save_path, save_name],'-dpng')
end

%% C1 vs C2, corrected vs uncorrected

for ii = 1:2
    if ii == 1
        x = input_data.diameter_C1;
        y = input_data.diameter_C2;
    elseif ii == 2
        x = C1;
        y = C2;
    end

bins_colormap = logspace(1,3,100);
dp_midpoint_colormap = bins_colormap(1:end-1)+diff(bins_colormap)./2;
[N,Xedges,Yedges] = histcounts2(x,y,bins_colormap,bins_colormap);

f = figure(78);
    subplot(1,2,ii)
pcolor(dp_midpoint_colormap,dp_midpoint_colormap,N')
    hold on
plot(bins_colormap,bins_colormap,'Color','white', 'LineWidth', 1)
    hold off
colormap(jet)
shading interp
colorbar
set(gca, 'XScale','log'), set(gca, 'YScale','log')
xticks ([10, 20, 50, 100, 200, 500, 1000]), yticks ([10, 20, 50, 100, 200, 500, 1000])
xlabel ('Diameter C1'), ylabel ('Diameter C2')
    if ii==1
        title ([campaign, ', ', particle_type, ', Unorrected'])
    elseif ii==2
        title ([campaign, ', ', particle_type, ', Corrected'])
    end
end

set(f, 'Units', 'normalized', 'Position', [0.3, 0.1, 0.5, 0.7]); %size 

save_name = ['\C1_vs_C2_',particle_type ,'_corrected_', campaign];
if save_status == true
   print([save_path, save_name],'-dpng')
end

%% return corrected value

input_data.diameter_C1 = C1;
input_data.diameter_C2 = C2;


end % end of function
