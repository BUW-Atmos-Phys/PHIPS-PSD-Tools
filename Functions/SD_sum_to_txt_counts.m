function DATA = SD_sum_to_txt_counts(DATA_in, campaign, flight, tstep_sec, a, b, bins)
% Makes a txt file from sum with header

% Create cell, that contains all data, DATA_out
DATA_out = {};
row1 = ['PHIPS Size Distribution for ',num2str(campaign),' Flight ',num2str(flight),'.'];

row2 = ['Preliminary data. Please contact M. Schnaiter (martin.schnaiter@kit.edu) before publishing.'];
%row2 = ['Filename for PHIPS raw data: ',filename];
%row3 = ['Calibration coefficients: a = ',num2str(a),', b = ',num2str(b),' and c = ',num2str(c)];
row3 = ['Data produced with ',num2str(tstep_sec),' s averaging'];
row4 = ['Calibration from integrated SPF data (Dp = a * I^b), Calibration Parameters:  a = ', num2str(a), ', b  = ', num2str(b)]; 

row5 = ['Size bin edges: ', num2str(bins)];
row6 = ['Size bin mean diameters:      ',num2str(DATA_in(1,5:end))];
%row72 = ['Use bin lower diameters for plotting in Matlab.'];

row7 = ['1st Column: Time bin midpoints (yyyy-mm-dd HH:MM:SS)'];
row8 = ['2nd Column: Qualityflag (based on the ratio of large 2DS particles, 1 means good data quality, 0 means potential shattering)'];
row9 = ['3rd Column: Total Counts'];
% row10 = ['4th Column: dNtot, Uncertainty in total concentration [L-1]'];
row10 = ['4th Column onwards: Counts per bin'];

DATA_out{1,1} = {};
DATA_out{2,1} = row1;
DATA_out{3,1} = {};
DATA_out{4,1} = row2;
DATA_out{5,1} = {};
DATA_out{7,1} = row3;
DATA_out{8,1} = row4;
DATA_out{9,1} = {};
DATA_out{10,1} = row5;
DATA_out{11,1} = row6;
DATA_out{12,1} = {};
DATA_out{13,1} = row7;
DATA_out{14,1} = row8;
DATA_out{15,1} = row9;
DATA_out{16,1} = row10;
%DATA_out{17,1} = row11;
DATA_out{17,1} = {};

num_data = DATA_in(2:end,2:end);
num = num2str(num_data);
timestr = datestr(DATA_in(2:end,1),'yyyy-mm-dd HH:MM:SS');

for i=1:size(DATA_in,1)-1
    DATA_out{17+i,1} = timestr(i,:);
    DATA_out{17+i,2} = num(i,:);
end
DATA = cell2dataset(DATA_out);

