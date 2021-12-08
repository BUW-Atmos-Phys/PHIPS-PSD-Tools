
function save_SD(DATA_in,DATA_in_counts, save_status, savepath, campaign, flight, tstep, a, b, bins)
% save SD as .sum and .txt

% savepath = 'C:\Users\wa9929\Desktop';

if save_status == 1
    % Check if quicklooks folder exists
    if ~exist(savepath,'dir')
        mkdir(savepath)
    end

    
    input_name = inputname(1); %get the name of the variable
    input_name = strsplit(input_name,'SD'); %delete the 'SD_' in the name since it appears double
    input_name = input_name{2};
    

    
    save_filename = ['SD_', flight , '_', num2str(tstep), 's_', input_name, '.sum'];
    save([savepath,filesep,save_filename], 'DATA_in' ,'-ascii', '-double');

    
    %produce .txt file
    DATA = SD_sum_to_txt(DATA_in, campaign, flight, tstep, a, b, bins); %function that introduces the header

    save_filename = ['SD_', flight , '_', num2str(tstep), 's_', input_name, '.txt'];
    export(DATA,'file',[savepath,filesep,save_filename],'WriteVarNames',false); %save

    %% same but for counts
    
    save_filename = ['SD_', flight , '_', num2str(tstep), 's_', input_name, '_counts.sum'];
    save([savepath,filesep,save_filename], 'DATA_in_counts' ,'-ascii', '-double');

    %produce .txt file
    DATA = SD_sum_to_txt_counts(DATA_in_counts, campaign, flight, tstep, a, b, bins); %function that introduces the header

    save_filename = ['SD_', flight , '_', num2str(tstep), 's_', input_name, '_counts.txt'];
    export(DATA,'file',[savepath,filesep,save_filename],'WriteVarNames',false); %save

end