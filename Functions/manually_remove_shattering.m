function [ShatteringFlag,quality_string] = manually_remove_shattering(time_PHIPS,ShatteringFlag,campaign,flight)
% Manually remove shattering influenced periods

quality_string = 'No manual quality control.';

% ACLOUD
if strcmp(campaign,'ACLOUD')
    if strcmp(flight,'Flight 170527')
        ind = find(time_PHIPS>=datenum('27-May-2017 14:54:33') & ...
            time_PHIPS<=datenum('27-May-2017 15:30:04'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 14:54:33 and 15:30:04 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170602')
        ind = find(time_PHIPS>=datenum('02-Jun-2017 10:39:44') & ...
            time_PHIPS<=datenum('02-Jun-2017 11:36:50'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 10:39:44 and 11:36:50 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170604')
        ind = find(time_PHIPS>=datenum('04-Jun-2017 10:53:39') & ...
            time_PHIPS<=datenum('04-Jun-2017 11:23:10'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 10:53:39 and 11:23:10 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170605')
        ind = find(time_PHIPS>=datenum('05-Jun-2017 12:12:41') & ...
            time_PHIPS<=datenum('05-Jun-2017 12:48:10'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 12:12:41 and 12:48:10 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170616')
        ind = find(time_PHIPS>=datenum('2017-06-16 07:33:00') & ...
            time_PHIPS<=datenum('2017-06-16 08:42:00'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 07:33:00 and 08:42:00 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170617')
        ind = find(time_PHIPS>=datenum('2017-06-16 11:05:00') & ...
            time_PHIPS<=datenum('17-Jun-2017 11:27:20'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 11:05:00 and 11:27:20 UTC was manually set to 0 (shattering influenced period).';
    elseif strcmp(flight,'Flight 170618')
        ind = find(time_PHIPS>=datenum('2017-06-18 16:12:30') & ...
            time_PHIPS<=datenum('2017-06-18 17:20:00'));
        ShatteringFlag(ind) = 0;
        quality_string = 'Quality flag between 16:12:30 and 17:20:00 UTC was manually set to 0 (shattering influenced period).';
    end
elseif strcmp(campaign,'CIRRUS-HL')
    quality_string = 'During CIRRUS-HL noise triggers were observed outside clouds. Noise triggers were removed based on the shape of their light scattering function using a correlation-based pattern matching approach, where known noise patterns from RF03 particles #8–26 were used to compute cosine similarity scores for each row in the dataset, and rows exceeding a predefined threshold were identified as noise and removed.';
end

end