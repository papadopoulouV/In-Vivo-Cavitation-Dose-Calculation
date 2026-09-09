
% Passive cavitation data analysis for one day of biofilm wound treatment
% Last updated: 8-01-02025 by Kelly VanTreeck
% DOI: 

% Before running this script, please review the readME, ensure you have required
% function, and download matfiles from: X

% REQUIRED FUNCTIONS
%     pcdFFT4
%     natsortfiles: https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

%Fully updated code using pcdFFT4 to filter and generate the following cavitation values:   
%       harmonicAUC
%       broadbandAUC
%       totalcavAUC (broadband + harmonics)
%This script will generate an export with the following:
%       PBSData: meanFFT, cavData
%       expData: net cavData (PBS subtracted), meanFFT

%Run this code for each day of example data. It will prompt you to save an export file, which can be loaded into
%   cavitationAnalysis_Step2

clear all; close all; clc

%% Load example data from .mat file
[matfile, matpath] = uigetfile('*.mat', 'Select example .mat file');
if isequal(matfile, 0)
    error('No file selected. Exiting script.');
end

load(fullfile(matpath, matfile), 'PBSData', 'expData');

%% Set parameters
cf = 1.1;
ALINES = 200;

Fs = 200; 
T = 1/Fs;
L = 40000;
T = (0:L-1)*T;
f = Fs*(0:(L/2))/L;
freq = f;

freqOptions.Fc = cf;
freqOptions.Fs = Fs;
freqOptions.scPassband = [3 5];
freqOptions.icPassband = [3 5];

%% Process PBS Data
for i = 1:length(PBSData)
    voltageData = PBSData(i).voltData;
    time = PBSData(i).time;

    [cavData, ~, fftdata] = pcdFFT4(voltageData, time, [], freqOptions);

    PBSData(i).freq          = fftdata.freq;
    PBSData(i).FFT           = fftdata.raw_fft;
    PBSData(i).FFTharmonic   = fftdata.fft_s;
    PBSData(i).FFTbroadband  = fftdata.fft_i;
    PBSData(i).FFTtotalcav   = fftdata.fft_total;
    PBSData(i).FFTmean       = mean(fftdata.raw_fft, 2);

    PBSData(i).harmonicAUC   = cavData.scd_AUC;
    PBSData(i).broadband     = cavData.icd;
    PBSData(i).totalcav      = cavData.totalcav;

    % Transmit OFF Detection
    [~, tx_row] = max(PBSData(i).FFTmean);
    if ~isempty(tx_row)
        max_value = max(PBSData(i).FFT(tx_row, :));
        threshold = 0.1 * max_value;

        transmitOff = PBSData(i).FFT(tx_row, :) < threshold;
        numOff = sum(transmitOff);

        if numOff > 0
            fprintf('⚠️ Transmit OFF detected in %d A-lines of %s\n', ...
                numOff, PBSData(i).filename);

            PBSData(i).FFT(:, transmitOff) = NaN;
            PBSData(i).FFTmean = mean(PBSData(i).FFT, 2, 'omitnan');

            metrics = {'harmonicAUC', 'broadband', 'totalcav'};
            for m = 1:numel(metrics)
                temp = PBSData(i).(metrics{m});
                if length(temp) == size(PBSData(i).FFT, 2)
                    temp(transmitOff) = NaN;
                    PBSData(i).(metrics{m}) = temp;
                end
            end
        end
    end
end

%% Process Experimental Data & Detect Transmit OFF
errorFiles = {};  % Initialize error tracking

for i = 1:length(expData)
    try
        % Match PBS data by animal number
        pbsIdx = find([PBSData.animalnum] == expData(i).animalnum, 1);
        if ~isempty(pbsIdx)
            PBS_stack.FFTmean = PBSData(pbsIdx).FFTmean;
        else
            warning(['No matching PBS pressure found for: ', expData(i).filename]);
            PBS_stack.FFTmean = [];
        end

        voltageData = expData(i).voltData;
        time = expData(i).time;

        % FFT and cavitation processing
        [cavData, ampdata, fftdata] = pcdFFT4(voltageData, time, [], freqOptions, PBS_stack, true);

        % Store results
        expData(i).freq       = fftdata.freq;
        expData(i).FFTraw     = fftdata.raw_fft;
        expData(i).FFTnet     = fftdata.FFTnet;
        expData(i).FFTmean    = mean(fftdata.raw_fft, 2, 'omitnan');

        expData(i).harmonicAUC = cavData.scd_AUC;
        expData(i).broadband   = cavData.icd;
        expData(i).totalcav    = cavData.totalcav;

        % Transmit-OFF detection
        if isfield(expData(i), 'timepoint') && expData(i).timepoint == 0
            [~, tx_row] = max(expData(i).FFTmean);
            if ~isempty(tx_row)
                max_val = max(expData(i).FFTraw(tx_row, :));
                threshold = 0.1 * max_val;

                transmitOff = expData(i).FFTraw(tx_row, :) < threshold;
                numOff = sum(transmitOff);

                if numOff > 0
                    fprintf('⚠️ Transmit OFF in %d A-lines: %s\n', ...
                            numOff, expData(i).filename);

                    expData(i).FFTraw(:, transmitOff) = NaN;
                    expData(i).FFTnet(:, transmitOff) = NaN;

                    metrics = {'harmonicAUC', 'broadband', 'totalcav'};
                    for m = 1:numel(metrics)
                        vec = expData(i).(metrics{m});
                        if length(vec) == size(expData(i).FFTraw, 2)
                            vec(transmitOff) = NaN;
                            expData(i).(metrics{m}) = vec;
                        end
                    end
                end
            end
        end

        % Validate filename field
        if ~ischar(expData(i).filename)
            error('Filename is not a character vector');
        end

    catch ME
        warning(['Failed to process file: ', expData(i).filename]);
        disp(['Error: ', ME.message]);
        for s = 1:length(ME.stack)
            fprintf('  In %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
        end
        errorFiles{end+1} = expData(i).filename;
        continue
    end
end

%% Net cavitation calculation (PBS baseline subtracted)
% Calculate mean of PBS baseline signal
for i = 1:length(PBSData)
    PBSData(i).harmonicAUC_mean   = mean(PBSData(i).harmonicAUC, 'omitnan');
    PBSData(i).harmonicAUC_SD     = std(PBSData(i).harmonicAUC, 'omitnan');
    
    PBSData(i).broadband_mean     = mean(PBSData(i).broadband, 'omitnan');
    PBSData(i).broadband_SD       = std(PBSData(i).broadband, 'omitnan');

    PBSData(i).totalcav_mean      = mean(PBSData(i).totalcav, 'omitnan');
    PBSData(i).totalcav_SD        = std(PBSData(i).totalcav, 'omitnan');

end

for q = 1:length(expData)
    expData(q).FFTmean = mean(expData(q).FFTraw, 2,'omitnan');
    matchIdx = find([PBSData.animalnum] == expData(q).animalnum, 1);
    
    if ~isempty(matchIdx)
        pbs = PBSData(matchIdx);
        
        % Subtract PBS mean from each 1x200 vector
        expData(q).harmonicAUCnet   = expData(q).harmonicAUC   - pbs.harmonicAUC_mean;
        expData(q).broadbandnet     = expData(q).broadband     - pbs.broadband_mean;
        expData(q).totalcavnet      = expData(q).totalcav      - pbs.totalcav_mean;

        % Clip negative values to zero
        fieldsToClip = {'harmonicAUCnet', 'broadbandnet', 'totalcavnet'};
        
        for f = 1:numel(fieldsToClip)
            val = expData(q).(fieldsToClip{f});
            if ~isempty(val) && isnumeric(val)
                % Only clip finite (non-NaN) values
                clipped = val; % preserve NaNs
                finiteMask = isfinite(val);  % true for non-NaN, non-Inf
                clipped(finiteMask & val < 0) = 0;
                expData(q).(fieldsToClip{f}) = clipped;
            end
        end
    else
        warning('No matching PBSData found for animal %d', expData(q).animalnum);
    end
end

disp('Net values calculated');

%% Cav data export: All data, no averaging

% List of fields to export from expData
fieldsToExport = {
    'filename', ...       
    'pressure', ...         
    'animalnum', ...        
    'treatmentnum', ...     
    'timepoint', ...        
    'freq', ...             
    'FFTmean', ...          
    'harmonicAUCnet', ...      
    'broadbandnet', ...       
    'totalcavnet', ...        
};

% Initialize export structure
exportCavData = struct();

% Loop through each element of expData and extract specified fields
for i = 1:numel(expData)
    for f = 1:numel(fieldsToExport)
        field = fieldsToExport{f};
        if isfield(expData(i), field)
            exportCavData(i).(field) = expData(i).(field);
        else
            warning('Field "%s" not found in expData.', field);
            exportCavData(i).(field) = [];
        end
    end
end

% Prompt user to select a folder
savepath = uigetdir(matpath, 'Select folder to save data');

% Check if user canceled the dialog
if savepath ~= 0
    % Remove .mat extension from original filename
    [~, baseName, ~] = fileparts(matfile);

    % Create new filename
    exportFilename = fullfile(savepath, [baseName, '_cavDataExport.mat']);

    % Save the variables
    save(exportFilename, 'exportCavData', 'PBSData');

    fprintf('✅ Data saved successfully to:\n%s\n', exportFilename);
else
    disp('❌ Save operation canceled by user.');
end
