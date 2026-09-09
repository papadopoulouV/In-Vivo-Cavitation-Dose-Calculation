function [cavData, ampdata, fftdata] = pcdFFT4(voltageData, time, saveDir, freqOptions, PBS_stack, subtractPBS)
% Passive Cavitation Detection (PCD) FFT
% Created by Kelly VanTreeck
% Last updated by Kelly VanTreeck on 7/28/2025
% DOI:
%  
% Outputs:
%   cavData  – struct of cavitation metrics for each A‐line
%   ampdata  – struct of raw and filtered time‐domain signals
%   fftdata  – struct of raw and band‐limited FFT results

    % Handle optional inputs for PBS subtraction to generate netFFT
    if nargin < 5
        PBS_stack.FFTmean = [];
    end
    if nargin < 6
        subtractPBS = false;
    end

    % Unpack frequency options
    Fc = freqOptions.Fc;
    Fs = freqOptions.Fs;
    scPassband = freqOptions.scPassband;
    icPassband = Fc * freqOptions.icPassband;

    numalines = size(voltageData, 1);

    % Preallocate cavData fields (1 × numalines)
    cavData.scd_AUC   = nan(1, numalines);
    cavData.scd_peaks = nan(1, numalines);
    cavData.icd       = nan(1, numalines);
    cavData.totalcav  = nan(1, numalines);
    cavData.ucd_AUC   = nan(1, numalines);
    cavData.ucd_peaks = nan(1, numalines);
    cavData.shd_AUC   = nan(1, numalines);
    cavData.shd_peaks = nan(1, numalines);

    % Compute FFT length and frequency axis
    L0 = size(voltageData, 2);
    L  = 2^(nextpow2(L0) + 1);
    Fn = Fs / 2;
    freq = Fs * (0:(L/2)) / L;   % [L/2+1 × 1]
    fftdata.freq = freq;

    % Precompute PBS baseline if requested
    if subtractPBS && ~isempty(PBS_stack.FFTmean)
        FFT_baseline = mean(PBS_stack.FFTmean, 2, 'omitnan');  % [L/2+1 × 1]
    end

    % Preallocate ampdata & fftdata arrays
    ampdata.raw_data  = zeros(L0, numalines);
    ampdata.amp       = zeros(L0, numalines);
    ampdata.amp_i     = zeros(L0, numalines);
    ampdata.amp_total = zeros(L0, numalines);

    fftdata.raw_fft   = zeros(L/2+1, numalines);
    fftdata.FFTnet    = zeros(L/2+1, numalines);
    fftdata.fft_s     = zeros(L/2+1, numalines);
    fftdata.fft_i     = zeros(L/2+1, numalines);
    fftdata.fft_total = zeros(L/2+1, numalines);
    fftdata.fft_u     = zeros(L/2+1, numalines);
    fftdata.fft_sh    = zeros(L/2+1, numalines);

%% Filters
    % Stable Cavitation (harmonic bands)
    sc_harmonics = scPassband(1) : scPassband(2);
    filters.SC = cell(1, length(sc_harmonics));
    for i = 1:length(sc_harmonics)
        harm = sc_harmonics(i);
        Wp = [ (harm - 0.05)*Fc, (harm + 0.05)*Fc ] / Fn;
        Ws = [ (harm - 0.10)*Fc, (harm + 0.10)*Fc ] / Fn;
        [n, Ws_adj] = cheb2ord(Wp, Ws, 1, 60);
        [z, p, k] = cheby2(n, 60, Ws_adj, 'bandpass');
        [sos, g] = zp2sos(z, p, k);
        filters.SC{i} = {sos, g, Ws_adj};
    end

    % Inertial Cavitation (broadband around icPassband)
    Wp_ic = icPassband / Fn;
    Ws_ic = [ (icPassband(1) - 0.1), (icPassband(2) + 0.1) ] / Fn;
    [n_ic, Ws_ic_adj] = cheb2ord(Wp_ic, Ws_ic, 1, 60);
    [z_ic, p_ic, k_ic] = cheby2(n_ic, 60, Ws_ic_adj);
    [sos_ic, g_ic] = zp2sos(z_ic, p_ic, k_ic);
    filters.IC = {sos_ic, g_ic, Ws_ic_adj};

    % Total Cavitation
    Wp_total = icPassband / Fn;
    Ws_total = [ (icPassband(1) - 0.15), (icPassband(2) + 0.15) ] / Fn;
    [n_total, Ws_total_adj] = cheb2ord(Wp_total, Ws_total, 1, 60);
    [z_total, p_total, k_total] = cheby2(n_total, 60, Ws_total_adj);
    [sos_total, g_total] = zp2sos(z_total, p_total, k_total);
    filters.TOTAL = {sos_total, g_total, Ws_total_adj};

%% Compute raw FFT & all cavitation metrics for each A-line
    for ii = 1:numalines
        amp0 = voltageData(ii, :)';
        amp  = amp0 - mean(amp0);   % remove DC offset

        ampdata.time       = time;
        ampdata.raw_data(:, ii) = amp0;
        ampdata.amp(:, ii)      = amp;

        % Raw FFT
        raw_fft = fft(amp, L);
        raw_fft = abs(raw_fft / L0);
        raw_fft = raw_fft(1 : L/2 + 1);
        raw_fft(2:end-1) = 2 * raw_fft(2:end-1);
        fftdata.raw_fft(:, ii) = raw_fft;

        % Subtract PBS baseline if expData
        if subtractPBS && exist('FFT_baseline', 'var')
            fftdata.FFTnet(:, ii) = raw_fft - FFT_baseline;
        else
            fftdata.FFTnet(:, ii) = raw_fft;
        end

        % --- Stable Cavitation (harmonics) ---
        fft_s_combined = zeros(length(freq), 1);
        peak_sum = 0;
        for i = 1:length(filters.SC)
            sos = filters.SC{i}{1};
            g   = filters.SC{i}{2};
            Ws  = filters.SC{i}{3};

            amp_s = filtfilt(sos, g, amp);
            fft_s = fft(amp_s, L);
            fft_s = abs(fft_s / L0);
            fft_s = fft_s(1 : L/2 + 1);
            fft_s(2:end-1) = 2 * fft_s(2:end-1);

            mask = (freq >= Ws(1)*Fn) & (freq <= Ws(2)*Fn);
            fft_s(~mask) = 0;

            fft_s_combined = fft_s_combined + fft_s;
            peak_sum       = peak_sum + max(fft_s(mask));
        end
        fftdata.fft_s(:, ii)    = fft_s_combined;
        cavData.scd_AUC(ii)     = trapz(freq * 1e6, fft_s_combined);
        cavData.scd_peaks(ii)   = peak_sum;

        % --- Inertial Cavitation (broadband) ---
        sos = filters.IC{1};
        g   = filters.IC{2};
        Ws  = filters.IC{3};
        amp_i = filtfilt(sos, g, amp);
        ampdata.amp_i(:, ii) = amp_i;

        fft_i = fft(amp_i, L);
        fft_i = abs(fft_i / L0);
        fft_i = fft_i(1 : L/2 + 1);
        fft_i(2:end-1) = 2 * fft_i(2:end-1);
        fft_i(freq < Ws(1)*Fn | freq > Ws(2)*Fn) = 0;

        % Zero out any SC harmonics in the IC band
        for icomb = freqOptions.icPassband(1) : freqOptions.icPassband(2)
            idx = (freq < Fc*(icomb + 0.1)) & (freq > Fc*(icomb - 0.1));
            fft_i(idx) = 0;
        end
        fftdata.fft_i(:, ii) = fft_i;
        cavData.icd(ii)     = trapz(freq * 1e6, fft_i);

        % --- Total Cavitation (harmonics + broadband) ---
        sos   = filters.TOTAL{1};
        g     = filters.TOTAL{2};
        Ws    = filters.TOTAL{3};
        amp_total = filtfilt(sos, g, amp);
        ampdata.amp_total(:, ii) = amp_total;

        fft_total = fft(amp_total, L);
        fft_total = abs(fft_total / L0);
        fft_total = fft_total(1 : L/2 + 1);
        fft_total(2:end-1) = 2 * fft_total(2:end-1);
        fft_total(freq < Ws(1)*Fn | freq > Ws(2)*Fn) = 0;

        fftdata.fft_total(:, ii) = fft_total;
        cavData.totalcav(ii)     = trapz(freq * 1e6, fft_total);

    end 

    fprintf('%s: Analysis done.\n', string(datetime('now','Format','HH:mm:ss')));
end
