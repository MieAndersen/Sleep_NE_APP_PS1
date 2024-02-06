%% 1) Define mouse data

clear all
close all

%1: FP file
%2: EEG file
%3: sleep score file
%4: rig 1 blue channel
%5: rig 1 violet channel
%6: TTL pulse channel
%7: start time for fit 
%8: stop time for fit 
%9: time correction for sleep scoring

example_mouseID = {'C:\Users\username\data\FP_data_folder' 'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' 'ch465' 'ch405' 'PtC0' 1000 10000 0};

mouse = example_mouseID;

%% 2) Load FP data (batch I)
  
data = TDTbin2mat(mouse{1});

signal_fs = data.streams.(mouse{4}).fs;

signal_465 = data.streams.(mouse{4}).data;

control_405 = data.streams.(mouse{5}).data;
 
% removing FP trace prior to first TTL pulse
TTL_FP = data.epocs.(mouse{6}).onset;
TTL_gap = diff(TTL_FP) > 5 + 1;
if isempty(find(TTL_gap == 1, 1))
    TTL_onset = TTL_FP(1);  % when TTL pulse train is only started once
else 
    TTL_onset = TTL_FP(find(TTL_gap==1)+1); % when TTL pulse train is started more than once
end

first_TTL = TTL_onset(1)*signal_fs;

onset_FP = first_TTL;

signal_465 = signal_465(round(onset_FP):end);
control_405 = control_405(round(onset_FP):end);

%% 3) Normalize and plot 

MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:1:length(signal_465);
sec_signal = fs_signal/signal_fs;

reg_465 = polyfit(control_405(round(mouse{7}*signal_fs):round(mouse{8}*signal_fs)), signal_465(round(mouse{7}*signal_fs):round(mouse{8}*signal_fs)), 1);
a = reg_465(1);
b = reg_465(2);
controlFit = a.*control_405 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat_465 = (signal_465 - controlFit)./controlFit;
delta_465 = normDat_465 * 100;

figure
a = subplot(4,1,1);
plot(sec_signal(1000:end), control_405(1000:end));
title('raw control');
b = subplot(4,1,2);
plot(sec_signal(1000:end), signal_465(1000:end));
title('raw signal');
c = subplot(4,1,3);
plot(sec_signal(1000:end), signal_465(1000:end));
hold on
plot(sec_signal(1000:end), controlFit(1000:end));
title('fitted control');
d = subplot(4,1,4);
plot(sec_signal(1000:end), delta_465(1000:end));
title('normalized signal 465');
linkaxes([a,b,c,d],'x');

% smoothing traces
delta465_filt = filtfilt(MeanFilter,1,double(delta_465));
delta465_filt = detrend(delta465_filt);

% downsampling traces for plotting
ds_factor_FP = 100; % also used for plotting later (section 9b)
ds_delta465_filt = downsample(delta465_filt, ds_factor_FP);

ds_sec_signal = downsample(sec_signal, ds_factor_FP); % for plotting

figure
plot( ds_sec_signal(1000:end),ds_delta465_filt(1000:end))
title('NE');
 
%% 4) loading and plotting EEG and EMG raw data

% Make sure the "ExpToolbox" is added to the matlab path - can be obtained from Viewpoint ( https://www.viewpoint.fr/)

% Import EEG raw data to matlab
Info=loadEXP([mouse{2}],'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=inf; %inf to include all data (until last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;


% Plot of EEG and EMG traces
figure
a = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');

%% 4b) Interpolate gaps in EEG/EMG 

% some recordsings were automatically backed up at 18.00 and therefore .exp
% were converted. These recordigns have a gap in EEG/EMG traces that must be
% removed to run the spectogram function. Here they are interpolated (the
% interpolated part should not be used for analysis - make sure they don't
% occur in vicinity of laser stim)

EEG_NaNs = find(isnan(EEG_rawtrace));

if ~isempty(EEG_NaNs)
    EEG_nans = isnan(EEG_rawtrace); % identify NaNs
    EEG_nans_n = 1:numel(EEG_rawtrace); % vector of indices, used for indexing below
    EEG_intpl = EEG_rawtrace; % copy of velocity used for interpolation
    EEG_intpl(EEG_nans) = interp1(EEG_nans_n(~EEG_nans), EEG_rawtrace(~EEG_nans), EEG_nans_n(EEG_nans));
    EEG_rawtrace = EEG_intpl;
    
    EMG_nans = isnan(EMG_rawtrace); % identify NaNs
    EMG_nans_n = 1:numel(EMG_rawtrace); % vector of indices, used for indexing below
    EMG_intpl = EMG_rawtrace; % copy of velocity used for interpolation
    EMG_intpl(EMG_nans) = interp1(EMG_nans_n(~EMG_nans), EMG_rawtrace(~EMG_nans), EMG_nans_n(EMG_nans));
    EMG_rawtrace = EMG_intpl;
end

% interpolated - therefore remove this one
if isnan(EEG_rawtrace(end))
    EEG_rawtrace = EEG_rawtrace(1:end-1);
    EMG_rawtrace = EMG_rawtrace(1:end-1);
    EEG_time = EEG_time(1:end-1);
end

%% 5)  Alingment of EEG recording and FP recording

% TTL pulse from FP
TTL_pulse = Data(3,1:end);
onset_EEG = find(diff(TTL_pulse>1*10^-3));
onset_EEG_time = onset_EEG/sampling_freq;
onset_EEG_time_diff = diff(onset_EEG_time);

TTL_gap_EEG = onset_EEG_time_diff > 10;
if isempty(find(TTL_gap_EEG==1, 1))
    onset_EEG = onset_EEG(1);                                              
else 
    onset_EEG = onset_EEG(find(onset_EEG_time_diff>10)+1);
end

TTL_EEG_onset = onset_EEG./sampling_freq;

%Cutting EEG/EMG traces leading up to first TTL 
% Removing first seconds of EEG and EMG raw traces to align with FP trace
EMG_rawtrace_cut = EMG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_rawtrace_cut = EEG_rawtrace(round(TTL_EEG_onset*sampling_freq):end);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;
ds_EEG_time = downsample(EEG_time_cut, 10);


figure
a = subplot(3,1,1);
plot( ds_sec_signal,ds_delta465_filt)
title('NE');
b = subplot(3,1,2);
    plot(EEG_time_cut, EEG_rawtrace_cut); 
    xlabel('time (s)');
    ylabel('EEG (V)');
c = subplot(3,1,3);
    plot(EEG_time_cut, EMG_rawtrace_cut); 
    xlabel('time (s)');
    ylabel('EMG (V)');
linkaxes([a,b,c],'x');


%% 6) open EEG scoring

time_correction = mouse{9};
EEG_sleepscore = xlsread(mouse{3});

% Create binary vectors for sleep stages
%Awake
wake_onset = rmmissing(EEG_sleepscore(:, 2)); % onset of each wake bout in sec (NaNs are removed)
wake_duration = rmmissing(EEG_sleepscore(:, 3)); % duration of each wake bout in sec (NaNs are removed)

%Slow-wave sleep
sws_onset = rmmissing(EEG_sleepscore(:, 6)); % onset of each SWS bout in sec (NaNs are removed)
duration_sws = rmmissing(EEG_sleepscore(:, 7)); % duration of each SWS bout in sec (NaNs are removed)

%REM
REM_onset = rmmissing(EEG_sleepscore(:, 10)); % onset of each REM bout in sec (NaNs are removed)
REM_duration = rmmissing(EEG_sleepscore(:, 11)); % duration of each REM bout in sec (NaNs are removed)


% Many EEG scorings don't start at time 0 - which shifts the timeline of the
% scoring relative to the EEG/EMG traces - this is corrected for below
if min([wake_onset(1), sws_onset(1), REM_onset(1)]) ~= 0
    EEG_scoring_onset = min([wake_onset(1), sws_onset(1), REM_onset(1)]); % determines the number of seconds to be subtracted
    wake_onset = wake_onset - EEG_scoring_onset;
    sws_onset = sws_onset - EEG_scoring_onset;
    REM_onset = REM_onset - EEG_scoring_onset;
end

wake_onset = wake_onset+time_correction;
sws_onset = sws_onset+time_correction;
REM_onset = REM_onset+time_correction;

hypno_duration = max([(wake_onset(end)+wake_duration(end)), (sws_onset(end)+duration_sws(end)), (REM_onset(end)+REM_duration(end))]);

% Assumption: For binary vectors index 1 = time 0-1s, index 2= time 1-2 sec, and so forth
wake_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(wake_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(sws_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(REM_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

figure;
a = subplot(1,1,1);
    plot_sleep(EEG_time, EMG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot_sleep(EEG_time, EEG_rawtrace, sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');
    
% 2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];

%% Percent time spent in states

sum_total_time= sum(wake_duration)+ sum(REM_duration)+sum(duration_sws);

% percentage
sum_NREM_perc = sum(duration_sws)/sum_total_time*100;
sum_REM_perc = sum(REM_duration)/sum_total_time*100;
sum_awake_perc = sum(wake_duration)/sum_total_time*100;


%% 6b) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(MA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    MA_binary_vector(t:t+d) = 1;
end

% remove micrarrousal from wake vectors
wake_woMA_onset = wake_onset;
wake_woMA_onset(MA_idx) = [];
wake_woMA_duration = wake_duration;
wake_woMA_duration(MA_idx) = [];
wake_woMA_binary_vector = zeros([1, hypno_duration+time_correction]); % vector of zeros matching the length of recording in seconds (+1 for last time interval).
for i=1:length(wake_woMA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];



%% 7) Aligning binary sleep score vectors and on/offsets

% Remove first seconds of EEG score to align with FP trace
wake_binary_vector_cut = wake_binary_vector(round(TTL_EEG_onset+1):end);
sws_binary_vector_cut = sws_binary_vector(round(TTL_EEG_onset+1):end);
REM_binary_vector_cut = REM_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[sws_onset_cut, sws_offset_cut] = binary_to_OnOff(sws_binary_vector_cut);
sws_duration_cut = sws_offset_cut - sws_onset_cut;

[REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
REM_duration_cut = REM_offset_cut - REM_onset_cut;


% Align period arrays according to TTL
wake_periods_cut = [wake_onset_cut wake_offset_cut];
sws_periods_cut = [sws_onset_cut sws_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];


%% 7b)  Alingment of MA vectors

% Remove first seconds of EEG score to align with FP trace
MA_binary_vector_cut = MA_binary_vector(round(TTL_EEG_onset+1):end);
wake_woMA_binary_vector_cut = wake_woMA_binary_vector(round(TTL_EEG_onset+1):end);

% Align onset, offset, and duration vectors based on TTL
[MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
MA_duration_cut = MA_offset_cut - MA_onset_cut;

[wake_woMA_onset_cut, wake_woMA_offset_cut] = binary_to_OnOff(wake_woMA_binary_vector_cut);
wake_woMA_duration_cut = wake_woMA_offset_cut - wake_woMA_onset_cut;

MA_periods_cut = [MA_onset_cut MA_offset_cut];
wake_woMA_periods_cut = [wake_woMA_onset_cut wake_woMA_offset_cut];


%% 8) Re-classify MA as NREM using boutscore_vector
%State transitions (uncut vectors)
% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, sum(hypno_duration+time_correction+6)]);

% Here using the unaligned "uncut" vectors
for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    boutscore_vector(t:t+d) = 1; % wake=1
end
for i=1:length(sws_onset)
    t = sws_onset(i)+1;
    d = duration_sws(i)-1;
    boutscore_vector(t:t+d) = 4; % sws=4
end
if ~isnan(REM_onset)
    for i=1:length(REM_onset)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        boutscore_vector(t:t+d) = 9; %REM=9
    end
end
for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    boutscore_vector(t:t+d) = 15; %MA=15
end
% re-classify MA as NREM
NREMinclMA_binary_vector = boutscore_vector==4 | boutscore_vector==15;
NREMinclMA_binary_vector_cut = NREMinclMA_binary_vector(round(TTL_EEG_onset+1):end);
[NREMinclMA_onset_cut, NREMinclMA_offset_cut] = binary_to_OnOff(NREMinclMA_binary_vector_cut);
NREMinclMA_duration_cut = NREMinclMA_offset_cut-NREMinclMA_onset_cut;
NREMinclMA_periods_cut = [NREMinclMA_onset_cut NREMinclMA_offset_cut];
%% 9) PSD on NE trace
% power spectral densities
t1 = NREMinclMA_periods_cut(:,1);
t2 = NREMinclMA_periods_cut(:,2);

tsamp1 = floor(t1*signal_fs); %eeg start time 
tsamp2 = floor(t2*signal_fs); %eeg end time
NREM_data = cell(1, numel(tsamp1));

PXX = [];
PXXlog = [];
PXX_pk_f = [];
PXX_pk = [];

NREM_data_collect = [];
period_duration = [];

for i=1:numel(tsamp1)
    period_length_i = tsamp2(i)-tsamp1(i);
    if period_length_i < 120*signal_fs % periods shorter than 120 s are excluded from analysis
        continue
    end
    period_duration = [period_duration period_length_i/signal_fs];
    NREM_data{i} = delta465_filt(tsamp1(i):tsamp2(i));
    timetrace_i = sec_signal(tsamp1(i):tsamp2(i));
    
    %detrend (and center around 0)
    [p,s,mu] = polyfit((1:numel(NREM_data{i}))',NREM_data{i},5);
    f_y = polyval(p,(1:numel(NREM_data{i}))',[],mu);
    detrend_data = NREM_data{i} - f_y';        % Detrend data
    
    [pxx, f] = pwelch(detrend_data, [], [],[0:0.002:0.1], signal_fs); %
    logpxx = 10*log10(pxx);
    FX{i} = f;
    [pxx_pk_psd, max_idx] = max(pxx);
    PXX_pk = [PXX_pk pxx_pk_psd];
    pxx_pk_f = f(max_idx);
    PXX_pk_f = [PXX_pk_f pxx_pk_f];
    PXX = [PXX pxx'];

    figure
    set(gcf, 'Position',  [100, 300, 1500, 250])
    a = subplot(1,2,1);
        plot(timetrace_i,NREM_data{i});
        hold on
        plot(timetrace_i,detrend_data);
        legend({'raw','fitted'})
        hold off
    b = subplot(1,2,2);
        plot(f,pxx);
end

weighted_mean_PXX_NE = sum(period_duration.*PXX,2)/sum(period_duration); %period duration is used as weights

[pxx_pk_power, pxx_pk_idx] = max(weighted_mean_PXX_NE); % mean peak power
pxx_pk_f = f(pxx_pk_idx); % mean peak frequency

prism_psd_NE = weighted_mean_PXX_NE;
prism_freq_NE = f;

% power spectral density plot
figure
    plot(prism_freq_NE,prism_psd_NE)
    xlabel('frequency (Hz)');
    ylabel('PSD');

%% 10) NE osciallation rate (method: findpeaks, number per time spent in NREM and inter-peak interval)


%save generated plots
path_save= 'C:\Users\username\data\MATLAB\24hrec_peaks';
my_folder= path_save;
folder = mkdir(path_save, my_folder);
path  = [path_save,filesep,my_folder];

pks_n = [];
pkwdths_n = [];
pkproms_n = [];
period_duration = [];
interpeak_interval = [];

min_pkDist = 10; % use 10
min_pkProm = 0.5; % decide for each animal 
min_pkWdth = 8; % use 8.

if NREMinclMA_periods_cut(1,1) == 0 % if first period starts at time=0, set to 1 in order to index further down
    NREMinclMA_periods_cut(1,1) = 1;
end

for i = 1:length(NREMinclMA_periods_cut)% -1% exclude last bout if it is outside NE recording
    period_duration_i = NREMinclMA_periods_cut(i,2)-NREMinclMA_periods_cut(i,1);
    if period_duration_i < 60 % periods shorter than 60 s are excluded from analysis
        continue
    end
    period_duration = [period_duration period_duration_i];
    period_i = NREMinclMA_periods_cut(i,:)*signal_fs;
    NEtrace_i = delta465_filt(round((period_i(1):period_i(2))));    
    timetrace_i = sec_signal(round((period_i(1):period_i(2))));
    [pks, pklocs, w, p] = findpeaks(NEtrace_i, timetrace_i, 'MinPeakDistance', min_pkDist, 'MinPeakWidth', min_pkWdth, 'MinPeakProminence',min_pkProm); 
    %[pks, pklocs, w, p] = findpeaks(NEtrace_i, timetrace_i, 'MinPeakDistance', min_pkDist, 'MinPeakProminence',min_pkProm); 
    pkwdths_n = [pkwdths_n w];
    pkproms_n = [pkproms_n p];
    pks_n = [pks_n length(pks)];
    interpeak_interval = [interpeak_interval diff(pklocs)];
    figure
    set(gcf, 'Position',  [100, 300, 1500, 250])
    findpeaks(NEtrace_i, timetrace_i, 'MinPeakDistance', min_pkDist,'MinPeakWidth', min_pkWdth, 'MinPeakProminence',min_pkProm, 'Annotate','extents'); % here you check that the detection settings are alright. If not, you adjust the settings in findpeaks (both for analysis and figure) – distance and prominence.
end

total_period_duration = sum(period_duration); % cummulative duration of NREM (s)
total_pks_n = sum(pks_n); % cummulative number of peaks (estimate of NE events)
oscillation_frq_per_NREM_s = total_pks_n/total_period_duration;
oscillation_interpk_intvl = mean(interpeak_interval);


%% calculate NE amplitude (90-10 percentile of each NREM bout)

% select sws periods
period = sws_periods_cut; % NB! select how many periods to inlcude in analysis
if period(end)>length(delta465_filt)/signal_fs || period(end)>length(EEG_rawtrace_cut)/sampling_freq
    period = period(1:end-1,:);
end
period_fs = period*signal_fs;

period_ampl = [];
trace_NE_stitch = [];
for i = 1:length(period_fs)
   period_i = period_fs(i,:);
   traceNE = delta465_filt(period_i(1):period_i(2));
   perc90 = prctile(traceNE,90);
   perc10 = prctile(traceNE,10);
   perc_amplitude = perc90-perc10;
   period_ampl = [period_ampl perc_amplitude];
   trace_NE_stitch = [trace_NE_stitch traceNE];

end

% weighted average based on period duration (short NREM periods weigh less than long ones)
period_dur = period(:,2)-period(:,1);
prism_mean_amplitude_weighted = sum(period_ampl*period_dur)/sum(period_dur);

%% 8) Plotting all traces and scorings together
% plot EEG sleep scoring bouts together with FP data

% Time vector for sleep scoring (1 Hz)
sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % should be same length for wake/sws/REM

fig = figure;
a = subplot(3,1,1);
    plot_sleep(ds_sec_signal, ds_delta465_filt, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    title('NE');
b = subplot(3,1,2);
    ds_EEG_time = downsample(EEG_time_cut, 10);
    ds_EMG_rawtrace = downsample(EMG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EMG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EMG (V)');
c = subplot(3,1,3);
    ds_EEG_rawtrace = downsample(EEG_rawtrace_cut, 10);
    plot_sleep(ds_EEG_time, ds_EEG_rawtrace, sleepscore_time_cut, wake_binary_vector_cut, sws_binary_vector_cut, REM_binary_vector_cut, MA_binary_vector_cut);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b,c],'x');
