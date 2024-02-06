%% 1) Define mouse data

clear all
close all

% data structure:
    % 1) EEG raw data
    % 2) EEG sleep score
    % 3) EEG sleepscore time correction

example_mouseID  = {'C:\Users\username\data\EEG_data.exp' 'C:\Users\username\data\sleep_score.xlsx' 0};

mouse = example_mouseID; 
%% 2) loading and plotting EEG and EMG raw data

% Make sure the "ExpToolbox" is added to the matlab path - can be obtained from Viewpoint ( https://www.viewpoint.fr/)

% Import EEG raw data to matlab
Info=loadEXP([mouse{1}],'no');

TimeReldebSec=0; %start extract data from the beginning (first bin)
TimeRelEndSec=inf; %inf to include all data (until last bin)

[Data,Time]=ExtractContinuousData([],Info,[],TimeReldebSec, TimeRelEndSec,[],1);

EMG_rawtrace = Data(1,1:end);
EEG_rawtrace = Data(2,1:end);

%time vector using sampling frequency
sampling_freq = Info.Fs;
EEG_time = (1:length(EEG_rawtrace))/sampling_freq;
Duration = EEG_time(end)+100;

% Plot of EEG and EMG traces
figure;
h(1) = subplot(2,1,1);
    plot(EEG_time, EMG_rawtrace); 
    xlabel('time (s)');
    ylabel('EMG (V)');
h(2) = subplot(2,1,2);
    plot(EEG_time, EEG_rawtrace); 
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([h(1),h(2)],'x');

%% 3) Interpolate gaps in EEG/EMG 
% used in case of missing samples in the EEG recording at backup time

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

if isnan(EEG_rawtrace(end))
    EEG_rawtrace = EEG_rawtrace(1:end-1);
    EMG_rawtrace = EMG_rawtrace(1:end-1);
    EEG_time = EEG_time(1:end-1);
end


%% 4) open EEG scoring

time_correction = mouse{3};
EEG_sleepscore = xlsread(mouse{2});

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

% NB! all EEG/EMG traces are not aligned properly with sleep score (~4 s delay)
wake_onset = wake_onset+time_correction;
sws_onset = sws_onset+time_correction;
REM_onset = REM_onset+time_correction;


% Assumption: For binary vectors index 1 = time 0-1s, index 2= time 1-2 sec, and so forth
wake_binary_vector = zeros([1, round(Duration)]); % vector of zeros matching the length of recording in seconds (+1 for last time interval). Sum is used because of instances where exp was backed up at 18 has two durations
for i=1:length(wake_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_onset(i)+1; % +1 to put time 0 as index 1
    d = wake_duration(i)-1; % -1 compensates for adding 1
    wake_binary_vector(t:t+d) = 1;
end

sws_binary_vector = zeros([1, round(Duration)]); % vector of zeros matching the length of recording in seconds  (+1 for last time interval)
for i=1:length(sws_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = sws_onset(i)+1; 
    d = duration_sws(i)-1;
    sws_binary_vector(t:t+d) = 1;
end

REM_binary_vector = zeros([1, round(Duration)]); % vector of zeros matching the length of recording in seconds (+1 for last time interval)
for i=1:length(REM_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = REM_onset(i)+1;
    d = REM_duration(i)-1;
    REM_binary_vector(t:t+d) = 1;
end

% Time vector for sleep scoring (1 Hz)
sleepscore_time = 0:length(wake_binary_vector)-1; % Should be same length for wake/sws/REM binary vectors

% 2-column vectors with on- and offsets for each state
wake_periods = [wake_onset wake_onset+wake_duration];
sws_periods = [sws_onset sws_onset+duration_sws];
REM_periods = [REM_onset REM_onset+REM_duration];

fig = figure;
a = subplot(2,1,1);
    plot_sleep(downsample(EEG_time, 10), downsample(EMG_rawtrace,10), sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EMG (V)');
b = subplot(2,1,2);
    plot_sleep(downsample(EEG_time, 10), downsample(EEG_rawtrace,10), sleepscore_time, wake_binary_vector, sws_binary_vector, REM_binary_vector);
    xlabel('time (s)');
    ylabel('EEG (V)');
linkaxes([a,b],'x');
   


%% 5) Dividing wake bouts into microarousals (MA) and wake w/o MA

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector = zeros([1, round(Duration)]);
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
wake_woMA_binary_vector = zeros([1, round(Duration)]);
for i=1:length(wake_woMA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

wake_binary_vector = wake_woMA_binary_vector;

% 2-column vectors with on- and offsets for each state
MA_periods = [MA_onset MA_onset+MA_duration];
wake_woMA_periods = [wake_woMA_onset wake_woMA_onset+wake_woMA_duration];

%% 6) Re-classify MA as NREM using boutscore_vector
% Here you can pool MAs with NREM sleep which can be beneficial for some
% analyses related to infraslow oscillations (eg. PSD analysis), where you
% don't want to divide your traces into short/pure NREM bouts

%State transitions (uncut vectors)
% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, round(Duration)]);

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


%% 7) make standard hypnogram for plotting

collect_ = NaN(1, length(wake_binary_vector));
for i = 1:length(wake_binary_vector)
    if wake_binary_vector(i) == 1
        collect_(i) = 3;
    end

     if sws_binary_vector(i) == 1
        collect_(i) = 1;
     end

     if REM_binary_vector(i) == 1
        collect_(i) = 2;
     end
     if MA_binary_vector(i) == 1
        collect_(i) = 4;
     end
end

 figure
a = subplot(3,1,1);
       plot(EEG_time, EEG_rawtrace);
b = subplot(3,1,2);
       plot(EEG_time, EMG_rawtrace);
c = subplot(3,1,3);
       plot(sleepscore_time, collect_);
linkaxes([a,b,c],'x');


%% 8) only for APP

sleep_scoring_interval_start = sws_onset(1); % sec
sleep_scoring_interval = sleep_scoring_interval_start + 3600; %one hour in sec

wake_binary_vector_cut = wake_woMA_binary_vector(sleep_scoring_interval_start:sleep_scoring_interval);
REM_binary_vector_cut = REM_binary_vector(sleep_scoring_interval_start:sleep_scoring_interval);
sws_binary_vector_cut = sws_binary_vector(sleep_scoring_interval_start:sleep_scoring_interval);
MA_binary_vector_cut = MA_binary_vector(sleep_scoring_interval_start:sleep_scoring_interval);
sleepscore_time_cut = 0:length(wake_binary_vector_cut)-1; % Should be same length for wake/sws/REM binary vectors
NREMinclMA_binary_vector_cut = NREMinclMA_binary_vector(sleep_scoring_interval_start:sleep_scoring_interval);
collect_cut = collect_(sleep_scoring_interval_start:sleep_scoring_interval);

EMG_rawtrace_cut = EMG_rawtrace(sleep_scoring_interval_start*sampling_freq:sleep_scoring_interval*sampling_freq);
EEG_rawtrace_cut = EEG_rawtrace(sleep_scoring_interval_start*sampling_freq:sleep_scoring_interval*sampling_freq);
EEG_time_cut = (1:length(EEG_rawtrace_cut))/sampling_freq;

 figure
a = subplot(3,1,1);
       plot(EEG_time_cut, EEG_rawtrace_cut);
b = subplot(3,1,2);
       plot(EEG_time_cut, EMG_rawtrace_cut);
c = subplot(3,1,3);
       plot(sleepscore_time_cut, collect_cut);
linkaxes([a,b,c],'x');



%% 9) EEG power spectrum analysis for normal trace

 Data_EEG = EEG_rawtrace_cut; 

analysis_window = 5; %sec. 1 for 30 sec

power_bands = {[0.5, 1.75], [1.75 4] [4, 8], [8, 15], [15, 30] [30, 100]}; % delta was 0.2 before
total_power_band = [0, 100];
frw = 0:0.2:100;


frq = sampling_freq;
  
%% 10) Align scored states to FP recording

% Align onset, offset, and duration vectors based on TTL
[wake_onset_cut, wake_offset_cut] = binary_to_OnOff(wake_binary_vector_cut);
wake_duration_cut = wake_offset_cut - wake_onset_cut;

[sws_onset_cut, sws_offset_cut] = binary_to_OnOff(sws_binary_vector_cut);
sws_duration_cut = sws_offset_cut - sws_onset_cut;

[REM_onset_cut, REM_offset_cut] = binary_to_OnOff(REM_binary_vector_cut);
REM_duration_cut = REM_offset_cut - REM_onset_cut;

[MA_onset_cut, MA_offset_cut] = binary_to_OnOff(MA_binary_vector_cut);
MA_duration_cut = MA_offset_cut - MA_onset_cut;

[NREMinclMA_onset, NREMinclMA_offset] = binary_to_OnOff(NREMinclMA_binary_vector_cut);
NREMinclMA_duration = NREMinclMA_offset-NREMinclMA_onset;


% Align period arrays according to TTL
wake_periods_cut = [wake_onset_cut wake_offset_cut];
sws_periods_cut = [sws_onset_cut sws_offset_cut];
REM_periods_cut = [REM_onset_cut REM_offset_cut];
MA_periods_cut = [MA_onset_cut MA_offset_cut];
NREMinclMA_periods = [NREMinclMA_onset NREMinclMA_offset];
   
%% 11) State transitions

% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1, round(Duration)]);

% Here using the aligned "cut" vectors
if isnan(wake_onset_cut)~=1
    for i=1:length(wake_onset_cut)
        t = wake_onset_cut(i)+1;
        d = wake_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 1; % wake=1
    end
end


if isnan(sws_onset_cut)~=1
    for i=1:length(sws_onset_cut)
        t = sws_onset_cut(i)+1;
        d = sws_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 4; % sws=4
    end
end

if isnan(REM_onset_cut)~=1
    for i=1:length(REM_onset_cut)
        t = REM_onset_cut(i)+1;
        d = REM_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 9; %REM=9
    end
end

if isnan(MA_onset_cut)~=1
    for i=1:length(MA_onset_cut)
        t = MA_onset_cut(i)+1;
        d = MA_duration_cut(i)-1;
        boutscore_vector(t:t+d) = 15; %MA=15
    end
end

% Vectors indicate time of transitions in seconds
transition_sws_wake =  find(diff(boutscore_vector)== -3);
transition_REM_wake =  find(diff(boutscore_vector)== -8);
transition_sws_MA =  find(diff(boutscore_vector)== 11);
transition_REM_sws =  find(diff(boutscore_vector)== -5);
transition_sws_REM =  find(diff(boutscore_vector)== 5);
transition_REM_MA =  find(diff(boutscore_vector)== 6);



%% 12) Get the spectogram our for different sleep phases

t1 = sws_onset_cut;
%t1 = REM_onset_cut;
%t1 = wake_onset_cut;
t2 = sws_offset_cut;
%t2 = REM_offset_cut;
%t2 = wake_offset_cut;

tsamp1 = floor(t1*sampling_freq); %eeg start time 
tsamp2 = floor(t2*sampling_freq); %eeg end time
NREM_data = cell(1, numel(tsamp1));

PXX = [];
NREM_data_collect = [];

for i=1:numel(tsamp1)
    if tsamp1(i) == 0
        tsamp1(i) = 1;
    else
    end
    NREM_data{i} = EEG_rawtrace_cut(tsamp1(i):tsamp2(i));
    NREM_data_cut = EEG_rawtrace_cut(tsamp1(i):tsamp2(i));
    NREM_data_collect = [NREM_data_collect NREM_data_cut];
    [pxx, f] = pwelch(NREM_data{i}, [], [],[0:0.2:100], sampling_freq);
    logpxx = 10*log10(pxx);
    FX{i} = f;
    PXX(:,i) = logpxx;
end

mean_PXX = mean(PXX,2);
f = f';

mean_delta_slow_power_density = mean(mean_PXX(f>0.5 & f<1.75));
mean_delta_fast_power_density = mean(mean_PXX(f>1.75 & f<4));
mean_delta_power_density = mean(mean_PXX(f>1 & f<4));
mean_sigma_power_density = mean(mean_PXX(f>8 & f<15));
mean_theta_power_density = mean(mean_PXX(f>4 & f<8));
mean_beta_power_density = mean(mean_PXX(f>15 & f<30));

prism_band_collect = [mean_delta_slow_power_density mean_delta_fast_power_density mean_delta_power_density mean_theta_power_density mean_sigma_power_density mean_beta_power_density]';

