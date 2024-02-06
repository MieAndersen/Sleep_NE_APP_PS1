# Sleep_NE_APP_PS1

Code to reproduce figures in paper by Cankar et al. (2024).

The script "sleep_APP_PS1" was used to analyze 24 h EEG/EMG recordings, while "APP_interval_spectral" was used to analyze combined EEG and fiber photometry recordings.

You can find the customary functions binary_to_OnOff and plot_sleep [here](https://github.com/MieAndersen/NE-oscillations/tree/main/functions). Additionally, to load fiber photometry data we use a custom function provided by Tucker Davis Technologies called "TDTbin2mat" that you can find in their [repository](https://github.com/tdtneuro/TDTMatlabSDK). To load EEG and EMG data obtained using Sleepscore (ViewPoint, France) we use a custom function provided by ViewPoint Behavior Technology called "loadEXP" - the toolbox needed for this is called "ExpToolbox" and can be obtained upon request (www.viewpoint.fr).
