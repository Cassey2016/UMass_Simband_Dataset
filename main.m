% 
% Main Function: loading Simband data for PPG peak detection
% Author: Dong Han (dong.han@uconn.edu)
% Date: 05/14/2020
%
% Please cite our paper if you use our data:
% A PPG Peak Detection Method for Accurate Determination of Heart Rate 
% during Sinus Rhythm and Cardiac Arrhythmia 
% 
% Last modification: Dong Han, 05/14/2020.
%


clear all;
close all;
clc;

fs = 128; % Hz, Simband original sampling frequency. Reference ECG was already downsampled from 180 Hz to 128Hz.
fs_ACC = 30; % Hz, downsampled ACC frequency.
fs_PPG = 50; % Hz, downsampled PPG frequency.

aaa = 1; % this is the beginning of main function.
output_struct = my_step_01_know_sub_seg();
output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,output_struct);
load_Simband_subject_name = output_data.load_Simband_subject_name;
for aaa = 1:size(load_Simband_subject_name,1)
    if aaa > 1
        % if we want to run more than one subject.
        output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,output_struct);
    end
    countlen = output_data.countlen;
    while output_data.iiii_PPG_end <= countlen
        %% 30-sec segment for 
        % PPG_buffer - Simband PPG, 
        % Simband_ECG_buffer - Simband ECG (not used in this study, only
        %                      provided here),
        % ACC_buffer - Simband ACC,
        % Ref_ECG_buffer - Holter monitor reference ECG.
        [PPG_buffer,Simband_ECG_buffer,ACC_buffer,Ref_ECG_buffer] = my_step_03_prepare_buffer(output_data);

        %% ECG beat annotation: if not labeled manually, 
        %      then the ECG beats from Pan-Tompkins' algorithm is correct.
        [refECG_pkloc,wbwrefHR] = my_ECG_Peak_Detection_concise(output_data,fs,Ref_ECG_buffer);
        
        %% Cardiac Arrhythmia Adjudication:
        disease_label = output_data.disease_label;
        ii = floor(output_data.iiii_PPG_end / (30 * fs_PPG)); % current segment index.
        ground_truth_of_this_seg = disease_label(ii,2); % the adjudicated ground truth of current segment.
        
        % update buffer pointing index.
        output_data.iiii_PPG_start = output_data.iiii_PPG_end + 1;
        output_data.iiii_PPG_end = output_data.iiii_PPG_start + 30 * fs_PPG - 1;
        output_data.iiii_ECG_start = output_data.iiii_ECG_end + 1;
        output_data.iiii_ECG_end = output_data.iiii_ECG_start + 30 * fs - 1;
        output_data.iiii_ACC_start = output_data.iiii_ACC_end + 1;
        output_data.iiii_ACC_end = output_data.iiii_ACC_start + 30 * fs_ACC - 1;
    end
end
