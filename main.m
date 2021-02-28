% 
% Main Function: loading Simband data for PPG peak detection
% Author: Dong Han (dong.han@uconn.edu)
% Create Date: 05/14/2020
%
% Please cite our paper if you use our data:
% A PPG Peak Detection Method for Accurate Determination of Heart Rate 
% during Sinus Rhythm and Cardiac Arrhythmia 
% 
% Last modification: Dong Han, 02/28/2021.
%


clear all;
close all;
clc;

fs = 128; % Hz, Simband original sampling frequency. Reference ECG was already downsampled from 180 Hz to 128Hz.
fs_ACC = 30; % Hz, downsampled ACC frequency.
fs_PPG = 50; % Hz, downsampled PPG frequency.

%% User Interact: Waiting for Command Window User input: subject & window
output_struct = my_step_01_know_sub_seg();
you_want_to_plot_sig = input('Do you want to plot the signal traces? y/Y = Yes, otherwise = no. \n','s');
if strcmp(you_want_to_plot_sig,'y') || strcmp(you_want_to_plot_sig,'Y')
    disp('--- Plotting signal traces. ---');
    you_want_to_plot_sig_flag = true;
else
    disp('--- Not plotting signal traces. ---');
    you_want_to_plot_sig_flag = false;
end
%% Start to run the main func
aaa = 1; % Inital aaa, must keep for func `my_step_02_load_data`.
output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,output_struct);
load_Simband_subject_name = output_data.load_Simband_subject_name;
for aaa = 1:size(load_Simband_subject_name,1)
    fprintf('--- Start to run on subject %d ---\n',load_Simband_subject_name(aaa));
    if aaa > 1
        % If we want to run more than one subject:
        output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,output_struct);
    end
    countlen = output_data.countlen; % The sample length of this subject's data.
    while output_data.iiii_PPG_end <= countlen
        %% 30-sec Segment of: 
        % PPG_buffer - Simband PPG, 
        % Simband_ECG_buffer - Simband ECG (not used in this study, only
        %                      provided here),
        % ACC_buffer - Simband ACC,
        % Ref_ECG_buffer - Holter monitor reference ECG.
        [PPG_buffer,Simband_ECG_buffer,ACC_buffer,Ref_ECG_buffer] = my_step_03_prepare_buffer(output_data);

        %% ECG Beat Annotation: if not corrected manually, 
        %      then the ECG beats from Pan-Tompkins' algorithm is correct.
        [refECG_pkloc,wbwrefHR] = my_ECG_Peak_Detection_concise(output_data,fs,Ref_ECG_buffer);
        
        %% Cardiac Arrhythmia Adjudication:
        % 0 - NSR:
        % 1 - AF;
        % 2 - PAC/PVC;
        % 3 - not sure if it is PAC/PVC or NSR. We considered them as NSR.
        % 5 - PPG noisy;
        % NaN - not enough reference ECG.
        disease_label = output_data.disease_label;
        ii = floor(output_data.iiii_PPG_end / (30 * fs_PPG)); % current segment index.
        ground_truth_of_this_seg = disease_label(ii,2); % the adjudicated ground truth of current segment.
        
        %% Display Info and Plot this 30-sec Segment
        fprintf('Loaded Subject %s, Segment %03d to Work Space!\n',output_data.Simband_Subject,ii);
        
        if you_want_to_plot_sig_flag
            switch ground_truth_of_this_seg
                case 0
                    plot_GT = 'NSR';
                case 1
                    plot_GT = 'AF';
                case 2
                    plot_GT = 'PAC/PVC';
                case 3
                    plot_GT = 'not sure if it is PAC/PVC or NSR. We considered them as NSR.';
                case 5
                    plot_GT = 'PPG noisy';
                otherwise
                    if isnan(ground_truth_of_this_seg)
                        plot_GT = 'not enough reference ECG';
                    else
                        disp('Unknown ground truth, please contact author (dong.han@uconn.edu) for this error!');
                        keyboard;
                    end
            end

            figure;
            sgtitle({['Subject ',output_data.Simband_Subject,', Segment ',num2str(ii,'%03d')],['Ground Truth: ',plot_GT]}); % Title of the figure.
            ax(1) = subplot(4,1,1);
            t_ECG = (output_data.iiii_ECG_start:output_data.iiii_ECG_end)./fs;
            plot(t_ECG,Ref_ECG_buffer);
            hold on;
            plot(t_ECG(refECG_pkloc),Ref_ECG_buffer(refECG_pkloc),'r.');
            ylabel('Reference ECG (a.u.)');

            ax(2) = subplot(4,1,2);
            t_PPG = (output_data.iiii_PPG_start:output_data.iiii_PPG_end)./fs_PPG;
            plot(t_PPG,PPG_buffer);
            ylabel('Simband PPG (a.u.)');

            ax(3) = subplot(4,1,3);
            plot(t_ECG(refECG_pkloc(2:end)),wbwrefHR,'r.-');
            ylabel('Reference HR (BPM)');

            ax(4) = subplot(4,1,4);
            t_ACC = (output_data.iiii_ACC_start:output_data.iiii_ACC_end)./fs_ACC;
            plot(t_ACC,ACC_buffer);
            ylabel('Simband ACC (a.u.)');

            linkaxes(ax,'x');
            xlim([t_PPG(1) t_PPG(end)]);
            xlabel('Time (sec)');

            close all;
        end
        %% Update 30-sec Buffer Pointer for Next Segment:
        output_data.iiii_PPG_start = output_data.iiii_PPG_end + 1;
        output_data.iiii_PPG_end = output_data.iiii_PPG_start + 30 * fs_PPG - 1;
        output_data.iiii_ECG_start = output_data.iiii_ECG_end + 1;
        output_data.iiii_ECG_end = output_data.iiii_ECG_start + 30 * fs - 1;
        output_data.iiii_ACC_start = output_data.iiii_ACC_end + 1;
        output_data.iiii_ACC_end = output_data.iiii_ACC_start + 30 * fs_ACC - 1;
    end
    fprintf('--- Finished subject %d, closing all figures now ---\n',load_Simband_subject_name(aaa));
    close all;
end
