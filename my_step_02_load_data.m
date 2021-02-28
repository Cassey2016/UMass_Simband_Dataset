% 
% function output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,...
%                                             output_struct)
% 
% Author: Dong Han (dong.han@uconn.edu)
% Date: 05/14/2020
% 
% Description: load Simband and reference ECG data for selected subjects.
%
% Return: PPG - 0: all subjects; 1-41: specific subject.
%         Simband_ECG - the subject name(s) user selected.
%         ACC - the data structure of Simband.
%         aligned_ChestECG - 0: all segments; 1-unknown: segment index.
%         new_ECG_peak_loc - starting segment index.
%         new_ECG - ending segment index.
%         seg_num_ECG - if user input running 'all segments.'
%         seg_num_peak - the path where Simband data stored.
%         disease_label - ground truth adjudicated by UConn and UMass
%                         group.
%         my_Simband_datafilename - current Simband subject name.
%         countseg - maximum total number of 30-sec seg for this subject.
%         countlen - the smallest recording sample size among PPG & ACC.
%         iiii_PPG_start - sample index of PPG starting point (30-sec).
%         iiii_PPG_end - sample index of PPG ending point (30-sec).
%         iiii_ECG_start - sample index of ECG starting point (30-sec).
%         iiii_ECG_end - sample index of ECG ending point (30-sec).
%         iiii_ACC_start - sample index of ACC starting point (30-sec).
%         iiii_ACC_end - sample index of ACC ending point (30-sec).
%         count_ii - count how many 30-sec has finished.
%         load_Simband_subject_struct - the data structure of Simband for
%                                       all subjects that user selected.
%         aaa - current subject index among user selected subjects.
%         Simband_Subject - same as 'my_Simband_datafilename'.
%         load_Simband_subject_name - all subjects user selected to run.
%         subject_win_idx - the interval array from 1 to the countseg.
%
% Last Modification: Dong, 02/28/2021.
function output_data = my_step_02_load_data(fs,fs_ACC,fs_PPG,aaa,...
                                             output_struct)
                                         
    user_input_subject = output_struct.user_input_subject;
    load_Simband_subject_name = output_struct.load_Simband_subject_name;
    load_Simband_subject_struct = output_struct.load_Simband_subject_struct;
    user_input_win = output_struct.user_input_win;
    start_win_idx = output_struct.start_win_idx;
    end_win_idx = output_struct.end_win_idx;
    all_win_flag = output_struct.all_win_flag;
    
    Simband_data_folder = output_struct.Simband_data_folder;

    [data,aligned_ChestECG,disease_label,my_Simband_datafilename] = my_load_data(load_Simband_subject_name,aaa,Simband_data_folder);


    % subfunction 3
    [PPG,Simband_ECG,ACC,aligned_ChestECG] = my_load_Simband_data(data,aligned_ChestECG,load_Simband_subject_struct(aaa),fs,fs_ACC,fs_PPG,my_Simband_datafilename);
    % subfunction 4
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_load_ref_ECG_peak(my_Simband_datafilename,Simband_data_folder,fs);

    shift_sec = 0;
    % subfunction 5
    [countseg,countlen] = my_count_total_seg(shift_sec,length(PPG),length(Simband_ECG),length(ACC),fs,fs_ACC,fs_PPG);

    % subfunction 6
    [start_win_idx,end_win_idx,subject_win_idx] = my_set_start_end_win_idx(user_input_subject,user_input_win,all_win_flag,aaa,countseg,start_win_idx,end_win_idx);
    iiii_PPG_start = (start_win_idx-1) * 30 * fs_PPG + 1;
    iiii_PPG_end = iiii_PPG_start + 30 * fs_PPG - 1;
    iiii_ECG_start = (start_win_idx-1) * 30 * fs + 1;
    iiii_ECG_end = iiii_ECG_start + 30 * fs - 1;
    iiii_ACC_start = (start_win_idx-1) * 30 * fs_ACC + 1;
    iiii_ACC_end = iiii_ACC_start + 30 * fs_ACC - 1;
    count_ii = 0;

    output_data = struct('PPG',PPG,...
        'Simband_ECG',Simband_ECG,...
        'ACC',ACC,...
        'aligned_ChestECG',aligned_ChestECG,...
        'new_ECG_peak_loc',new_ECG_peak_loc,...
        'new_ECG',new_ECG,...
        'seg_num_ECG',seg_num_ECG,...
        'seg_num_peak',seg_num_peak,...
        'disease_label',disease_label,...
        'my_Simband_datafilename',my_Simband_datafilename,...
        'countseg',countseg,...
        'countlen',countlen,...
        'iiii_PPG_start',iiii_PPG_start,...
        'iiii_PPG_end',iiii_PPG_end,...
        'iiii_ECG_start',iiii_ECG_start,...
        'iiii_ECG_end',iiii_ECG_end,...
        'iiii_ACC_start',iiii_ACC_start,...
        'iiii_ACC_end',iiii_ACC_end,...
        'count_ii',count_ii,...
        'load_Simband_subject_struct',load_Simband_subject_struct,...
        'aaa',aaa,...
        'Simband_Subject',my_Simband_datafilename,...
        'load_Simband_subject_name',load_Simband_subject_name,...
        'subject_win_idx',subject_win_idx);

    
end

function [data,aligned_ChestECG,disease_label,my_Simband_datafilename] = my_load_data(load_Simband_subject_name,aaa,Simband_data_folder)


    my_Simband_datafilename = num2str(load_Simband_subject_name(aaa));


    myfilename_new = [Simband_data_folder,'\',my_Simband_datafilename];
    load(myfilename_new); % loaded name is 'data'

    myfilename = [myfilename_new,'RefECG'];

    load(myfilename); % variable name will be: aligned_ChestECG

    myfilename = [myfilename_new,'_ground_truth']; % % this is the ground truth of AF, NSR and PAC/PVC. it is segment-wise. 0-NSR, 1-AF, 2-PAC/PVC, 3-not sure if is PAC/PVC or NSR, 5-noisy, NaN-not enough reference ECG.
    load(myfilename); % loaded name: disease_label

    
    if ((strcmp(my_Simband_datafilename,'4019') || (strcmp(my_Simband_datafilename,'4026')) || strcmp(my_Simband_datafilename,'4030')) || (strcmp(my_Simband_datafilename,'4034')))
        aligned_ChestECG = Aligned_ECG_Ref_ch2;
    elseif (strcmp(my_Simband_datafilename,'4021') || strcmp(my_Simband_datafilename,'4031'))
        aligned_ChestECG = Aligned_ECG_Ref_ch3;
    elseif (strcmp(my_Simband_datafilename,'4018'))
    %     aligned_ChestECG; % 4018 do not have Holter ECG data.
    else
        aligned_ChestECG = Aligned_ECG_Ref_ch1;
    end
end

function [PPG,Simband_ECG,ACC,aligned_ChestECG] = my_load_Simband_data(data,aligned_ChestECG,myoldstruct,fs,fs_ACC,fs_PPG,my_Simband_datafilename)
%---Introduction---
% This function just split the entire Simband data that we need to use into
% 30 seconds chunk.
% -Input:
%   --data: It is the Simband data we downloaded from ARTIK cloud. It is a
%   struct data type.
%   --aligned_ChestECG: It is the Reference ECG we collected from our chest
%   band.
%   --myoldstruct: Indicating if the input data has reference ECG, or which
%   kind of data structure it is:
%       ---'0': it is latest data structure of Simband data.
%       ---'1': it is the version of data structure downloaded from UMass
%       Hospital toolbox.
% -Output:
%   --PPG: just original PPG with its original length and size.
%   --PPG_30sec: it is a matrix, with each row represent each 30 second PPG data.
%   --PPG_HR_30sec: initial a matrix to store each 30 second PPG HR with
%   the same generating method in PPG_30sec.
%   --ECG_30sec: it is a matrix with each row represent each 30 second PPG data.
%   --ECG_HR_30sec: initial a matrix to store each 30 second PPG HR with
%   the same generating method in ECG_30sec.
%   --Simband_ECG_30sec: same kind of matrix storing 30 second Simband ECG
%   data.
%   --ACC_30sec: same kind of matrix storing 30 seconds Simband ACC data.
%   --estimated_HR_30sec:  initial a matrix to store each 30 second estimated HR 
%   for usable portion in the function 'usabilityQI_estimatedHR'.
%   --usability_30sec: initial a matrix to store the usability detection
%   results in the function 'usabilityQI_estimatedHR'.
if myoldstruct == 4
    [PPG,Simband_ECG,ACC] = myfillgap(data,fs);
    PPG_down = resample(PPG,fs_PPG,fs);
    PPG = PPG_down;
    ACC_down = resample(ACC,fs_ACC,fs);
    ACC = ACC_down;
elseif myoldstruct == 5
    PPG = data.physiosignal.ppg.e.visual.signal;
    PPG_down = resample(PPG,fs_PPG,fs);
    PPG = PPG_down;
    Simband_ECG = data.physiosignal.ecg.visual.signal;
    ACC = data.accel.magnitude.signal;
    ACC_down = resample(ACC,fs_ACC,fs);
    ACC = ACC_down;
elseif myoldstruct == 6
    % only for 4018
    PPG = data.band.physiosignal.ppg4.visual.signal;
    PPG_down = resample(PPG,fs_PPG,fs);
    PPG = PPG_down;
    Simband_ECG = data.band.physiosignal.ecg.visual.signal;
    ACC = data.band.accelerometer.magnitude.signal;
    ACC_down = resample(ACC,fs_ACC,fs);
    ACC = ACC_down;
else
    [PPG,Simband_ECG,ACC] = myfillgap(data,fs);
    PPG_down = resample(PPG,fs_PPG,fs);
    PPG = PPG_down;
    ACC_down = resample(ACC,fs_ACC,fs);
    ACC = ACC_down;
end

if strcmp(my_Simband_datafilename,'4006')
    sub_point = 596;
    Simband_ECG = Simband_ECG(sub_point:end);
    aligned_ChestECG = aligned_ChestECG(sub_point:end);
    sub_point_ACC = 72;
    ACC = ACC(sub_point_ACC:end);
elseif strcmp(my_Simband_datafilename,'4007')
    sub_point = 661;
    Simband_ECG = Simband_ECG(sub_point:end);
    aligned_ChestECG = aligned_ChestECG(sub_point:end);
        sub_point_ACC = 125;
    ACC = ACC(sub_point_ACC:end);
elseif strcmp(my_Simband_datafilename,'4009')
    sub_point = 34;
    Simband_ECG = Simband_ECG(sub_point:end);
    aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4008')
    add_point_ACC = 28;
    ACC = [zeros(1,add_point_ACC),ACC];
elseif strcmp(my_Simband_datafilename,'4015')
    add_point = 44;
    Simband_ECG = [zeros(1,add_point),Simband_ECG];
    aligned_ChestECG = [zeros(add_point,1);aligned_ChestECG];
elseif strcmp(my_Simband_datafilename,'4016')
    add_point = 51;
    Simband_ECG = [zeros(1,add_point),Simband_ECG];
    aligned_ChestECG = [zeros(add_point,1);aligned_ChestECG];
elseif strcmp(my_Simband_datafilename,'4017')
    sub_point = 150;
    Simband_ECG = Simband_ECG(sub_point:end);
    aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4001')
    add_point = 49;
    Simband_ECG = [zeros(1,add_point),Simband_ECG];
    aligned_ChestECG = [zeros(add_point,1);aligned_ChestECG];
elseif strcmp(my_Simband_datafilename,'4002')
    add_point = 41;
    Simband_ECG = [zeros(1,add_point),Simband_ECG];
    aligned_ChestECG = [zeros(add_point,1);aligned_ChestECG];
elseif strcmp(my_Simband_datafilename,'4012') || strcmp(my_Simband_datafilename,'4005') || strcmp(my_Simband_datafilename,'4013')
    add_point = 48;
    Simband_ECG = [zeros(1,add_point),Simband_ECG];
    aligned_ChestECG = [zeros(add_point,1);aligned_ChestECG];
end

if strcmp(my_Simband_datafilename,'4043')
   sub_point = round(1.2 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4038')
   sub_point = round(0.7 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4027')
   sub_point = round(0.74 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4022')
   sub_point = round(0.2 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4021')
   sub_point = round(0.66 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4037')
   sub_point = round(1.17 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
elseif strcmp(my_Simband_datafilename,'4034')
   sub_point = round(5.729 * fs);
   Simband_ECG = Simband_ECG(sub_point:end);
   aligned_ChestECG = aligned_ChestECG(sub_point:end);
       sub_point_ACC = round(3.607 * fs_ACC);
    ACC = ACC(sub_point_ACC:end);
end

PPG = PPG(:); % make sure PPG is column vector

Simband_ECG = Simband_ECG(:);  % make sure is column vector
ACC = ACC(:); % make sure is column vector
aligned_ChestECG = aligned_ChestECG(:); % make sure is column vector
end

% subfunction 4
function [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_load_ref_ECG_peak(my_Simband_datafilename,Simband_data_folder,fs)
if (strcmp(my_Simband_datafilename,'4005'))
    load([Simband_data_folder,'\','4005_ECG_30sec']);
    load([Simband_data_folder,'\','4005_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4002'))
    load([Simband_data_folder,'\','4002_ECG_30sec']);
    load([Simband_data_folder,'\','4002_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4001'))
    load([Simband_data_folder,'\','4001_ECG_30sec']);
    load([Simband_data_folder,'\','4001_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4006'))
    load([Simband_data_folder,'\','4006_ECG_30sec']);
    load([Simband_data_folder,'\','4006_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4012'))
    load([Simband_data_folder,'\','4012_ECG_30sec']);
    load([Simband_data_folder,'\','4012_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4013'))
    load([Simband_data_folder,'\','4013_ECG_30sec']);
    load([Simband_data_folder,'\','4013_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4016'))
    load([Simband_data_folder,'\','4016_ECG_30sec']);
    load([Simband_data_folder,'\','4016_Peak_Ref_ECG']);
    [new_ECG_peak_loc,new_ECG,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs);
elseif (strcmp(my_Simband_datafilename,'4033') ...
        || strcmp(my_Simband_datafilename,'4045') ...
        || strcmp(my_Simband_datafilename,'4043') ...
        || strcmp(my_Simband_datafilename,'4042') ...
        || strcmp(my_Simband_datafilename,'4040') ...
        || strcmp(my_Simband_datafilename,'4038') ...
        || strcmp(my_Simband_datafilename,'4037') ...
        || strcmp(my_Simband_datafilename,'4027') ...
        || strcmp(my_Simband_datafilename,'4025') ...
        || strcmp(my_Simband_datafilename,'4024') ...
        || strcmp(my_Simband_datafilename,'4022') ...
        || strcmp(my_Simband_datafilename,'4021') ...
        || strcmp(my_Simband_datafilename,'4019') ...
        || strcmp(my_Simband_datafilename,'4015'))

    load([Simband_data_folder,'\',my_Simband_datafilename,'_Peak_Ref_ECG']);
    new_ECG_peak_loc = save_ECG_peak_loc;
    new_ECG = [];
    seg_num_ECG = [];
    seg_num_peak = [];
else
    new_ECG_peak_loc = [];
    new_ECG = [];
    seg_num_ECG = [];
    seg_num_peak = [];
end
end

% subfunction 5
function [countseg,countlen] = my_count_total_seg(shift_sec,len_PPG,len_Simband_ECG,len_ACC,fs,fs_ACC,fs_PPG)
split_sec = 30;

new_len_PPG = len_PPG - shift_sec*fs_PPG;
new_len_Simband_ECG = len_Simband_ECG - shift_sec*fs_ACC;
new_len_ACC = len_ACC - shift_sec*fs;

countseg_PPG = floor(new_len_PPG/split_sec/fs_PPG);
countseg_ECG = floor(new_len_Simband_ECG/split_sec/fs);
countseg_ACC = floor(new_len_ACC/split_sec/fs_ACC);
countseg = max([countseg_PPG,countseg_ECG,countseg_ACC]);
countlen = min([len_PPG,floor(len_Simband_ECG/fs*fs_PPG),floor(len_ACC/fs_ACC*fs_PPG)]);
end

% subfunction 6
function [start_win_idx,end_win_idx,subject_win_idx] = my_set_start_end_win_idx(user_input_subject,user_input_win,all_win_flag,aaa,countseg,start_win_idx,end_win_idx)
subject_win_idx = 1:countseg;
if user_input_subject == 0  % all subjects
% all windows or custom window?
    if all_win_flag == true
        end_win_idx = countseg; % It was NaN.
    end
else % user specific subject.
    if user_input_win == 0
        end_win_idx = countseg; % It was NaN.
    end
end
end

% subfunction 7
function [new_PPG_signal,new_ECG_signal,new_ACC_signal] = myfillgap(data,fs)
% clear all;close all;clc;
% warning off all

PPG = data.band.physiosignal.ppg4.visual.signal;
% PPG = data.band.physiosignal.ppg4.signal;
Simband_ECG = data.band.physiosignal.ecg.visual.signal;
ACC = data.band.accelerometer.magnitude.signal;

Simband_PPG_timestamps = data.band.physiosignal.ppg4.visual.timestamps;
% Simband_PPG_timestamps = data.band.physiosignal.ppg4.timestamps;
Simband_ECG_timestamps = data.band.physiosignal.ecg.visual.timestamps;
Simband_ACC_timestamps = data.band.accelerometer.magnitude.timestamps;

PPG_diff = diff(Simband_PPG_timestamps);
ECG_diff = diff(Simband_ECG_timestamps);
ACC_diff = diff(Simband_ACC_timestamps);

% figure;plot(PPG_diff);title('derivative of the PPG time stamps');ylabel('time (sec)');xlabel('index of array (integer)');
% figure;plot(ECG_diff);title('derivative of the ECG time stamps');ylabel('time (sec)');xlabel('index of array (integer)');
% figure;plot(ACC_diff);title('derivative of the ACC time stamps');ylabel('time (sec)');xlabel('index of array (integer)');
% 
% PPG_two_length_diff = length(Simband_PPG_timestamps) - length(PPG)
% ECG_two_length_diff = length(Simband_ECG_timestamps) - length(Simband_ECG)
% ACC_two_length_diff = length(Simband_ACC_timestamps) - length(ACC)

% ------------For PPG filling gaps----------------
% Simband_PPG_timestamps = Simband_PPG_timestamps';
PPG_gaps = find(PPG_diff ~= (1/fs));
PPG_gaps_locs = [0,PPG_gaps,length(Simband_PPG_timestamps)];
new_Simband_PPG_timestamps = 0;
new_PPG_signal = 0;
if isempty(PPG_gaps) ~= 1
    for ii = 1:(length(PPG_gaps_locs)-2)
        middle_time_start = (Simband_PPG_timestamps(PPG_gaps_locs(ii+1))+(1/fs)):(1/fs):(Simband_PPG_timestamps(PPG_gaps_locs(ii+1)+1)-(1/fs));
        new_Simband_PPG_timestamps = [new_Simband_PPG_timestamps,Simband_PPG_timestamps((PPG_gaps_locs(ii)+1):PPG_gaps_locs(ii+1)),middle_time_start];
        PPG_zero = zeros(size(middle_time_start));
        new_PPG_signal = [new_PPG_signal,PPG((PPG_gaps_locs(ii)+1):PPG_gaps_locs(ii+1)),PPG_zero];
    end
    new_Simband_PPG_timestamps = [new_Simband_PPG_timestamps,Simband_PPG_timestamps((PPG_gaps_locs(end-1)+1):PPG_gaps_locs(end))];
    new_Simband_PPG_timestamps = new_Simband_PPG_timestamps(2:end);
    new_PPG_signal = [new_PPG_signal,PPG((PPG_gaps_locs(end-1)+1):PPG_gaps_locs(end))];
    new_PPG_signal = new_PPG_signal(2:end);
else
    new_Simband_PPG_timestamps = Simband_PPG_timestamps;
    new_PPG_signal = PPG;
end

% figure;plot(diff(new_Simband_PPG_timestamps));title('the derivative of PPG stamps');
% figure;plot(new_PPG_signal,'linewidth',2);
% ------------For ECG filling gaps----------------
ECG_gaps = find(ECG_diff ~= (1/fs));
ECG_gaps_locs = [0,ECG_gaps,length(Simband_ECG_timestamps)];
new_Simband_ECG_timestamps = 0;
new_ECG_signal = 0;
if isempty(ECG_gaps) ~= 1
    for ii = 1:(length(ECG_gaps_locs)-2)
        middle_time_start = (Simband_ECG_timestamps(ECG_gaps_locs(ii+1))+(1/fs)):(1/fs):(Simband_ECG_timestamps(ECG_gaps_locs(ii+1)+1)-(1/fs));
        new_Simband_ECG_timestamps = [new_Simband_ECG_timestamps,Simband_ECG_timestamps((ECG_gaps_locs(ii)+1):ECG_gaps_locs(ii+1)),middle_time_start];
        ECG_zero = zeros(size(middle_time_start));
        new_ECG_signal = [new_ECG_signal,Simband_ECG((ECG_gaps_locs(ii)+1):ECG_gaps_locs(ii+1)),ECG_zero];
    end
    new_Simband_ECG_timestamps = [new_Simband_ECG_timestamps,Simband_ECG_timestamps((ECG_gaps_locs(end-1)+1):ECG_gaps_locs(end))];
    new_Simband_ECG_timestamps = new_Simband_ECG_timestamps(2:end);
    new_ECG_signal = [new_ECG_signal,Simband_ECG((ECG_gaps_locs(end-1)+1):ECG_gaps_locs(end))];
    new_ECG_signal = new_ECG_signal(2:end);
else
    new_Simband_ECG_timestamps = Simband_ECG_timestamps;
    new_ECG_signal = Simband_ECG;
end

% figure;plot(diff(new_Simband_ECG_timestamps));title('the derivative of ECG stamps');
% figure;plot(new_ECG_signal,'linewidth',2);
% ------------For ACC filling gaps--------------
ACC_gaps = find(ACC_diff ~= (1/fs));
ACC_gaps_locs = [0,ACC_gaps,length(Simband_ACC_timestamps)];
new_Simband_ACC_timestamps = 0;
new_ACC_signal = 0;
if isempty(ACC_gaps) ~= 1
    for ii = 1:(length(ACC_gaps_locs)-2)
        middle_time_start = (Simband_ACC_timestamps(ACC_gaps_locs(ii+1))+(1/fs)):(1/fs):(Simband_ACC_timestamps(ACC_gaps_locs(ii+1)+1)-(1/fs));
        new_Simband_ACC_timestamps = [new_Simband_ACC_timestamps,Simband_ACC_timestamps((ACC_gaps_locs(ii)+1):ACC_gaps_locs(ii+1)),middle_time_start];
        ACC_zero = zeros(size(middle_time_start));
        new_ACC_signal = [new_ACC_signal,ACC((ACC_gaps_locs(ii)+1):ACC_gaps_locs(ii+1)),ACC_zero];
    end
    new_Simband_ACC_timestamps = [new_Simband_ACC_timestamps,Simband_ACC_timestamps((ACC_gaps_locs(end-1)+1):ACC_gaps_locs(end))];
    new_Simband_ACC_timestamps = new_Simband_ACC_timestamps(2:end);
    new_ACC_signal = [new_ACC_signal,ACC((ACC_gaps_locs(end-1)+1):ACC_gaps_locs(end))];
    new_ACC_signal = new_ACC_signal(2:end);
else
    new_Simband_ACC_timestamps = Simband_ACC_timestamps;
    new_ACC_signal = ACC;
end

% figure;plot(diff(new_Simband_ACC_timestamps));title('the derivative of ACC time stamps');
% figure;plot(new_ACC_signal,'linewidth',2);


end

function [new_ECG_peak_loc,temp,seg_num_ECG,seg_num_peak] = my_ref_ECG_concatenate(Peak_Ref_ECG,ECG_30sec,fs)
new_ECG_peak_loc = [];
for ii = 1:size(Peak_Ref_ECG,2) % column number is the 30 sec index
    temp = Peak_Ref_ECG(:,ii);
    temp(isnan(temp)) = [];
    temp = temp + fs*(ii-1)*30; % add previous time into this ECG index.
    new_ECG_peak_loc = [new_ECG_peak_loc;temp];
end
temp = ECG_30sec(:);
temp = temp(:);
seg_num_ECG = size(ECG_30sec,1);
seg_num_peak = size(Peak_Ref_ECG,2);
end