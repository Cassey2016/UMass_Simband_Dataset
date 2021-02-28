% 
% function output_struct = my_step_01_know_sub_seg()
% 
% Author: Dong Han (dong.han@uconn.edu)
% Date: 05/14/2020
% 
% Description: Ask user loading which Simband subject and which segment 
%              of that subject.
% Return: user_input_subject - 0: all subjects; 1-41: specific subject.
%         load_Simband_subject_name - the subject name(s) user selected.
%         load_Simband_subject_struct - the data structure of Simband.
%         user_input_win - 0: all segments; 1-unknown: segment index.
%         start_win_idx - starting segment index.
%         end_win_idx - ending segment index.
%         all_win_flag - if user input running 'all segments.'
%         Simband_data_folder - the path where Simband data stored.
%
% Last Modification: Dong, 02/28/2021.

function output_struct = my_step_01_know_sub_seg()
    Simband_data_folder = '..\Data';
    load([Simband_data_folder,'\','UMass_SimbandInfo']);
    UMassSimbandInfo
    prompt = 'Please select which subject you would like to run: \n   0. All subjects;\n';
    user_input_subject = input(prompt);
    if user_input_subject == 0
       disp('--- Selected all users for running! ---');
       load_Simband_subject_name = UMassSimbandInfo(:,2);
       load_Simband_subject_struct = UMassSimbandInfo(:,3);

       prompt = 'Please select which 30-sec segment you want to run:\n   0. All segments;\n   2. Custom segment;\n';
       user_input_win = input(prompt);
       if user_input_win == 0 % run all windows
           start_win_idx = 1;
%            end_win_idx = true;
           end_win_idx = NaN; % Will be assigned in func 'my_step_02_load_data.m'.
           all_win_flag = true;
           fprintf('--- Selected all segments to run on all subjects! ---\n');
       else % custom windows
           prompt = 'Please input segment index:\n';
           user_input_win = input(prompt);
           start_win_idx = user_input_win;
           end_win_idx = start_win_idx;
           all_win_flag = false;
           fprintf('--- Selected segment %d in all subjects for running! ---\n',start_win_idx);
       end
    else % run specific subject.
       load_Simband_subject_name = UMassSimbandInfo(user_input_subject,2);
       load_Simband_subject_struct = UMassSimbandInfo(user_input_subject,3);
       fprintf('--- Selected subject %d for running! ---\n',load_Simband_subject_name);
       prompt = 'Please select which 30-sec segment you want to run:\n   0. All segments;\n   2. Custom segment;\n';
       user_input_win = input(prompt);
       if user_input_win == 0 % run all windows
           start_win_idx = 1;
%            end_win_idx = true;
            end_win_idx = NaN; % Will be assigned in func 'my_step_02_load_data.m'.
            all_win_flag = true;
            fprintf('--- Selected all segments to run on subject %d! ---\n',load_Simband_subject_name);
       else % custom windows
           prompt = 'Please input segment index:\n';
           user_input_win = input(prompt);
           start_win_idx = user_input_win;
           end_win_idx = start_win_idx;

           all_win_flag = false;
           fprintf('--- Selected segment %d in subject %d for running! ---\n',start_win_idx,load_Simband_subject_name);
       end
    end
    output_struct = struct('user_input_subject',user_input_subject,...
                        'load_Simband_subject_name',load_Simband_subject_name,...
                        'load_Simband_subject_struct',load_Simband_subject_struct,...
                        'user_input_win',user_input_win,...
                        'start_win_idx',start_win_idx,...
                        'end_win_idx',end_win_idx,...
                        'all_win_flag',all_win_flag,...
                        'Simband_data_folder',Simband_data_folder);
end % of function