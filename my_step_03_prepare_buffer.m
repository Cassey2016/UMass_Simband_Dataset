function [PPG_buffer,Simband_ECG_buffer,ACC_buffer,Ref_ECG_buffer,mytitle] = my_step_03_prepare_buffer(output_data)

% undock all data:
PPG = output_data.PPG;
Simband_ECG = output_data.Simband_ECG;
ACC = output_data.ACC;
aligned_ChestECG = output_data.aligned_ChestECG;
iiii_PPG_start = output_data.iiii_PPG_start;
iiii_PPG_end = output_data.iiii_PPG_end;
iiii_ECG_start = output_data.iiii_ECG_start;
iiii_ECG_end = output_data.iiii_ECG_end;
iiii_ACC_start = output_data.iiii_ACC_start;
iiii_ACC_end = output_data.iiii_ACC_end;
load_Simband_subject_struct = output_data.load_Simband_subject_struct;
my_Simband_datafilename = output_data.my_Simband_datafilename;
aaa = output_data.aaa;

PPG_buffer = PPG(iiii_PPG_start:iiii_PPG_end);
Simband_ECG_buffer = Simband_ECG(iiii_ECG_start:iiii_ECG_end);
ACC_buffer = ACC(iiii_ACC_start:iiii_ACC_end);
if iiii_ECG_start > length(aligned_ChestECG) 
    Ref_ECG_buffer = zeros(size(Simband_ECG_buffer)); % no enough reference ECG
elseif iiii_ECG_end > length(aligned_ChestECG)
    Ref_ECG_buffer = [aligned_ChestECG(iiii_ECG_start:end);zeros(length(Simband_ECG_buffer) - (length(aligned_ChestECG) - iiii_ECG_start) - 1,1)]; % some ref ECG with padded zeros.
else % enough ECG
    Ref_ECG_buffer = aligned_ChestECG(iiii_ECG_start:iiii_ECG_end);% Reference HR: peak detection for ECG
end

end