function [  ] = load_gaze_data_texas(  )
%function [  ] = load_data_oleg_fixed(  )
%   This function does preliminary processing of the data
%   collected by Oleg. We load the gaze data, and filter.
%

global ROOT_DIR SAMP_RATE B_FILT A_FILT NUM_TRIAL

ROOT_DIR = '~/Desktop/GitHub/eye_trace_analysis';

FILT_ORDER = 3;
FILT_CUTOFF = 90;
NUM_BLOCK = 4;
NUM_TRIAL = 100;

%Butterworth filter parameters
SAMP_RATE = 1017.3;
[B_FILT, A_FILT] = butter(FILT_ORDER, 2*FILT_CUTOFF/SAMP_RATE, 'low');

%load subject data
subj_texas = load_subject_data_texas_xlsx();
NUM_SUBJ = length(subj_texas.sleep);

for kk = 1:NUM_SUBJ
  if subj_texas.cut(kk); continue; end

  subj_id = num2str(kk);
  if (kk < 10)
    subj_id = ['00', subj_id];
  elseif (kk < 100)
    subj_id = ['0', subj_id];
  end

  if ~mod(kk, 25); fprintf('Loading data for Subject %s\n', subj_id); end
  
  tgt = new_struct({'time','x','y'}, 'dim',[1,NUM_BLOCK]);
  gaze = new_struct({'time','x','y','vx','vy','v'}, 'dim',[1,NUM_BLOCK]);
  
  [gaze(1), tgt(1)] = load_data_single_block([subj_id, '_S1_HSS.asc'], 'HSS');
  [gaze(2), tgt(2)] = load_data_single_block([subj_id, '_S2_HSS.asc'], 'HSS');
  [gaze(3), tgt(3)] = load_data_single_block([subj_id, '_S1_RAN.asc'], 'RAN');
  [gaze(4), tgt(4)] = load_data_single_block([subj_id, '_S2_RAN.asc'], 'RAN');
  
  save([ROOT_DIR, 'filt/data_', subj_id, '.mat'], 'gaze','tgt')

end%for:subjects(kk)

end%function:load_gaze_data_texas()


function [ gaze , tgt ] = load_data_single_block( file_id , blk_type )

global ROOT_DIR SAMP_RATE B_FILT A_FILT NUM_TRIAL

data = importdata([ROOT_DIR, 'raw/S_', file_id], '\t', 1);
data = data.data(:,[2,8,9,13,14]);

time = uint64(data(:,1));

idx_tgt = find((diff(data(:,4)) ~= 0) | (diff(data(:,5)) ~= 0)) + 1;

if strcmp(blk_type, 'HSS')
  idx_tgt(end) = [];
else
  idx_tgt = [1; idx_tgt];
end

%last block for subject S178 with only 99 trials -- fix this issue now
if strcmp(file_id, '178_S2_RAN.asc')
  idx_tgt = [idx_tgt; idx_tgt(end)];
end

if (length(idx_tgt) ~= NUM_TRIAL)
  warning('Number of trials for File %s is not %d\n', file_id, NUM_TRIAL)
end

tgt = struct('time',time(idx_tgt), 'x',data(idx_tgt, 4), 'y',data(idx_tgt, 5));

%filter gaze data
gaze_x = single(filtfilt(B_FILT, A_FILT, data(:,2)));
gaze_y = single(filtfilt(B_FILT, A_FILT, data(:,3)));

%get rid of NaN's based on gaze position
idx_nan = (abs(gaze_x) > 20) | (abs(gaze_y) > 13);
gaze_x(idx_nan) = NaN;
gaze_y(idx_nan) = NaN;

%differentiate gaze data
gaze_vx = diff(gaze_x) * SAMP_RATE;
gaze_vy = diff(gaze_y) * SAMP_RATE;
gaze_v = sqrt(gaze_vx.^2 + gaze_vy.^2);
gaze_x(end) = [];
gaze_y(end) = [];
time(end) = [];

gaze = struct('time',time, 'x',single(gaze_x), 'y',single(gaze_y), ...
'vx',single(gaze_vx), 'vy',single(gaze_vy), 'v',single(gaze_v));

end%function:load_data_single_block()

