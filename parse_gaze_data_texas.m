function [ ] = parse_gaze_data_texas( )
%[ ] = parse_gaze_data_texas( )

global ALLOT APPEND VEL_CUT MIN_IMI MIN_HOLD
global NUM_TRIAL NUM_SAMPLES_TRIAL
global MIN_DISP MAX_CURVE MAX_ERR_INIT MAX_PEAKVEL MAX_SKEW
global LIM_DISP LIM_DUR LIM_RT_RAN LIM_RT_HSS
global FIELDS_UINT16 FIELDS_SINGLE FIELDS_VECTOR
global ALLOC_MOVES

ROOT_DIR = '~/Documents/vigor_across_modality/texas/';

MIN_MOVE_HORZ = 50;
MIN_MOVE_VERT = 25;

ALLOT = 180;
APPEND = 30;

VEL_CUT = [15, 50, 80];

MIN_IMI = 20;
MIN_HOLD = 20;

NUM_BLOCK = 4;
NUM_TRIAL = 100;
NUM_SAMPLES_TRIAL = 1000;

%criteria: task-relevant
LIM_DISP = [0.6, 1.4];
MAX_ERR_INIT = 4.0;
% LIM_RT_RAN = [150, 400];
LIM_RT_RAN = [100, 350];
LIM_RT_HSS = [50, 350];

%criteria: saccade
MAX_CURVE = 3.0;
LIM_DUR = [20, 150];
MIN_DISP = 2.0;
MAX_PEAKVEL = 850;
MAX_SKEW = 0.75;

FIELDS_UINT16 = {'block','horz','quadrant','trial'};
FIELDS_SINGLE = {'init_x','init_y','init_t','fin_x','fin_y', ...
  'err_x','err_y','err_init','tgt_x','tgt_y','duration','rxntime','peakvel','skew', ...
  'amplitude','displacement','amp_err','amp_tgt','amp_pcent','dir_err'};
FIELDS_VECTOR = {'x','y','vx','vy','v'};

%load subject data
subjects = load_subject_data_texas_xlsx();
NUM_SUBJ = 1;%length(subjects.sleep);
subj_id = cell(1,NUM_SUBJ);

for kk = 1:NUM_SUBJ
  if subjects.cut(kk); continue; end
  
  moves_tr = new_struct([FIELDS_UINT16,FIELDS_SINGLE,FIELDS_VECTOR], 'dim',[1,NUM_BLOCK]);
  moves_tr = populate_struct(moves_tr, FIELDS_UINT16, uint16(zeros(1,NUM_TRIAL)));
  moves_tr = populate_struct(moves_tr, FIELDS_SINGLE, single(NaN(1,NUM_TRIAL)));
  moves_tr = populate_struct(moves_tr, FIELDS_VECTOR, single(NaN(ALLOT,NUM_TRIAL)));
  moves_tr = orderfields(moves_tr);
  
  moves_all = moves_tr;

  subj_id{kk} = num2str(kk);
  
  if (kk < 10)
   subj_id{kk} = ['00', subj_id{kk}];
  elseif (kk < 100)
   subj_id{kk} = ['0', subj_id{kk}];
  end
  
  load([ROOT_DIR, 'filt/data_', subj_id{kk}, '.mat'], 'gaze','tgt')
  
  for jj = 1:NUM_BLOCK
    
    ALLOC_MOVES = NUM_TRIAL;
    idx_all = 1;
    
    for tt = 1:NUM_TRIAL
      
      idx_trial = find(gaze(jj).time == tgt(jj).time(tt));
      idx_trial = (idx_trial : idx_trial+NUM_SAMPLES_TRIAL);
      
      idx_lim = identify_movement_candidates(gaze(jj), idx_trial);
      
      if isempty(idx_lim); continue; end
      
      NUM_CAND = size(idx_lim, 2);
      cands = new_struct([FIELDS_UINT16,FIELDS_SINGLE,FIELDS_VECTOR], 'dim',[1,NUM_CAND]);
      cands = populate_struct(cands, FIELDS_UINT16, uint16(0));
      cands = populate_struct(cands, FIELDS_SINGLE, single(NaN));
      cands = populate_struct(cands, FIELDS_VECTOR, single(NaN(ALLOT,1)));
      cands = orderfields(cands);
      
      for cc = NUM_CAND:-1:1
        
        [idx_lim(1,cc), idx_lim(2,cc), skip] = identify_movement_bounds(gaze(jj).v(idx_trial), idx_lim(1,cc), ...
          idx_lim(2,cc), 'v_cut',VEL_CUT(2), 'min_hold',MIN_HOLD);
        
        if (skip); cands(cc) = []; continue; end
        
        idx_init_cc = idx_trial(idx_lim(1,cc));
        idx_fin_cc = idx_trial(idx_lim(2,cc));
        
        if (idx_init_cc <= APPEND); cands(cc) = []; continue; end
        
        %provide precise account of reaction time
        offset_init = find(gaze(jj).v(idx_init_cc : -1 : idx_init_cc-APPEND) < VEL_CUT(1), 1, 'first');
        if isempty(offset_init)
          cands(cc) = []; continue
        else
          idx_init_cc = idx_init_cc - offset_init;
        end
        
        cands(cc) = save_candidate_parm(cands(cc), gaze(jj), tgt(jj), tt, jj, ...
          idx_init_cc, idx_fin_cc);
        
      end%for:candidates(cc)
      
      [moves_all(jj), cands, idx_all] = identify_all_movement(moves_all(jj), cands, idx_all);
      
      if (jj < 3)
        moves_tr(jj) = identify_taskrel_movement(moves_tr(jj), cands, tt, 'blk_type','HSS');
      else
        moves_tr(jj) = identify_taskrel_movement(moves_tr(jj), cands, tt, 'blk_type','RAN');
      end
      
    end%for:trials(tt)
    
    %remove extra memory from struct with all movements
    for ff = 1:length(FIELDS_UINT16)
      moves_all(jj).(FIELDS_UINT16{ff})(idx_all:ALLOC_MOVES) = [];
    end
    for ff = 1:length(FIELDS_SINGLE)
      moves_all(jj).(FIELDS_SINGLE{ff})(idx_all:ALLOC_MOVES) = [];
    end
    for ff = 1:length(FIELDS_VECTOR)
      moves_all(jj).(FIELDS_VECTOR{ff})(:,idx_all:ALLOC_MOVES) = [];
    end
    
  end%for:blocks(jj)
  
  taskrel = ~isnan([moves_tr.displacement]);
  move_horz = logical([moves_tr.horz]);
  task_ran = [ false(1,200) , true(1,200) ];
  
  idx_horz_HSS = (taskrel & ~task_ran);
  idx_horz_RAN = (taskrel & task_ran & move_horz);
  idx_vert = (taskrel & ~move_horz);
  
  if ( (sum(idx_horz_HSS) >= MIN_MOVE_HORZ) && (sum(idx_horz_RAN) >= MIN_MOVE_HORZ) && ...
      (sum(idx_vert) >= MIN_MOVE_VERT) )
    save([ROOT_DIR, 'moves/moves_parsed_', subj_id{kk}, '.mat'], 'moves_all', 'moves_tr')
  else
    fprintf('**Warning S%s -- Only %d HSS horz and %d RAN horz and %d vert moves\n', ...
      subj_id{kk}, sum(idx_horz_HSS), sum(idx_horz_RAN), sum(idx_vert))
    save([ROOT_DIR, 'moves/moves_parsed_', subj_id{kk}, '.mat'], 'moves_all', 'moves_tr')
  end
  
end%for:subjects(kk)


end%function:parse_gaze_data_texas()


function [ idx_lim ] = identify_movement_candidates( gaze , idx_trial )

global APPEND VEL_CUT MIN_IMI
global NUM_SAMPLES_TRIAL

idx_lim = [];

idx_prelim = find(gaze.v(idx_trial) > VEL_CUT(3)); %all samples with velocity greater than cutoff
if isempty(idx_prelim)
%   fprintf('**No data points with velocity above %d\n', VEL_CUT(3))
  return
end

idx_jump = find(diff(idx_prelim) > MIN_IMI);
idx_start = [idx_prelim(1); idx_prelim(idx_jump + 1)]';
idx_end   = [idx_prelim(idx_jump); idx_prelim(end)]';

idx_clipped = ((idx_start < APPEND) | (idx_end > (NUM_SAMPLES_TRIAL-100)));
idx_start(idx_clipped) = [];
idx_end(idx_clipped) = [];

%identify preliminary candidates that are blinks
num_prelim = length(idx_start);
for cand = num_prelim:-1:1
  vel_start = gaze.v( idx_trial(idx_start(cand))-APPEND+1 : idx_trial(idx_start(cand)) );
  vel_end = gaze.v( idx_trial(idx_end(cand)) : idx_trial(idx_end(cand))+APPEND-1 );
  if sum(isnan([vel_start;vel_end]))
    idx_start(cand) = [];
    idx_end(cand) = [];
  end
end

idx_lim = [ idx_start ; idx_end ];

end%function:identify_movement_candidates()

function [ cands ] = save_candidate_parm( cands , gaze , tgt , trial , block , idx_init , idx_fin )

global ALLOT APPEND

%get indexes for full plotting of traces
i_rec = idx_init - APPEND : idx_init - APPEND + ALLOT - 1;
i_rec(i_rec < 1) = 1; %if first few points lie outside trial limits

%save candidate kinematics
cands.x(:) = gaze.x(i_rec);
cands.y(:) = gaze.y(i_rec);
cands.vx(:) = gaze.vx(i_rec);
cands.vy(:) = gaze.vy(i_rec);
cands.v(:) =  gaze.v(i_rec);

%saccade timing and trial number
cands.init_t = gaze.time(idx_init);
cands.duration = gaze.time(idx_fin) - cands.init_t;
cands.rxntime = single(gaze.time(idx_init)) - single(tgt.time(trial));

cands.trial = trial;
cands.block = block;

%% Gaze shift parameters

%start- and end-point
cands.init_x = gaze.x(idx_init);
cands.init_y = gaze.y(idx_init);
cands.fin_x =  gaze.x(idx_fin);
cands.fin_y =  gaze.y(idx_fin);

if (trial == 1)
  err_init_x = gaze.x(idx_init);
  err_init_y = gaze.y(idx_init);
else
  err_init_x = gaze.x(idx_init) - tgt.x(trial - 1);
  err_init_y = gaze.y(idx_init) - tgt.y(trial - 1);
end
cands.err_init = sqrt(err_init_x^2 + err_init_y^2);

%amplitude
amp_x = cands.fin_x - cands.init_x;
amp_y = cands.fin_y - cands.init_y;
cands.amplitude =  sqrt(amp_x^2 + amp_y^2);

%displacement
disp_x = sum(abs(diff(gaze.x(idx_init:idx_fin)))); %saccades(kk).a_disp = a_disp;
disp_y = sum(abs(diff(gaze.y(idx_init:idx_fin)))); %saccades(kk).e_disp = e_disp;
cands.displacement = sqrt(disp_x^2 + disp_y^2);

%target location
cands.tgt_x = tgt.x(trial);
cands.tgt_y = tgt.y(trial);

%target jump amplitude and direction
if (trial > 1)
  amp_tgt_x = cands.tgt_x - tgt.x(trial - 1);
  amp_tgt_y = cands.tgt_y - tgt.y(trial - 1);
else
  amp_tgt_x = NaN;
  amp_tgt_y = NaN;
end
cands.amp_tgt = sqrt(amp_tgt_x^2 + amp_tgt_y^2);

%direction x target direction
vec_sacc = double([amp_x, amp_y]);
vec_tgt = [amp_tgt_x, amp_tgt_y];
cands.dir_err = acosd(dot(vec_sacc,vec_tgt)/(norm(vec_sacc)*norm(vec_tgt)));

%displacement x target displacement
cands.amp_pcent = cands.amplitude / cands.amp_tgt;

%endpoint x target location
cands.err_x = abs(cands.fin_x) - abs(cands.tgt_x);
cands.err_y = abs(cands.fin_y) - abs(cands.tgt_y);
cands.amp_err = sqrt(cands.err_x^2 + cands.err_y^2);

%saccade quadrant
move_dir = atan2d(amp_y, amp_x);
if (move_dir > -45 && move_dir <=  45)
  cands.quadrant = 1; %right
  cands.horz = 1;
elseif (move_dir >  45 && move_dir <=  135)
  cands.quadrant = 2; %top
elseif (move_dir > 135 || move_dir <= -135)
  cands.quadrant = 3; %left
  cands.horz = 1;
elseif (move_dir >-135 && move_dir <= -45)
  cands.quadrant = 4; %bottom
else
  error('Movement quadrant undefined.')
end

%velocity parameters
[cands.peakvel, i_pv] = max(gaze.v(idx_init:idx_fin));
cands.skew = i_pv / (idx_fin - idx_init);

end%function:save_movement_param

function [ moves , cands , index ] = identify_all_movement( moves , cands , index )

global MAX_CURVE MIN_DISP LIM_DUR MAX_PEAKVEL MAX_SKEW
global FIELDS_UINT16 FIELDS_SINGLE FIELDS_VECTOR
global NUM_TRIAL ALLOT ALLOC_MOVES

%% Identify saccades (both task-irrelevant and task-relevant)

icut_curve = ([cands.displacement]-[cands.amplitude] > MAX_CURVE);
icut_disp = ([cands.displacement] < MIN_DISP);
icut_dur = ([cands.duration] < LIM_DUR(1) | [cands.duration] > LIM_DUR(2));
icut_pv = ([cands.peakvel] > MAX_PEAKVEL);
icut_skew = ([cands.skew] > MAX_SKEW);

icut_saccade = (icut_curve | icut_disp | icut_dur | icut_pv | icut_skew);
idx_saccade = find(~icut_saccade);
num_saccade  = length(idx_saccade);

%check for memory re-alloc
if (index + num_saccade > ALLOC_MOVES)
  ALLOC_MOVES = ALLOC_MOVES + NUM_TRIAL;
  
  for ff = 1:length(FIELDS_UINT16)
    moves.(FIELDS_UINT16{ff}) = [moves.(FIELDS_UINT16{ff}), uint16(zeros(1,NUM_TRIAL))];
  end
  for ff = 1:length(FIELDS_SINGLE)
    moves.(FIELDS_SINGLE{ff}) = [moves.(FIELDS_SINGLE{ff}), single(NaN(1,NUM_TRIAL))];
  end
  for ff = 1:length(FIELDS_VECTOR)
    moves.(FIELDS_VECTOR{ff}) = [moves.(FIELDS_VECTOR{ff}), single(NaN(ALLOT,NUM_TRIAL))];
  end
end


%save all saccades for this trial
for ff = 1:length(FIELDS_UINT16)
  moves.(FIELDS_UINT16{ff})(index:index+num_saccade-1) = [cands(idx_saccade).(FIELDS_UINT16{ff})];
end
for ff = 1:length(FIELDS_SINGLE)
  moves.(FIELDS_SINGLE{ff})(index:index+num_saccade-1) = [cands(idx_saccade).(FIELDS_SINGLE{ff})];
end
for ff = 1:length(FIELDS_VECTOR)
  moves.(FIELDS_VECTOR{ff})(:,index:index+num_saccade-1) = [cands(idx_saccade).(FIELDS_VECTOR{ff})];
end

cands = cands(idx_saccade);

index = index + num_saccade;

end%function:identify_all_movement()

function [ moves ] = identify_taskrel_movement( moves , cands , trial , varargin )

args = getopt(varargin, {{'blk_type=','HSS'}});

global LIM_DISP LIM_RT_RAN LIM_RT_HSS MAX_ERR_INIT
global FIELDS_UINT16 FIELDS_SINGLE FIELDS_VECTOR

%% Identify task-relevant saccades

ratio_disp = [cands.amplitude] ./ [cands.amp_tgt];
icut_disp = ((ratio_disp < LIM_DISP(1)) | (ratio_disp > LIM_DISP(2)));

if strcmp(args.blk_type, 'RAN')
  icut_rt = (([cands.rxntime] < LIM_RT_RAN(1)) | ([cands.rxntime] > LIM_RT_RAN(2)));
else
  icut_rt = (([cands.rxntime] < LIM_RT_HSS(1)) | ([cands.rxntime] > LIM_RT_HSS(2)));
end

icut_xinit = ([cands.err_init] > MAX_ERR_INIT);

icut_taskrel = (icut_disp | icut_rt | icut_xinit);
idx_taskrel = find(~icut_taskrel);

num_taskrel = length(idx_taskrel);

if (num_taskrel == 1)
  
  for ff = 1:length(FIELDS_UINT16)
    moves.(FIELDS_UINT16{ff})(trial) = cands(idx_taskrel).(FIELDS_UINT16{ff});
  end
  for ff = 1:length(FIELDS_SINGLE)
    moves.(FIELDS_SINGLE{ff})(trial) = cands(idx_taskrel).(FIELDS_SINGLE{ff});
  end
  for ff = 1:length(FIELDS_VECTOR)
    moves.(FIELDS_VECTOR{ff})(:,trial) = cands(idx_taskrel).(FIELDS_VECTOR{ff});
  end
  
end

end%function:identify_taskrel_movement()


