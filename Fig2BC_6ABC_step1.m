
clc
clear
close all

% subject index
subs = [1:3,5,7,9:12,14:15,19:27,29:38];
filename = 'Exp1.csv';
t1 = [-3500,0]; % FW
t2 = [-350,-150]; % ap

% load behavioral data
load threshold.mat

% input path of preprocessed EEG
inputParentFolder = '.../Data/Exp1';

% channels of interest to analyze traveling waves
Chan = {'Oz','O1' 'PO3' 'POz' 'PO4' 'O2','Pz','CPz','Cz','FCz','Fz'};
channels = zeros(numel(Chan),1);

% prepare long data for R
[block,ratings,thres_all,FW,BW,AP] = deal(zeros(30*16*40,1)); % 30 subs, 16 blocks, 40 trials in each block

for sub_index = 1:numel(subs)
    
    % load behavioral data
    thres = threshold{sub_index};
    thres = reshape(thres',640,1);
    
    % load preprocessed eeg
    subName = ['sub' num2str(subs(sub_index))]; % number to stringr
    EEG = pop_loadset('filename',[subName '.set'],'filepath',inputParentFolder);

    % interpolate FCz
    EEG.chanlocs(13).labels = 'FCz';
    EEG.chanlocs(13).sph_radius = 1;
    EEG.chanlocs(13).sph_theta = 0;
    EEG.chanlocs(13).sph_phi = 67;
    EEG.chanlocs(13).theta = 0;
    EEG.chanlocs(13).radius = 0.127777777777778;
    EEG.chanlocs(13).X = 0.390731128489274;
    EEG.chanlocs(13).Y = 0;
    EEG.chanlocs(13).Z = 0.920504853452440;
    EEG.chanlocs(13).type = [];
    EEG.chanlocs(13).ref = 'TP9 TP10';
    EEG.chanlocs(13).urchan = [];
 
    EEG = pop_interp(EEG, [13], 'spherical');
 
    % find raw trial and block index corresponding to preprocessed trials
    total_trial = [];
    for i = 1:numel(EEG.urevent)
        if EEG.urevent(i).type(3)=='1'& EEG.urevent(i).type(4)=='3'
            total_trial = [total_trial;EEG.urevent(i).bvmknum];
        else
            total_trial = total_trial;
        end
        
    end
    
    total_trial(:,2) = [1:640];
    total_trial(:,3) = [repmat(1,40,1);
        repmat(2,40,1);
        repmat(3,40,1);
        repmat(4,40,1);
        repmat(5,40,1);
        repmat(6,40,1);
        repmat(7,40,1);
        repmat(8,40,1);
        repmat(9,40,1);
        repmat(10,40,1);
        repmat(11,40,1);
        repmat(12,40,1);
        repmat(13,40,1);
        repmat(14,40,1);
        repmat(15,40,1);
        repmat(16,40,1);
        ];
    
    remain_trial = [];
    for i = 1:numel(EEG.event)
        if EEG.event(i).type(3)=='1'& EEG.event(i).type(4)=='3'
            remain_trial = [remain_trial;EEG.event(i).bvmknum];
        else
            remain_trial = remain_trial;
        end
        
    end
    
    thres_remain = thres(ismember(total_trial(:,1)',remain_trial));

    % calculate BW and FW traveling waves
    
    [~,t1_idx(1)] = min(abs(t1(1)-EEG.times));
    [~,t1_idx(2)] = min(abs(t1(2)-EEG.times));
    [~,t2_idx(1)] = min(abs(t2(1)-EEG.times));
    [~,t2_idx(2)] = min(abs(t2(2)-EEG.times));

    EEG = pop_select( EEG, 'channel', Chan);

    data_trials = EEG.data; 
    [fw,bw,ap1,ap2] = deal(zeros(EEG.trials,1));

    for tt = 1:EEG.trials
        [fw(tt),bw(tt)] = quantifyingTW(data_trials([11,10,6,9,5,3,4],t1_idx(1):t1_idx(2),tt),EEG.srate);
    
    end    
    
    ap = tfpower_exp1(EEG,t2);
    
    % prepare long data for R analysis
    block(16*40*(sub_index-1)+1:16*40*sub_index) = total_trial(:,3);
    ratings(16*40*(sub_index-1)+1:16*40*sub_index) = fati_all;
    thres_all(16*40*(sub_index-1)+1:16*40*sub_index) = thres;
    [temp1, temp2, temp3,temp4, temp5, temp6] = deal(nan(640,1));
    temp1(ismember(total_trial(:,1)',remain_trial)) = fw;
    temp2(ismember(total_trial(:,1)',remain_trial)) = bw;
    temp3(ismember(total_trial(:,1)',remain_trial)) = ap;

    FW(16*40*(sub_index-1)+1:16*40*sub_index) = temp1;
    BW(16*40*(sub_index-1)+1:16*40*sub_index) = temp2;
    AP(16*40*(sub_index-1)+1:16*40*sub_index) = temp3;

    clear fw bw ap
     
end

% generate long data for Fig2BC_6ABC_step2.R
id = repmat([1:numel(subs)],640,1);
id = reshape(id,[],1);

dataTable = table(id,block,ratings,thres_all,FW,BW,AP);
writetable(dataTable,filename);



