
clc
clear
close all

t1= [-3500 0]; % traveling wave
t2 = [-3500,-150];
lap = 0;
filename = 'Exp2.csv';

% select subjects
subs = [1,3:13,15:16,19:37]; % Sub 14&38 had bad EEG signals. Sub 17 did not finish the experiment. Sub 18 had strange behavioral performance.
sti = [1,2,3]; % stimulation group: 1 denotes 45° TW-tACS, 2 denotes -45° TW-tACS, 3 denotes sham

% load stimulation sequence
load('sequence.mat');
sequence1 = sequence(subs',:);
[a,seqIndex] = sort(sequence1');

%% analyze behavioral data

% initiate visual threholds
threshold = zeros(numel(subs),3,50,9); % analyze blocks 2-10

% tired ratings
T = zeros(numel(subs),3,2,9); % tiredness rating: 3 stimulation groups, 2 denotes before and after each block, 9 blocks

for i = 1:numel(subs)
    for j = 1:numel(sti)
        
        % load behavioral data: visual threshold
        subName = num2str(subs(i));
        folderpath =  ['...Data/Exp2_behData/sub' subName '/' num2str(sti(j))];
        cd(folderpath)
        load(fullfile(folderpath,'behavioral_data.mat'))
        rsp = rsp(:,:,2);
        rsp = reshape(rsp',500,1);
        load(fullfile(folderpath,'contrast2.mat'))
        thr = reshape(contrast',500,1);
        
        % load behavioral data: tiredness rating
        load(fullfile(folderpath,'tiredScores1.mat'))
        load(fullfile(folderpath,'tiredScores2.mat'))
        T(i,j,1,:) = tiredScores1(2:end,1);
        T(i,j,2,:) = tiredScores2(2:end,1);
        
        % extract threshold in each block
        for block_idx = 2:10 % ten blocks in total, but analyze block 2-9
            threshold(i,j,:,block_idx-1) = thr(1+50*(block_idx-1):50*(block_idx));
        end
        
    end
end

% arrange data sequence according to stimulation sequence intead of arrival
% sequence
for i = 1:numel(subs)
    threshold(i,:,:,:) = threshold(i,seqIndex(:,i),:,:);
end


%% calculate traveling wave

% input folder
inputFolder = '.../Data/Exp2';

% initiate output variables
[FW,BW,AP] = deal(nan(numel(subs),3,450));% 3 stimulation groups, 50*9 trials

for sub_index = 1:numel(subs)
    for sti_index = 1:numel(sti)
        
        % load preprocessed data using eeglab
        subName = ['sub' num2str(subs(sub_index)) '_' num2str(sti_index)]; % number to stringr
        EEG = pop_loadset('filename',[subName '.set'],'filepath',inputFolder);
        
        if lap ==1
            EEG.data = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
        end

        % select channels to calculate traveling waves
        Channels = {'Oz','POz','Pz','CPz','Cz','FCz','Fz'}; % select middle line channels to analyze traveling waves
        for cc = 1:numel(Channels)
            channel_index(cc) = find(strcmpi(Channels{cc}, {EEG.chanlocs.labels}));
        end
        
        % select channels to calculate prestimulus alpha power
        Channels_occipital = {'Oz','POz','PO3','PO4','O1','O2'}; % select middle line channels to analyze traveling waves
        for cc = 1:numel(Channels_occipital)
            channel_occi_index(cc) = find(strcmpi(Channels_occipital{cc}, {EEG.chanlocs.labels}));
        end

        % select time interval to analyze
        toi = dsearchn(EEG.times',t1');
        
        % find trial index because some trials were removed due to noise
        total_trial = [];
        for i = 1:numel(EEG.urevent)
            if EEG.urevent(i).type(3)=='1'& EEG.urevent(i).type(4)=='3'
                total_trial = [total_trial;EEG.urevent(i).bvmknum];
            else
                total_trial = total_trial;
            end
            
        end
        
        if sub_index == 1 & sti_index==1
            total_trial = [total_trial;nan(50,1)];
        end
        
        total_trial(:,2) = [1:450];
        total_trial(:,3) = [
            repmat(1,50,1);
            repmat(2,50,1);
            repmat(3,50,1);
            repmat(4,50,1);
            repmat(5,50,1);
            repmat(6,50,1);
            repmat(7,50,1);
            repmat(8,50,1);
            repmat(9,50,1);
            ]
        
        remain_trial = [];
        for i = 1:numel(EEG.event)
            if EEG.event(i).type(3)=='1'& EEG.event(i).type(4)=='3'
                remain_trial = [remain_trial;EEG.event(i).bvmknum];
            else
                remain_trial = remain_trial;
            end
        end
        
        trial_index = total_trial(ismember(total_trial(:,1)',remain_trial),2);
        block_index = total_trial(ismember(total_trial(:,1)',remain_trial),3);
        
        % calculate traveling waves：BW, FW
        data_TW = EEG.data(channel_index,toi(1):toi(2),:);
        for tt = 1:EEG.trials
            [fw(tt),bw(tt)] = quantifyingTW(sub_index,data_TW(:,:,tt),EEG.srate);
        end

        % calculate prestimulus alpha power
        EEG = pop_select( EEG, 'channel', Channels_occipital);
        ap = tfpower_exp3(EEG,t2);

        % save variables of interest
        FW(sub_index,sti_index,trial_index) = fw;
        BW(sub_index,sti_index,trial_index) = bw;
        AP(sub_index,sti_index,trial_index) = ap;

        clear fw bw ap
        
    end
    
end


% arrange data sequence according to stimulation sequence intead of arrival
% sequence
for i = 1:numel(subs)
    FW(i,:,:,:) = FW(i,seqIndex(:,i),:,:);
    BW(i,:,:,:) = BW(i,seqIndex(:,i),:,:);
    AP(i,:,:,:) = AP(i,seqIndex(:,i),:,:);
end

% reshape data
FW = reshape(FW,numel(subs),3,50,9);
BW = reshape(BW,numel(subs),3,50,9);
AP = reshape(AP,numel(subs),3,50,9);

%% process data to statistically analyze in Fig4_5_step2.R

block_num = 9; % because the block just before the stimulation is considered the covariate, so 6 blocks are left.
trial_num = 50;
pre_data = zeros(numel(subs)*3*block_num*trial_num,4); % baseline: block3
data = zeros(numel(subs)*3*block_num*trial_num,4);

% preprocess visual threshold data
threshold_1 = squeeze(threshold(:,1,:,:));
threshold_2 = squeeze(threshold(:,2,:,:));
threshold_3 = squeeze(threshold(:,3,:,:));

threshold_1(threshold_1==0) = NaN;
threshold_2(threshold_2==0) = NaN;
threshold_3(threshold_3==0) = NaN;

pre_threshold_1 = mean(threshold_1(:,:,2),2,'omitnan');
pre_threshold_2 = mean(threshold_2(:,:,2),2,'omitnan');
pre_threshold_3 = mean(threshold_3(:,:,2),2,'omitnan');

% preprocess BW
BW_1 = squeeze(BW(:,1,:,:));
BW_2 = squeeze(BW(:,2,:,:));
BW_3 = squeeze(BW(:,3,:,:));

pre_BW_1 = mean(BW_1(:,:,2),2,'omitnan');
pre_BW_2 = mean(BW_2(:,:,2),2,'omitnan');
pre_BW_3 = mean(BW_3(:,:,2),2,'omitnan');


% preprocess FW
FW_1 = squeeze(FW(:,1,:,:));
FW_2 = squeeze(FW(:,2,:,:));
FW_3 = squeeze(FW(:,3,:,:));

pre_FW_1 = mean(FW_1(:,:,2),2,'omitnan');
pre_FW_2 = mean(FW_2(:,:,2),2,'omitnan');
pre_FW_3 = mean(FW_3(:,:,2),2,'omitnan');

% preprocess AP
AP_1 = squeeze(AP(:,1,:,:));
AP_2 = squeeze(AP(:,2,:,:));
AP_3 = squeeze(AP(:,3,:,:));

pre_AP_1 = mean(AP_1(:,:,2),2,'omitnan');
pre_AP_2 = mean(AP_2(:,:,2),2,'omitnan');
pre_AP_3 = mean(AP_3(:,:,2),2,'omitnan');

% initiate tiredness rating in baseline for covariate
pre_rating = zeros(numel(subs)*3*block_num*trial_num,1);


% calculate baseline tiredness ratings
rating_pre_1 = mean(mean(squeeze(T(:,1,:,1:2)),3),2);
rating_pre_2 = mean(mean(squeeze(T(:,2,:,1:2)),3),2);
rating_pre_3 = mean(mean(squeeze(T(:,3,:,1:2)),3),2);

% prepare covariates for linear mixed model analysis in R
for i = 1:numel(subs)
    for j = 1:3
        switch j
            case 1
            
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = pre_threshold_1(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,2) = pre_BW_1(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,3) = pre_FW_1(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,4) = pre_AP_1(i);
            pre_rating(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = rating_pre_1(i);
            
            case 2
            
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = pre_threshold_2(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,2) = pre_BW_2(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,3) = pre_FW_2(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,4) = pre_AP_2(i);
            pre_rating(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = rating_pre_2(i);
            
            otherwise
            
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = pre_threshold_3(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,2) = pre_BW_3(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,3) = pre_FW_3(i);
            pre_data(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,4) = pre_AP_3(i);
            pre_rating(3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+1:3*block_num*trial_num*(i-1)+trial_num*block_num*(j-1)+trial_num*block_num,1) = rating_pre_3(i);
            
        end
    end
end

% prepare fixed variables for linear mixed model analysis 

% initiate tiredness ratings in each block
block_rating = zeros(numel(subs)*3*block_num*trial_num,1);

for i = 1:numel(subs)
    for j = 1:3
        for k = 1:block_num
            switch j
                case 1
                
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,1) = threshold_1(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,2) = BW_1(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,3) = FW_1(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,4) = AP_1(i,:,k);
                block_rating(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num) = mean(squeeze(T(i,1,:,k)));
                
                case 2
                
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,1) = threshold_2(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,2) = BW_2(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,3) = FW_2(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,4) = AP_2(i,:,k);
                block_rating(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num) = mean(squeeze(T(i,2,:,k)));
                
                otherwise
                
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,1) = threshold_3(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,2) = BW_3(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,3) = FW_3(i,:,k);
                data(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num,4) = AP_3(i,:,k);
                block_rating(3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+1:3*block_num*trial_num*(i-1)+block_num*trial_num*(j-1)+trial_num*(k-1)+trial_num) = mean(squeeze(T(i,3,:,k)));
                
            end
        end
    end
end


%% prepare longData in Fig4_5_step2.R

% subject id
id = repmat([1:numel(subs)]',1,3*block_num*trial_num);
id = reshape(id',[],1);

% stimulation group
sti1 = repmat(1,block_num*trial_num,1);
sti2 = repmat(2,block_num*trial_num,1);
sti3 = repmat(3,block_num*trial_num,1);
group = repmat([sti1;sti2;sti3],numel(subs),1);

% block index
d=[];
for i = 1:9
    d(i,:) = repmat(i,1,trial_num);
end
d = reshape(d',[],1);
blockIndex = repmat(d,numel(subs)*3,1);

interval = zeros(numel(blockIndex),1);
for jj=1:numel(blockIndex)
    if blockIndex(jj)<3
        interval(jj) = 1;
    elseif blockIndex(jj)>5
        interval(jj) = 3;
    else
        interval(jj)=2;
    end
    
end
    
% trial ind
trialIndex = repmat([1:trial_num]',numel(id)/trial_num,1);

% covariates
pre_thres = pre_data(:,1);
pre_BW = pre_data(:,2);
pre_FW = pre_data(:,3);
pre_AP = pre_data(:,4);

% fixed effects
thres = data(:,1);
BW_TW = data(:,2);
FW_TW = data(:,3);
AlphaPower = data(:,4);

% data sequence
seq = [1:numel(id)]';

% generate data in .csv format
dataTable = table(id,group,blockIndex,interval,block_rating,pre_rating,pre_thres,pre_BW,pre_FW,pre_AP,...
    thres,BW_TW,FW_TW,AlphaPower,trialIndex,seq);
writetable(dataTable,filename);
