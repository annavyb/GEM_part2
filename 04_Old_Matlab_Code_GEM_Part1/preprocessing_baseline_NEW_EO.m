function preprocessing_baseline_NEW_EO

clear all
close all
c_dir=cd;


%% Future input arguments and global variables

w_dir1='E:\Data\GetReconnect\Patients\p121\rest_d1\Evening\Eyes_open';

sr=1024;
sr_new=128;
low_cf=1;
high_cf=40;
filt_order=8;
rule_ICA=20;
subject = {'p121'};
day = {'d1'};
moment = {'e'};
blocks={'p121_reo_d1_e_s11.bdf'};
n_blocks=length(blocks);
namefolder=['E:\Data\GetReconnect\Patients\p121\rest_d1\Evening\Eyes_open\Proc'];

%% Preprocessing

% Data loading
cd(w_dir1);
data_conc=[];
for ii=1:n_blocks
    disp(['Loading data... ' blocks{ii}]);
    data=pop_biosig(blocks{ii});
    data_conc=[data_conc data.data];
end
base_original=pop_importdata('setname',['Sub' subject ...
    '_resting_original'],'data',data_conc,'dataformat','array','srate',...
    data.srate,'nbchan',data.nbchan,'chanlocs',data.chanlocs);

% Add channel position
cd(c_dir)
load('chanlocs.mat')
n_channel_activity = 32;
for ch1 = 1:n_channel_activity
    for ch2 = 1:size(EEG.chanlocs,2)
        if strcmp(base_original.chanlocs(ch1,1).labels,EEG.chanlocs(1,ch2).labels)
            base_original.chanlocs(ch1,1).X = EEG.chanlocs(1,ch2).X;
            base_original.chanlocs(ch1,1).Y = EEG.chanlocs(1,ch2).Y;
            base_original.chanlocs(ch1,1).Z = EEG.chanlocs(1,ch2).Z;
            
            base_original.chanlocs(ch1,1).theta = EEG.chanlocs(1,ch2).theta;
            base_original.chanlocs(ch1,1).radius = EEG.chanlocs(1,ch2).radius;
            base_original.chanlocs(ch1,1).sph_theta = EEG.chanlocs(1,ch2).sph_theta;
            base_original.chanlocs(ch1,1).sph_phi = EEG.chanlocs(1,ch2).sph_phi;
            base_original.chanlocs(ch1,1).sph_phi_besa = EEG.chanlocs(1,ch2).sph_phi_besa;
            base_original.chanlocs(ch1,1).sph_theta_besa = EEG.chanlocs(1,ch2).sph_theta_besa;
            base_original.chanlocs(ch1,1).type = EEG.chanlocs(1,ch2).type;
        end
    end
end

mkdir(namefolder)
cd(namefolder);
save('base_original.mat','base_original','-v7.3');

% Channels removel
bad_chans=[41:64];
base_ch_removed1=base_original;
good_chans=setdiff(1:40,bad_chans);
n_channels=length(good_chans);
base_ch_removed1.data=base_ch_removed1.data(good_chans,:);
base_ch_removed1.nbchan=n_channels;
base_ch_removed1.chanlocs=base_ch_removed1.chanlocs(good_chans);
base_ch_removed1.setname=['Sub' subject '_resting_ch_removed1'];
save('base_ch_removed1.mat','base_ch_removed1');

% Data filtering
disp('Filtering data...');
base_filtered=base_ch_removed1;
base_filtered.data = double(base_filtered.data);
base_filtered.data=base_filtered.data-repmat(mean(base_filtered.data,2),...
    [1 size(base_filtered.data,2)]);
h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr);
d2=design(h2,'Butter');
freqz(d1)
freqz(d2)
for nn=1:n_channel_activity
    base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
       base_filtered.data(nn,:));
    base_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
        double(base_filtered.data(nn,:)));
end
for nn=n_channel_activity+1:size(base_ch_removed1.data,1)
    base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
       base_filtered.data(nn,:));
end
base_filtered.setname=['Sub' subject '_resting_filtered'];
save('base_filtered.mat','base_filtered','-v7.3');
clear base_ch_removed1

% Data resampling
disp('Resampling data...');
base_resampled=pop_resample(base_filtered,sr_new);        
base_resampled.setname=['Sub' subject '_resting_resampled'];          
save('base_resampled.mat','base_resampled','-v7.3');
clear base_filtered

% Channels visual removal
disp('Channels visual removal...');
base_ch_removed2=base_resampled;
eegplot(base_ch_removed2.data(1:n_channel_activity,:),'srate',sr_new);
bad_chans=input('Bad channels (vector): ');
good_chans=setdiff(1:n_channel_activity,bad_chans);
n_channels=length(good_chans);
base_ch_removed2.data=base_ch_removed2.data(good_chans,:);
base_ch_removed2.nbchan=n_channels;
base_ch_removed2.chanlocs=base_ch_removed2.chanlocs(good_chans);
base_ch_removed2.setname=['Sub' subject '_resting_ch_removed2'];
save('base_ch_removed2.mat','base_ch_removed2');

% Common average re-referencing
disp('Re-referencing data...');
base_rereferenced=base_ch_removed2;
base_rereferenced.data=base_rereferenced.data-repmat(mean...
    (base_rereferenced.data,1),[size(base_rereferenced.data,1) 1]);
base_rereferenced.setname=['Sub' subject '_resting_rereferenced'];
save('base_rereferenced.mat','base_rereferenced');

% Run ICA
disp('Running Infomax ICA...');
tmpdata=base_rereferenced.data-repmat(mean(base_rereferenced.data,2),...
    [1 size(base_rereferenced.data,2)]);
n_comp=min([rank(tmpdata),floor(sqrt((size(base_rereferenced.data,2)/...
    rule_ICA)))]);
base_ICA=pop_runica(base_rereferenced,'icatype','runica','extended',1,'pca',...
    n_comp,'stop',1e-7);
base_ICA.setname=['Sub' subject '_resting_ICA'];
save('base_ICA.mat','base_ICA');

cond_num=cond(base_ICA.icaweights);
disp(['Cond. number -> ' num2str(cond_num)]);

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
base_IC_visual_removal1=base_ICA;
pop_topoplot(base_IC_visual_removal1,0,1:n_comp,'Independent Components');
eegplot(eeg_getdatact(base_IC_visual_removal1,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=base_IC_visual_removal1.icawinv(:,bad_IC)*...
        eeg_getdatact(base_IC_visual_removal1,'component',bad_IC,...
        'reshape','2d');
    base_IC_visual_removal1.data=base_IC_visual_removal1.data-projection;
    good_IC=setdiff(1:size(base_IC_visual_removal1.icaweights,1),bad_IC);
    base_IC_visual_removal1.icawinv=...
        base_IC_visual_removal1.icawinv(:,good_IC);
    base_IC_visual_removal1.icaweights=...
        base_IC_visual_removal1.icaweights(good_IC,:);
end
base_IC_visual_removal1.reject=[];
base_IC_visual_removal1.setname=['Sub' subject '_resting_IC_visual_removal1'];
save('base_IC_visual_removal1.mat','base_IC_visual_removal1');

channels_eyes = [34 35 36 37 40];
base_IC_visual_removal1.data(end+1:end+5,:) = base_resampled.data(channels_eyes,:);
% Identify part where the patient opens the eye
TMPREJ = [];
eegplot(base_IC_visual_removal1.data,'srate',sr_new,'command','close');
rm_data=[];
for nn=1:size(TMPREJ,1)
    rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
end

% Region of interest selection
base_final=base_IC_visual_removal1;
base_final.data(:,rm_data)= [];
base_final.pnts=size(base_final.data,2);
base_final.xmax=(base_final.pnts-1)/sr_new;
base_final.times=linspace(base_final.xmin,base_final.xmax,...
    size(base_final.data,2));
eegplot(base_final.data,'srate',sr_new,'command','close');
save('base_final.mat','base_final');
save('rm_data', 'rm_data')

%% cacluate how many block of 1 minute I have for test and re-test
n_blk = floor(size(base_final.data,2)/(sr_new*60));

for blk = 1:n_blk
    base_cartool = base_final;
    base_cartool.data = base_final.data(:,(blk-1)*(sr_new*60)+1:blk*(sr_new*60));
    base_cartool.data(end-4:end,:)=[];
    eegplot(base_cartool.data,'srate',sr_new,'command','close');
    base_cartool.pnts=size(base_cartool.data,2);
    base_cartool.xmax=(base_cartool.pnts-1)/sr_new;
    base_cartool.times=linspace(base_cartool.xmin,base_cartool.xmax,...
        size(base_cartool.data,2));
    base_cartool.setname=['Sub' subject '_resting_scrolled'];
    save(['base_cartool_' num2str(blk) '.mat'],'base_cartool');
    pop_writeeeg(base_cartool,[subject{1} '_reo_' day{1} '_' moment{1} '_piece_' num2str(blk) '.bdf'],'TYPE','BDF');
    
    if blk ~= n_blk
        clear base_cartool
    end
end

% CARTOOL
output=fopen([subject{1} '_reo_' day{1} '_' moment{1} '_electrodes.xyz'],'w');
len=length(base_cartool.chanlocs);
fprintf(output,'%u\t%u',len,1);
for nn=1:len
    fprintf(output,'\n%f\t%f\t%f\t%s',-base_final.chanlocs(nn).X,...
        -base_final.chanlocs(nn).Y,base_final.chanlocs(nn).Z,...
        base_final.chanlocs(nn).labels);
end
fclose(output);
% LORETA

cd(c_dir);