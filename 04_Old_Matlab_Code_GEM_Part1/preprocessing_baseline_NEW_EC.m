function preprocessing_baseline_NEW_EC

clear all
close all
c_dir=cd;


%% Future input arguments and global variables

w_dir1='E:\Data\GetReconnect\Patients\p112\D2\Evening\Eyes_close';

sr=1024;
sr_new=128;
low_cf=1;
high_cf=40;
filt_order=8;
subject = {'p112'};
day = {'d2'};
moment = {'6'};
blocks={'p112_rec_d2_e_s26.bdf'};
n_blocks=length(blocks);
namefolder=['E:\Data\GetReconnect\Patients\p112\D2\Evening\Eyes_close\Proc'];

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
            base_original.chanlocs(ch1,1).sph_radius = EEG.chanlocs(1,ch2).sph_radius;
            base_original.chanlocs(ch1,1).type = EEG.chanlocs(1,ch2).type;
        end
    end
end

mkdir(namefolder)
cd(namefolder);
save('base_original.mat','base_original','-v7.3');

% Channels removel
bad_chans=[41:size(base_original.data,1)];
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
good_chans=setdiff(1:size(base_ch_removed2.data,1),bad_chans);
n_channels=length(good_chans);
base_ch_removed2.data=base_ch_removed2.data(good_chans,:);
base_ch_removed2.nbchan=n_channels;
base_ch_removed2.chanlocs=base_ch_removed2.chanlocs(good_chans);
base_ch_removed2.setname=['Sub' subject '_resting_ch_removed2'];
eegplot(base_ch_removed2.data(1:n_channel_activity-length(bad_chans),:),'srate',sr_new);
save('base_ch_removed2.mat','base_ch_removed2');

% Channel interpolation
base_ch_interp = base_resampled;
base_ch_interp = pop_interp(base_ch_interp,bad_chans,'spherical');
base_ch_interp.setname=['Sub' subject '_inter'];
eegplot(base_ch_interp.data(1:n_channel_activity,:),'srate',sr_new);
save('base_ch_interp.mat','base_ch_interp');

% Common average re-referencing
disp('Re-referencing data...');
base_rereferenced=base_ch_interp;
base_rereferenced.data(1:n_channel_activity,:)=base_rereferenced.data(1:n_channel_activity,:)-repmat(mean...
    (base_rereferenced.data(1:n_channel_activity,:),1),[size(base_rereferenced.data(1:n_channel_activity,:),1) 1]);
base_rereferenced.setname=['Sub' subject '_resting_rereferenced'];
eegplot(base_rereferenced.data(1:n_channel_activity,:),'srate',sr_new);
save('base_rereferenced.mat','base_rereferenced');
channels_eyes = [1:32 34 35 36 37 40];

% Identify part where the patient opens the eye
TMPREJ = [];
eegplot(base_rereferenced.data(channels_eyes,:),'srate',sr_new,'command','close');
rm_data=[];
for nn=1:size(TMPREJ,1)
    rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
end

% Region of interest selection
base_final=base_rereferenced;
base_final.data(:,rm_data)= [];
base_final.pnts=size(base_final.data,2);
base_final.xmax=(base_final.pnts-1)/sr_new;
base_final.times=linspace(base_final.xmin,base_final.xmax,...
    size(base_final.data,2));
eegplot(base_final.data(1:n_channel_activity,:),'srate',sr_new,'command','close');
save('base_final.mat','base_final');
save('rm_data', 'rm_data')

%% cacluate how many block of 1 minute I have for test and re-test
% n_blk = round(size(base_final.data,2)/(sr_new*60));
% 
% for blk = 1:n_blk
%     base_cartool = base_final;
%     base_cartool.data = base_final.data(:,(blk-1)*(sr_new*60)+1:blk*(sr_new*60));
%     base_cartool.data(end-4:end,:)=[];
%     eegplot(base_cartool.data,'srate',sr_new,'command','close');
%     base_cartool.pnts=size(base_cartool.data,2);
%     base_cartool.xmax=(base_cartool.pnts-1)/sr_new;
%     base_cartool.times=linspace(base_cartool.xmin,base_cartool.xmax,...
%         size(base_cartool.data,2));
%     base_cartool.setname=['Sub' subject '_resting_scrolled'];
%     save(['base_cartool_' num2str(blk) '.mat'],'base_cartool');
%     pop_writeeeg(base_cartool,[subject{1} '_rec_' day{1} '_' moment{1} '_piece_' num2str(blk) '.bdf'],'TYPE','BDF');
%     
%     if blk ~= n_blk
%         clear base_cartool
%     end
% end

% CARTOOL
% output=fopen([subject{1} '_rec_' day{1} '_' moment{1} '_electrodes.xyz'],'w');
% len=length(base_cartool.chanlocs);
% fprintf(output,'%u\t%u',len,1);
% for nn=1:len
%     fprintf(output,'\n%f\t%f\t%f\t%s',-base_final.chanlocs(nn).X,...
%         -base_final.chanlocs(nn).Y,base_final.chanlocs(nn).Z,...
%         base_final.chanlocs(nn).labels);
% end
% fclose(output);
% LORETA

cd(c_dir);