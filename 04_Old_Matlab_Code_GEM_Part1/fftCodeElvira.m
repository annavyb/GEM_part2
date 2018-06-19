%% Compute FFT
data = ALLEEG(1,2).data;
winsize = ALLEEG(1,2).srate;
windows = [1:2*winsize:size(data,2)];
for ch2 = 1:length(all_channel)
    for wn = 1:length(windows)-1
        data_cut = data(ch2,windows(wn):windows(wn+1));
        f =  winsize*(0:(size(data_cut,2)/2))/size(data_cut,2);
        Y = fft(data_cut);
        P2 = abs(Y/size(data_cut,2));
        P1 = P2(:,1:size(data_cut,2)/2+1);
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        Y_all(:,wn) = P1;
    end
    Y_mean(ch2,:) = mean(Y_all,2);
end
 
f_plot = find(f>100,1);
figure
hold on
for ch2 = 1:length(all_channel)
    ind = intersect(ch2,ind_gel);
    if isempty(ind)
        plot(f(1:f_plot),Y_mean(ch2,1:f_plot), 'k')
    elseif ind == 2
        plot(f(1:f_plot),Y_mean(ch2,1:f_plot), 'r')
    elseif ind == 7
        plot(f(1:f_plot),Y_mean(ch2,1:f_plot), 'c')
    elseif ind == 8
        plot(f(1:f_plot),Y_mean(ch2,1:f_plot), 'b')
    end
    xlim([-2 f(f_plot+1)])
end