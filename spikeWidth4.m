function waveforms = spikeWidth4(wv,dur_Threshold)
%SPIKEWIDTH2 Determine spike duration
% Geneology: spikeWidth.m,spikeduration.m

%clear;
APmaxnumber=2000;  % largest number of APs to consider
APmaxnumber=5000;  % largest number of APs to consider
%dur_Threshold=380e-6;
testmode=0;

%nargin=0;  % TEST
if nargin == 0
dur_Threshold=300e-6;
elseif nargin==1
dur_Threshold=300e-6;
end

clear spike_dur; clear thresh1; clear thresh2;  %CHECK

%  Limit number of waveforms for calculation to prevent huge clusters from
%  requiring huge CPU time
%size(wv,1)
for i=1:min(size(wv,1),APmaxnumber)
    wv1=interp(wv(i,:),5);  % upsample 5x
    wvpeak=min(wv1);  % waveform peak
    wvpeakincr=find(wv1==wvpeak);
    wvpeakincr = min(wvpeakincr);  % prevent case when more than 1 peak
    wv_interp(i,:)=wv1;
    
    if wvpeakincr<0.4*length(wv1)  % otherwise its noise   
        wvbaseline=max(wv1(1:wvpeakincr));
        wvbaselineincr=find(wv1(1:wvpeakincr)==wvbaseline);
        wvbaselineincr=min(wvbaselineincr);  % prevent case when max of wv1 occurs more than once
        
        wvdelta=wvpeak-wvbaseline;
        wvdelta=wvbaseline-wvpeak;
        wv75=wvpeak+0.75*wvdelta;

        z1=abs(wv1(wvbaselineincr:wvpeakincr)-wv75);
        z1=abs(wv1(wvbaselineincr:wvpeakincr)-wvbaseline);
        nearestx1=find(min(z1)==z1);
        nearestx1=min(nearestx1);  % prevent case when min of z1 occurs more than once
        thresh1(i)=nearestx1+wvbaselineincr-1;
        thresh1Z(i)=min(z1);
        
        z2=abs(wv1(wvpeakincr:min(wvpeakincr+2*(wvpeakincr-wvbaselineincr),length(wv1)))-wvbaseline); % assumes waveform ~ symmetrical
        nearestx2=find(min(z2)==z2)+wvpeakincr-1;
        nearestx2=min(nearestx2); % prevent case when min of z2 occurs more than once
        %nearestx2
        thresh2(i)=nearestx2; 
        thresh2Z(i)=min(z2);
        
     %   nearestx2
     %   nearestx1
        spike_dur(i)=(nearestx2-nearestx1)/(40000*5);
        
        % Find trough (Peyrache et al 2015)
        %[wvtrough_val,wvtrough_indx]=max(wv1);  % waveform trough
        [wvtrough_val,wvtrough_indx]=max(wv1(wvpeakincr(1):end));  % waveform trough
        %spike_peak_trough_dur(i) = (wvtrough_indx(1) - wvpeakincr(1))/(40000*5);
        %spike_peak_trough_dur(i) = (wvpeakincr(1)+wvtrough_indx(1) - wvpeakincr(1))/(40000*5);
        spike_peak_trough_dur(i) = (wvtrough_indx(1))/(40000*5);
        
        % Find time from peak to where spike crosses threshold after peak
        %peak2thresh = find(wv1(wvpeakincr(1):end)>wv1(1),1);
        %nearestx2
        peak2thresh = find(wv1(wvpeakincr(1):end)>wv1(nearestx2),1);
        if isempty(peak2thresh)
            peak2thresh = 0;
            %i
            %pause
        end
        spike_peak_threshold2_dur(i) = (peak2thresh)/(40000*5);
        
        % Find threshold to peak excursion
%         wv(1)
%         wvpeak
        threshold2peak_ampl(i)=wv1(1)-wvpeak;
%        pause
        
        if spike_peak_trough_dur(i) < 0
            stop
        end
        
    else
        spike_dur(i)=0;
        thresh1(i)=1; % dummy to prevent error in plot
        thresh2(i)=1;
        spike_peak_trough_dur(i)=0;
        spike_peak_threshold2_dur(i)=0;
        threshold2peak_ampl(i)=0;
    end
end
%spike_dur(1:10)
Percent_less_than_380=100*length(find(spike_dur<350e-6))/length(spike_dur);

spike_dur_non0=spike_dur(find(spike_dur));
spike_dur_mean=mean(spike_dur_non0);
spike_dur_std=std(spike_dur_non0);

% spike_peak_trough_dur(1:10)
% mean(spike_peak_trough_dur(1:10))
% median(spike_peak_trough_dur(1:10))
% mode(spike_peak_trough_dur(1:10))
spike_peak_trough_dur_non0=spike_peak_trough_dur(find(spike_peak_trough_dur));
spike_peak_trough_dur_mean=mean(spike_peak_trough_dur_non0);
spike_peak_trough_dur_std=std(spike_peak_trough_dur_non0);

Percent_less_than_threshold=100*length(find(spike_dur_non0<dur_Threshold))/length(spike_dur_non0);

if Percent_less_than_threshold < 50
	cellType='Pyramidal';
else
	cellType='Interneuron';
end


%function [spike_dur_mean,spike_dur_std,cellType,Percent_less_than_threshold,wv_interp,spike_peak_trough_dur_mean,spike_peak_trough_dur_std,spike_peak_trough_dur,spike_peak_threshold2_dur] = spikeWidth4(wv,dur_Threshold)
waveforms = struct();
waveforms.spike_dur_mean = spike_dur_mean;
waveforms.spike_dur_std = spike_dur_std;
waveforms.cellType = cellType;
waveforms.Percent_less_than_threshold = Percent_less_than_threshold;
waveforms.wv_interp = wv_interp;
waveforms.spike_peak_trough_dur_mean = spike_peak_trough_dur_mean;

waveforms.spike_peak_trough_dur_std = spike_peak_trough_dur_std;
waveforms.spike_peak_trough_dur = spike_peak_trough_dur;
waveforms.spike_peak_threshold2_dur = spike_peak_threshold2_dur;
waveforms.threshold2peak_ampl = threshold2peak_ampl;

 
% Only analyze good waveforms

%mean(wvfms.spike_peak_threshold2_dur(find(wvfms.threshold_peak_dur>mean(wvfms.threshold_peak_dur)-std(wvfms.threshold_peak_dur))))
% find units with threshold to peak value > 1 std less than mean value
%un1=find(wvfms.threshold_peak_dur>mean(wvfms.threshold_peak_dur)-std(wvfms.threshold_peak_dur))

%%try again
% m=mean(threshold2peak_ampl)
% s=std(threshold2peak_ampl)
% threshold2peak_ampl(1:12)
% mean(threshold2peak_ampl)
% 0.5*mean(threshold2peak_ampl)
% % Find waveforms that have significant peak response
% mm=find(threshold2peak_ampl>0.65*mean(threshold2peak_ampl));
% mm(1:10)
% mean(spike_peak_threshold2_dur(mm))
% std(spike_peak_threshold2_dur(mm))
% stop

%%
% Display histogram of spike times    
%figure(1);clf
time_range=[200:10:500];  % microseconds
% spiketime_threshold=380;
% spike_dur_below_thresh=hist(1e6*spike_dur,time_range(time_range<spiketime_threshold))
% spike_dur_above_thresh=hist(1e6*spike_dur,time_range(time_range>=spiketime_threshold))
% spike_barplot=[[spike_dur_below_thresh,zeros(1,(length(time_range)-length(spike_dur_below_thresh)))]'...
%     [zeros(1,(length(time_range)-length(spike_dur_above_thresh))),spike_dur_above_thresh]'];
% %bar(time_range,spike_barplot);

if testmode
	fprintf('channel is %d\n',channel)
	fprintf('unit is %d\n',unit)

[spike_dur_hist,spike_dur_hist_ind]=hist(1e6*spike_dur,time_range);
%spike_dur_hist
spike_dur_below_thresh=spike_dur_hist.*(spike_dur_hist_ind<380);
spike_dur_above_thresh=spike_dur_hist.*(spike_dur_hist_ind>=380);

figure(1);clf
bar(time_range,[spike_dur_below_thresh' spike_dur_above_thresh']);
title('Histogram of spike widths');
xlabel('Spike Width (us)');
ylabel('Frequency');

xt_interp=[1:size(wv_interp,2)]/(40*5);  % time increment, in ms
figure(2);clf
wv_interp_plot=wv_interp(1:5,:);
for j=1:size(wv_interp,1)
    wv_interp_y1(j)=wv_interp(j,thresh1(j));
    wv_interp_y2(j)=wv_interp(j,thresh2(j));
end
thresh1(1:5)/(40*5);
thresh2(1:5)/(40*5);
[(thresh1(1:5)/(40*5))' (thresh2(1:5)/(40*5))'];
[wv_interp_y1' wv_interp_y2'];

plot(xt_interp,wv_interp(1:5,:),'.');
hold on
plot([(thresh1(1:5)/(40*5))' (thresh2(1:5)/(40*5))']',[wv_interp_y1(1:5)' wv_interp_y2(1:5)']');
hold off
grid

figure(3);clf
for i=1:20
    subplot(5,4,i);plot(xt_interp,wv_interp(i,:),'.');
    hold on
    plot([(thresh1(i)/(40*5))' (thresh2(i)/(40*5))']',[wv_interp_y1(i)' wv_interp_y2(i)']');
    hold off
    title(strcat('wv ',num2str(i)));
    L1=legend(strcat(num2str(1e6*spike_dur(i)),'us'));
    set(L1,'FontSize',10,'Location','SouthEast')
    legend('boxoff')
end

end

