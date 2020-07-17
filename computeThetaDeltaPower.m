function res=computeThetaDeltaPower(lfp,ts,fs,windowLength,thetaDeltaThreshold)

%function res=plotThetaDeltaPower(lfp,ts,fs,windowLength,LEDts,speed,showPlot)

% Compute theta/delta power ratio
% in 2s windows, which is
% theta/delta power ratio < 4 ~ quiet restfulness (Witton 2014)
%windowLength = 1;  % s, should be 2 in real life
stepIncr = windowLength*fs;  % lfp increments

for i=1:floor((ts(end)-ts(1))/windowLength)
    startIncr=(i-1)*stepIncr+1;
    stopIncr=i*stepIncr;
    ts_segment_means(i) = mean(ts(startIncr:stopIncr));
    lfp_segment_ = lfp(startIncr:stopIncr);
    theta_delta_ratios(i) = computeThetaDeltaPowerRatio(lfp_segment_,fs);
end

% low_theta_delta_power_fraction(i,j)=sum(thetaDeltaPower{i,j}.theta_delta_ratios < theta_delta_power_threshold) /...
%             length(thetaDeltaPower{i,j}.theta_delta_ratios);

time_theta_delta_power_below_threshold = sum(theta_delta_ratios < thetaDeltaThreshold) * windowLength;
low_theta_delta_power_fraction = sum(theta_delta_ratios < thetaDeltaThreshold) /...
            length(theta_delta_ratios);        

res = struct();
res.ts_segment_means = ts_segment_means;
res.theta_delta_ratios = theta_delta_ratios;
res.low_theta_delta_power_fraction = low_theta_delta_power_fraction;
res.duration_theta_delta_power_below_threshold = time_theta_delta_power_below_threshold/60;  % minutes
