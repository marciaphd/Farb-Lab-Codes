function [theta_delta_ratio] = computeThetaDeltaPowerRatio(lfp,fs)

% Compute theta-delta power ratio of LFP

% INPUTS
% lfp -- lfp segment
% fs -- lfp sample rate

if nargin == 0
    fs = 1000;  % sample rate
    fc=3;  % delta frequency 0.0211
    fc=7;  % theta frequency 1703
    t=(1:1000)/1000;  % 1ms samples
    noise = rand(1, length(t));
    lfp=cos(2*pi*fc*t);
    lfp=1*cos(2*pi*3*t)+5*cos(2*pi*7*t)+noise;  % at least 7 cycles, ~ 1 s.
end

n=2^nextpow2(length(lfp));  % next larger power of 2
Y=fft(lfp,n);
f = (0:n-1)*(fs/n);
power = Y.*conj(Y)/n;

% determine power ratio
delta_first_incr = find(f > 1.5,1);
delta_last_incr = find(f > 4,1)-1;
theta_first_incr = find(f > 4,1);
theta_last_incr = find(f > 10,1)-1;
delta_power = sum(power(delta_first_incr:delta_last_incr));
theta_power = sum(power(theta_first_incr:theta_last_incr));
theta_delta_ratio = theta_power / delta_power;
%power
%stop
if nargin == 0
    % plot lfp wrt time
    figure(1);clf
    plot(t,lfp)
    xlabel('Time (s)')
    ylabel('Amplitude')
    grid
    
    % Plot lft segment power spectrum, colors delineate frequency bands
    figure(2);clf
    plot(f(1:delta_first_incr-1),power(1:delta_first_incr-1), '.-c')
    hold on
    plot(f(delta_first_incr:delta_last_incr),power(delta_first_incr:delta_last_incr),'.-')
    plot(f(theta_first_incr:theta_last_incr),power(theta_first_incr:theta_last_incr),'.-')
    %plot(f(theta_last_incr+1:500),power(theta_last_incr+1:500),'k')
    fcTimes10 = find(f > 10*fc,1)-1;  % plot spectrum to 10 times highest freq
    plot(f(theta_last_incr+1:fcTimes10),power(theta_last_incr+1:fcTimes10),'.-k')
    hold off
    grid
    xlabel('Frequency (Hz)')
    ylabel('Power')
    legend('low-freq','delta-freq','theta-freq','high-freq')
    legend boxoff
    text(0.6,0.5,['theta/delta ratio=' num2str(theta_delta_ratio,'%3.1f')],'Units','Normalized','Interpreter','None')
    title('\bf Periodogram for lfp')
end


