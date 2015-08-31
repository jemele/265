% Joshua Emele <jemele@acm.org>
% WES265, Final, Problem 3

% clean up environment before starting
clear all;
close all;

% filter design parameters
% passband peak-to-peak ripple (dB)
rp = .1
% stopband attenuation (dB)
rs = 100
% convert from dB to volt
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)];
% the transition and stop band, in kHz
f = [2.5 17.5]
a = [1 0]
% the sample rate, in kHz
fs = 1280
% the number of fft bins
bins = 2048

% Estimate the filter order using firpmord and the Harris approximation.
% Use the the worst-case value and generate the filter.
[n,fo,ao,w]=firpmord(f,a,dev,fs)
n = max(ceil((fs/(f(2)-f(1)))*(rs/22)),n)
h = firpm(n,fo,ao,w);
h = h/max(h);
frequency_response = 20*log10(abs(fftshift(fft(h,bins))));

% 3a. Show the time series and spectrum of the prototype filter. Show the
% filter specifications on the spectrum.
figure(1);
subplot(3,1,1);
hold on;
plot(h);
axis([0 length(h) -0.5 1.5]);
grid on;
title('Prototype Filter, Time Series, Real Part');
xlabel('Time Index');
ylabel('Normalized Amplitude');

frequency_response = fftshift(20*log10(abs(fft(h/sum(h),bins))));
subplot(3,1,2);
hold on;
plot(-0.5:1/bins:0.5-1/bins,frequency_response);
grid on;
axis([-0.1-fo(3)/2 0.1+fo(3)/2 -110 1]);
title('Prototype Filter, Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% plot the filter specifications
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[+rp/2 +rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[-rp/2 -rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(1)/fs -2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(1)/fs +2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(2)/fs -2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(2)/fs +2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')

subplot(3,2,5);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([-3*fo(2)/2 3*fo(2)/2 -rp/2 rp/2]);
title('Zoom to Passband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[+rp/2 +rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[-rp/2 -rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(1)/fs -2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(1)/fs +2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(2)/fs -2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(2)/fs +2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')

subplot(3,2,6);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([fo(3)/2-0.01 fo(3)/2+0.05 -rs-10 -rs+10]);
title('Zoom to Stopband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[+rp/2 +rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[-rp/2 -rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(1)/fs -2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(1)/fs +2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(2)/fs -2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(2)/fs +2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')

% 3b. Show the time series and spectrum pf the derivative filter dh(n) formed
% by convolution with [1 0 -1]/2 phase aligned with h(n) by discarding end
% points.
dh = conv([1 0 -1]/2,h);
dh = dh(2:length(dh)-1);

figure(2);
subplot(2,1,1);
hold on;
plot(dh);
grid on;
title('Prototype Derivative Filter, Time Series, Real Part');
xlabel('Time Index');
ylabel('Amplitude');

frequency_response = fftshift(20*log10(abs(fft(dh,bins))));
subplot(2,1,2);
hold on;
plot(-0.5:1/bins:0.5-1/bins,frequency_response);
grid on;
title('Prototype Derivative Filter, Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

