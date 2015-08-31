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
bins = 4096

% Estimate the filter order using firpmord and the Harris approximation.
% Use the the worst-case value and generate the filter.
% We ignore the (large) w value created, and use our own empirically derived
% value.
[n,fo,ao,w]=firpmord(f,a,dev,fs)
n = max(ceil((fs/(f(2)-f(1)))*(rs/22)),n)
h = firpm(447,fo,ao,[1 4]);
%h = h/max(h);

% 3a. Show the time series and spectrum of the prototype filter. Show the
% filter specifications on the spectrum.
figure(1);
subplot(3,2,1:2);
hold on;
plot(impz(h));
grid on;
title('Remez Filter, Impulse Response');
xlabel('Time Index');
ylabel('Normalized Amplitude');

frequency_response = 20*log10(abs(fftshift(fft(h/sum(h),bins))));
subplot(3,2,3:4);
hold on;
plot(-0.5:1/bins:0.5-1/bins, frequency_response);
grid on;
axis([-0.5 0.5 -150 1]);
title('Remez Filter, Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% plot the filter specifications
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(1)/fs -2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(1)/fs +2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(2)/fs -2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(2)/fs +2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')

subplot(3,2,5);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([-2*fo(2)/3 2*fo(2)/3 -rp/128 rp/128]);
title('Zoom to Passband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% These aren't relevent, since there's almost no passband ripple.
line('XData',[-0.5 0.5],'YData',[+rp/2 +rp/2], 'LineWidth', 1, 'Color','r')
line('XData',[-0.5 0.5],'YData',[-rp/2 -rp/2], 'LineWidth', 1, 'Color','r')

subplot(3,2,6);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([fo(3)/2-0.01 fo(3)/2+0.05 -rs-5 -rs+5]);
title('Zoom to Stopband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(1)/fs -2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(1)/fs +2*f(1)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[-2*f(2)/fs -2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')
line('XData',[+2*f(2)/fs +2*f(2)/fs],'YData',ylim, 'LineWidth', 1, 'Color','r')

% 3b. Show the time series and spectrum pf the derivative filter dh(n) formed
% by convolution with [1 0 -1]/2 phase aligned with h(n) by discarding end
% points.
dh = conv([1 0 -1]/2,h);
dh = dh(2:end-1);

figure(2);
subplot(2,1,1);
hold on;
plot(dh);
grid on;
title('Prototype Derivative Filter, Time Series, Real Part');
xlabel('Time Index');
ylabel('Amplitude');

subplot(2,1,2);
hold on;
plot(-0.5:1/bins:0.5-1/bins,fftshift(20*log10(abs(fft(dh*64,bins)))));
grid on;
title('Derivative Remez Filter, Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% 3c. Partition the time series into a pair of 64-path polyphase filters
% h2(m,n) and dh2(m,n). Write the Matlab script to program their use for a
% 1-to-64 up-sampling filter. Apply 200 samples of a complex sinewave of unit
% amplitude and of input frequency 2.4 KHz, sampled at 20 kHz. Plot the time
% response and the spectrum of the time response and note the frequency of the
% spectrum and the amplitude of any spectral artifacts.
extend = 64-mod(length(h),64);
x = exp(j*2*pi*(2.4/20)*(1:200));
h2 = reshape(h, 64, 7);
dh2 = reshape(dh, 64, 7);

index = 0;
reg = zeros(1:7);
x2 = zeros(1,length(x));
for n=1:length(x)
    reg = [x(n) reg(1:6)];
    for m=1:64
        x2(index+m) = reg*h2(m,:)';
    end
    index = index + 64;
end

figure(3);
subplot(2,1,1);
hold on;
plot(real(x2), 'r');
plot(imag(x2), 'b');
grid on;
axis([10000 inf -inf inf]);
title('Polyphase Filter, Impulse Response, Real and Imag');

w=kaiser(length(x2),20)';
w=20*w/sum(w);
subplot(2,1,2);
hold on;
plot(-0.5:1/(4*bins):0.5-1/(4*bins),20*log10(abs(fftshift(fft(x2.*w,4*bins)))))
grid on;
title('Polyphase Filter, Windowed Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% 3d.
