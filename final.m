% Joshua Emele <jemele@acm.org>
% WES265, Final, Problem 1

% clean up environment before starting
clear all;
close all;

% design parameters
alpha = 0.3
sps = 4
delay = 10
symbol_count = 2000
bins = 1024

% shaping and match filter
shape = sqrt_nyq_y2(sps, alpha, delay, 0);
shape = shape/max(shape);
shape_frequency_response = fftshift(20*log10(abs(fft(shape/sum(shape),bins))));
match = conv(shape,shape)/(shape*shape');
match_frequency_response = fftshift(20*log10(abs(fft(match/sum(match),bins))));

% 1a. Impulse and frequency response of shape filter.
% Impulse and frequency response of cascade shape and match filter.
figure(11);
subplot(4,1,1);
hold on;
grid on;
plot(-delay:1:delay, shape(1:sps:length(shape)), 'ro');
plot(-delay:1/sps:delay, shape, 'b');
title(['Shape Filter, Impulse Response'])
xlabel('Time Index');
ylabel('Normalized Amplitude');

subplot(4,1,2);
hold on;
grid on;
plot((-0.5:1/bins:0.5-1/bins)*sps,shape_frequency_response);
axis([-2 2 -150 0.1]);
title(['Shape Filter, Frequency Resonse'])
xlabel('Frequency');
ylabel('Normalized Log Magnitude (dB)');

subplot(4,1,3);
hold on;
grid on;
plot(-delay:1:delay, match(1:2*sps:length(match)), 'ro');
plot(-delay:1/(2*sps):delay, match, 'b');
title(['Match Filter, Impulse Resonse'])
xlabel('Time Index');
ylabel('Normalized Amplitude');

subplot(4,1,4);
hold on;
grid on;
plot((-0.5:1/bins:0.5-1/bins)*sps,match_frequency_response);
axis([-2 2 -150 0.1]);
title(['Match Filter, Frequency Resonse'])
xlabel('Normalized Frequency');
ylabel('Normalized Log Magnitude (dB)');

% generate symbols
symbols = [ -1-j -1+j 1-j 1+j ];
data = arrayfun(@(i) symbols(i),randi(length(symbols),symbol_count,1));
upsampled_data = upsample(data, sps);
sample_count = length(upsampled_data);

% shape data, ideal channel, match data
shaped_data = filter(shape, 1, upsampled_data);
matched_data = conv(shaped_data,shape)/(shape*shape');

domain_samples = [ 4*delay+1:1:4*delay+1+400 ];
domain_symbols = [ 4*delay+1:sps:4*delay+1+400 ];

% 1b. Show 400 samples of the time response (real part), the eye diagram, and
% the constellation formed by the shape filter.
figure(12);
subplot(3,1,1);
hold on;
plot(domain_samples, real(shaped_data(domain_samples)));
plot(domain_symbols, real(shaped_data(domain_symbols)),'ro');
axis([min(domain_samples) max(domain_samples) -2 2]);
grid on;
title(['Shaped Data, Real part']);
xlabel('Sample');
ylabel('Amplitude');

subplot(3,1,2);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/4:1,real(shaped_data(n:n+2*sps)),'b');
end
hold off;
grid on;
title(['Eye Diagram, Shaped Data , Real Part']);
ylabel('Amplitude');

subplot(3,1,3);
plot(shaped_data(81:sps:4000),'r.');
grid on;
axis('equal');
axis([-1.5 1.5 -1.5 1.5]);
title(['Constellation Diagram']);

domain_samples = [ 8*delay+1:1:8*delay+1+400 ];
domain_symbols = [ 8*delay+1:sps:8*delay+1+400 ];

% 1c. Show 400 samples of the time response (real part) the eye diagram and the
% constellation of the signal formed at the matched filter output.
figure(13);
subplot(3,1,1);
hold on;
plot(domain_samples, real(matched_data(domain_samples)));
plot(domain_symbols, real(matched_data(domain_symbols)),'ro');
axis([min(domain_samples) max(domain_samples) -2 2]);
grid on;
title(['Matched Data, Real part']);
xlabel('Sample');
ylabel('Amplitude');

subplot(3,1,2);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/sps:1,real(matched_data(n:n+2*sps)),'b');
end
grid on;
title(['Eye Diagram, Matched Data, Real Part']);
ylabel('Amplitude');

subplot(3,1,3);
hold on;
plot(matched_data(domain_symbols),'r.');
axis('equal');
axis([-1.5 1.5 -1.5 1.5])
grid on;
title(['Constellation Diagram, Matched Data'])

domain_samples = [ 4*delay+1:1:4*delay+1+400 ];
domain_symbols = [ 4*delay+1:sps:4*delay+1+400 ];

% 1d. Show 400 samples of the time response (real part), the eye diagram and the
% constellation of the signal received at the channel output.
channel = [1 0 0 0 0 0.3 0 0 j*0.1];
received_data = conv(shaped_data,channel);
matched_data = conv(received_data,shape)/(shape*shape');

figure(14);
subplot(3,1,1);
hold on;
plot(domain_samples, real(received_data(domain_samples)));
plot(domain_symbols, real(received_data(domain_symbols)),'ro');
axis([min(domain_samples) max(domain_samples) -2 2]);
grid on;
title(['Channel Output, Real part']);
xlabel('Sample');
ylabel('Amplitude');

subplot(3,1,2);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/sps:1,real(received_data(n:n+2*sps)),'b');
end
hold off;
grid on;
title(['Eye Diagram, Channel Output, Real Part']);
ylabel('Amplitude');

subplot(3,1,3);
hold on;
plot(received_data(domain_symbols),'r.');
axis('square');
axis([-1.5 1.5 -1.5 1.5]);
grid on;
title(['Constellation Diagram']);

domain_samples = [ 8*delay+1:1:8*delay+1+400 ];
domain_symbols = [ 8*delay+1:sps:8*delay+1+400 ];

% 1e. Show 400 samples of the time response (real part) the eye diagram and the
% constellation of the signal formed at the matched filter output.
figure(15);
subplot(3,1,1);
hold on;
plot(domain_samples, real(matched_data(domain_samples)));
plot(domain_symbols, real(matched_data(domain_symbols)),'ro');
axis([min(domain_samples) max(domain_samples) -2 2]);
grid on;
title(['Matched Output, Real part']);
xlabel('Sample');
ylabel('Amplitude');

subplot(3,1,2);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/sps:1,real(received_data(n:n+2*sps)),'b');
end
hold off;
grid on;
title(['Eye Diagram, Matched Output, Real Part']);
ylabel('Amplitude');

subplot(3,1,3);
hold on;
plot(matched_data(domain_symbols),'r.');
axis('square');
axis([-1.5 1.5 -1.5 1.5]);
grid on;
title(['Constellation Diagram']);

% 1f. Show the three learning curves for the three initial conditions.
figure(16);
taps = 65;
initial_taps = [33 33-24 33+24];
for n=1:length(initial_taps)

    tap=initial_taps(n);
    [output, err] = equalizer(sps,matched_data,taps,tap);

    subplot(length(initial_taps),1,n);
    hold on;
    plot(20*log10(abs(err)));
    plot(20*log10(filter(0.05,[1 -0.95], abs(err))),'r');
    axis([0 symbol_count -60 5]);
    grid on;
    title(['Learning Curve (initial tap ',num2str(tap), ')']);
    xlabel('Symbol');
    ylabel('Error Magnitude (dB)');
end

% 1g. Show 400 samples of the time response (real part) the eye diagram and the
% constellation of the signal formed at the equalizer filter output after
% reaching steady state

% Based on the learning curve, the equalizer is in steady state after 400 samples.
domain_samples = [400+1:1:400+1+400];
domain_symbols = [400+1:sps:400+1+400];

figure(17);
[output, err] = equalizer(sps,matched_data,taps,33);

subplot(3,1,1);
hold on;
plot(domain_samples, real(output(domain_samples)));
plot(domain_symbols, real(output(domain_symbols)),'ro');
axis([min(domain_samples) max(domain_samples) -2 2]);
grid on;
title(['Equalizer Output, Real part']);
xlabel('Sample');
ylabel('Amplitude');

subplot(3,1,2);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/sps:1,real(output(n:n+2*sps)),'b');
end
hold off;
grid on;
title(['Eye Diagram, Equalizer Output, Real Part']);
ylabel('Amplitude');

subplot(3,1,3);
hold on;
plot(output(domain_symbols),'r.');
axis('square');
axis([-1.5 1.5 -1.5 1.5]);
grid on;
title(['Constellation Diagram']);


% Joshua Emele <jemele@acm.org>
% WES265, Final, Problem 2
 
% clean up environment before starting
clear all;
close all;

% design parameters
alpha = 0.5
sps = 4
delay = 8
bins = 2000
symbol_count = 1000
theta_0=2*pi/600

eta = sqrt(2)/2; eta=1.0*eta;
denom = (1+2*eta*theta_0+theta_0*theta_0);
k_i = (4*theta_0*theta_0)/denom;
k_p = (4*eta*theta_0)/denom;

% shape filters
shape = rcosine(1,sps,'sqrt',alpha,delay);
shape = shape/max(shape);
match = conv(shape,shape)/(shape*shape');

% 2a. Show the impulse response and frequency response of the shaping filter
% and the cascade of the two filters.
figure(21);
subplot(2,2,1);
hold on;
plot(1:length(impz(shape)), shape);
grid on;
axis([1 length(shape) -inf inf]);
title(['Shape Filter, Impulse Response'])
xlabel('Time index');
ylabel('Normalized Amplitude');

subplot(2,2,3);
hold on;
grid on;
plot((-0.5:1/bins:0.5-1/bins)*sps,fftshift(20*log10(abs(fft(shape/sum(shape),bins)))));
axis([-inf inf -100 10]);
title(['Shape Filter, Frequency Response'])
xlabel('Frequency');
ylabel('Normalized Log Magnitude (dB)');

subplot(2,2,2);
hold on;
plot(1:length(impz(match)), match);
grid on;
axis([1 length(match) -inf inf]);
title(['Match Filter, Impulse Response'])
xlabel('Time index');
ylabel('Normalized Amplitude');

subplot(2,2,4);
hold on;
plot((-0.5:1/bins:0.5-1/bins)*sps,fftshift(20*log10(abs(fft(match/sum(match),bins)))));
grid on;
axis([-inf inf -150 10]);
title(['Match Filter, Frequency Resonse'])
xlabel('Normalized Frequency');
ylabel('Normalized Log Magnitude (dB)');

% 2b. Form the pair of band edge filters and show on a single page the impulse
% response and frequency response of the filters.
[left,right] = band_edge_harris(sps,alpha,delay);

figure(22);
subplot(3,2,1);
hold on;
plot(1:length(impz(real(left))), real(left), 'b');
grid on;
axis([1 length(left) -inf inf]);
title(['Left Band Edge Filter, Impulse Response, Real Part'])
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,2,3);
hold on;
plot(1:length(impz(imag(left))), imag(left), 'b');
grid on;
axis([1 length(left) -inf inf]);
title(['Left Band Edge Filter, Impulse Response, Imag Part'])
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,2,2);
hold on;
plot(1:length(impz(real(right))), real(right), 'r');
grid on;
axis([1 length(right) -inf inf]);
title(['Right Band Edge Filter, Impulse Response, Real Part'])
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,2,4);
hold on;
plot(1:length(impz(real(left))), imag(right), 'r');
grid on;
axis([1 length(right) -inf inf]);
title(['Left Band Edge Filter, Impulse Response, Imag Part'])
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,2,5:6);
hold on;
plot((-0.5:1/bins:0.5-1/bins)*sps,fftshift(20*log10(abs(fft(left,bins)))),'b');
plot((-0.5:1/bins:0.5-1/bins)*sps,fftshift(20*log10(abs(fft(right,bins)))),'r');
grid on;
axis([-inf inf -100 10]);
title(['Band Edge Filters, Frequency Response'])
xlabel('Normalized Frequency');
ylabel('Normalized Log Magnitude (dB)');

% generate symbols
symbols = [ -1-j -1+j 1-j 1+j ];
data = arrayfun(@(i) symbols(i),randi(length(symbols),symbol_count,1)).';
upsampled_data = upsample(data, sps);
sample_count = length(upsampled_data);
shaped_data = filter(shape, 1, upsampled_data);

% 2c. Form a 1000 symbols of the modulator signal and connect it though the
% multiplier to the demodulator shaped filter and with no input spinning and
% decoupled loop filter form the constellation and eyediagram at the matched
% filter output.
acc = 0;
shaped_mult_data = shaped_data.*exp(-2*pi*j*acc*(1:sample_count));
matched_data = conv(shaped_mult_data,shape)/(shape*shape');

domain_samples = [ 8*delay+1:1:sample_count-8*delay ];
domain_symbols = [ 8*delay+1:sps:sample_count-8*delay ];

figure(23);
subplot(1,4,1);
hold on;
plot(matched_data(domain_symbols),'r.');
axis('equal');
axis([-2 2 -2 2])
grid on;
title(['Constellation Diagram, Decoupled Matched Data'])

subplot(1,4,2:4);
plot(0,0);
hold on;
for n=min(domain_samples):2*sps:max(domain_samples)
plot(-1:1/sps:1,real(matched_data(n:n+2*sps)),'b');
end
grid on;
axis([-1 1 -2.0 2.0]);
title(['Eye Diagram, Decoupled Matched Data, Real Part']);
ylabel('Amplitude');

% 2d. Form and plot the time series at the output of the band edge filters
% detectors and of the output of the summing junction. Plot the spectrum at the
% filter outputs.
left_output = conv(left,shaped_mult_data);
left_detect = left_output.*conj(left_output);
right_output = conv(right,shaped_mult_data);
right_detect = right_output.*conj(right_output);
detect = right_detect - left_detect;

domain_samples = [ 1:sample_count ];

figure(24);
subplot(4,1,1);
hold on;
plot(left_detect(domain_samples),'r');
grid on;
axis([1 200 -inf inf]);
title('Time Response, Left Band Edge Detector');
xlabel('Time Index');
ylabel('Amplitude');

subplot(4,1,2);
hold on;
plot(right_detect(domain_samples),'b');
grid on;
axis([1 200 -inf inf]);
title('Time Response, Right Band Edge Detector');
xlabel('Time Index');
ylabel('Amplitude');

subplot(4,1,3);
hold on;
plot(detect(domain_samples));
grid on;
axis([1 200 -inf inf]);
title('Time Response, Band Edge Detector Sum');
xlabel('Time Index');
ylabel('Amplitude');

w = kaiser(bins,10);
w = 20*w.'/sum(w);
subplot(4,1,4);
hold on;
plot((-0.5:1/bins:0.5-1/bins),fftshift(20*log10(abs(fft(left_output(1:bins).*w,bins)))),'r');
plot((-0.5:1/bins:0.5-1/bins),fftshift(20*log10(abs(fft(right_output(1:bins).*w,bins)))),'b');
grid on;
axis([-inf inf -100 30]);
title('Windowed Frequency Response, Band Edge Filters');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% 2e. Offset the received signal by a spinner of frequency 0.01 per sample.
% With the decoupled loop filter repeat part d.
shaped_data = shaped_data.*exp(2*pi*j*.01*(0:sample_count-1));
shaped_mult_data = shaped_data.*exp(-2*pi*j*acc*(1:sample_count));

left_output = conv(left,shaped_mult_data);
left_detect = left_output.*conj(left_output);
right_output = conv(right,shaped_mult_data);
right_detect = right_output.*conj(right_output);
detect = right_detect - left_detect;

domain_samples = [ 1:sample_count ];

figure(25);
subplot(3,1,1);
hold on;
plot(left_detect(domain_samples),'r');
plot(right_detect(domain_samples),'b');
grid on;
axis([1 200 -inf inf]);
title('Time Response, Band Edge Detectors');
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,1,2);
hold on;
plot(detect(domain_samples));
grid on;
axis([1 200 -inf inf]);
title('Time Response, Band Edge Detector Sum');
xlabel('Time Index');
ylabel('Amplitude');

subplot(3,1,3);
hold on;
plot((-0.5:1/bins:0.5-1/bins),fftshift(20*log10(abs(fft(left_output(1:bins).*w,bins)))),'r');
plot((-0.5:1/bins:0.5-1/bins),fftshift(20*log10(abs(fft(right_output(1:bins).*w,bins)))),'b');
grid on;
axis([-inf inf -100 20]);
title('Band Edge Filters, Windowed Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% 2f. Close the loop of the frequency locked loop. Plot the output time series
% of the two detectors, of the summing junction and of the phase profiles of
% the input spinner phase accumulator and of the frequency lock loop phase
% accumulator.
right_reg = zeros(1,length(right));
left_reg = zeros(1,length(left));
acc = 0;
int = 0;
dphi = 0;

output = zeros(1,sample_count);
left_output = zeros(1,sample_count);
right_output = zeros(1,sample_count);
detect = zeros(1,sample_count);
control = zeros(1,sample_count);

for n = 1:sample_count

    output(n) = shaped_data(n)*exp(-j*2*pi*acc);
    left_reg = [output(n) left_reg(1:end-1)];
    right_reg = [output(n) right_reg(1:end-1)];

    left_output(n) = left_reg * left.';
    left_detect(n) = left_output(n).*conj(left_output(n));
    right_output(n) = right_reg * right.';
    right_detect(n) = right_output(n).*conj(right_output(n));
    detect(n) = right_detect(n)-left_detect(n);

    % loop filter
    int = int + k_i*detect(n);
    dphi = k_p*detect(n)+int;

    % dds control
    acc = dphi + acc;
    control(n) = acc;
end

figure(26);
subplot(4,1,1);
hold on;
plot(1:sample_count, left_detect(1:sample_count), 'b');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response, Left Band Edge Filter with Spin');

subplot(4,1,2);
hold on;
plot(1:sample_count, right_detect(1:sample_count), 'r');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response, Right Band Edge Filter with Spin');

subplot(4,1,3);
hold on;
plot(1:sample_count, detect(1:sample_count), 'b');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response, Band Edge Detector Sum with Spin');

subplot(4,1,4);
hold on;
plot(1:sample_count, .01*(1:sample_count), 'b');
plot(1:sample_count, control(1:sample_count), 'r');
grid on;
axis([1 sample_count -inf inf])
title('DDS Control, FLL Phase Accumulator');

% 2g. Plot the constellations at the input and output of the matched filter for
% the time interval beyond the loops transient response.
matched_output = conv(shape,output)/(shape*shape');

figure(27);
subplot(1,2,1);
hold on;
plot(output(3001:sps:sample_count),'r.');
grid on;
axis('equal');
axis([-2 2 -2 2]);
title('FLL Output, Constellation Diagram');

subplot(1,2,2);
hold on;
plot(matched_output(3001:sps:sample_count),'r.');
grid on;
axis('equal');
axis([-2 2 -2 2]);
title('Matched FLL Output, Constellation Diagram');

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

% 3a. Show the time series and spectrum of the prototype filter. Show the
% filter specifications on the spectrum.
figure(31);
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
axis([-2*fo(2)/3 2*fo(2)/3 -rp/2 rp/2]);
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

figure(32);
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

figure(33);
subplot(2,1,1);
hold on;
plot(real(x2), 'r');
plot(imag(x2), 'b');
grid on;
axis([10000 inf -inf inf]);
title('Polyphase Filter, Impulse Response, Real and Imag');

w = kaiser(length(x2),20)';
w = 20*w/sum(w);
subplot(2,1,2);
hold on;
plot(-0.5:1/(4*bins):0.5-1/(4*bins),20*log10(abs(fftshift(fft(x2.*w,4*bins)))));
grid on;
title('Polyphase Filter, Windowed Frequency Response');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

% 3d. Use an address accumulator initialized at 1 and use an increment of 10
% (up 64 and down 10) to skip output filter ports to obtain a sample rate
% increase of 1-to-6.4, from 20 to 128 kHz . Apply 200 samples of a complex
% sinewave of unit amplitude and of input frequency 2.4 KHz. Plot the time
% response and the spectrum of the time response and note the frequency of the
% spectrum and the amplitude of any spectral artifacts.
% The new sample rate is 1/RR...
RR=20/128;
del=64*RR;
reg=zeros(1,7);
accum=1;
mm=1;                                   % output clock
nn=1;                                   % input clock
while nn<200
    reg=[x(nn) reg(1:6)];               % deliver new input when accum overflows
    while accum<65                      % compute outputs till accum overflows
         pntr=floor(accum);             % integer part of polyphase pointer
         y=reg*h2(pntr,:)';             % compute amplitude estimate
         dy=reg*dh2(pntr,:)';           % compute derivative estimate
         x3(mm)=y+dy*(accum-pntr)*1;    % compute derivative corrected output
         mm=mm+1;                       % Increment output clock
         accum=accum+del;               % increment accumulator
    end
    while accum>=65                     % detect accumulator overflow
         accum=accum-64;                % decrement accumulator 
         if accum>65                    % determine if need new input sample
             nn=nn+1;                   % increment input clock
             reg_2=[x(nn) reg(1:6)];    % deliver new input sample
         end
    end
    nn=nn+1;                            % increment input clock to deliver input
end

figure(34);
subplot(2,1,1);
hold on;
plot(real(x3),'r');
plot(imag(x3),'b');
grid on;
axis([1000 inf -inf inf]);
title('Polyphase Filter, Impulse Response, Real and Imag');
xlabel('Time Index');
ylabel('Amplitude');

w = kaiser(length(x3),20)';
w = 20*w/sum(w);
subplot(2,1,2);
hold on;
plot(2*64*(-0.5:1/(4*bins):0.5-1/(4*bins)),20*log10(abs(fftshift(fft(x3.*w,4*bins)))));
grid on;
axis([-64 64 -180 50]);
title('Polyphase Filter, Windowed Frequency Response');

% 3e. Use the address accumulator initialized at 1 and use an increment of
% 58.1818 (up 64 and down 58.1818) with integer address to use nearest neighbor
% skipping of output filter ports to obtain a sample rate increase of
% 1-to-64/58.1818) or 20*63/58.1818 = 22 kHz, a 10% increase in sample rate.
% Apply 200 samples of a complex sinewave of unit amplitude and of input
% frequency 2.4 KHz. Plot the time response and the spectrum of the time
% response and note the frequency of the spectrum and the amplitude of any
% spectral artifacts.
RR=20/22;
del=64*RR;
reg=zeros(1,7);
accum=1;
mm=1;                                   % output clock
nn=1;                                   % input clock
while nn<200
    reg=[x(nn) reg(1:6)];               % deliver new input when accum overflows
    while accum<65                      % compute outputs till accum overflows
         pntr=floor(accum);             % integer part of polyphase pointer
         y=reg*h2(pntr,:)';             % compute amplitude estimate
         dy=reg*dh2(pntr,:)';           % compute derivative estimate
         x4(mm)=y+dy*(accum-pntr)*1;    % compute derivative corrected output
         mm=mm+1;                       % Increment output clock
         accum=accum+del;               % increment accumulator
    end
    while accum>=65                     % detect accumulator overflow
         accum=accum-64;                % decrement accumulator 
         if accum>65                    % determine if need new input sample
             nn=nn+1;                   % increment input clock
             reg_2=[x(nn) reg(1:6)];    % deliver new input sample
         end
    end
    nn=nn+1;                            % increment input clock to deliver input
end

figure(35);
subplot(2,1,1);
hold on;
plot(real(x4),'r');
plot(imag(x4),'b');
grid on;
axis([100 inf -inf inf]);
title('Polyphase Filter, Impulse Response, Real and Imag');
xlabel('Time Index');
ylabel('Amplitude');

w = kaiser(length(x4),20)';
w = 20*w/sum(w);
subplot(2,1,2);
hold on;
plot(2*64*(-0.5:1/(4*bins):0.5-1/(4*bins)),20*log10(abs(fftshift(fft(x4.*w,4*bins)))));
grid on;
axis([-64 64 -180 50]);
title('Polyphase Filter, Windowed Frequency Response');

% 3f. Repeat part e except do not use the derivative filter (multiply  by 0).
% Plot the time response and the spectrum of the time response and note and
% comment on the frequency of the spectrum and the amplitude of any spectral
% artifacts.
RR=20/22;
del=64*RR;
reg=zeros(1,7);
accum=1;
mm=1;                                   % output clock
nn=1;                                   % input clock
while nn<200
    reg=[x(nn) reg(1:6)];               % deliver new input when accum overflows
    while accum<65                      % compute outputs till accum overflows
         pntr=floor(accum);             % integer part of polyphase pointer
         y=reg*h2(pntr,:)';             % compute amplitude estimate
         x5(mm)=y;                      % omit derivative corrected output
         mm=mm+1;                       % Increment output clock
         accum=accum+del;               % increment accumulator
    end
    while accum>=65                     % detect accumulator overflow
         accum=accum-64;                % decrement accumulator 
         if accum>65                    % determine if need new input sample
             nn=nn+1;                   % increment input clock
             reg_2=[x(nn) reg(1:6)];    % deliver new input sample
         end
    end
    nn=nn+1;                            % increment input clock to deliver input
end

figure(36);
subplot(2,1,1);
hold on;
plot(real(x5),'r');
plot(imag(x5),'b');
grid on;
axis([100 inf -inf inf]);
title('Polyphase Filter without Derivative Corrected Output, Impulse Response, Real and Imag');
xlabel('Time Index');
ylabel('Amplitude');

w = kaiser(length(x5),20)';
w = 20*w/sum(w);
subplot(2,1,2);
hold on;
plot(2*64*(-0.5:1/(4*bins):0.5-1/(4*bins)),20*log10(abs(fftshift(fft(x5.*w,4*bins)))));
grid on;
axis([-64 64 -180 50]);
title('Polyphase Filter without Derivative Corrected Output, Windowed Frequency Response');
ylabel('Log Magnitude (dB)');


