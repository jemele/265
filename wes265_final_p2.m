% Joshua Emele <jemele@acm.org>
% WES265, Final, Problem 2
 
% clean up environment before starting
clear all;
close all;

% design parameters
sps = 4
alpha = 0.5
delay = 8
bins = 1024
symbol_count = 1000
theta_0=2*pi/600

% shape and match filters
shape = rcosine(1,sps,'sqrt',alpha,delay);
shape = shape/max(shape);
shape_frequency_response = fftshift(20*log10(abs(fft(shape/sum(shape),bins))));
match = conv(shape,shape)/(shape*shape');
match_frequency_response = fftshift(20*log10(abs(fft(match/sum(match),bins))));

% 2a. Show the impulse response and frequency response of the shaping filter
% and the cascade of the two filters.
figure(1);
subplot(4,1,1);
hold on;
grid on;
plot(-delay:1:delay, shape(1:sps:length(shape)), 'ro');
plot(-delay:1/sps:delay, shape, 'b');
title(['Shape Filter, Impulse Resonse'])
xlabel('Time index');
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
xlabel('Time index');
ylabel('Normalized Amplitude');

subplot(4,1,4);
hold on;
grid on;
plot((-0.5:1/bins:0.5-1/bins)*sps,match_frequency_response);
axis([-2 2 -150 0.1]);
title(['Match Filter, Frequency Resonse'])
xlabel('Normalized Frequency');
ylabel('Normalized Log Magnitude (dB)');

% 2b. Form the pair of band edge filters and show on a single page the impulse
% response and frequency response of the filters.
[left,right] = band_edge_harris(sps,alpha,delay);
left_frequency_response = fftshift(20*log10(abs(fft(left,bins))));
right_frequency_response = fftshift(20*log10(abs(fft(right,bins))));
shape_frequency_response = fftshift(20*log10(abs(fft(shape/sum(shape),bins))));

figure(2);
hold on;
plot((-0.5:1/bins:0.5-1/bins)*sps,left_frequency_response,'r');
plot((-0.5:1/bins:0.5-1/bins)*sps,right_frequency_response,'r');
plot((-0.5:1/bins:0.5-1/bins)*sps,shape_frequency_response,'b');
axis([-2 2 -60 10]);
grid on;
title(['Band Edge Filter, Frequency Response'])
xlabel('Frequency');
ylabel('Normalized Log Magnitude (dB)');

% generate symbols
symbols = [ -1-j -1+j 1-j 1+j ];
data = arrayfun(@(i) symbols(i),randi(length(symbols),symbol_count,1));
upsampled_data = upsample(data, sps);
sample_count = length(upsampled_data);
shaped_data = filter(shape, 1, upsampled_data);
matched_data = conv(shaped_data,shape)/(shape*shape');

domain_samples = [ 8*delay+1:1:length(matched_data)-8*delay ];
domain_symbols = [ 8*delay+1:sps:length(matched_data)-8*delay ];

% 2c. Form a 1000 symbols of the modulator signal and connect it though the
% multiplier to the demodulator shaped filter and with no input spinning and
% decoupled loop filter form the constellation and eyediagram at the matched
% filter output.
figure(3);
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

% 2d. Form and plot the time series at the output of the band edge filters
% detectors and of the output of the summing junction. Plot the spectrum at the
% filter outputs.
left_output = filter(left,1,shaped_data);
right_output = filter(right,1,shaped_data);


% fll loop filter parameters
eta = sqrt(2)/2; eta=1.0*eta;
denom = (1+2*eta*theta_0+theta_0*theta_0);
k_i = (4*theta_0*theta_0)/denom;
k_p = (4*eta*theta_0)/denom;


