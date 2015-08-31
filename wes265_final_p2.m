% Joshua Emele <jemele@acm.org>
% WES265, Final, Problem 2
 
% clean up environment before starting
clear all;
%close all;

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

% 2a. Show the impulse response and frequency response of the shaping filter
% and the cascade of the two filters.
figure(1);
subplot(2,1,1);
hold on;
plot(1:length(impz(shape)), shape);
grid on;
axis([1 length(shape) -inf inf]);
title(['Shape Filter, Impulse Response'])
xlabel('Time index');
ylabel('Normalized Amplitude');

subplot(2,1,2);
hold on;
grid on;
plot((-0.5:1/bins:0.5-1/bins)*sps,fftshift(20*log10(abs(fft(shape/sum(shape),bins)))));
axis([-inf inf -100 10]);
title(['Shape Filter, Frequency Response'])
xlabel('Frequency');
ylabel('Normalized Log Magnitude (dB)');

match = conv(shape,shape)/(shape*shape');

figure(2);
subplot(2,1,1);
hold on;
plot(1:length(impz(match)), match);
grid on;
axis([1 length(match) -inf inf]);
title(['Match Filter, Impulse Response'])
xlabel('Time index');
ylabel('Normalized Amplitude');

subplot(2,1,2);
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

figure(3);
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

figure(4);
subplot(2,1,1);
hold on;
plot(matched_data(domain_symbols),'r.');
axis('equal');
axis([-2 2 -2 2])
grid on;
title(['Constellation Diagram, Decoupled Matched Data'])

subplot(2,1,2);
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

figure(5);
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

% 2e. 
shaped_data = shaped_data.*exp(2*pi*j*.01*(0:sample_count-1));
shaped_mult_data = shaped_data.*exp(-2*pi*j*acc*(1:sample_count));

left_output = conv(left,shaped_mult_data);
left_detect = left_output.*conj(left_output);
right_output = conv(right,shaped_mult_data);
right_detect = right_output.*conj(right_output);
detect = right_detect - left_detect;

domain_samples = [ 1:sample_count ];

figure(6);
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

% 2f.
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

figure(7);
subplot(4,1,1);
hold on;
plot(1:sample_count, left_detect(1:sample_count), 'b');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response of Band Edge Filters with Spin')

subplot(4,1,2);
hold on;
plot(1:sample_count, right_detect(1:sample_count), 'r');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response of Band Edge Filters with Spin')

subplot(4,1,3);
hold on;
plot(1:sample_count, detect(1:sample_count), 'b');
grid on;
axis([1 sample_count -inf inf]);
title('Time Response of Band Edge Sum Filters with Spin')

subplot(4,1,4);
hold on
plot(1:sample_count, .01*(1:sample_count), 'b')
plot(1:sample_count, control(1:sample_count), 'r')
axis([1 sample_count -inf inf])
hold off
title('Phase Profile for input Spin and Output of FLL')

% 2g.
