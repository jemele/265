% Joshua Emele <jemele@acm.org>
% WES265, Final

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
figure(2);
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
figure(3);
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
title(['Eye Diagram, Matched Data , Real Part']);
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

figure(4);
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
figure(5);
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
figure(6);
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

figure(7);
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


