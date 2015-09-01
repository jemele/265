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
fo = [0 2.5 17.5 640]/640
bins = 4096

% The Harris approximation yields N=389, but since we know we are going to
% upsample by 64, we choose the next largest filter with order divisible by 64. 
% The w value was empirically determined and appears to be the smallest
% admissisble value.
N = 447
ao = {'myfrf', [1 1 0 0]}
w = [1 2]
h = firpm(N,fo,ao,w);

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
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'Color','r')
line('XData',[-fo(2) -fo(2)],'YData',ylim, 'Color','r')
line('XData',[+fo(2) +fo(2)],'YData',ylim, 'Color','r')
line('XData',[-fo(3) -fo(3)],'YData',ylim, 'Color','r')
line('XData',[+fo(3) +fo(3)],'YData',ylim, 'Color','r')

subplot(3,2,5);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([-2*fo(2)/3 2*fo(2)/3 -rp/2 rp/2]);
title('Zoom to Passband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');

subplot(3,2,6);
hold on;
plot((-0.5:1/bins:0.5-1/bins), frequency_response);
grid on;
axis([fo(3)/2-0.01 fo(3)/2+0.05 -rs-5 -rs+5]);
title('Zoom to Stopband Ripple');
xlabel('Normalized Frequency');
ylabel('Log Magnitude (dB)');
line('XData',[-0.5 0.5],'YData',[-rs -rs], 'Color','r')
line('XData',[-fo(2) -fo(2)],'YData',ylim, 'Color','r')
line('XData',[+fo(2) +fo(2)],'YData',ylim, 'Color','r')
line('XData',[-fo(3) -fo(3)],'YData',ylim, 'Color','r')
line('XData',[+fo(3) +fo(3)],'YData',ylim, 'Color','r')

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

figure(4);
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

figure(5);
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

figure(6);
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


