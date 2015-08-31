function h4=sqrt_nyq_x(f_smpl,alpha,n_sym,flg_plot)
%function hh=sqrt_nyq_x(f_smpl,alpha,n_sym,flg_plot)
% if flg_plt=1 plots impulse response and spectrum
%
%  hh=sqrt_nyq_x(8,0.25,6,0) 
%  sqrt-nyq filter, alpha=0.25,8-samples/symbol, 
%  6-symbols delay from start to center of filter
%  filter coefficients normalized to unity amplitude for modulator

% compute impulse response
nn=(-n_sym:1/f_smpl:n_sym)+0.0000000001;
h4=(4*alpha*nn).*cos(pi*(1+alpha)*nn)+sin(pi*(1-alpha)*nn);
h4=h4./(pi*nn.*(1-(4*alpha*nn).^2));
h4=h4/max(h4);
h4=(h4+fliplr(h4))/2;

% plot response if required
if flg_plot==1
figure(1)
subplot(2,1,1)
plot(nn,h4,'-o','linewidth',2)
grid on
axis([-n_sym-0.5 n_sym+0.5 -0.3 1.2])
title('Impulse Response: SQRT Nyquist Filter')

text(1.0,0.9,['         \alpha = ',num2str(alpha)])
text(1.0,0.7,['Samples/Symbol = ',num2str(f_smpl)])

subplot(2,1,2)

n_fft=1024;
if length(nn)>n_fft;
    n_fft=2*2^(ceil(log10(length(nn))/log10(2)));
end

fh4=fftshift(abs(fft(h4,n_fft)));
q1=max(fh4(n_fft/2:floor((n_fft/2)*(1+0.3/f_smpl))));
q2=min(fh4(n_fft/2:floor((n_fft/2)*(1+0.3/f_smpl))));
scl=(q1+q2)/2;
plot((-0.5:1/n_fft:0.5-1/n_fft)*f_smpl,(20*log10(fh4/scl)))
hold on

plot([-0.5 -0.5 0.5 0.5],[-60 0 0 -60],'r')
hold off
grid on
axis([-1 1 -60 5])
title('Scaled Spectrum')
xlabel('Normalized Frequency: f/f_S_Y_M_B_O_L')
ylabel('Log Magnitude (dB)')
N=2*n_sym*f_smpl+1;
figure(2)
subplot(2,1,1)
plot(0:N-1,h4,'-o','linewidth',2)
grid on
axis([-2 N+2 -0.3 1.2])
title('SQRT Nyquist Filter Impulse Response')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/1000:0.5-1/1000)*f_smpl,fftshift(20*log10(abs(fft(h4/sum(h4),1000)))),'linewidth',2)
grid on
axis([-f_smpl/2 f_smpl/2 -80 10])
title('SQRT Nyquist Spectrum')
xlabel('Frequency')
ylabel('Magnitude')

figure(3)
hh=conv(h4,h4)/(h4*h4');
subplot(2,1,1)
plot(-(N-1):(N-1),hh,'linewidth',2)
hold on
plot(-(N-1):f_smpl:(N-1),hh(1:f_smpl:2*N),'ro','linewidth',2)
hold off
grid on
axis([-N N -0.3 1.2])
title('SQRT Nyquist Filter Matched Filter Response')
xlabel('Time Index')
ylabel('Amplitude')

hh=conv(h4,h4)/(h4*h4');
subplot(2,1,2)
plot(-(N-1):(N-1),hh,'linewidth',2)
hold on
plot(-(N-1):f_smpl:(N-1),hh(1:f_smpl:2*N),'ro','linewidth',2)
hold off
grid on
axis([-N N -0.001 0.001])
title('Zoom to Zero Crossings Matched Filter Response')
xlabel('Time Index')
ylabel('Amplitude')
end
