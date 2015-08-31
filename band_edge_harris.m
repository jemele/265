function [left,right] = band_edge_harris(sps,alpha,delay)

% band_edge_harris(sps,alpha,delay), 
% positive freq band edge filter
% baseband band edge filter
% sps samples/symbol, alpha excess bw
% delay is delay in symbols to filter center
% n_len=2*delay+1
%h=rcosine(1,sps,'sqrt',alpha,delay);
%h=h/max(h);
%n=length(h)
n=2*delay*sps+1;

t = (-delay*sps:1:delay*sps);
base = sinc(2*alpha/sps*t-0.5)+sinc(2*alpha/sps*t+0.5);
base = base/sum(base);
m=n/2;

% now translate center to (1+alpha)*fs/2
phi=(1+alpha)*(-m:0.5:m)/(2*sps);
phi=phi(2:2:length(phi));
right=exp(j*2*pi*phi).*base;
left=conj(right);
end
