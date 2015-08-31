function [output, err] = equalizer(sps,samples,taps,initial_tap)
mu=0.005;
reg=zeros(1,taps);
wts=zeros(1,taps);
wts(initial_tap)=1;
m=1;
for n=1:sps:length(samples)-sps
    reg=[samples(n) reg(1:taps-1)];
    output(n)=reg*wts';
    
    err(m)=0.01;
    if n>=20
        x_det=sign(real(output(n)))+j*sign(imag(output(n)));
        err(m)=x_det-output(n);
        wts=wts+mu*conj(err(m))*reg;
    end

    m=m+1;
    for k=1:1:sps-1
        reg=[samples(n+k) reg(1:taps-1)];
        output(n+k)=reg*wts';
    end
end
end

