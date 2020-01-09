

dk = 1.0/((npts)*dx);
knyq = 1.0/(2.0*dx);
k_np = (npts/2);
k_ax = zeros(1,k_np);

for ii=1:k_np
    k_ax(1,ii) = dk*(ii-1);
end

fft_sig = (fft(imag(rf_ey)));

plot(k_ax,abs(fft_sig(1:npts/2))/npts)