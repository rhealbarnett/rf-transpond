fslow = 1.0e6;
ffast = 1.0e6;

Afast = 1.0;
Aslow = 1.0;

Tfast = 1.0/ffast;
Tslow = 1.0/fslow;
tmax = 7*Tslow;
dt = Tfast/20;

NP = floor(tmax/dt);

tseries_fast = zeros(1,NP);
tseries_slow = zeros(1,NP);

df = 1.0/(tmax);
fmax = 1.0/(2.0*dt);
freq_arr = linspace(1,fmax,fmax/df);


for ii=1:NP
    
    tseries_fast(1,ii) = Afast*cos(2.0*pi*ffast*dt*ii);
    tseries_slow(1,ii) = Aslow*cos(2.0*pi*fslow*dt*ii);
    
end

tseries_comb = tseries_fast.*tseries_slow;

fft_fast = fft(tseries_fast)/NP;
fft_slow = fft(tseries_slow)/NP;
fft_comb = fft(tseries_comb)/NP;


