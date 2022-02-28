
%output = surrogates(x) returns a 
%phase-randomized copy of x

function output = surrogates(x)

%length of x
N = length(x);

%apply DFT
s = fft(x); 
%plot(abs(s),'g'), length(s)
%construct (conjugate symmetric) array of random phases
phase_rnd = zeros(1,N);
%define first phase
phase_rnd(1) = 0;
if (odd(N) == 1);
    ph = 2*pi.*rand(1,(N-1)/2) - pi; 
    phase_rnd(2:N) = [ph,-flipdim(ph,2)];
end
if (odd(N) == 0);
    ph(1:(N-2)/2) = 2*pi.*rand(1,(N-2)/2) - pi;
    phase_rnd(2:N) = [ph,0,-flipdim(ph,2)];
end
%randomize phases
s_rnd = zeros(1,N); %initialization
s_rnd = abs(s).*exp(i.*(angle(s) + phase_rnd));
s_rnd(1) = s(1); %to make sure that s_rnd(1) is real 

%apply inverse DFT
x_rnd = ifft(s_rnd,'symmetric'); %use "symmetric" because of round-off errors

%define output
output = x_rnd;

%---------------------------------------------------------
% o = odd(n) equals 1 if n is odd, 0 if n is even
function outp = odd(n);

for i = 1:round(n/2);
    if (2*i == n);
        outp = 0;
    else 
        outp = 1;
    end
end

if (n == 0);
    outp = 0;
end
end
%---------------------------------------------------------
%plots
% subplot(2,2,1);
% plot(x,'k'), title('ORIGINAL SIGNAL');
% subplot(2,2,2);
% plot(x_rnd,'k'), title('PHASE RANDOMIZED SIGNAL');
% subplot(2,2,3);
% plot(angle(s)), title('ORIGINAL PHASES');
% subplot(2,2,4);
% plot(angle(s_rnd)), title('RANDOMIZED PHASES');

end

