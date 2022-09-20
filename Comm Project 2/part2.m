clc
close all
clear all

bits = randi([0,1],1,10000);
p_bits = [];
for i = 1:length(bits)
    if bits(i) == 0
        p_bits = [p_bits -1];
    elseif bits(i) == 1
        p_bits = [p_bits 1];
    end
end
p = [5 4 3 2 1]/sqrt(55);
p_bits_sampled = upsample(p_bits,length(p));
yn = conv(p_bits_sampled,p);%output of the transmitter
yn = yn(1:50000);
MF = fliplr(p);
corr = [5 5 5 5 5]/sqrt(125);

snr_db = -2:1:5;
snr = 10.^(snr_db./10);
N0 = 1./snr;
 
BER_MF = [];
BER_corr = [];

for i =1:length(N0)
Nn = randn(1,length(yn)).*sqrt(N0(i)/2);
Vn = Nn + yn;

%at the receiver
MF_ = conv(Vn,MF);
MF_ = MF_(1:50000);

corr_ = conv(Vn,corr);
corr_ = corr_(1:50000);

xn = [];
for i = 1:10000
    xn = [xn MF_(i*5)];
end
yn_hat_MF =[];
for i = 1:length(xn)
    if xn(i) < 0
        yn_hat_MF = [yn_hat_MF 0];
    elseif xn(i) > 0
        yn_hat_MF = [yn_hat_MF 1];
    end
end

count = 0;
for i = 1:length(yn_hat_MF)
    if bits(i) ~= yn_hat_MF(i)
        count = count+1;
    end
end
BER_MF = [BER_MF count/length(bits)];


xn_corr = [];
for i = 1:10000
    xn_corr = [xn_corr corr_(i*5)];
end

%decision
yn_hat_corr =[];
for i = 1:length(xn_corr)
    if xn_corr(i) < 0
        yn_hat_corr = [yn_hat_corr 0];
    elseif xn_corr(i) > 0
        yn_hat_corr = [yn_hat_corr 1];
    end
end

count = 0;
for i = 1:length(yn_hat_corr)
    if bits(i) ~= yn_hat_corr(i)
        count = count+1;
    end
end
BER_corr = [BER_corr count/length(bits)];

end

BER_theoretical = (1/2).*erfc(sqrt(1./N0));

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------
figure(1)
semilogy(snr_db,BER_MF)
hold on 
semilogy(snr_db,BER_corr)
hold on 
semilogy(snr_db,BER_theoretical)
legend('BER MF','BER Rect Filter','BER theoretical')
ylabel('BER')
xlabel('SNR (dB)')
grid on
