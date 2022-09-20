clc
close all
clear all

bits = randi([0,1],1,10);
p_bits = [];
for i = 1:length(bits)
    if bits(i) == 0
        p_bits = [p_bits -1];
    elseif bits(i) == 1
        p_bits = [p_bits 1];
    end
end

p = [5 4 3 2 1]/sqrt(55);
p_bits_sampled = upsample(p_bits,5);
yn = conv(p_bits_sampled,p);%output of the transmitter
yn=yn(1:50);
%at the receiver
%--------------------------------------------------------------------------
%Matched filter
%--------------------------------------------------------------------------

MF = fliplr(p);
MF_ = conv(MF,yn);
MF_ = MF_(1:50);
xn_MF =[];
for i = 1:10
    xn_MF = [xn_MF MF_(i*5)];
end

%decision
yn_hat_MF =[];
for i = 1:length(xn_MF)
    if xn_MF(i) <= 0
        yn_hat_MF = [yn_hat_MF -1];
    elseif xn_MF(i) > 0
        yn_hat_MF = [yn_hat_MF 1];
    end
end

%--------------------------------------------------------------------------
%normalized filter
%--------------------------------------------------------------------------

nor = [1 1 1 1 1]/sqrt(5);
nor_ = conv(nor,yn);
nor_=nor_(1:50);
xn_nor =[];
for i = 1:10
    xn_nor = [xn_nor nor_(i*5)];
end

%decision
yn_hat_nor =[];
for i = 1:length(xn_nor)
    if xn_nor(i) <= 0
        yn_hat_nor = [yn_hat_nor -1];
    elseif xn_nor(i) > 0
        yn_hat_nor = [yn_hat_nor 1];
    end
end


%--------------------------------------------------------------------------
%correlator
%--------------------------------------------------------------------------

% corr_ = conv(yn,p,'same');
p_ = repmat(p,1,10); 
corr = p_.*yn;

corr_ = [];
sum = 0;
for i = 1:50
    sum = sum + corr(i);
    corr_ = [corr_ sum];
    if mod(i,5)==0
        sum=0;
    end
end

xn_corr =[];
for i = 1:10
    xn_corr = [xn_corr corr_(i*5)];
end

% decision
yn_hat_corr =[];
for i = 1:length(xn_corr)
    if xn_corr(i) <= 0
        yn_hat_corr = [yn_hat_corr -1];
    elseif xn_corr(i) > 0
        yn_hat_corr = [yn_hat_corr 1];
    end
end

%--------------------------------------------------------------------------
%Figures
%--------------------------------------------------------------------------
t1=0:0.2:50*0.2-0.2;
t2 = 0:10-1;

figure(1)
subplot(2,1,1)
stem(t1,MF_,'blue')
hold on 
plot(t1,MF_,'blue')
hold off
title('Matched Filter output')

subplot(2,1,2)
stem(t1,nor_,'red')
hold on 
plot(t1,nor_,'red')
hold off
title('Normalized Filter output')

figure(2)
subplot(3,1,1)
stem(t2,p_bits,'black')
title('Transmitted bits')

subplot(3,1,2)
stem(t2,yn_hat_MF,'blue')
title('Matched Filter output at each Ts')

subplot(3,1,3)
stem(t2,yn_hat_nor,'red')
title('Normalized Filter output at each Ts')

figure(3)
plot(t1,MF_,'blue')
hold on
plot(t1,corr_,'red')
hold off
legend('Matched Filter convloution','Correlator convloution')

figure(4)
stem(t1,MF_,'blue')
hold on
stem(t1,corr_,'red')
hold off
legend('Matched Filter convloution','Correlator convloution')


figure(5)
subplot(3,1,1)
stem(t2,p_bits,'black')
title('Transmitted bits')

subplot(3,1,2)
stem(t2,yn_hat_MF,'blue')
title('Matched Filter output')

subplot(3,1,3)
stem(t2,yn_hat_corr,'red')
title('correlator output')

