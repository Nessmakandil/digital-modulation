clc
close all
clear all

bits = [0,1,1,0,0,0,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,1,0,0,1,1,0,1,0,1,0,1];
p_bits = [];
for i = 1:length(bits)
    if bits(i) == 0
        p_bits = [p_bits -1];
    elseif bits(i) == 1
        p_bits = [p_bits 1];
    end
end
p_bits_sampled = upsample(p_bits,5);

sps = 5;
R=[0 0 1 1];
D = [2 8 2 8];
for i = 1:4
[NUM, DEN] = rcosine(1, sps, 'sqrt', R(i), D(i));
filterA = filter(NUM,DEN,p_bits_sampled);
filterB = filter(NUM,DEN,filterA);
eyediagram(filterA,2*sps)
title(['eye diagram of A R =', num2str(R(i)), ' D = ', num2str(D(i))])
grid on
eyediagram(filterB,2*sps)
title(['eye diagram of B R =', num2str(R(i)), ' D = ', num2str(D(i))])
grid on
end