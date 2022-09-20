clc  
close all 
clear all 
%---------------------------------------------------------------- 
%% BPSK 
%---------------------------------------------------------------- 
bits1 = randi([0,1],1,100000); 
BPSK_constellation = [1, -1]; 
BPSK_constellation_points = [1, 0]; 
S_BPSK = []; 
for i = 1:length(bits1) 
    if bits1(i) == 0 
        S_BPSK = [S_BPSK -1]; 
    elseif bits1(i) == 1 
        S_BPSK = [S_BPSK 1]; 
    end 

end 
  
snr_db = -2:0.5:5; 
snr = 10.^(snr_db./10); 
L =length(S_BPSK); 
Eb = sum(BPSK_constellation.^2)/length(BPSK_constellation); 
N0 = Eb./snr; 
BER_BPSK = []; 
  
for i =1:length(N0) 
Nv = randn(1,L).*sqrt(N0(i)/2); 
X = S_BPSK + Nv; 
  
X_demap_BPSK = []; 
for k = 1:L 
    [x, z] = min(abs(BPSK_constellation-X(k)));  
    X_hat(k) =  BPSK_constellation(z); 
    X_demap_BPSK = [X_demap_BPSK BPSK_constellation_points(z)]; 
end 
  
count = 0; 
for i = 1:L 
    if bits1(i) ~= X_demap_BPSK(i) 
        count = count+1; 
    end 
end 
BER_BPSK = [BER_BPSK count/L]; 
  
end 
  
BER_theoretical = (1/2).*erfc(sqrt(Eb./N0)); 
  
figure(1) 
semilogy(snr_db,BER_BPSK,'Color','B','LineWidth',1.0) 
hold on  
semilogy(snr_db,BER_theoretical,'Color','R','LineWidth',1.0) 
legend('BPSK BER','BER theoretical') 
ylabel('BER') 
xlabel('SNR (dB)') 
grid on 
  
%---------------------------------------------------------------- 
%% QPSK 
%---------------------------------------------------------------- 
QPSK_constellation = [-1-1i, -1+1i, 1+1i,  1-1i]; 
QPSK_constellation_points = [0 0; 0 1; 1 1; 1 0]; 
  
S_QPSK = []; 
for i = 1:2:length(bits1) 
    if bits1(i) == 0 && bits1(i+1) == 0 
        S_QPSK = [S_QPSK -1-1i]; 
    elseif bits1(i) == 0 && bits1(i+1) == 1  
        S_QPSK = [S_QPSK -1+1i]; 
    elseif bits1(i) == 1 && bits1(i+1) == 0  
        S_QPSK = [S_QPSK 1-1i]; 
    elseif bits1(i) == 1 && bits1(i+1) == 1  
        S_QPSK = [S_QPSK 1+1i]; 
    end 
end 
  
L =length(S_QPSK);  
Eb = sum(abs(QPSK_constellation).^2)/(2*length(QPSK_constellation)); 
N0 = Eb./snr; 
BER_QPSK = []; 
  
for i =1:length(N0) 
Nv = randn(1,L).*sqrt(N0(i)/2)+randn(1,L).*sqrt(N0(i)/2)*1i; 
X = S_QPSK + Nv; 
  
X_demap_QPSK = []; 
for k = 1:L     
    [x, z] = min((real(QPSK_constellation-X(k))).^2+(imag(QPSK_constellation-X(k))).^2);  
    X_hat(k) =  QPSK_constellation(z); 
    X_demap_QPSK = [X_demap_QPSK QPSK_constellation_points(z,:)]; 
end 
  
count = 0; 
for u = 1:length(bits1) 
    if bits1(u) ~= X_demap_QPSK(u) 
        count = count+1; 
    end 
end 
BER_QPSK = [BER_QPSK count/length(bits1)]; 
  
end 
  
BER_theoretical = (1/2).*erfc(sqrt(Eb./N0)); 
  
figure(2) 
semilogy(snr_db,BER_QPSK,'Color','B','LineWidth',1.0) 
hold on  
semilogy(snr_db,BER_theoretical,'Color','R','LineWidth',1.0) 
legend('QPSK BER','BER theoretical') 
ylabel('BER') 
xlabel('SNR (dB)') 
grid on 
  
%---------------------------------------------------------------- 
%% 8-PSK 
%---------------------------------------------------------------- 
bits2 = bits1(1:99999); 

 

S_MPSK = []; 
M = 8; 
MPSK_constellation = [exp((2*pi/M)*0*1i), exp((2*pi/M)*1*1i),exp((2*pi/M)*2*1i),... 
                      exp((2*pi/M)*3*1i), exp((2*pi/M)*4*1i),exp((2*pi/M)*5*1i),... 
                      exp((2*pi/M)*6*1i), exp((2*pi/M)*7*1i)]; 
MPSK_constellation_points = [0 0 0; 0 0 1; 0 1 1;... 
                             0 1 0; 1 1 0; 1 1 1;... 
                             1 0 1; 1 0 0]; 
for i = 1:3:length(bits2) 
    if bits2(i) == 0 && bits2(i+1) == 0 && bits2(i+2) == 0 
        S_MPSK = [S_MPSK exp((2*pi/M)*0*1i)]; 
    elseif bits2(i) == 0 && bits2(i+1) == 0 && bits2(i+2) == 1 
        S_MPSK = [S_MPSK exp((2*pi/M)*1*1i)]; 
    elseif bits2(i) == 0 && bits2(i+1) == 1 && bits2(i+2) == 1 
        S_MPSK = [S_MPSK exp((2*pi/M)*2*1i)]; 
    elseif bits2(i) == 0 && bits2(i+1) == 1 && bits2(i+2) == 0 
        S_MPSK = [S_MPSK exp((2*pi/M)*3*1i)]; 
    elseif bits2(i) == 1 && bits2(i+1) == 1 && bits2(i+2) == 0 
        S_MPSK = [S_MPSK exp((2*pi/M)*4*1i)]; 
    elseif bits2(i) == 1 && bits2(i+1) == 1 && bits2(i+2) == 1 
        S_MPSK = [S_MPSK exp((2*pi/M)*5*1i)]; 
    elseif bits2(i) == 1 && bits2(i+1) == 0 && bits2(i+2) == 1 
        S_MPSK = [S_MPSK exp((2*pi/M)*6*1i)]; 
    elseif bits2(i) == 1 && bits2(i+1) == 0 && bits2(i+2) == 0 
        S_MPSK = [S_MPSK exp((2*pi/M)*7*1i)]; 
    end 
end 
  
L =length(S_MPSK);  
Eb = sum(abs(MPSK_constellation).^2)/(3*length(MPSK_constellation)); 
N0 = Eb./snr; 
BER_MPSK = []; 
  
  
for i =1:length(N0) 
Nv = randn(1,L).*sqrt(N0(i)/2)+randn(1,L).*sqrt(N0(i)/2)*1i; 
X = S_MPSK + Nv; 
  
X_hat =[]; 
X_demap_MPSK = []; 
for k = 1:L 
    [x, z] = min((real(MPSK_constellation-X(k))).^2+(imag(MPSK_constellation-X(k))).^2);  
    X_hat(k) =  MPSK_constellation(z); 
    X_demap_MPSK = [X_demap_MPSK MPSK_constellation_points(z,:)]; 
  
end 
count = 0; 
for u = 1:length(bits2)

    if bits2(u) ~= X_demap_MPSK(u) 
        count = count+1; 
    end 
end 
BER_MPSK = [BER_MPSK count/length(bits2)]; 
  
end 
  
BER_theoretical = (1/3)*erfc(sqrt(3*Eb./N0).*sin(pi/M)); 
  
figure(3) 
semilogy(snr_db,BER_MPSK,'Color','B','LineWidth',1.0) 
hold on  
semilogy(snr_db,BER_theoretical,'Color','R','LineWidth',1.0) 
legend('8-PSK BER','BER theoretical') 
ylabel('BER') 
xlabel('SNR (dB)') 
grid on 
  
%---------------------------------------------------------------- 
%% 16-QAM 
%---------------------------------------------------------------- 
S_QAM = []; 
QAM_constellation = [-3+3*1i; -3+1i; -3-1i; -3-3*1i;... 
                     -1+3*1i; -1+1i; -1-1i; -1-3*1i;... 
                      1+3*1i;  1+1i;  1-1i;  1-3*1i;... 
                      3+3*1i;  3+1i;  3-1i;  3-3*1i]; 
constellation_points = [0 0 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0;... 
                        0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0;... 
                        1 1 1 0; 1 1 1 1; 1 1 0 1; 1 1 0 0;... 
                        1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0]; 
for i = 1:4:length(bits1) 
    if bits1(i:i+3) == constellation_points(1,:) 
        S_QAM = [S_QAM QAM_constellation(1)]; 
    elseif bits1(i:i+3) == constellation_points(2,:)  
        S_QAM = [S_QAM QAM_constellation(2)]; 
    elseif bits1(i:i+3) == constellation_points(3,:)  
        S_QAM = [S_QAM QAM_constellation(3)]; 
    elseif bits1(i:i+3) == constellation_points(4,:)  
        S_QAM = [S_QAM QAM_constellation(4)]; 
    elseif bits1(i:i+3) == constellation_points(5,:)  
        S_QAM = [S_QAM QAM_constellation(5)]; 
    elseif bits1(i:i+3) == constellation_points(6,:)   
        S_QAM = [S_QAM QAM_constellation(6)]; 
    elseif bits1(i:i+3) == constellation_points(7,:)   
        S_QAM = [S_QAM QAM_constellation(7)]; 
    elseif bits1(i:i+3) == constellation_points(8,:)   
        S_QAM = [S_QAM QAM_constellation(8)]; 
    elseif bits1(i:i+3) == constellation_points(9,:)  
        S_QAM = [S_QAM QAM_constellation(9)]; 
    elseif bits1(i:i+3) == constellation_points(10,:)   
        S_QAM = [S_QAM QAM_constellation(10)]; 
    elseif bits1(i:i+3) == constellation_points(11,:)   
        S_QAM = [S_QAM QAM_constellation(11)]; 
    elseif bits1(i:i+3) == constellation_points(12,:)   
        S_QAM = [S_QAM QAM_constellation(12)]; 
    elseif bits1(i:i+3) == constellation_points(13,:)   
        S_QAM = [S_QAM QAM_constellation(13)]; 
    elseif bits1(i:i+3) == constellation_points(14,:)   
        S_QAM = [S_QAM QAM_constellation(14)]; 
    elseif bits1(i:i+3) == constellation_points(15,:)   
        S_QAM = [S_QAM QAM_constellation(15)]; 
    elseif bits1(i:i+3) == constellation_points(16,:)    
        S_QAM = [S_QAM QAM_constellation(16)]; 
    end 
end 
  
L =length(S_QAM);  
Eb = sum(abs(QAM_constellation).^2)/(4*length(QAM_constellation)); 
N0 = Eb./snr; 
BER_QAM = []; 
  
for i =1:length(N0) 
Nv = randn(1,L).*sqrt(N0(i)/2)+randn(1,L).*sqrt(N0(i)/2)*1i; 
X = S_QAM + Nv; 
  
X_demap =[]; 
for k = 1:L 
   [x, z] = min((real(QAM_constellation-X(k))).^2+(imag(QAM_constellation-X(k))).^2);  
    X_hat(k) =  QAM_constellation(z); 
    X_demap = [X_demap constellation_points(z,:)]; 
end 
  
count = 0; 
for u = 1:length(bits1) 
    if bits1(u) ~= X_demap(u) 
        count = count+1; 
    end 
end 
BER_QAM = [BER_QAM count/length(bits1)]; 
end 
  
BER_theoretical = (1.5/4)*erfc(sqrt(1./(N0))); 
  
figure(4) 
semilogy(snr_db,BER_QAM,'Color','B','LineWidth',1.0) 
hold on  
semilogy(snr_db,BER_theoretical,'Color','R','LineWidth',1.0) 
legend('QAM BER','BER theoretical') 
ylabel('BER') 
xlabel('SNR (dB)') 
grid on 
  
figure(5) 
semilogy(snr_db,BER_BPSK,'Color','R','LineWidth',1.0) 
hold on 
semilogy(snr_db,BER_QPSK,'Color','B','LineWidth',1.0) 
hold on 
semilogy(snr_db,BER_MPSK,'Color','K','LineWidth',1.0) 
hold on 
semilogy(snr_db,BER_QAM,'Color','M','LineWidth',1.0) 
hold on  
legend('BPSK BER','QPSK BER','8-PSK BER','16-QAM BER') 
ylabel('BER') 
xlabel('SNR (dB)') 
grid on 