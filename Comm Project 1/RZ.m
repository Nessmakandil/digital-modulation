clc
close all
clear all
%---------------------------------------------------------
%Data && activating the DAC for 70 ms
%----------------------------------------------------------
Data = randi([0,1],500,100);
Data =transpose(Data );
data=reshape(Data,1,[]);
A=4;
Tx=((2*data)-1)*A; % maping for 0 to be –A, 1 to be A
Tx2=repmat(Tx,7,1);
Tx_out=reshape(Tx2,size(Tx2,1)* size(Tx2,2),1);
Tx_out1=reshape(Tx_out,[],500);
Tx_out2=transpose(Tx_out1);
%----------------------------------------------------------
%Data with Delay
%---------------------------------------------------------
no_of_delay=randi([0,6],500,1);
[W L] = size(no_of_delay);
rand_delay= randi([0,1],500,1);
Tx3=((2*rand_delay)-1)*A;
delay=[];
for i=1:W
r=no_of_delay(i,1);
Tx_out2(i,1:r)=Tx3(i,1);
end
data_with_delay=Tx_out2;
%------------------------------------------------------------
%RZ
%------------------------------------------------------------
data2=reshape(data_with_delay,[],1);
bitrate=200;
T = length(data2)/bitrate; % full time of bit sequence
n = 1;
N = n*length(data2);
dt = T/N;
t= 0:dt:T-dt;
[W2 L2]=size(data_with_delay);
for i = 1:7:L2
  data_with_delay(:,i+4:i+6) = 0;
end
d=data_with_delay';
c= reshape(d,1,[]);
%d=reshape(data_with_delay,[],1);
D= stairs(t,c)
ylim([-5 5])

%-------------------------------------------------------------
%statistical mean
%-------------------------------------------------------------
stat_mean_mat = [];
[W1 L1] = size(data_with_delay);
for i = 1:L1
    stat_mean_mat = [stat_mean_mat sum(data_with_delay(:,i))/W1];
end
stat_mean = sum(stat_mean_mat)/L1;
%--------------------------------------------------------------
%time mean
%---------------------------------------------------------------
time_mean = sum(data_with_delay(1,:))/L1;
%---------------------------------------------------------------
%the ensemble autocorrelation function
%--------------------------------------------------------------
tau=0;
q=zeros(L1);
 while tau<L1
     k=[];
for i=1:L1-tau
k=[k  sum(data_with_delay(:,i).*data_with_delay(:,i+tau))./W1];
len=length(k);
q(i,1:len)=k(1,1:len);
q(i,len+1:end)=0;
end
tau=tau+1;
 end
B = flipud(q);
[P O]=size(B);
w=[];
 for j=0:P-1
     w=[w sum(B(j+1,:))./(O-j)];
 end
rho1 = [fliplr(w) w(1:end)];
figure(1)
tau1=-L1:1:L1-1
plot(tau1,rho1)
%----------------------------------------------------------------
%time autocorrelation function
%------------------------------------------------------------------
R = data_with_delay(1,:);
rm=R;
for j = 1:7:L1
    rand_shift = randi([0,1],1,1);
    for i = 1:7
     if  rand_shift == 1
         if i <=3
       rm = [rm; zeros(1,i) rm(j,1:L1-i);];
         elseif i > 3
             rm = [rm; repmat(4,1,i-3) zeros(1,3) rm(j,1:L1-i);];
         end
     elseif rand_shift == 0
        if i <=3
       rm = [rm; zeros(1,i) rm(j,1:L1-i);];
         elseif i > 3
             rm = [rm; repmat(-4,1,i-3) zeros(1,3) rm(j,1:L1-i);];
         end
     end     
    end
    
end
[N V]=size(rm);
time_acf = [];
for i=1:N
    time_acf=[time_acf sum(R.*rm(i,:))./V];
end

rho = [fliplr(time_acf) time_acf(2:end)];
figure(2)
tau2=-length(rho)/2:length(rho)/2-1;
plot(tau2,rho)
%---------------------------------------------------------------
%B.W
%-----------------------------------------------------------------
B=abs(fftshift(fft(rho)));
figure (3)
plot(B)
