% To see which conditions lead to stable states (No stimulation)

clear all;
clc;

%Randomly generating initial condition for each binary string using the threshold 
%Threshold for binary state: m to 0.07 is state '0' and 0.07 to M is state'1'

N=10;
f=fopen('Timeseries.txt','r');
IC=textscan(f,'%12.8f\t %12.8f', 'delimiter', '\n');
TP=size(IC{1});
fclose(f);

m=min(IC{2});
M=max(IC{2});

%Creating a matrix of all the binary strings of length N
N=10;
IP=zeros(2^N, N, 3);
B=dec2bin([0:2^N-1]);
for i=1:2^N
    for j=1:N
        IP(i,j,1)=str2num(B(i,j));
    end
end  

for k=1:16  
for i=1:N
    if IP(k,i,1)==0 
        IP(k,i,2)=(0.07-m)*rand(1,1) + m;
    else
        IP(k,i,2)=(M-0.07)*rand(1,1) + 0.07;
    end
    E=find(abs(IP(k,i,2)-IC{2}(:)) < 0.001);
    s=size(E);
    j=randi([1,s(1)]);

    IP(k,i,2)=IC{2}(E(j,1));
    IP(k,i,3)=IC{1}(E(j,1));
end
end


% Solving the system

h=0.01; 
tfinal=2000;
NT=ceil(tfinal/h);
E=zeros(1,1);

for k=1:16
k;
clearvars NU NV;

NU=IP(k,:,3);
NV=IP(k,:,2);

for T=1:NT
    t=h*T;
    [dU,dV]=RK(t,T,NU,NV,N,h,E);
    NU=[NU; NU(T,:) + dU(1,:)];
    NV=[NV; NV(T,:) + dV(1,:)];
end


%Output State
for i=1:N
      if NV(NT+1,i) <= 0.07 
          OP(k,i,1)=0;
      else
          OP(k,i,1)=1;
      end  
end  

% OPU(:,:,k)=NU;
% OPV(:,:,k)=NV;

OP(k,:,2)=NV(NT+1,:);
OP(k,:,3)=NU(NT+1,:);

end