clear all;
clc;

% Find the range of u and v values in a single cycle

N=1;
h=0.1;
tfinal=5000;
NT=ceil(tfinal/h);

WU=rand(1,1);
WV=rand(1,1);
for T=1:NT
        t=h*T;
        [dU,dV]=RK(t,T,WU,WV,N,h);
        WU=[WU; WU(T,1) + dU(1,1)];
        WV=[WV; WV(T,1) + dV(1,1)];
end

M=max(WU);
I=find(abs(M-WU(1:end))<0.00001)
M=max(WU(I(1,1)+1:end));
J=find(abs(M-WU(I(1,1)+1:end))<0.00001);
for i=1:size(J)
    if J(i+1)-J(i)>1 
        TP=J(i+1);
        break;
    end
end    
IC=[WU(I(1,1):I(1,1)+TP),WV(I(1,1):I(1,1)+TP)];

%Writing Initial Conditions onto a file
f=fopen('Timeseries.txt','w');
for i=1:TP
fprintf(f,'%12.8f\t %12.8f\n',IC(i,1),IC(i,2));
end
fclose(f);
