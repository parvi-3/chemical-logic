function stim_str_pn(N, stim, nt_stim, PN, strnum)

h=0.1;
tfinal=10000;
NT=ceil(tfinal/h);

fname1=sprintf('FSUBucp_%d.dat',N);
f1=fopen(fname1,'r');
FS=importdata(fname1,'\t');
fclose(f1);
FS=sortrows(FS,N+3);

icv=FS(strnum,1:N+1);
icu=FS(strnum,1:N+1);

%Substituting peak values in place of 0's and 1's
for j=1:N
    if icv(j)==0
        icv(j)=0.03;
        icu(j)=-0.1;
    else
        icv(j)=0.093;
        icu(j)=0.85;
    end;
end

%Finding the stimulation protocol
fname1=sprintf('Stimulationsites_%d.dat',N);
f1=fopen(fname1,'r');
SP=importdata(fname1,'\t');
fclose(f1);

N_sim=SP(PN,1:N);
N_stim=zeros(1,N);
fl=0;
for j=1:N
    if N_sim(j)==1
        fl=fl+1;
        N_stim(fl)=j;
    end
end

NU=zeros(NT+1,N);
NV=zeros(NT+1,N);
fsu=zeros(1,N+3);
fsv=zeros(1,N+3);
flag1=zeros(1,1);
flag2=zeros(1,1);
flag3=zeros(1,1);
flag4=zeros(1,1);
thresv=0.06;
thresu=0.4;

%Solving the system of coupled oscillators using R-K Method
    
E=zeros(1,1);
NU(1,:)=icu(1:N);
NV(1,:)=icv(1:N);

for T=1:ceil(NT/2)-1
    [dU,dV]=RKstim(NU(T,:),NV(T,:),N,h,E,E);
    NU(T+1,:) = NU(T,:) + dU(1,:);
    NV(T+1,:) = NV(T,:) + dV(1,:);
end
for T=ceil(NT/2):ceil(NT/2) + nt_stim
    [dU,dV]=RKstim(NU(T,:),NV(T,:),N,h,N_stim,stim);
    NU(T+1,:) = NU(T,:) + dU(1,:);
    NV(T+1,:) = NV(T,:) + dV(1,:);
end
for T=ceil(NT/2) + nt_stim + 1: NT
    [dU,dV]=RKstim(NU(T,:),NV(T,:),N,h,E,E);
    NU(T+1,:) = NU(T,:) + dU(1,:);
    NV(T+1,:) = NV(T,:) + dV(1,:);
end

%Determining whether it reaches a SPOD state
NU_var=var(NU(ceil(4*NT/5):end,:));
N_no=0;
for j=1:N
    if (NU_var(j) < 1e-3)
        N_no=N_no+1;
    end
end

if N_no==N; flag1=1; end;

NV_var=var(NV(ceil(4*NT/5):end,:));
N_no=0;
for j=1:N
    if (NV_var(j) < 1e-3)
        N_no=N_no+1;
    end
end

if N_no==N; flag2=1; end;

opu=NU(NT,:);
opv=NV(NT,:);

%Converting to binary
for j=1:N
    if opv(j) <= thresv
        fsv(j)=0;
    else
        fsv(j)=1;
    end;
    if opu(j) <= thresu
        fsu(j)=0;
    else
        fsu(j)=1;
    end;
end

clear vars A1;
clear vars str1;
temp=fsu(1:N);
A1=mat2cell(temp');
str1=num2str(A1{1});
str1=str1';
fsu(N+1)=bin2dec(str1);
clear vars A1;
clear vars str1;
temp=fsv(1:N);
A1=mat2cell(temp');
str1=num2str(A1{1});
str1=str1';
fsv(N+1)=bin2dec(str1);


% Checking whether input and output states are the same
if fsu(N+1)==icu(N+1); flag3=1; end;
if fsv(N+1)==icv(N+1); flag4=1; end;

%Writing O/P onto file
fname1 = sprintf('fsu_%d_%d_%5.4f_%d.dat',strnum, PN, stim, nt_stim);
fname2 = sprintf('fsv_%d_%d_%5.4f_%d.dat',strnum, PN, stim, nt_stim);

f1 = fopen(fname1,'w');
f2 = fopen(fname2,'w');
for j=1:N+1
    fprintf(f1,'%12.8f\t',fsu(j));
    fprintf(f2,'%12.8f\t',fsv(j));
end
fprintf(f1,'%12.8f\t %12.8f',flag1, flag3);
fprintf(f2,'%12.8f\t %12.8f',flag2, flag4);
fclose(f1);
fclose(f2);

end