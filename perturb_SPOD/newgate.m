%Obtain Output for different stimulation protocol

function newgate(N, stim, nt_stim, batchnum)

if ischar(N);           N=str2double(N);                end;
if ischar(stim);        stim=str2double(stim);          end;
if ischar(nt_stim);     nt_stim=str2double(nt_stim);    end;
if ischar(batchnum);    batchnum=str2double(batchnum);  end;
%if ischar(b);           b=str2double(b);                end;
%if ischar(D);           D=str2double(D);                end;

h=0.1;
tfinal=10000;
NT=ceil(tfinal/h);

fname1=sprintf('FSUBucp_%d.dat',N);
f1=fopen(fname1,'r');
FS=importdata(fname1,'\t');
%FS=cell2mat(textscan(f1,'%12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f', 'delimiter', '\n'));
fclose(f1);
FS=sortrows(FS,N+3);

icv=FS(:,1:N+1);
icu=FS(:,1:N+1);
[R,~]=size(icv);


%Substituting peak values in place of 0's and 1's
for i=1:R
    for j=1:N
        if icv(i,j)==0
            icv(i,j)=0.03;
            icu(i,j)=-0.1;
        else
            icv(i,j)=0.093;
            icu(i,j)=0.85;
        end;
    end
end

for PN=(batchnum-1)*100+1:batchnum*100

%Finding the stimulation protocol
fname1=sprintf('Stimulationsites_%d.dat',N);
f1=fopen(fname1,'r');
SP=importdata(fname1,'\t');
%SP=cell2mat(textscan(f1,'%12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f\t %12.8f', 'delimiter', '\n'));
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
fsu=zeros(R,N+3);
fsv=zeros(R,N+3);
flag1=zeros(R,1);
flag2=zeros(R,1);
flag3=zeros(R,1);
flag4=zeros(R,1);
thresv=0.06;
thresu=0.4;

%Solving the system of coupled oscillators using R-K Method
for i=1:R
    
    E=zeros(1,1);
    NU(1,:)=icu(i,1:N);
    NV(1,:)=icv(i,1:N);
    
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
    %flag1=0;
    N_no=0;
    for j=1:N
        if (NU_var(j) < 1e-3)
            N_no=N_no+1;
        end
    end
    
    if N_no==N; flag1(i)=1; end;
    
    NV_var=var(NV(ceil(4*NT/5):end,:));
    %flag2=0;
    N_no=0;
    for j=1:N
        if (NV_var(j) < 1e-3)
            N_no=N_no+1;
        end
    end
    
    if N_no==N; flag2(i)=1; end;
    
    opu=NU(NT,:);
    opv=NV(NT,:);
    
    %Converting to binary
    for j=1:N
        if opv(j) <= thresv
            fsv(i,j)=0;
        else
            fsv(i,j)=1;
        end;
        if opu(j) <= thresu
            fsu(i,j)=0;
        else
            fsu(i,j)=1;
        end;
    end
    
    clear vars A1;
    clear vars str1;
    temp=fsu(i,1:N);
    A1=mat2cell(temp');
    str1=num2str(A1{1});
    str1=str1';
    fsu(i,N+1)=bin2dec(str1);
    clear vars A1;
    clear vars str1;
    temp=fsv(i,1:N);
    A1=mat2cell(temp');
    str1=num2str(A1{1});
    str1=str1';
    fsv(i,N+1)=bin2dec(str1);
    
    
    % Checking whether input and output states are the same
    %flag3=0;
    %flag4=0;
    if fsu(i,N+1)==icu(i,N+1); flag3(i)=1; end;
    if fsv(i,N+1)==icv(i,N+1); flag4(i)=1; end;
    
    %Writing O/P onto file
    fname1 = sprintf('fsu_%d_%d_%5.4f_%d.dat',i, PN, stim, nt_stim);
    fname2 = sprintf('fsv_%d_%d_%5.4f_%d.dat',i, PN, stim, nt_stim);
    
    f1 = fopen(fname1,'w');
    f2 = fopen(fname2,'w');
    for j=1:N+1
        fprintf(f1,'%12.8f\t',fsu(i,j));
        fprintf(f2,'%12.8f\t',fsv(i,j));
    end
    fprintf(f1,'%12.8f\t %12.8f',flag1(i), flag3(i));
    fprintf(f2,'%12.8f\t %12.8f',flag2(i), flag4(i));
    fclose(f1);
    fclose(f2);
end

%Writing O/P onto file
fname1 = sprintf('fsu_%d_%5.4f_%d.dat', PN, stim, nt_stim);
fname2 = sprintf('fsv_%d_%5.4f_%d.dat', PN, stim, nt_stim);

f1 = fopen(fname1,'w');
f2 = fopen(fname2,'w');
for i=1:R
    for j=1:N+1
        fprintf(f1,'%12.8f\t',fsu(i,j));
        fprintf(f2,'%12.8f\t',fsv(i,j));
    end
    fprintf(f1,'%12.8f\t %12.8f \n',flag1(i), flag3(i));
    fprintf(f2,'%12.8f\t %12.8f \n',flag2(i), flag4(i));
end
fclose(f1);
fclose(f2);

end

end
