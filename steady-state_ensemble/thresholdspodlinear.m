function thresholdspodlinear(N, batchnum, b, D)

if ischar(N);        N = str2double(N);              end;
if ischar(batchnum); batchnum = str2double(batchnum); end;
if ischar(b);        b = str2double(b);              end;
if ischar(D);        D = str2double(D);              end;

%Replace random stream with one based on SEED
RandStream.setDefaultStream(RandStream('mt19937ar','Seed',batchnum));

h=0.1;
tfinal=10000;
NT=ceil(tfinal/h);

%Reading timeseries data from file
f=fopen('TimeSeriesEqArL.dat','r');
IC=cell2mat(textscan(f,'%12.8f\t %12.8f', 'delimiter', '\n'));
TP=size(IC);
fclose(f);

NU=zeros(NT,N+2);
NV=zeros(NT,N+2);

for trial=(batchnum-1)*200 + 1 : batchnum*200
    
    % Getting IC from Time Series
    P=randi([1,TP(1,1)],N+2,1);
    U0 = zeros(1,N+2);
    V0 = zeros(1,N+2);
    for i=1:N+2
        U0(i)=IC(P(i),1);
        V0(i)=IC(P(i),2);
    end
        
    %Solving the system of coupled oscillators using R-K Method
    NU(1,:)=U0;
    NV(1,:)=V0;
    for T=1:NT-1
        [dU,dV]=RKlinear(NU(T,:),NV(T,:),N,h, b, D);
        NU(T+1,:) = NU(T,:) + dU(1,:);
        NV(T+1,:) = NV(T,:) + dV(1,:);
    end
    
    %Determining whether it reaches a SPOD state
    NU_var=var(NU(ceil(NT/2):end,:));
    flag1=0;
    N_no=0;
    for j=2:N+1
        if (abs(NU_var(j)) < 1e-3)
            N_no=N_no+1;
        end
    end
    
    if N_no==N; flag1=1; end;
    
    NV_var=var(NV(ceil(NT/2):end,:));
    flag2=0;
    N_no=0;
    for j=2:N+1
        if (abs(NV_var(j)) < 1e-3)
            N_no=N_no+1;
        end
    end
    
    if N_no==N; flag2=1; end;
    
    %Storing the I/P and O/P
    IP(:,1)=U0;
    IP(:,2)=V0;
    OP(:,1)=NU(NT,:);
    OP(:,2)=NV(NT,:);
    
    %Writing I/P and O/P onto file
    fname1 = sprintf('ICU_%d_%2.3f_%2.4f.dat',trial, b, D);
    fname2 = sprintf('OPU_%d_%2.3f_%2.4f.dat',trial, b, D);
    fname3 = sprintf('ICV_%d_%2.3f_%2.4f.dat',trial, b, D);
    fname4 = sprintf('OPV_%d_%2.3f_%2.4f.dat',trial, b, D);
    
    f1 = fopen(fname1,'w');
    f2 = fopen(fname2,'w');
    f3 = fopen(fname3,'w');
    f4 = fopen(fname4,'w');
    
    for j=1:N+2
        fprintf(f1,'%12.8f\t',IP(j,1));
        fprintf(f2,'%12.8f\t',OP(j,1));
        fprintf(f3,'%12.8f\t',IP(j,2));
        fprintf(f4,'%12.8f\t',OP(j,2));
    end
    fprintf(f2,'%d',flag1);
    fprintf(f4,'%d',flag2);
    
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    
end

end
