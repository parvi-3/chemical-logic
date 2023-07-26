%Use peak values for the unique binary SPOD states and check their stability

function [fsu,fsv,opu,opv]=spodstable(icu,icv,N)

if ischar(N);           N=str2double(N);                 end;
h=0.1;
tfinal=5000;
NT=ceil(tfinal/h);

% icv=FSVBucp(:,1:10);
% icu=FSUBucp(:,1:10);
[R,C]=size(icv);

%Substituting peak values in place of 0's and 1's
for i=1:R
    for j=1:C
        if icv(i,j)==0
            icv(i,j)=0.03; 
            icu(i,j)=-0.1;
        else
            icv(i,j)=0.093;
            icu(i,j)=0.85;
        end;
    end
end

opu=zeros(R,C);
opv=zeros(R,C);

%Solving the system of coupled oscillators using R-K Method
for i=1:R
    
    NU=icu(i,:);
    NV=icv(i,:);
    for T=1:NT
        [dU,dV]=RK(NU,NV,N,h);
        NU = NU + dU(1,:);
        NV = NV + dV(1,:);
    end
    
    opu(i,:)=NU;
    opv(i,:)=NV;
    
end    

fsu=zeros(R,C);
fsv=zeros(R,C);

%Converting to binary
for i=1:R
    for j=1:C
        if opv(i,j) <= 0.06
            fsv(i,j)=0;
        else
            fsv(i,j)=1;
        end;
        if opu(i,j) <= 0.4
            fsu(i,j)=0;
        else
            fsu(i,j)=1;
        end;
    end
end

end