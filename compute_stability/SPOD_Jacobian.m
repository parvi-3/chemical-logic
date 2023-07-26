N=10;
a=0.139;
k=0.6;
e=0.001;
D=0.002;

fname1=sprintf('FSUBucp_%d.dat',N);
f1=fopen(fname1,'r');
FSU=importdata(fname1,'\t');
fclose(f1);
[R1,~]=size(FSU);
FSU=sortrows(FSU,N+3);

fname2=sprintf('FSUBucp_sorted_num_%d.dat',N);
f2=fopen(fname2,'r');
num=importdata(fname2,'\t');
fclose(f2);

min_ev=zeros(1,R1);
max_ev=zeros(1,R1);

for j=1:R1
    u=zeros(1,N);
    
    for l=1:N
        
        if FSU(j,l)==0
            u(l)=num(j,1);
        end;
        if FSU(j,l)==1
            u(l)=num(j,2);
        end
    end
    
    J=zeros(2*N,2*N);
    for i=1:N
        J(i,i)=-3*u(i)^2+2*(1-a)*u(i)-a;
        J(i,N+i)=-1;
        J(N+i,i)=e*k;
        J(N+i,N+i)=-e-2*D;
        %J(N+i,N+mod(i+1,N))=D;
        J(N+i,N+1+mod(i,N))=D;
        J(N+i,2*N-mod(N-i+1,N))=D;
    end
    
    [s1,s2]=eig(J);
    
%     fname1 = sprintf('Evector_SPOD_N%d_%d.dat',N,j);
%     f1 = fopen(fname1,'w');
%     for i=1:2*N
%         for l=1:2*N
%             fprintf(f1,'%10.8f\t',s1(i,l));
%         end
%         fprintf(f1,'\n');
%     end
%     fclose(f1);
%     
%     fname1 = sprintf('Evalue_SPOD_N%d_%d.dat',N,j);
%     f1 = fopen(fname1,'w');
%     for i=1:2*N
%         fprintf(f1,'%10.8f\n',s2(i,i));
%     end
%     fclose(f1);

    max_ev(j)=max(real(eig(J)));
    min_ev(j)=min(real(eig(J)));
end

fname1 = sprintf('Min_Evalue_SPOD_N_%d.dat',N);
f1 = fopen(fname1,'w');
for i=1:R1
   fprintf(f1,'%10.8f\n',min_ev(i));
end
fclose(f1);

fname1 = sprintf('Max_Evalue_SPOD_N_%d.dat',N);
f1 = fopen(fname1,'w');
for i=1:R1
   fprintf(f1,'%10.8f\n',max_ev(i));
end
fclose(f1);


%plot eigen value
% 
% plot(max_ev, 'r-*');
% set(gca,'xtick',1:R1);
% set(gca,'ticklength',[0.005,0.005]);
% set(gca,'xTicklabel',Labels, 'Fontsize',14);
% rotateXLabels(gca(),45);
% xlabel('Binary String','Fontsize',14);
% ylabel('Max. Eigen Value', 'Fontsize',14);

