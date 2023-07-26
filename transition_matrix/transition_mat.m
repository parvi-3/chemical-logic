function transition_mat(R, stim, nt_stim)

MV=zeros(R,R);

for i=1:R
    fname1=sprintf('PrtclxStr_%d_%5.4f_%d.dat',i,stim,nt_stim);
    f1=fopen(fname1,'r');
    FS=importdata(fname1,'\t');
    fclose(f1);
    for j=1:R        
        MV(i,j)=sum(FS(:,j));
    end
end

%Writing the matrix onto file
fname1 = sprintf('StrxStr_%5.4f_%d.dat',stim,nt_stim);
f1 = fopen(fname1,'w');
for i=1:R
    for j=1:R
        fprintf(f1,'%d\t',MV(i,j));
    end
    fprintf(f1,'\n');
end
fclose(f1);

%Writing the transition matrix onto file
fname1 = sprintf('StrxStr_Transition_%5.4f_%d.dat',stim,nt_stim);
f1 = fopen(fname1,'w');
for i=1:R
    for j=1:R
        fprintf(f1,'%d\t',MV(i,j)/2^N-1);
    end
    fprintf(f1,'\n');
end
fclose(f1);

end