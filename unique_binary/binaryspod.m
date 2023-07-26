%Converting the  strings to binary and finding the unique binary strings

function binaryspod(V,N,batchnum,var)

if ischar(N);           N=str2double(N);                 end;
if ischar(batchnum);    batchnum=str2double(batchnum);   end;

numtrials=10000;
B=zeros(N+2,numtrials);
B(1:N,:) = V(1:N,(batchnum-1)*numtrials+1:batchnum*numtrials);
B(N+3,:)=ones(1,numtrials);
B(N+2,:)=[(batchnum-1)*numtrials+1 : batchnum*numtrials];
B(N+1,:)=zeros(1,numtrials);

if strcmp(var,'u')==1;      thres=0.4;      end;
if strcmp(var,'v')==1;      thres=0.06;     end;

%Converting the frozen states to binary using threshold
for j=1:numtrials
for i=1:N
      if B(i,j) <= thres 
          B(i,j)=0;
      else
          B(i,j)=1;
      end  
end
end
%Writing onto file
fname=sprintf('FSB_%d.dat',batchnum);
f1 = fopen(fname,'w');
for j=1:numtrials
    for i=1:N+1
        fprintf(f1,'%d\t',B(i,j)); 
    end
    fprintf(f1,'\n');  
end
fclose(f1);

%Finding the strings of all low and high states
a=zeros(1,N);
Z=find(sum(abs(B(1:10,:)'-ones(numtrials,1)*a),2)==0);
a=ones(1,N);
O=find(sum(abs(B(1:10,:)'-ones(numtrials,1)*a),2)==0);

%Check for number of occuarances of each of the unique strings
k=0;
for j=1:numtrials
    if B(N+3,j)==1;
        k=k+1;
        clear vars A1;
        clear vars str1;
        A1=mat2cell(B(1:N,j));
        str1=num2str(A1{1});
        str1=str1';
        B(11,j)=bin2dec(str1);
        for i=j+1:numtrials
            if B(N+3,i)==1
                A1=mat2cell(B(1:N,i));
                str2=num2str(A1{1});
                str2=str2';
                if strcmp(str1,str2) ==1; 
                    B(N+3,i)=0; 
                    B(N+3,j)=B(N+3,j)+1; 
                end;
            end
        end
        Bu(:,k)=B(:,j);
    end            
end

C=sortrows(Bu',11);
Bu=C';

%writing onto file
fname=sprintf('Unique_%d.dat',batchnum);
f1 = fopen(fname,'w');
for j=1:k
    for i=1:N+3
        fprintf(f1,'%d\t',Bu(i,j)); 
    end
    fprintf(f1,'\n');  
end
fclose(f1);


%Check for strings which are circular permutations of each other
for j=1:k
    if Bu(N+3,j)>=1;
        clear vars str1;
        clear vars A1;
        A1=mat2cell(Bu(1:N,j));
        str1=num2str(A1{1});
        str1=str1';
        str1=strcat(str1,str1);
        for i=j+1:k
            if Bu(N+3,i)>=1
                clear vars str2;
                clear vars A1;
                A1=mat2cell(Bu(1:N,i));
                str2=num2str(A1{1});
                str2=str2';
                if strfind(str1,str2) > 0 ; 
                    Bu(N+3,j)=Bu(N+3,j)+Bu(N+3,i);
                    Bu(N+3,i)=0; 
                end;
            end
        end    
    end            
end    

%Writing onto file
fname = sprintf('Uniquecp_%d.dat',batchnum);
f1 = fopen(fname,'w');
for j=1:k
    if Bu(N+3,j)>=1
        for i=1:N+3
            fprintf(f1,'%d\t',Bu(i,j)); 
        end
        fprintf(f1,'\n');
    end    
end
fclose(f1);

end