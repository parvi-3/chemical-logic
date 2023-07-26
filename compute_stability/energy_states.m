function energy_states(N,B,Labels)

    fname=sprintf('FSUBucp.dat');
    f=fopen(fname,'r');
    V=importdata(fname,'\t');
    fclose(f);
    V=sortrows(V,N+3);
    [R,~]=size(V);
    E=zeros(R,1);
    
    for i=1:R
        for j=1:N
            if V(i,j)==0
                V(i,j)=-1;
            end    
        end
    end     
     
%     for i=1:R
%         for j=1:N
%             if j < N-1
%                 E(i)= E(i) + (V(i,j+1)+V(i,j))^2 + B*(V(i,j+2)-V(i,j))^2;
%             else
%                 if j==N-1
%                     E(i)= E(i) + (V(i,N)+V(i,N-1))^2 + B*(V(i,1)-V(i,N-1))^2;
%                 else
%                     E(i)= E(i) + (V(i,1)+V(i,N))^2 + B*(V(i,2)-V(i,N))^2;
%                 end    
%             end 
%         end
%     end
    
    %Product rule
    for i=1:R
        for j=1:N
            if j>1 && j < N
                E(i)= E(i) + V(i,j+1)*V(i,j) - B*V(i,j+1)*V(i,j-1);
            else
                if j==N
                    E(i)= E(i) + V(i,N)*V(i,1) - B*V(i,1)*V(i,N-1);
                else
                    E(i)= E(i) + V(i,2)*V(i,1) - B*V(i,2)*V(i,N);
                end    
            end
        end
    end
    
    %Writing
    fname = sprintf('State_energy_pdt_%d_%2.1f.dat',N,B);
    f=fopen(fname,'w');
    for i=1:R
        for j=1:N
            fprintf(f,'%d\t',V(i,j));
        end
        fprintf(f,'%3.2f\n',E(i));
    end
    
    %Plot Energy
    figure;
    bar(E,0.5);
    set(gca,'xtick',1:R);
    set(gca,'ticklength', [0 0]);
    set(gca(),'Xticklabel', Labels);
    rotateXLabels( gca, 45 );
    title('Energy Bar Diagram of the Stable Binary Strings','fontsize',16);
    xlabel('Binary String','fontsize',16);
    ylabel('Energy','fontsize',16);
end