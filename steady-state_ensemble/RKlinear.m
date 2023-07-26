function [dU,dV] = RKlinear(U,V,N,h,b,D)
if N==1
    k1u=f(U(1),V(1));
    k1v=g(U(1),V(1),V(1),V(1),b,D);
    k2u=f(U(1)+h*k1u/2, V(1)+h*k1v/2);
    k2v=g(U(1)+h*k1u/2, V(1)+h*k1v/2, V(1)+h*k1v/2, V(1)+h*k1v/2,b,D);
    k3u=f(U(1)+h*k2u/2, V(1)+h*k2v/2);
    k3v=g(U(1)+h*k2u/2, V(1)+h*k2v/2, V(1)+h*k2v/2, V(1)+h*k2v/2,b,D);
    k4u=f(U(1)+h*k3u, V(1)+h*k3v);
    k4v=g(U(1)+h*k3u, V(1)+h*k3v, V(1)+h*k3v, V(1)+h*k3v,b,D);
end
if N>1
for i=1:N+2 
    k1u(i)=f( U(i), V(i) );
    if i==1 
        k1v(i)=D*(V(i+1)-V(i));
    elseif i==N+2    
        k1v(i)=D*(V(i-1)-V(i));
    else
        k1v(i)=g( U(i), V(i), V(i-1), V(i+1), b, D);
    end
 end
for i=1:N+2 
    k2u(i)=f( U(i)+h*k1u(i)/2, V(i)+h*k1v(i)/2);
    if i==1 
        k2v(i)=D*( (V(i+1)+h*k1v(i+1)/2) - (V(i)+h*k1v(i)/2) );
    elseif i==N+2    
        k2v(i)=D*( (V(i-1)+h*k1v(i-1)/2) - (V(i)+h*k1v(i)/2));
    else
        k2v(i)=g( U(i)+h*k1u(i)/2, V(i)+h*k1v(i)/2, V(i-1)+h*k1v(i-1)/2, V(i+1)+h*k1v(i+1)/2, b, D);  
    end
end
for i=1:N+2
    k3u(i)=f( U(i)+h*k2u(i)/2, V(i)+h*k2v(i)/2);
    if i==1 
        k3v(i)=D*( (V(i+1)+h*k2v(i+1)/2) - (V(i)+h*k2v(i)/2));
    elseif i==N+2    
        k3v(i)=D*( (V(i-1)+h*k2v(i-1)/2) - (V(i)+h*k2v(i)/2)); 
    else
        k3v(i)=g( U(i)+h*k2u(i)/2, V(i)+h*k2v(i)/2, V(i-1)+h*k2v(i-1)/2, V(i+1)+h*k2v(i+1)/2, b, D);
    end    
end
for i=1:N+2
    k4u(i)=f( U(i)+h*k3u(i), V(i)+h*k3v(i));
    if i==1 
        k4v(i)=D* ( (V(i+1)+h*k3v(i+1)) - (V(i)+h*k3v(i)) );
    elseif i==N+2    
        k4v(i)=D* ( (V(i-1)+h*k3v(i-1)) - (V(i)+h*k3v(i)) );
    else
        k4v(i)=g( U(i)+h*k3u(i), V(i)+h*k3v(i), V(i-1)+h*k3v(i-1), V(i+1)+h*k3v(i+1), b, D);
    end
end
end
dU=zeros(1,N+2);
dV=zeros(1,N+2);
for i=1:N+2
    dU(1,i)=h*(k1u(i)+2*k2u(i)+2*k3u(i)+k4u(i))/6;
    dV(1,i)=h*(k1v(i)+2*k2v(i)+2*k3v(i)+k4v(i))/6;
end    
end