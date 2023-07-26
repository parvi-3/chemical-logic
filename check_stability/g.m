function [gv] = g(u_i,v_i,v_h,v_j,b,D)
epsilon=0.001;
k=0.6;
%b=0.155;
%D=0.0018;
    gv = epsilon*( k*u_i - v_i - b)+D*(v_h + v_j - 2*v_i);
end

