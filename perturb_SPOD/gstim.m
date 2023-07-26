function [gv] = gstim(u_i,v_i,v_h,v_j,stim)

if ischar(u_i);      u_i=str2double(u_i);           end;
if ischar(v_i);      v_i=str2double(v_i);           end;
if ischar(v_h);      v_h=str2double(v_h);           end;
if ischar(v_j);      v_j=str2double(v_j);           end;
if ischar(stim);     stim=str2double(stim);         end;

epsilon=0.001;
k=0.6;
b=0.16;
D=0.002;

% if fl==1    
    gv = epsilon*(k*u_i-v_i-b)+D*(v_h + v_j - 2*v_i) + stim; 
% else
%     gv = epsilon*(k*u_i-v_i-b)+D*(v_h + v_j - 2*v_i) ; 
% end    
end