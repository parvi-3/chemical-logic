function [fu] = f(u_i, v_i)
alpha=0.139;
fu=u_i*(1-u_i)*(u_i-alpha)-v_i;
end

