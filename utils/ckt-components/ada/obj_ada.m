function obj_sim_nega = obj_ada(x_delta,func_source,Xt,yt)
[nt,dim] = size(Xt);
Xt_ada = Xt - repmat(x_delta,nt,1);
yt_ada = yt;
ranks_target = zeros(nt,1);
objs_val = zeros(nt,1);
ranks_val = zeros(nt,1);
for i = 1:nt
    ranks_target(i) = sum(yt_ada<yt_ada(i))+1;
    objs_val(i) = func_source(Xt_ada(i,:));
end
for i = 1:nt
    ranks_val(i) = sum(objs_val<objs_val(i))+1;
end
uncertainty = 1-max(abs(x_delta));
obj_sim_nega = -uncertainty*mySpearman(ranks_target,ranks_val);