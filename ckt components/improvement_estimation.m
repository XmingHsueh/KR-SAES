function imp = improvement_estimation(objs_source,objs_val,objs_target,sim)

no_samples_source = length(objs_source);
no_samples_target = length(objs_target);
coe_source  = ([ones(no_samples_source,1), transpose([1:no_samples_source])])\log(objs_source);
ps = [exp(coe_source(1)) -coe_source(2)];
coe_target = ([ones(no_samples_target,1), transpose([1:no_samples_target])])\log(objs_target);
pt = [exp(coe_target(1)) -coe_target(2)];
tv = (-log(min(objs_val)/ps(1))/ps(2)>0)*(-log(min(objs_val)/ps(1))/ps(2))+(-log(min(objs_val)/ps(1))/ps(2)<0)*0;
decay_target = @(x,a)a(1)*exp(-a(2)*x);

if sim>=0
    imp = sim*(min(objs_target)-decay_target(sim*(no_samples_target+no_samples_source-tv),pt));
else
    imp = 0;
end