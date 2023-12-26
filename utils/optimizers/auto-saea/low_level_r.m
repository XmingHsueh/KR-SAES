function reward = low_level_r(Ydata, Ynew)

% Update the low-level arm reward
N = length(Ydata);
Yall = [reshape(Ydata, 1, N), Ynew];
[~, idx] = sort(Yall);
in = find(idx == (N+1));
reward = -1/N*in + (N+1)/N;