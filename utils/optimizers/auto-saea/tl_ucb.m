function u_value = tl_ucb(sum_reward, id, q_value)

alpha = 2.5;
u_value = q_value + sqrt((alpha*log(id))/length(sum_reward));