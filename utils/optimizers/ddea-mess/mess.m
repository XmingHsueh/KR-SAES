function strategy_id = mess(FEsUsed,FEsMax)

beta = (1-(FEsUsed/FEsMax)^3)^2;
alpha = abs(beta*sin(3*pi/2+sin(beta*3*pi/2)));
P3 = (1-alpha)*(alpha>2/3)+1/3*(alpha<=2/3);
P1 = (1-P3)/2;
r = rand;
if r <= P1
    strategy_id = 1;
elseif r > P1 && r <= 2*P1
    strategy_id = 2;
else
    strategy_id = 3;
end