function [coeff,Xrank,Yrank] = mySpearman(X , Y)  

if length(X) ~= length(Y)  
    error('different dimensions!');  
    return;  
end  

N = length(X); 
Xrank = zeros(1 , N);
Yrank = zeros(1 , N);
for i = 1 : N  
    cont1 = 1;
    cont2 = -1;
    for j = 1 : N  
        if X(i) < X(j)  
            cont1 = cont1 + 1;  
        elseif X(i) == X(j)  
            cont2 = cont2 + 1;  
        end  
    end  
    Xrank(i) = cont1 + mean([0 : cont2]);  
end  
for i = 1 : N  
    cont1 = 1;
    cont2 = -1;
    for j = 1 : N  
        if Y(i) < Y(j)  
            cont1 = cont1 + 1;  
        elseif Y(i) == Y(j)  
            cont2 = cont2 + 1;  
        end  
    end  
    Yrank(i) = cont1 + mean([0 : cont2]);  
end  

fenzi = 6 * sum((Xrank - Yrank).^2);  
fenmu = N * (N^2 - 1);  
coeff = 1 - fenzi / fenmu;  
  
end