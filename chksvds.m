%%checking speed of svds

clear
clc

n = 1000;
t_max = 5000;

r_1 = 100;
r_0 = 1;
t_0 = 0;
t_1 = 0;

k = 500;

for ii = 1 : 100
    DataTemp = randn(n, k);
    Data = DataTemp * DataTemp';
    
    tt0 = tic;
    X_0 = svds(Data, r_0);
    t_0 = t_0 + toc(tt0);
    
    tt1 = tic;
    X_1 = svds(Data, r_1);
    t_1 = t_1 + toc(tt1);
end


t_0/100
t_1/100