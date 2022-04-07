function [ Uxt ] = Q1_InitialConditions_kChange (x)
%在此函数中定义初值条件
if x==1
    Uxt=75;
else
    Uxt=37;
end
end