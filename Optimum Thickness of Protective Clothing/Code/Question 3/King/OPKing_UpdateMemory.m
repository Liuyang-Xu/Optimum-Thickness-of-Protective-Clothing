function [Pbest] = OPKing_UpdateMemory(G,Pbest,human_number,dimension,city_number)
%EVAL 此处显示有关此函数的摘要
%   此处显示详细说明
%   输入G，计算每个个体的适应度函数，得到每个个体的历史最优位置和当前人群的最优位置

[P_Eval] = OP_evalCal(G,human_number,dimension,city_number);                  %计算当前人群的适应度
[PB_EvalValue] = OP_evalCal(Pbest,human_number,dimension,city_number);        %计算人群最佳位置的适应度
for i=1:human_number
    for j=1:city_number
        if PB_EvalValue(i,j)<P_Eval(i,j)
            Pbest(i,dimension*j-1:dimension*j)=G(i,dimension*j-1:dimension*j);             %更新单个粒子记忆
        end
    end
end
%得到新的种群最佳位置




end

