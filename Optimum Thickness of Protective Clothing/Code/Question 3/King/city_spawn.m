function [G] = city_spawn(human_number,dimension,city_number,range,ChangeMatrix)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%   粒子群算法，输入粒子维数及范围限制，生成初始粒子群(个数为m)
%   range是一个n*2的矩阵，第一列表示最小值，第二列表示最大值

max=range(:,2)';                                  %各个维度的最大范围限制
min=range(:,1)';                                  %各个维度的最小范围限制

G=repmat(min*ChangeMatrix,human_number,1)+repmat(((max-min)*ChangeMatrix),human_number,1).*rand(human_number,dimension*city_number);%生成初始种群

G=ceil(G*10)/10;

end

