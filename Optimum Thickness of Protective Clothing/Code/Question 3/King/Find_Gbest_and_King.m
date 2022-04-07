function [Gbest,Gbest_eval,King,King_eval,King_number]=Find_Gbest_and_King(G,List,city_number,dimension,EvalValue)
%FIND_GBEST 此处显示有关此函数的摘要
%   此处显示详细说明
Gbest=zeros(1,dimension*city_number);                             %初始化诸侯
Gbest_eval=zeros(1,city_number);                         %诸侯适应值
for i=1:city_number
    Gbest(1,dimension*i-1:dimension*i)=G(List(1,i),dimension*i-1:dimension*i);
    Gbest_eval(1,i)=EvalValue(List(1,i),i);
end
[~,GbestList]=sort(Gbest_eval,'descend');
King_number=GbestList(1);
King=Gbest(1,dimension*King_number-1:dimension*King_number);
King_eval=Gbest_eval(1,King_number);
end
