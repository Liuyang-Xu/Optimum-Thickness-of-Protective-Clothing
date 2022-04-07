clc
clear
human_number=10;                                    %人群人数
city_number=3;                                      %城邦个数
dimension=2;                                        %位置维数
w=[0.2,0.2,0.2];        %惯性系数
c1=[0.25,0.3,0.3];       %自我学习系数
c2=[0.55,0.6,0.6];       %群体学习系数
c3=[0.55,0.65,0.75];   %追随系数
c4=[0.15,0.2,0.2];       %远离系数

g=10;                                              %学习代数
ChangeMatrix=repmat(eye(2),1,city_number);          %转化矩阵

range=[0.6,25;0.6,6.4];                                  %位置限制
Vmax=[3.5,1]*ChangeMatrix;     %最大速度限制
Vmin=[0.1,0.1]*ChangeMatrix;    %最小速度限制
c5=0.1*(range(:,2)-range(:,1))';                     %散步系数
radius=0.2*(range(:,2)-range(:,1))'*ChangeMatrix;       %势力范围

[G] = city_spawn(human_number,dimension,city_number,range,ChangeMatrix);  %生成初始解集
%G 横向分为数个城邦，各城邦中依次排序其个体各维度

fprintf('初始化个体最优位置\n');
Pbest=G;                                                                  %初始化个体最优位置
[EvalValue] = OP_evalCal(G,human_number,dimension,city_number);           %计算初始种群适应度
[~,List]=sort(EvalValue,'descend');
fprintf('计算初始种群适应度结束\n');

[Gbest,Gbest_eval,King,King_eval,King_number]=Find_Gbest_and_King(G,List,city_number,dimension,EvalValue);          %得到诸侯及国王
V=rand(human_number,dimension*city_number);                               %初始化速度
rememberSL=cell(g,1);                                 %记录
remember=zeros(g,1);
remembertheKing=cell(g,1); 
for k=1:g
     fprintf('开始第%d次学习\n',k);
     [G,V] = people_learn(G,Pbest,Gbest,V,w,c1,c2,c3,c4,c5,range,Vmax,Vmin,human_number,dimension,city_number,radius,King_number,King,List);    %粒子群学习并更新位置
     
     [Pbest] = OPKing_UpdateMemory(G,Pbest,human_number,dimension,city_number);      %更新粒子群记忆
     [EvalValue] = OP_evalCal(G,human_number,dimension,city_number);
     [~,List]=sort(EvalValue,'descend');
     [Gbest,Gbest_eval,King,King_eval,King_number]=Find_Gbest_and_King(G,List,city_number,dimension,EvalValue);
     rememberSL{k,1}=King;   %记录
     [rememberValue] = OP_evalCal(Pbest,human_number,dimension,city_number);
     [~,rememberList]=sort(rememberValue,'descend');
     [~,~,rememberKing,rememberKing_eval,~]=Find_Gbest_and_King(Pbest,rememberList,city_number,dimension,rememberValue);
     remember(k)=rememberKing_eval;
     remembertheKing{k,1}=rememberKing;
end
x=1:g;
plot(x,remember,'r');
    