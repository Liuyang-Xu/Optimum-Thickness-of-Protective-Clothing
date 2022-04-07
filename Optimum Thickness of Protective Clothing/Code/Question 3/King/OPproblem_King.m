clc
clear
human_number=10;                                    %��Ⱥ����
city_number=3;                                      %�ǰ����
dimension=2;                                        %λ��ά��
w=[0.2,0.2,0.2];        %����ϵ��
c1=[0.25,0.3,0.3];       %����ѧϰϵ��
c2=[0.55,0.6,0.6];       %Ⱥ��ѧϰϵ��
c3=[0.55,0.65,0.75];   %׷��ϵ��
c4=[0.15,0.2,0.2];       %Զ��ϵ��

g=10;                                              %ѧϰ����
ChangeMatrix=repmat(eye(2),1,city_number);          %ת������

range=[0.6,25;0.6,6.4];                                  %λ������
Vmax=[3.5,1]*ChangeMatrix;     %����ٶ�����
Vmin=[0.1,0.1]*ChangeMatrix;    %��С�ٶ�����
c5=0.1*(range(:,2)-range(:,1))';                     %ɢ��ϵ��
radius=0.2*(range(:,2)-range(:,1))'*ChangeMatrix;       %������Χ

[G] = city_spawn(human_number,dimension,city_number,range,ChangeMatrix);  %���ɳ�ʼ�⼯
%G �����Ϊ�����ǰ���ǰ�����������������ά��

fprintf('��ʼ����������λ��\n');
Pbest=G;                                                                  %��ʼ����������λ��
[EvalValue] = OP_evalCal(G,human_number,dimension,city_number);           %�����ʼ��Ⱥ��Ӧ��
[~,List]=sort(EvalValue,'descend');
fprintf('�����ʼ��Ⱥ��Ӧ�Ƚ���\n');

[Gbest,Gbest_eval,King,King_eval,King_number]=Find_Gbest_and_King(G,List,city_number,dimension,EvalValue);          %�õ�������
V=rand(human_number,dimension*city_number);                               %��ʼ���ٶ�
rememberSL=cell(g,1);                                 %��¼
remember=zeros(g,1);
remembertheKing=cell(g,1); 
for k=1:g
     fprintf('��ʼ��%d��ѧϰ\n',k);
     [G,V] = people_learn(G,Pbest,Gbest,V,w,c1,c2,c3,c4,c5,range,Vmax,Vmin,human_number,dimension,city_number,radius,King_number,King,List);    %����Ⱥѧϰ������λ��
     
     [Pbest] = OPKing_UpdateMemory(G,Pbest,human_number,dimension,city_number);      %��������Ⱥ����
     [EvalValue] = OP_evalCal(G,human_number,dimension,city_number);
     [~,List]=sort(EvalValue,'descend');
     [Gbest,Gbest_eval,King,King_eval,King_number]=Find_Gbest_and_King(G,List,city_number,dimension,EvalValue);
     rememberSL{k,1}=King;   %��¼
     [rememberValue] = OP_evalCal(Pbest,human_number,dimension,city_number);
     [~,rememberList]=sort(rememberValue,'descend');
     [~,~,rememberKing,rememberKing_eval,~]=Find_Gbest_and_King(Pbest,rememberList,city_number,dimension,rememberValue);
     remember(k)=rememberKing_eval;
     remembertheKing{k,1}=rememberKing;
end
x=1:g;
plot(x,remember,'r');
    