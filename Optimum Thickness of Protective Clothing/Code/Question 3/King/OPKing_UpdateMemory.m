function [Pbest] = OPKing_UpdateMemory(G,Pbest,human_number,dimension,city_number)
%EVAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   ����G������ÿ���������Ӧ�Ⱥ������õ�ÿ���������ʷ����λ�ú͵�ǰ��Ⱥ������λ��

[P_Eval] = OP_evalCal(G,human_number,dimension,city_number);                  %���㵱ǰ��Ⱥ����Ӧ��
[PB_EvalValue] = OP_evalCal(Pbest,human_number,dimension,city_number);        %������Ⱥ���λ�õ���Ӧ��
for i=1:human_number
    for j=1:city_number
        if PB_EvalValue(i,j)<P_Eval(i,j)
            Pbest(i,dimension*j-1:dimension*j)=G(i,dimension*j-1:dimension*j);             %���µ������Ӽ���
        end
    end
end
%�õ��µ���Ⱥ���λ��




end

