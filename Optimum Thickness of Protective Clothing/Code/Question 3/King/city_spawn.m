function [G] = city_spawn(human_number,dimension,city_number,range,ChangeMatrix)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   ����Ⱥ�㷨����������ά������Χ���ƣ����ɳ�ʼ����Ⱥ(����Ϊm)
%   range��һ��n*2�ľ��󣬵�һ�б�ʾ��Сֵ���ڶ��б�ʾ���ֵ

max=range(:,2)';                                  %����ά�ȵ����Χ����
min=range(:,1)';                                  %����ά�ȵ���С��Χ����

G=repmat(min*ChangeMatrix,human_number,1)+repmat(((max-min)*ChangeMatrix),human_number,1).*rand(human_number,dimension*city_number);%���ɳ�ʼ��Ⱥ

G=ceil(G*10)/10;

end

