function [ Uxt ] = Q1_BoundaryConditions_kChange( LeftorRight,n,LeftBoundry)
%�ڴ˺����ж����ֵ����
%LeftOrRight��ʶ��߽����������ұ߽����� 0������1������

if LeftorRight==0
        Uxt=75;
elseif LeftorRight==1
    kk=floor(n/100)+1;
    Uxt=LeftBoundry(kk);
else
    %��Ӧ�ó��ִ���������쳣
    error('�߽������������');
end

end