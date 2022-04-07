function [ Uxt ] = Q1_BoundaryConditions_kChange( LeftorRight,n,LeftBoundry)
%在此函数中定义边值条件
%LeftOrRight标识左边界条件还是右边界条件 0代表左，1代表右

if LeftorRight==0
        Uxt=75;
elseif LeftorRight==1
    kk=floor(n/100)+1;
    Uxt=LeftBoundry(kk);
else
    %不应该出现此情况，抛异常
    error('边界条件输入错误');
end

end