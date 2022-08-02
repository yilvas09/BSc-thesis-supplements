function x_new = Homotopy(fun,grad_fun,x0)
% Homotopy ͬ�׷��������Է�����
%   fun - ��������ߺ�����n*1��������
%   grad_fun - ��������ߺ�����Jacobi����n*n
%   x0 - ��ʼֵ��n*1����

d_lambda = 1e-2; % �����޸���ֵODE�Ĳ���
epsilon = 1e-3;

n = length(x0);
lambda_vec = 0:d_lambda:1;
x_vec = [x0,zeros(n,length(lambda_vec)-1)];
for i = 2:length(lambda_vec)
   x_vec(:,i) = x_vec(:,i-1) - d_lambda * (grad_fun(x_vec(:,i-1)) \ fun(x0));
   fval = round(norm(fun(x_vec(:,i))),4);
   disp(['iter: ',num2str(i),', fvalue: ',num2str(fval)])
end
disp(['fval = ',num2str(fval)])
if(fval>epsilon)
    disp('stopped prematurely!')
end
x_new = x_vec(:,end);

end

