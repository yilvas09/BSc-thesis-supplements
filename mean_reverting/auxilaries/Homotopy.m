function x_new = Homotopy(fun,grad_fun,x0)
% Homotopy 同伦法求解非线性方程组
%   fun - 方程组左边函数：n*1向量函数
%   grad_fun - 方程组左边函数的Jacobi矩阵：n*n
%   x0 - 初始值：n*1向量

d_lambda = 1e-2; % 可以修改数值ODE的步长
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

