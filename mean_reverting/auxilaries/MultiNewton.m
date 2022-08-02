function x_new = MultiNewton(fun,grad_fun,x0)
%MultiNewton 多元牛顿法求解n元非线性方程组
%   fun - 方程组左边函数：n*1向量函数
%   grad_fun - 方程组左边函数的Jacobi矩阵：n*n
%   x0 - 初始值：n*1向量
iter_max = 1e03;
epsilon = 1e-8;

iter = 0;
while (iter<=iter_max)
    x_new = x0 - grad_fun(x0)\fun(x0);
    fval = norm(fun(x_new));
    if(fval<epsilon)
        break
    end
    x0 = x_new;
    iter = iter + 1;
    disp(['newton_iter: ',num2str(iter),', fvalue: ',num2str(round(fval,4))]);
end

end

