function x_new = MultiNewton(fun,grad_fun,x0)
%MultiNewton ��Ԫţ�ٷ����nԪ�����Է�����
%   fun - ��������ߺ�����n*1��������
%   grad_fun - ��������ߺ�����Jacobi����n*n
%   x0 - ��ʼֵ��n*1����
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

