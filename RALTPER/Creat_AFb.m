function [A,F,b] = Creat_AFb(x_pre_set, xP_pre_set, alpha, radius, N)
    
    % 计算置信度
    c = -2*log(alpha);

    A = zeros(2,2,N);
    F = zeros(2,2,N);
    b = zeros(2,N);

    for k = 1:N
        % 计算椭圆的特征值和特征向量
        Sigma = xP_pre_set(:,:,k);
        [eigenVec, eigenVal] = eig(Sigma);
        
        % 计算椭圆的半长轴和半短轴
        a2 = sqrt(c * eigenVal(1, 1)) + radius;
        b2 = sqrt(c * eigenVal(2, 2)) + radius;

        A(:,:,k) = eigenVec;
        % F(:,:,k) = [a2,0;0,b2];
        F(:,:,k) = [1/(a2*a2),0;0,1/(b2*b2)];
        b(:,k) = x_pre_set(1:2,k);
        
    end
end