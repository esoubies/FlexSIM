function [C, label, J] = kmeans(I, k) 
[m, n, p] = size(I);%图片的大小m*n，p代表RGB三层
X = reshape(double(I), m*n, p);
rng('default');
C = X(randperm(m*n, k), :);%随机选三个聚类中心
J_prev = inf; iter = 0; J = []; tol = 1e-11;%容忍度tol，inf为无穷大
 
while true
    iter = iter + 1;
    dist = sum(X.^2, 2)*ones(1, k) + (sum(C.^2, 2)*ones(1, m*n))' - 2*X*C';%计算图片中各个点到K个聚类中心的距离
    [~,label] = min(dist, [], 2) ;  %label记录最小值的行数
    for i = 1:k
       C(i, :) = mean(X(label == i , :)); %取新的k个聚类中心
    end
    J_cur = sum(sum((X - C(label, :)).^2, 2));%距离之和
    J = [J, J_cur];

    if norm(J_cur-J_prev, 'fro') < tol% A和A‘的积的对角线和的平方根，即sqrt(sum(diag(A'*A)))，本次与上次距离之差
        break;
    end
    if (iter==10)% A和A‘的积的对角线和的平方根，即sqrt(sum(diag(A'*A)))，本次与上次距离之差
        break;
    end
    J_prev = J_cur;
end
