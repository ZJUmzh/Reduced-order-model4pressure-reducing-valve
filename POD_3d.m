function [An,Ds,phi,U0x] = POD_3d(data_pod);

nx = size(data_pod,1);
ny = size(data_pod,2);
nz = size(data_pod,3);
nxyz = nx*ny*nz;
n = size(data_pod,4);
ALL = zeros(n,nxyz);
for i=1:1:n
    ALL(i,:) = reshape(data_pod(:,:,:,i),1,nxyz);
end
U0x = mean(ALL,1);%均值
U_m = ALL-U0x;%去除均值

[U,S,phi] = svd(U_m,'econ');%svd decomposition

An = U*S;%得到一个n阶的数
Ds = diag(S).^2/n;%特征值

end