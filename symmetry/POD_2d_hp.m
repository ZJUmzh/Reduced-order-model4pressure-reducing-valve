function [An,Ds,phi,U0x] = POD_2d_hp(samp_data);

nx = size(samp_data,1);
ny = size(samp_data,2);
nxy = nx*ny;
nh = size(samp_data,3);
np = size(samp_data,4);
n = nh*np;
ALL = zeros(n,nxy);
for ih = 1:nh
    for ip = 1:np
        ALL(ih*np+ip-np,:) = reshape(samp_data(:,:,ih,ip),1,nxy);
    end
end
%ALL = ALL([17,19,20],:);
U0x = mean(ALL,1);%均值
U_m = ALL-U0x;%去除均值

[U,S,phi] = svd(U_m,'econ');%svd decomposition

An = U*S;%得到一个n阶的数
Ds = diag(S).^2/n;%特征值

end