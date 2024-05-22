function [Homo_flag] = Judge_TwoEllipses_Dist(elp1, elp2, Homo_Dist_Tolerance, coeff1, coeff2)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
dist = dist_two__Elps(elp1, elp2, coeff1, coeff2);
%% 判断0.5*(dist1+dist2)距离, 离散度判断，等价于相似度判断
if dist >= Homo_Dist_Tolerance
    Homo_flag = true;
else
    Homo_flag = false;
end

end

