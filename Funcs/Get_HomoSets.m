function [HomoSets] = Get_HomoSets(ellipses, Mahala_Dist_Tolerance, coeff)
%找到候选椭圆中的同源子集
%   输入：
%   ellipse(i,5)：候选椭圆集合
%   Mahala_Dist_Tolerance：距离阈值
%   输出：
%   HomoSets{n}：n个同源集合

elps_num = size(ellipses, 1);
candidates_labels = false(elps_num,1);
num = 1;
for i=1:elps_num
    % 判断椭圆是否使用过了
    if candidates_labels(i) == true
        continue;
    end
    if i==elps_num
        HomoSets{num} = i;
        num = num + 1;
        return;
    end
    
    % 定义与 ith 椭圆同源的label
    Homo_labels = false(elps_num,1);
    
    % 初始化ith参数
    ci = ellipses(i,1:2);
    axesi = ellipses(i,3:4);
    thetai = ellipses(i,5);
    Homo_labels(i) = true;
    for j=i+1:elps_num
        cj = ellipses(j,1:2);
        axesj = ellipses(j,3:4);
        thetaj = ellipses(j,5);
        if norm(ci-cj) > 5 || norm(axesi-axesj) > 5 || norm(thetai-thetaj) > 5*(pi/180) % 加速计算，跳过明显不可能是同源的目标5 5 5
            continue
        end
        
        % 判断椭圆是否使用过了
        if candidates_labels(j) == true
            continue;
        end
        
        % 判断ith与jth椭圆的距离，如果距离合格将jth椭圆加入Homolabels
        if Judge_TwoEllipses_Dist(ellipses(i,:), ellipses(j,:), Mahala_Dist_Tolerance, coeff(i,:), coeff(j,:))
             Homo_labels(j) = true;
        end
    end
    
    % 将用过的elp label设为true，后续循环不再使用
    candidates_labels(Homo_labels) = true;
    
    % 将同源elps的索引输出
    index = find(Homo_labels);
    HomoSets{num} = index;
    num = num + 1;
end


end

