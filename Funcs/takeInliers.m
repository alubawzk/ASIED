%函数6
%连通性分析，对圆周上的内点进行提纯
%输入
%x：椭圆周上的内点(x,y),设为inlier_n x 2 
%center：一个椭圆的中心(x,y) 1x2
%tbins: 分区 = min( 180 , pi*(1.5*(a+b)-sqrt(a*b)) ) 
%输出
%idx：为与x一样长的，inlier_n x 1的logical向量，返回有效的满足一定连通长度的内点，对应位置有效则为1，否则为0
function idx = takeInliers(x, center, tbins)
   [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%得到[-pi,pi]的方位角，等价于 theta = atan2(x(:, 2) - center(2) , x(:, 1) - center(1)); 
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%将内点分区到[1 tbins]
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);%h为直方图[1 tbins]的统计结果
    mark = zeros(tbins, 1);
    compSize = zeros(tbins, 1);
    nComps = 0;
    queue = zeros(tbins, 1);
    du = [-1, 1];
    for i = 1 : tbins
        if (h(i) > 0 && mark(i) == 0)%如果落在第i个分区内的值大于0，且mark(i)为0
            nComps = nComps + 1;
            mark(i) = nComps;%标记第nComps个连通区域
            front = 1; rear = 1;
            queue(front) = i;%将该分区加入队列，并以此开始任务
            while (front <= rear)
                u = queue(front);
                front = front + 1;
                for j = 1 : 2
                    v = u + du(j);
                    if (v == 0)
                        v = tbins;
                    end
                    if (v > tbins)
                        v = 1;
                    end
                    if (mark(v) == 0 && h(v) > 0)
                        rear = rear + 1;
                        queue(rear) = v;
                        mark(v) = nComps;%标记第nComps个连通区域
                    end
                end
            end
            compSize(nComps) = sum(ismember(tt, find(mark == nComps)));%得到构成连通域为nComps的内点数量
        end
    end
    compSize(nComps + 1 : end) = [];
    maxCompSize = max(compSize);
    validComps = find(compSize >= maxCompSize * 0.1 & compSize > 10);%大于等于最大连通长度的0.1倍的连通区域是有效的
    validBins = find(ismember(mark, validComps));%有效的分区
    idx = ismember(tt, validBins);%有效的内点
end