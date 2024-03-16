function flag = Judge_Ellipse(ellipse_i, p_j)
    flag = false;
    N = size(p_j,1);
    s = Score(ellipse_i,p_j);
    H = s/N;
    if H >= 0.95
        flag = true;
    end
end