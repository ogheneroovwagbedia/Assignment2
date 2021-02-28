function I = totalI(nx, ny, V1, sigma_conduct, sigma_insulate, Wb, Lb)
%TOTALI Calculate the total current for given paramaters

    cMap = sigma_conduct*ones(nx, ny);
    cMap(1:Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
    cMap((1:Wb)+nx-Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
    V = numericSolution(nx, ny, cMap, Inf, Inf, 0, V1);
    [Ex, Ey] = gradient(V);
    Ex = -Ex;
    Ey = -Ey;
    Jx = cMap.*Ex;
    I = (abs(sum(Jx(1,:))) + abs(sum(Jx(nx,:))))/2;
end