function carbonformdecom2atmo  = decomposition(functionname,carboninit,internow,carbonformplant2atmo,UB,LB,Dim,Ps)
beta0 = 2;
beta = beta0*exp(-internow);
a = 2;
carbontemp = carboninit - beta*(carboninit-carbonformplant2atmo)+a*(rand()-0.5);
carbontemp = (carbontemp<=repmat(UB,Ps,1)).*carbontemp+(carbontemp>repmat(UB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));
carbontemp = (carbontemp>=repmat(LB,Ps,1)).*carbontemp+(carbontemp>repmat(LB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));
fitness = Evaluation(functionname,carbontemp);
fitnessinit = Evaluation(functionname,carboninit);
[len,~] = size(fitnessinit);
for i = 1:len
    if fitness(1,:)<fitnessinit(i,:)
        carbonformdecom2atmo(i,:)  = carbontemp(i,:);
    else
        carbonformdecom2atmo(i,:)  = carboninit(i,:);
    end
end
end