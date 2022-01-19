function [carbonformplant2animal,fitnessformplant2animal] = predation(functionname,carboninit,carbonformplant2atmo,LB,UB,Ps,Dim,F_num,fitnessinit)


carbontemp = carboninit - rand().*(carboninit-carbonformplant2atmo).*0.5;
carbontemp = (carbontemp<=repmat(UB,Ps,1)).*carbontemp+(carbontemp>repmat(UB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));
carbontemp = (carbontemp>=repmat(LB,Ps,1)).*carbontemp+(carbontemp<repmat(LB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));
fitness = Evaluation(functionname,carbontemp,F_num);
% fitnessinit = Evaluation(functionname,carboninit,F_num);
carbonnew = [carbontemp;carboninit];
fitnessnew = [fitness;fitnessinit];

Responsematrix = [carbonnew,fitnessnew];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
carbonformplant2animal = Sortmatrix(1:Ps,1:end-1);
fitnessformplant2animal = Sortmatrix(1:Ps,end);
end
