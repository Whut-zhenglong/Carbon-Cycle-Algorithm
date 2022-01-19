function [carbonformdecom2atmo,fitnessformdecom2atmo]  = decomposition(functionname,carboninit,internow,carbonformplant2atmo,LB,UB,Ps,Dim,F_num,fitnessinit)
beta0 = 1.5;
beta = beta0*exp(-internow);
a = 1.5;
carbontemp = carboninit - beta*(carboninit-carbonformplant2atmo)+a*(rand()-0.5);
carbontemp = (carbontemp<=repmat(UB,Ps,1)).*carbontemp+(carbontemp>repmat(UB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));
carbontemp = (carbontemp>=repmat(LB,Ps,1)).*carbontemp+(carbontemp<repmat(LB,Ps,1)).*(repmat(LB,Ps,1)+rand(Ps,Dim).*repmat(UB-LB,Ps,1));

% for i=1:size(carbontemp,1)
%     Flag4ub=carbontemp(i,:)>ub;
%     Flag4lb=carbontemp(i,:)<lb;
%     carbontemp(i,:)=(carbontemp(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% end
fitness = Evaluation(functionname,carbontemp,F_num);
% fitnessinit = Evaluation(functionname,carboninit,F_num);
carbonnew = [carbontemp;carboninit];
fitnessnew = [fitness;fitnessinit];

Responsematrix = [carbonnew,fitnessnew];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
carbonformdecom2atmo = Sortmatrix(1:Ps,1:end-1);
fitnessformdecom2atmo = Sortmatrix(1:Ps,end);
end