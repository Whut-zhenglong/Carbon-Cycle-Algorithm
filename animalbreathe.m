function [carbonformanimal2atmo,fitnessformanimal2atmo]  = animalbreathe(functionname,organicsinit,LB,UB,Ps,Dim,F_num,fitnessinit)

switchrate = ceil(Ps*(rand*0.5+0.5))*2;
beta = 1.5;
theta = ((gamma(1+beta)*sin(beta*pi()/2))/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = normrnd(0,theta);
v = normrnd(0,1);
S = u./(abs(v).^(1/beta));
carbonresult = [];
[len,~] = size(organicsinit);
for i=1:len
    carbontemp = organicsinit(i,:)+ 0.1.*(2.*rand(switchrate,Dim)-1).*S.*(UB-LB);
    carbontemp = (carbontemp<=repmat(UB,switchrate,1)).*carbontemp+(carbontemp>repmat(UB,switchrate,1)).*(repmat(LB,switchrate,1)+rand(switchrate,Dim).*repmat(UB-LB,switchrate,1));
    carbontemp = (carbontemp>=repmat(LB,switchrate,1)).*carbontemp+(carbontemp<repmat(LB,switchrate,1)).*(repmat(LB,switchrate,1)+rand(switchrate,Dim).*repmat(UB-LB,switchrate,1));
    carbonresult = [carbonresult;carbontemp];
end


% for i=1:size(carbontemp,1)
%     Flag4ub=carbontemp(i,:)>ub;
%     Flag4lb=carbontemp(i,:)<lb;
%     carbontemp(i,:)=(carbontemp(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% end
fitness = Evaluation(functionname,carbonresult,F_num);
% fitnessinit = Evaluation(functionname,organicsinit,F_num);
carbonnew = [carbonresult;organicsinit];
fitnessnew = [fitness;fitnessinit];

Responsematrix = [carbonnew,fitnessnew];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
carbonformanimal2atmo = Sortmatrix(1:len,1:end-1);
fitnessformanimal2atmo = Sortmatrix(1:len,end);
end