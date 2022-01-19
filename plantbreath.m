function [carbonformplant2atmo,fitnessformplant2atmo] = plantbreath(functionname,organicsinit,LB,UB,Ps,Dim,F_num,fitnessinit)
%PLANTRESPIRATION 此处显示有关此函数的摘要
%   此处显示详细说明


switchrate = ceil(Ps*(rand*0.5+0.5))*5;
beta = 1.5;
theta = ((gamma(1+beta)*sin(beta*pi()/2))/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = normrnd(0,theta);
v = normrnd(0,1);
S = u./(abs(v).^(1/beta));
carbontemp = organicsinit(1,:)+ (0.1).*rand(switchrate,Dim).*S.*(UB-LB);
carbontemp = (carbontemp<=repmat(UB,switchrate,1)).*carbontemp+(carbontemp>repmat(UB,switchrate,1)).*(repmat(LB,switchrate,1)+rand(switchrate,Dim).*repmat(UB-LB,switchrate,1));
carbontemp = (carbontemp>=repmat(LB,switchrate,1)).*carbontemp+(carbontemp<repmat(LB,switchrate,1)).*(repmat(LB,switchrate,1)+rand(switchrate,Dim).*repmat(UB-LB,switchrate,1));


% for i=1:size(carbontemp,1)
%     Flag4ub=carbontemp(i,:)>ub;
%     Flag4lb=carbontemp(i,:)<lb;
%     carbontemp(i,:)=(carbontemp(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% end
fitness = Evaluation(functionname,carbontemp,F_num);
% fitnessinit = Evaluation(functionname,organicsinit,F_num);
carbonnew = [carbontemp;organicsinit];
fitnessnew = [fitness;fitnessinit];

Responsematrix = [carbonnew,fitnessnew];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
carbonformplant2atmo = Sortmatrix(1,1:end-1);
fitnessformplant2atmo= Sortmatrix(1,end);
end

