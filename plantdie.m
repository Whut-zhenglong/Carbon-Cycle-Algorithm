function [carbonfromplant2decom,fitnessfromplant2decom] = plantdie(functionname,carboninit,internow,maxinter,LB,UB,Ps,Dim,F_num,fitnessinit)
T = (maxinter-internow)/maxinter;
a = 1.5;
f = (a/(1-maxinter))*internow+a-a/(1-maxinter);
Beta = (f-f*internow/maxinter)^2;
for i=1:size(carboninit,1)
    temp(i,:) = (1-(1+T)*rand()).*carboninit(1,:);
    carbontemp(i,:) = temp(i,:)+(carboninit(i,:)-temp(i,:))*exp(-Beta*T);
end
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
carbonfromplant2decom = Sortmatrix(1:Ps,1:end-1);
fitnessfromplant2decom = Sortmatrix(1:Ps,end);
end