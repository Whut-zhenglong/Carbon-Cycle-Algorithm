function [bast,myneed] = CCA()
Psize = 50;
dim = 38;
% lb = ones(1,dim).*-10;
% ub = ones(1,dim).*10;
% funcobj = @(x)sum(abs(x))+prod(abs(x));
F_num =  15;
FuncName = @BenchmarkFunction;
[~,lb,ub] = FuncName(zeros(1,dim),F_num);
[organicsinit,fitnessinit] = Initialize(Psize,dim,lb,ub,FuncName,F_num);
t=1;
maxiter = 1000;
fitness = fitnessinit;
organics = organicsinit;
organicstoatmo = [];
organicstoanimal = [];
organicstodecom = [];
% f = figure;
while t<=maxiter

     
    
    
    eatsize = 0.6;
    Responsematrix = [organics,fitness];
    [~,w] = size(Responsematrix);
    Sortmatrix = sortrows(Responsematrix,w);
    organicstoatmo = Sortmatrix(1,1:end-1);
    fitnesstoatmo = Sortmatrix(1,end);
    organicstoanimal = Sortmatrix(2:Psize*eatsize+1,1:end-1);
    fitnesstoanimal = Sortmatrix(2:Psize*eatsize+1,end);
    organicstodecom = Sortmatrix(Psize*eatsize+2:Psize,1:end-1);
    fitnesstodecom = Sortmatrix(Psize*eatsize+2:Psize,end);
%     fitnesstemp = sort(fitness);
%     for i = 2:Psize*eatsize+1
%         indextemp = find(fitness == fitnesstemp(i,:));
%         a = carbon(indextemp,:);
%         carbontoanimal(i-1,:) = carbon(indextemp,:);
%     end
    
%     for i = Psize*eatsize+2:Psize
%         indextemp = find(fitness == fitnesstemp(i,:));
%         carbontodecom(i-Psize*eatsize-1,:) = carbon(indextemp,:);
%     end
%     
    
    %进行植物呼吸作用
    [carbonformplant2atmo,fitnessformplant2atmo] = plantbreath(FuncName,organicstoatmo,lb,ub,Psize,dim,F_num,fitnesstoatmo);
    %进行植物死亡
    [carbonfromplant2decom,fitnessfromplant2decom] =  plantdie(FuncName,organicstodecom,t,maxiter,lb,ub,size(organicstodecom,1),dim,F_num,fitnesstodecom);
    %植物被动物捕食
    [carbonformplant2animal,fitnessformplant2animal]  = predation(FuncName,organicstoanimal,organicstoatmo,lb,ub,size(organicstoanimal,1),dim,F_num,fitnesstoanimal);
    [len,~]=size(carbonformplant2animal);
    carbonforanimalbeath = carbonformplant2animal(1:round(len/2),:);
    fitnessforanimalbeath = fitnessformplant2animal(1:round(len/2),:);
    carbonforanimaldei = carbonformplant2animal(round(len/2)+1:len,:);
    fitnessforanimaldei  = fitnessformplant2animal(round(len/2)+1:len,:);
    %进行动物呼吸
    [carbonformanimal2atmo,fitnessformanimal2atmo] = animalbreathe(FuncName,carbonforanimalbeath,lb,ub,Psize,dim,F_num,fitnessforanimalbeath);
    
    %进行动物死亡
    [carbonfromanimal2decom,fitnessfromanimal2decom] = animaldie(FuncName,carbonforanimaldei,t,maxiter,lb,ub,size(carbonforanimaldei,1),dim,F_num,fitnessforanimaldei);
    
    carbonindecom = [carbonfromplant2decom;carbonfromanimal2decom];
    fitnessindecom = [fitnessfromplant2decom;fitnessfromanimal2decom];
    %进行分解作用
%     carbonformdecom2atmo  = decomposition(FuncName,carbonindecom,t,carbontoatmo,lb,ub,dim,size(carbonindecom,1));
    [carbonformdecom2atmo,fitnessformdecom2atmo]  = decomposition(FuncName,carbonindecom,t,organicstoatmo,lb,ub,size(carbonindecom,1),dim,F_num,fitnessindecom);
    carboninatmo = [carbonformplant2atmo;carbonformanimal2atmo;carbonformdecom2atmo];
    fitnessininatmo = [fitnessformplant2atmo;fitnessformanimal2atmo;fitnessformdecom2atmo];
    organics = carboninatmo;
    fitness = fitnessininatmo;
    for i=1:size(organics,2)

        if all(~diff(organics(:,1)'))
            temp(i) = 1;
        else
            temp(i) = 0;
        end
    end
    if all(temp==1)
        break
    end
    Responsematrixtemp = [organics,fitness];
    [~,w] = size(Responsematrixtemp);
    Sortmatrixtemp = sortrows(Responsematrixtemp,w);
    myneed(t) = Sortmatrixtemp(1,end);
%     figure(f)
%     scatter(t,myneed(t))
%     hold on
    t=t+1;
   
end

Responsematrix = [organics,fitness];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
bast.carbon = Sortmatrix(1,1:end-1);
bast.fitness = Sortmatrix(1,end);

end