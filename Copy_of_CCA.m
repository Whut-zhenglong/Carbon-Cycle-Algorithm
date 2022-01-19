function bast = CCA()
Psize = 30;
dim = 30;
% lb = ones(1,dim).*-10;
% ub = ones(1,dim).*10;
% funcobj = @(x)sum(abs(x))+prod(abs(x));
F_num = 10;
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

while t<maxiter

     
    

    eatsize = 0.6;
    Responsematrix = [organics,fitness];
    [~,w] = size(Responsematrix);
    Sortmatrix = sortrows(Responsematrix,w);
    organicstoatmo = Sortmatrix(1,1:end-1);
    organicstoanimal = Sortmatrix(2:Psize*eatsize+1,1:end-1);
    organicstodecom = Sortmatrix(Psize*eatsize+2:Psize,1:end-1);
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
    carbonformplant2atmo = plantbreath(FuncName,organicstoatmo,lb,ub,Psize,dim,F_num);
    %进行植物死亡
    carbonfromplant2decom =  plantdie(FuncName,organicstodecom,t,maxiter,lb,ub,size(organicstodecom,1),dim,F_num);
    %植物被动物捕食
    carbonformplant2animal  = predation(FuncName,organicstoanimal,organicstoatmo,lb,ub,size(organicstoanimal,1),dim,F_num);
    [len,~]=size(carbonformplant2animal);
    carbonforanimalbeath = carbonformplant2animal(1:round(len/2),:);
    carbonforanimaldei = carbonformplant2animal(round(len/2)+1:len,:);
    
    %进行动物呼吸
    carbonformanimal2atmo  = animalbreathe(FuncName,carbonforanimalbeath,lb,ub,Psize,dim,F_num);
    
    %进行动物死亡
    carbonfromanimal2decom = animaldie(FuncName,carbonforanimaldei,t,maxiter,lb,ub,size(carbonforanimaldei,1),dim,F_num);
    
    carbonindecom = [carbonfromplant2decom;carbonfromanimal2decom];
    
    %进行分解作用
%     carbonformdecom2atmo  = decomposition(FuncName,carbonindecom,t,carbontoatmo,lb,ub,dim,size(carbonindecom,1));
    carbonformdecom2atmo  = decomposition(FuncName,carbonindecom,t,organicstoatmo,lb,ub,size(carbonindecom,1),dim,F_num);
    carboninatmo = [carbonformplant2atmo;carbonformanimal2atmo;carbonformdecom2atmo];
    
    organics = carboninatmo;
    fitness = Evaluation(FuncName,organics,F_num);
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
    t=t+1;      
end
Responsematrix = [organics,fitness];
[~,w] = size(Responsematrix);
Sortmatrix = sortrows(Responsematrix,w);
bast.carbon = Sortmatrix(1,1:end-1);
bast.fitness = Sortmatrix(1,end);

end