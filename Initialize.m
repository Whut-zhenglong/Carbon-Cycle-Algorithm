function [organics,fitness] = Initialize(Psize,Dim,LB,UB,FuncName,F_num)
    
    organics = repmat(LB,Psize,1)+rand(Psize,Dim).*repmat((UB-LB),Psize,1);

    fitness = Evaluation(FuncName,organics,F_num);
end