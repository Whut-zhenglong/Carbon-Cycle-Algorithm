function fitness = Evaluation(functionname,variable,F_num)
[Ps, ~] = size(variable);

fitness = zeros(Ps,1);

for i = 1:Ps
    fitness(i,:) = functionname(variable(i,:),F_num);
end

 end