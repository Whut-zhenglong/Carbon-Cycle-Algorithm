function [f,LB,UB] = BenchmarkFunction(x,F_num)

Dim = size(x,2);

switch F_num
    case 1
        f = sum((x).^2);
        LB = ones(1,Dim).*-100;
        UB = ones(1,Dim).*100;
    case 2
        f = sum(abs(x))+prod(abs(x));
        LB = ones(1,Dim).*-10;
        UB = ones(1,Dim).*10;
    case 3
        fun = zeros(Dim,1);
        for i = 1 : Dim
            fun(i,1)= sum(x(1:i))^2;
        end
        f = sum(fun);
        LB = ones(1,Dim).*-100;
        UB = ones(1,Dim).*100;
    case 4
        f = max(abs(x));
        LB = ones(1,Dim).*-100;
        UB = ones(1,Dim).*100;
    case 5
        fun = zeros(Dim,1);
        for i = 1 : (Dim-1)
            fun(i,1) = 100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
        end
        f = sum(fun);
        LB = ones(1,Dim).*-30;
        UB = ones(1,Dim).*30;
    case 6
        f = sum((x+0.5).^2);
        LB = ones(1,Dim).*-100;
        UB = ones(1,Dim).*100;
    case 7
        fun = zeros(Dim,1);
        for i = 1 : Dim
            fun(i,1) = i*x(i)^4;
        end
        f = sum(fun)+rand();
        LB = ones(1,Dim).*-1.28;
        UB = ones(1,Dim).*1.28;
    case 8
        f = sum(-x.*sin(abs(x).^0.5));
        LB = ones(1,Dim).*-500;
        UB = ones(1,Dim).*500;
    case 9
        f = sum(x.^2-10*cos(2*pi*x)+10);
        LB = ones(1,Dim).*-5.12;
        UB = ones(1,Dim).*5.12;
    case 10
        f = -20*exp(-0.2*(1/Dim*sum(x.^2)).^0.5)-exp(1/Dim*sum(cos(2*pi*x)))+20+exp(1);
%         f = 20+exp(1)-exp(1/Dim*sum(cos(2*pi*x)))-20*exp(-0.2*(1/Dim*sum(x.^2)).^0.5);
        LB = ones(1,Dim).*-32;
        UB = ones(1,Dim).*32;
    case 11
        fun = zeros(1,Dim);
        for i = 1 : Dim
            fun(1,i) = cos(x(i)/(i^0.5));
        end
        f = 1/4000*sum(x.^2)-prod(fun)+1;
        LB = ones(1,Dim).*-600;
        UB = ones(1,Dim).*600;
    case 12
        a = 10;
        k = 100;
        m = 4;
        y = zeros(1,Dim);
        u = zeros(1,Dim);
        for i = 1 : Dim
            y(1,i) = 1+(x(i)+1)/4;
            if x(i)>a
                u(1,i) = k*(x(i)-a).^m;
            elseif x(i)>-a
                u(1,i) = 0;
            else
                u(1,i) = k*(-x(i)-a).^m;
            end
        end
        f = (pi/Dim)*(10*sin(pi*y(1,1))^2+sum((y(1,1:Dim-1)-1).^2.*...
            (1+10*sin(pi*y(1,2:Dim)).^2))+(y(1,Dim)-1).^2)+sum(u);
        LB = ones(1,Dim).*-50;
        UB = ones(1,Dim).*50;
    case 13
        a = 10;
        k = 100;
        m = 4;
        u = zeros(1,Dim);
        for i = 1 : Dim
            if x(i)>a
                u(1,i) = k*(x(i)-a).^m;
            elseif x(i)>-a
                u(1,i) = 0;
            else
                u(1,i) = k*(-x(i)-a).^m;
            end
        end
        f = 0.1*(sin(3*pi*x(1,1)).^2+sum((x-1).^2.*(1+sin(3*pi*x+1).^2))+...
            (x(1,Dim)-1).^2*(1+sin(2*pi*x(1,Dim)).^2))+sum(u);
        LB = ones(1,Dim).*-50;
        UB = ones(1,Dim).*50;
    case 14
            J = [0.2,0.3,0.4,0.5,0.6];
    uinangle = x(1);	 
    uinshift = x(2);	
    uoutangle = x(3);	
    uoutshift = x(4);	
    linangle = x(5);	
    linshift = x(6);		
    loutangle = x(7);	
    loutshift  = x(8);
    if linshift ==0
        linshift = 0.00001;
    end
    f =  -(1.18287757075455*J(5)*sin(cos((loutangle - uinshift)*cos(cos(J(5))) + sin(J(5)) - sin(uoutangle)))...
            + 1.18287757075455*sin(cos(tan(J(5) - sin(uinangle)))) - 0.979720804709793+0.216016981165576*cos(J(5)*uoutshift*cos(loutshift)...
            - sin(linangle + linshift + loutshift*sin(linshift*sin(J(5)) + cos(loutshift/linshift) - tan(1))))...
            + 0.216016981165576*cos(tan(loutshift)) - 0.413293125924208);
        LB = [0,0,-0.11, 0, 0, 0,-0.11, 0.1];                   % Lower bound
        UB = [0.5,0.5,0.5, 0.5,0.5,0.5 ,0.5, 0.5];                     % Upper bound
    case 15
           %   判断线段与多边形某一条边是否相交，相交返回路径距离为100000，不相交返回两点欧式距离

%   输出 distance ：整条路径距离

%   输入 postion1  ：路径第i个坐标点的坐标[x11,y11]
%   输入 postion2  ：路径第i+1个坐标点的坐标点[x12,y12]
%   输入 postion1  ：多边形某边第i个坐标点的坐标[x21,y21]
%   输入 postion2  ：多边形某边第i+1个坐标点的坐标点[x22,y22]
% x=[0,-0.23,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,5,-5,-5,-5,0];
% x=[0,-3,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,5,-5,-5,-5,0];
pointx=[0,0.887200000000000,1.77430000000000,2.66150000000000,3.54860000000000,4.43580000000000,5.32290000000000,6.21010000000000,7.09720000000000,7.98440000000000,8.87150000000000,9.75870000000000,10.6458000000000,11.5330000000000,12.4202000000000,13.3073000000000,14.1945000000000,15.0816000000000,15.9688000000000,16.8559000000000,17.7431000000000,18.6302000000000,19.5174000000000,20.4045000000000,21.2917000000000,22.1788000000000,23.0660000000000,23.9532000000000,24.8403000000000,25.7275000000000,26.6146000000000,27.5018000000000,28.3889000000000,29.2761000000000,30.1632000000000,31.0504000000000,31.9375000000000,32.8247000000000,33.7118000000000,34.5990000000000];
pointy = [0,x,0];
% postion3,postion4
load obstacle_new
for  k= 1:length(pointx)-1
    postion1 = [pointx(k),pointy(k)];
    postion2 = [pointx(k+1),pointy(k+1)];

% postion1 = [0,4];%碰撞
% % postion1 = [0,10];%不碰撞
% postion2 = [3.84,8];

% %相交测试算例
% postion1 = [12.04,9.67];
% postion2 = [46.59,7.83];
% postion3 = [30.44,8.95];
% postion4 = [31.53,8.41];
% %不相交测试算例
% postion1 = [12.04,9.67];
% postion2 = [46.59,7.83];
% postion3 = [15.92,10.49];
% postion4 = [18.34,11.03];
    for i = 1:length(obstacle_new)
        obstacle_temp{i,1} = [obstacle_new{i,1},obstacle_new{i,1}(:,1)];
        for j =1:length(obstacle_new{i,1})
            postion3 = [obstacle_temp{i,1}(1,j),obstacle_temp{i,1}(2,j)];
            postion4 = [obstacle_temp{i,1}(1,j+1),obstacle_temp{i,1}(2,j+1)];
            %% 将两条线段平移，将AB左端点平移至原点
            x11 = postion1(1,1);
            y11 = postion1(1,2);
            x12 = postion2(1,1);
            y12 = postion2(1,2);
            x21 = postion3(1,1);
            y21 = postion3(1,2);
            x22 = postion4(1,1);
            y22 = postion4(1,2);

            %得到平移后的新的坐标点
            X11 = 0;
            Y11 = 0;
            X12 = x12 - x11;
            Y12 = y12 - y11;
            X21 = x21 - x11;
            Y21 = y21 - y11;
            X22 = x22 - x11;
            Y22 = y22 - y11; 


            %% 将两条线段旋转，逆时针旋转为正；顺时针旋转为正

            %旋转角度求解
            tha = atan2(Y12,X12);
            %坐标旋转矩阵
            A = [cos(tha),sin(tha);-sin(tha),cos(tha)];
            %旋转后的新坐标
            postion2_new = A*[X12;Y12];
            postion3_new = A*[X21;Y21];
            postion4_new = A*[X22;Y22];
            X12 = postion2_new(1,1);
            Y12 = postion2_new(2,1);
            X21 = postion3_new(1,1);
            Y21 = postion3_new(2,1);
            X22 = postion4_new(1,1);
            Y22 = postion4_new(2,1);

            %% 判断是否两条线相交 相交条件：1.线段2两个端点纵坐标异号；2.直线交点的位置
            collision(i,j) = 0;
    %         distance = sqrt(X12^2 + Y12 ^2);
            if Y21 * Y22 < 0  % 与X轴平行时  坐标为0判断
                PX = X22 + (X21 - X22) * Y22 / (Y22 - Y21);
                if PX > 0 && PX < X12   
    %                 distance = 100000;
                    collision(i,j) = 1;
                end
            end
        end
    end

%     if  sum(sum(collision)) == 0
%         fprintf("不碰撞\n") ;
%         distance(k) = sqrt(X12^2 + Y12 ^2); 
%     else
%         fprintf("碰撞\n");
%         distance(k) = 100000;
% %         break
            
   %% 计算路线与障碍物相交的数量，给一个惩罚值
   safe_value(k) = sum(sum(collision));
            
   %% 计算该条路线的实际距离
   distance(k) = sqrt(X12^2 + Y12 ^2);         
end
distance_all = sum(distance);
safe_all = sum(safe_value)*1000;
f=distance_all+safe_all;
        
        LB = ones(1,37).*-6;
        UB = ones(1,37).*1;
        LB = [-0.3,LB];
        UB = [-0.2,UB];
end
















