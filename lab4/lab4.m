N = 8;
Nx = N;
Ny = N;
Nz = N;
hx = 1/(Nx-1);
hy = 1/(Ny-1);
hz = 1/(Nz-1);

f = zeros(Nx,Ny,Nz);

% точное решение - для сравнения
f_exact = zeros(Nx,Ny,Nz);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            f_exact(i,j,k)  = func(i,j,k,hx,hy,hz);
        end
    end
end


% граничные условия по Z
for i = 1:Nx
    for j = 1:Ny
        f(i,j,1)  = func(i,j,1,hx,hy,hz);
        f(i,j,Nz) = func(i,j,Nz,hx,hy,hz);
    end
end
X = 1:Nx;% набор точек по Х
Y = 1:Ny;% набор точек по Y
%z = f(:,:,Nz/2); 
figure
imagesc(X,Y,reshape(f(Nx/2,:,:),Ny,Nz))
colorbar('eastoutside')
title('Граничные условия по Z для сечения Nx/2 ')
xlabel('X')
ylabel('Y')

% граничные условия по Y
for i = 1:Nx
    for k = 1:Nz
        f(i,1,k)  = func(i,1,k,hx,hy,hz);
        f(i,Ny,k) = func(i,Ny,k,hx,hy,hz);
    end
end
X = 1:Nx;% набор точек по Х
Z = 1:Nz;% набор точек по Z
%z = f(:,:,Nz/2); 
figure
imagesc(X,Z,reshape(f(Nx/2,:,:),Nx,Nz))
colorbar('eastoutside')
title('Граничные условия по Y для сечения Nx/2 ')
xlabel('X')
ylabel('Z')

% граничные условия по X
for j = 1:Ny
    for k = 1:Nz
        f(1,j,k)  = func(1,j,k,hx,hy,hz);
        f(Nx,j,k) = func(Nx,j,k,hx,hy,hz);
    end
end
Y = 1:Ny;% набор точек по Y
Z = 1:Nz;% набор точек по Z
%z = f(:,:,Nz/2); 
figure
imagesc(Y,Z,reshape(f(Nx/2,:,:),Nx,Nz))
colorbar('eastoutside')
title('Граничные условия по X для сечения Nx/2 ')
xlabel('Y')
ylabel('Z')

figure
plot(f_exact(:,Nx/2,Nz/2))
legend_strings = ["exact solution"];
legend(legend_strings);
%legend;
hold on;

%вычисления по формуле (2) из описания лабораторной работы
% здесь a == 0 и rho == 6

residual = 1e6; % невязка
eps      = 1e-4; % значение невязки для выхода из цикла
f1 = zeros(Nx,Ny,Nz) % массив для хранения следующей итерации
n = 0; % номер итерации
while residual > eps
    for i=2:Nx-1
      for j = 2:Ny-1
        for k = 2:Nz-1
            f1(i,j,k) = 1/(2/hx^2+2/hy^2+2/hz^2)*(...
                           (f(i+1,j,k)+f(i-1,j,k))/hx^2 + ...
                           (f(i,j-1,k)+f(i,j+1,k))/hy^2+...
                           (f(i,j,k-1)+f(i-1,j,k+1))/hz^2 - 6);
        end
      end
    end
    % вычисление поэлементного максимального модуля разности двух матриц
    % только по внутренней части  расчетной области
    residual = 0;
    for i=2:Nx-1
      for j = 2:Ny-1
        for k = 2:Nz-1
            t = abs(f1(i,j,k) - f(i,j,k));
            if t > residual
               residual = t;
            end
        end
      end
    end
    
    plot(f(:,Nx/2,Nz/2))
    str = sprintf('iter %d',n);
    legend_strings(end+1) = str;
    %leg_str = strcat(leg_str,num2str(n));
    legend(legend_strings);
    hold on;
    
    
    n = n+1;
    for i=2:Nx-1
      for j = 2:Ny-1
        for k = 2:Nz-1
            f(i,j,k) = f1(i,j,k)
        end
      end
    end
    n,residual
end
  disp("convergence reached!!!!")
  
% вычисление погрешности как 
% поэлементного максимального модуля разности решения
% и эталонной функции func
    err = 0;
    for i=2:Nx-1
      for j = 2:Ny-1
        for k = 2:Nz-1
            t = abs(f_exact(i,j,k) - f(i,j,k));
            
            if t > err
               err = t;
               max_i = i
               max_j = j
               max_k = k
            end
        end
      end
    end
    err

