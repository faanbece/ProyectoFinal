function [med,des] = perform_cross_validation(X,Y,C,pc,m,c)

    %med: error medio    des: desviación estandar del error
    %C: # de veces que hago validación cruzada
    %m: cuantos puntos voy a usar en la validación cruzada
    %c: puntos de corte
    
    [M,N] = size(X);

    for k = 1:C

        P = randperm(M);

        X = X(P,:);
        Y = Y(P,:);
        
        %Para entrenamiento
        Xr = X(1:m,:);
        Yr = Y(1:m,:);

        W = zeros(m,pc);

        for i = 1:m
            xi = Xr(i,:);
            W(i,:) = fsubs(xi,c);
        end

        B = solve_system(W,Yr);

        Xv = X(m+1:M,:);
        Yv = Y(m+1:M,:);

        mv =  M-m-1+1;
        Wv = zeros(mv,pc);

        for i = 1:mv
            xi = Xv(i,:);
            Wv(i,:) = fsubs(xi,c);
        end


        E(k) = norm(Yv-Wv*B);

    end

    med = mean(E);
    des = std(E);

end