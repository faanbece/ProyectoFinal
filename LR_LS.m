function [ solGi, meGi , deGi,pcGi,posGi, time ] = LR_LS(X,Y,m,Q,I,pc,C)
    %LR_LS Summary of this function goes here
    %   Detailed explanation goes here
    %   X:  Data
    %   Y:  Results
    %   m:  puntos para entrenamiento
    %   Q:  Cambio en el punto de corte
    %   pc: puntos de corte que definen la base
    %   I:  Iteraciones
    %   C:  Cantidad de validaciones cruzadas
    fprintf('\n *Búsqueda Tabú.');
    %   M:  Muestras
	%   N:  Atributos
    tic;
    [M,N]=size(X);
    m = round(M*m);
    
    %   sol: Vector solucion
    %   c: ¿?
    %   pos: ¿?
    sol = 1:N;
    c = sort(randperm(N,pc));
    pos = zeros(1,N);
    pos(c) = 1;
    
    %   meG: ¿?
    %   solG: ¿?
    %   pcG: ¿?
    %   posG: ¿?
    meG = Inf;
    solG = sol;
    pcG = c;
    posG = pos;
    
    %   meGi:
    %   deGi:
    %   solGi:
    %   pcGi:
    %   posGi:
    meGi=zeros(1,Q);
    deGi=zeros(1,Q);
    solGi=cell(Q,1);
    pcGi =cell(Q,1);
    posGi=cell(Q,1);
    
    fprintf('\n     -Iteración: ');
    for iq = 1:Q
        fprintf(num2str(iq));
        fprintf(',');
        %fprintf('\n');
        %iq
        for it = 1:I
             X = X(:,sol);
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);
            if(me<meG)
                meG = me;
                deG = de;
                solG = sol;
                pcG = c;
                posG = pos;
            end
            
            p1 = randperm(N,1);
            p2 = randperm(N,1);
            sol([p2 p1]) = solG([p1 p2]);
        end
        
        meGi(iq) = meG; 
        deGi(iq) = deG; 
        solGi(iq)={solG};
        pcGi(iq) ={pcG};
        posGi(iq)={posG};

        pos = posG;
        p3 = randperm(N,1);
        pos(p3) = ~pos(p3);

        c = find(pos>0);
        pc = length(c);
    end
    fprintf('\b');
    
    %fig = figure;
    %plot(meGi,'--*');
    %title('Buqueda Tabú - Media del error');
    %fig = figure;
    %plot(deGi,'--*');
    %title('Buqueda Tabú - desviación del error');
    
    time=toc;
end

