function [ solGi, meGi , deGi,pcGi,posGi,time  ]  = LR_LS_GA(X,Y,m,Q,I,pc,C, mutacion,numerOfGenerations)
    %LR_LS_AC Summary of this function goes here
    %   Detailed explanation goes here
    %   X:  Data
    %   Y:  Results
    %   m:  puntos para entrenamiento
    %   Q:  Cambio en el punto de corte
    %   pc: puntos de corte que definen la base
    %   I:  Iteraciones
    %   C:  Cantidad de validaciones cruzadas
    % ------Variables Genetico------
    %   mutacion:  Porcentaje de mutacion de los genes en el cruzamiento.
    %   numerOfGenerations: Cantidad de generaciones o iteraciones para hacer el mejoramiento de la especie
    fprintf('\n *Algoritmo Genetico.');
    
    tic;
    %   M:  Muestras
	%   N:  Atributos
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

    % Variables ALGORITMO GENETICO:
    % Poblacion tours, filas: individuos (cromosomas), columnas: genes
    PobTour=zeros(I,N);
    PobMe=zeros(I,1);
    PobDe=zeros(I,1);
    
    fprintf('\n     -Iteración: ');
    for iq = 1:Q
        fprintf(num2str(iq));
        fprintf(',');
        %iq
        PobTour=zeros(I,N);
        PobMe=zeros(I,1);
        PobDe=zeros(I,1);
         for it = 1:I
            X = X(:,sol);         
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);        
            % Obtiene la poblacion inicial de tours en que se generaron las
            % combinaciones de atributos para ese corte. y guarda el cromosoma
            % y la evaluacion de este (media y desviacion).
            PobTour(it,:)=sol;
            PobMe(it)=me;
            PobDe(it)=de;

            p1 = randperm(N,1);
            p2 = randperm(N,1);
            sol([p2 p1]) = solG([p1 p2]);
         end

        %ALGORITMO GENETICO:
         for generation=1:numerOfGenerations

            %SELECCIONAR -> ¿probabilidad de elecccion de una ciudad con 
            %mejores genes al azar o los mejores pobladores o uno a alzar y 
            %otro que sea un poblador con buena genetica?

            %azar: Se elemina el elemento para eviatar que sea escogido el
            %elemento tambien como elemento mejor. Se eliminan tambien los ratros
            %en los vectores que guardan las evaluaciones de su genetica.
            ind1=randi(N);
            tour1= PobTour(ind1,:);
            PobTour(ind1,:)=[];
            tme=PobMe(ind1,:);
            PobMe(ind1,:)=[];
            tde=PobDe(ind1,:);
            PobDe(ind1,:)=[];
            %El mejor poblador: una vez elegido se reingresa el elmento
            %escogido al azar.
            [~,ind2]=max(PobMe);
            tour2= PobTour(ind2,:);
            PobTour(end+1,:)=tour1;
            PobMe(end+1,:)=tme;
            PobDe(end+1,:)=tde;

            %PobTour(ind2,:)=[];
            %N=N-2;

            %CRUZAR -> Para el cruzamiento se copia una de las mitades de un 
            %cromosoma en el otro a la vez que son reorganizados (mutados) el
            %resto de elementos de la otra mitad para que no pierda sentido el
            %tour

            %Aqui copia la primera mitad del cromosoma del tour 1 en el tour2 y se
            %reorganizar la segunda mitad del tour2 para mutar y mentener la
            %consistencia del tour
            ttour=tour2;
            for gen=1:round(N/2)
                tour2(tour2==tour1(gen))=tour2(gen);
                tour2(gen)=tour1(gen);
            end
            %acá hace lo mismo pero copiando la segunda mitad del tour2 Original
            %en el tour 1.
            for gen=round(N/2+1):N
                tour1(tour1==ttour(gen))=tour1(gen);
                tour1(gen)=ttour(gen);
            end

            %MUTAR
            %Muta cierto porcentaje de los genes, ciudades, en los dos
            %cromosomas involucrados
            mut1=randi(N,1,int64(N*mutacion));
            mut2=randi(N,1,int64(N*mutacion));
            repo=tour1(mut1);
            tour1(mut1)=tour1(mut2);
            tour1(mut2)=repo;

            mut1=randi(N,1,int64(N*mutacion));
            mut2=randi(N,1,int64(N*mutacion));
            repo=tour2(mut1);
            tour2(mut1)=tour2(mut2);
            tour2(mut2)=repo;
            %Evaluar -> si los nuevos pobladores son mejor que algun ya
            %existente entonces desplazar a este y ocupar su lugar.
            X = X(:,tour1);         
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);
            [minv,ind]=min(PobMe);
            if (minv<=me)
                PobMe(ind)=me;
                PobDe(ind)=de;
                PobTour(ind,:)=tour1;
            end
            X = X(:,tour2);         
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);   
            [minv,ind]=min(PobMe);
            if (minv<=me)
                PobMe(ind)=me;
                PobDe(ind)=de;
                PobTour(ind,:)=tour1;
            end

            %Esta parte donde iria ? "iq" hace parte del ciclo más externo, por lo que 
            %ubicarlo aquí no me parece la mejor opcion, sin embargo en los
            %anteriores algoritmos el profesor me ha dicho que lo deje allí,
            %así que bueno, allí se queda hasta nueva orden.

            if(me<meG)
                [~,ind]=max(PobMe);
                meG = PobMe(ind); 
                deG = PobDe(ind);
                solG = PobTour(ind,:);
                pcG = c;
                posG = pos;
            end
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
    
    fig = figure;
    %plot(meGi,'--*');
    title('Algoritmo Genético - Media del error');

    fig = figure;
    %plot(deGi,'--*');
    title('Algoritmo Genético  - desviación del error');
    
    time=toc;
end

