function [ solGi, meGi , deGi,pcGi,posGi,time  ]  = LR_LS_AC(X,Y,m,Q,I,pc,C,numberOfAnts,coefPheromone,beta,alfa)
    %LR_LS_AC Summary of this function goes here
    %   Detailed explanation goes here
    %La matriz de feromonas se crea inicialmente (Hormigas que vagan 
    %aleatorieamente en busca de aliemntos) con los recorridos aleatorios.
    %Estos luego van a ser recorridos por hormigas con capacidad de desicion,
    %las cuales van tienen se paran en un ciudad aletoria y evaluan todos los
    %posibles CiudadesEnTour-1 caminos teniendo en cuenta la probabilidad y
    %escogiendo la de mayor atractivo, guardandola en el tour, luego calcular
    %con validacaion cruzada la nueva cantidad de feromonas a aplicar y
    %aplicarla al tour; luego evaporar las feromonas.
    %   X:  Data
    %   Y:  Results
    %   m:  puntos para entrenamiento
    %   Q:  Cambio en el punto de corte
    %   pc: puntos de corte que definen la base
    %   I:  Iteraciones
    %   C:  Cantidad de validaciones cruzadas
    % ------Variables ACO------
    %   numberOfIterations: Numero de veces que se manda el batallon de hormigas.
    %   numberOfAnts: Numero de hormigas para recorrer los caminos.
    %   coefPheromone: Coeficiente de evaporacion de la feromona.
    %   beta: Influencia de Distancia.
    %   alfa: Influencia de Feromona.
    fprintf('\n *Colonia de Hormigas.');
    
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
    
    % ------Variables ACO------
    %   distMat:Matriz de distancias.
    %   feroMat:Matriz de feromonas.
    %   rounds: veces que un camino es recorrido.
    distMat=zeros(N);
    feroMat=zeros(N);
    rounds=zeros(N);

    fprintf('\n     -Iteración: ');
    for iq = 1:Q  
        fprintf(num2str(iq));
        fprintf(',');
        %iq
        for it = 1:I
            X = X(:,sol);         
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);        
            % Obtiene la matriz de distancias (Error medio), en el orden en que se generaron las
            % combinaciones de atributos para ese corte.
            for dis=1:(N-1)
                distMat(sol(dis),sol(dis+1))=distMat(sol(dis),sol(dis+1))+me;
                feroMat(sol(dis),sol(dis+1))=feroMat(sol(dis),sol(dis+1))+1/me;
                rounds(sol(dis),sol(dis+1))=rounds(sol(dis),sol(dis+1))+1;
            end 
            distMat(sol(N),sol(1))=distMat(sol(N),sol(1))+me;        
            feroMat(sol(N),sol(1))=feroMat(sol(N),sol(1))+1/me;        
            rounds(sol(N),sol(1))=rounds(sol(N),sol(1))+1;  

            p1 = randperm(N,1);
            p2 = randperm(N,1);
            sol([p2 p1]) = solG([p1 p2]);
        end
        %Obtiene la media del peso de pasar de un estado a otro (distancia)
        distMat=distMat./rounds;
        %asigna el resto de valores como la media de conjunto generado
        distMat(isnan(distMat))=0;
        distMat(distMat==0)=sum(sum(distMat))/nnz(rounds);

        % ACO    
        %Cantidad de Feromonas a depositar
        Qu=sum(sum(1./distMat));
        for ant=1:numberOfAnts-1
            tour=zeros(1,N);
            tourInit=randperm(N);
            %Elige la ciudad de Origen
            antSteep=tourInit(randi(N));
            %Coloca la primera ciudad en el tour
            tour(1)=antSteep;
            %Elimina de la ciudad del conjunto disponible de ciudades
            tourInit(tourInit==antSteep)=[];        
            %Valor ocumulado inicial del tour
            tourFero=0.00000001;
            tourDist=0.00000001;
            %Creando el tour de la hormiga ant
            for city=1:N-1
                maxProb=-1;
                opDist=9999;
                opFero=0;
                opCity=-1;
                %Elige entre todos los posibles caminos segun distancias y
                %feromonas.
                for cityAnt=1:N-nnz(tour)
                    dist=distMat(antSteep,tourInit(cityAnt));
                    fero=feroMat(antSteep,tourInit(cityAnt));
                    prob=(1/dist^beta)*fero^alfa / (tourDist^beta*tourFero^alfa);
                    if (prob>maxProb)
                        maxProb=prob;
                        opDist=dist;
                        opFero=fero;
                        opCity=tourInit(cityAnt);
                    end
                end 
                %Actualiza la foromna el tour y la distancia del mismo
                tourDist= tourDist+ opDist;
                tourFero=tourFero+opFero;
                tour(city+1)=opCity; 
                %Elminar la ciudad elegida de la lista de ciudades
                tourInit(tourInit==opCity)=[];
            end        
            %Cancula el peso real elegir ese camino
            X = X(:,tour);
            [me,de] = perform_cross_validation(X,Y,C,pc,m,c);      

            pheromoneToAdd =Qu/tourDist;
            %Actualizar la matriz de feromonas ... y la distancias en caso que
            %que se hayan seleccionado ciudades nunca antes visitadas.
            for city=1:N-1
                %revisar la agregacion de feromonas
                feroMat(tour(city),tour(city+1)) = feroMat(tour(city),tour(city+1))+pheromoneToAdd;
                %Cuenta el paso sobre un camino para recalcular la media de las
                %dintacias de ese camino
                rounds(tour(city),tour(city+1))=rounds(tour(city),tour(city+1))+1;
                %Actualiza la media de las distancias de los caminos.
                distMat(tour(city),tour(city+1))=( distMat(tour(city),tour(city+1))+ me ) / rounds(tour(city),tour(city+1));
            end        
            feroMat(tour(N),tour(1)) = feroMat(tour(N),tour(1))+pheromoneToAdd;
            rounds(tour(N),tour(1))=rounds(tour(N),tour(1))+1;
            distMat(tour(N),tour(1))=( distMat(tour(N),tour(1))+ me ) / rounds(tour(N),tour(1));

            %Evaporación de la Feromona 
            feroMat = (1-coefPheromone).*feroMat;

            %Guarda las variables de interés 
             if(me<meG)
                meG = me;
                deG = de;
                solG = tour;
                pcG = c;
                posG = pos;
            end
        end
        meGi(iq) = meG; 
        deGi(iq) = deG; 
        solGi(iq)={tour};
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
    title('Colonia de Hormigas - Media del error');

    fig = figure;
    %plot(deGi,'--*');
    title('Colonia de Hormigas - desviación del error');
    
    time=toc;
end

