clc;clear all;
%csvwrite('AntiVirus.csv', inputs)

%filename = 'DataSets/forestfires.csv';
%filename = 'DataSets/bank.csv';
%filename = 'DataSets/AntiVirus.csv';

%lines=25;
%filename = 'DataSets/movement_libras.csv';
%[numbers, TEXT, everything] =xlsread(filename);
%[M,N]=size(numbers);
%if(M>lines)
%   M=lines;
%end
%X = numbers(:,1:N-1);     Y = numbers(:,N);

rng(10);
N = 8;     M = 30;
X = randn(M,N);     Y = randn(M,1);

%LS:
%X = numbers(:,1:N-1);     Y = numbers(:,N);
m=0.4;
Q = 30;
I = 100; 
pc = 5;
C = 10;
[LSsolGi,LSmeGi,LSeGi,LSpcGi,LSposGi,LStime]=LR_LS(X,Y,m,Q,I,pc,C);
fprintf(num2str(LStime));

%ACO:
%X = numbers(:,1:N-1);     Y = numbers(:,N);
m=0.4;  
Q = 30;
I = 100; 
pc = 5;
C = 10;
numberOfAnts=200;
coefPheromone=0.1;
beta=1;
alfa=1;
[ACOsolGi,ACOmeGi,ACOdeGi,ACOpcGi,ACOposGi,ACOtime]=LR_LS_AC(X,Y,m,Q,I,pc,C,numberOfAnts,coefPheromone,beta,alfa);
fprintf(num2str(ACOtime));

%GA:
%X = numbers(:,1:N-1);     Y = numbers(:,N);
m=0.4;
Q = 30;
I = 100; 
pc = 5;
C = 10;
mutacion=0.2;
numerOfGenerations=200;
[GAsolGi,GAmeGi,GAdeGi,GApcGi,GAposGi,GAtime]=LR_LS_GA(X,Y,m,Q,I,pc,C, mutacion,numerOfGenerations);
fprintf(num2str(GAtime));

%SA:
%X = numbers(:,1:N-1);     Y = numbers(:,N);
m=0.4;
Q = 30;
I = 100; 
pc = 5;
C = 10;
IniTemp=500;
coolingRate=0.80;
minTemp=10^-2;
I = 20; 
[SAsolGi,SAmeGi ,SAdeGi,SApcGi,SAposGi,SAtime]=LR_LS_SA(X,Y,m,Q,I,pc,C,IniTemp,coolingRate,minTemp);
fprintf(num2str(SAtime));




