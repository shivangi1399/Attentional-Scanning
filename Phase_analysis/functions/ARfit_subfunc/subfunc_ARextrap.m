function [exdat] = subfunc_ARextrap(dat,extrams,fsamprat,A,C)

origsamps = numel(dat); % Number of samples in the to-be-extrapolated signal

arord = numel(A); % Order of AR

exdat = [dat nan(1,round(extrams/fsamprat))];

rng('shuffle','multFibonacci')
%rng('shuffle','twister')

% Extrapolate

for es = 1 :round(extrams/fsamprat)
    
    currsamp = origsamps+es; % Location of new sample in the vector
    
    % For a n order AR model a with noise variance c, value x at time t is given by the
    % following equation : x(t) = a(1)*x(t-1) + a(2)*x(t-2) + ... +
    % a(n-1)*x(t-n+1) + a(n)*x(t-n) + sqrt(c)*randnoise
    
    exdat(1,currsamp) = sum(A.*fliplr(exdat((currsamp-arord):(currsamp-1)))) ...
        +sqrt(C)*randn(1,1); %normrnd(0,0.1);%(-1 + (1+1)*rand(1,1));
        
end