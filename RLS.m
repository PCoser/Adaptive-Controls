
%Parameters estimates by Recursive Least Square Algorithm
clc
clear all
close all

tic                 %measure the processing time

%Input signal
n = 400;            %number of iterations
u = idinput(n);     %input signal with binary random generator

for i = 100:300
    u(i) = 0;     %Set the input to 0 between 100 and 300 
end

Theta = transpose([1.5 -0.7 1]);     %input of parameters
Y = zeros([n 1]);                    %Initialize the system
PHI = zeros([3 n]);                  %Initialize the system 

%Initialize RLS Matrix
ThetaRLS = zeros([3,n]);       %Initialize matrix Theta for Recursive Least Squares
ThetaRLS(:,1) = [10,10,10];    %First Guess
Lambda = 0.95;                 %Forgetting coefficient
p = 100;                       %First value for Covariance matrix
P = p*eye(3);                  %Initial Covariance matrix with order = 3

for k = 3:n

    Change of parameters at n=250
    if k == 250
        Theta = transpose([1.7 -0.9 5]);
    end
 
    %Buiding PHI signal
    phi = [Y(k-1) Y(k-2) u(k-1)];
    for i = 1:length(phi)
        %build the PHI dataset
        PHI(i,k) = phi(i);
    end
    %Update the Y signal
    Y(k) = phi*Theta + sqrt(4)*randn(); %adding noise()
    
    %RLS Algorithm
    %Update P
    P = (1/Lambda)*(P - ((P*PHI(:,k)*PHI(:,k)'*P)/(Lambda + PHI(:,k)'*P*PHI(:,k))));
    L = (P*PHI(:,k))/(Lambda + PHI(:,k)'*P*PHI(:,k));
    %Estimate Theta
    ThetaRLS(:,k) = ThetaRLS(:,k-1) + L*(Y(k) - PHI(:,k)'*ThetaRLS(:,k-1));
    %Measure uncertanty by computing the trace of covariance matrix
    Trace(k) = sum(diag(P));
end

plot([1:n],ThetaRLS(:,1:n),[1:n],Trace)
title('Parameter Estimates') ;
xlabel('Time Steps')
ylabel('Coefficients')
legend('Param a1','Param a2','Param b','Trace of Cov.');
toc