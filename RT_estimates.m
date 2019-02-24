%Parameters estimates by Least Mean Square Algorithm
clc
clear all
close all
tic
% the input X(n) is white noise with variance 1 and mean=0
N=1000;             %Number of iteration
x=randn(N,1);       %System input X(n) with white noise variance 1 and mean=0

%Matrix initialization
d=zeros(N,1);           %Output of the system
d_hat=zeros(N,1);       %Predicted output
H_vector(:,1) = [0;0];  %Parameters at each iteration
R = zeros(2,2);         %Information Matrix 
p = zeros(2,1);         %Cross Correlation Matrix
alg = 1                 %Choose the algorithm to use 1 for RLS and 0 for SD/LMS

%Initializing LMS/SD Parameters
H = [0;0];  % The parameter estimate is initialised to zero
gamma = 0.05;  % step size

%Initializing RLS Parameters
Lambda = 0.95;                 %Forgetting coefficient
P = 100*eye(2);                %Initial Covariance matrix with order = 2

%Iterative system with N iterations
for n=2:N
    
   % Data Set Generation
   d(n)=0.2*x(n-1)+0.8*d(n-1);      %Build the system decribed as d/x = 0.2/z - 0.8
   d_hat(n)=[x(n-1) d(n-1)]*H;      %Predicting the output d^ to compare with the real system
   
   if alg == 0
       %Start the Steppest Descent Algorithm
       X=[x(n-1);d(n-1)];               %Computing X(n)
       H = H + 2*gamma*X*(d(n) - X'*H); %Estimate the parameters
       H_vector(:,n) = H;               %Build the H vector to plot the results
   end
   
   if alg == 1
       %RLS Algorithm
       X=[x(n-1);d(n-1)];               %Computing X(n)
       P = (1/Lambda)*(P - ((P*X*X'*P)/(Lambda + X'*P*X)));    %Update P
       L = (P*X)/(Lambda + X'*P*X);
       H = H + L*(d(n) - X'*H);         %ThetaRLS(:,n) = ThetaRLS(:,n-1) + L*(d(n) - X'*ThetaRLS(:,n-1));
       H_vector(:,n) = H;               %Build the H vector to plot the results
   end
      
   %Sum up R and p matrix; 
   R = R + X*X';        %Information Matrix
   p = p + d(n)*X;      %Cross correlation Matrix
end

%Computing the matrix R and p matrix
R = R/N;        %Using the sum made on the loop and dividing by the number of points
p = p/N;        %Using the sum made on the loop and dividing by the number of points

%Ploting the MSE Based on the predicted R and p
step = 100;
h1 = linspace(-0.5,1,step)';    %defining parameter a
h2 = linspace(-0.5,1,step)';    %defining parameter b
for i = 1:step
    for j = 1:step
        MESH(i,j) = 1 -2*[h1(j) h2(i)]*p + [h1(j) h2(i)]*R*[h1(j) h2(i)]'; %Computing cost function
    end
end

%Plot the Contours of MSE
figure                  
contour(h1,h2,MESH)     
hold on
plot(H_vector(1,:),H_vector(2,:))
title('Steepest Descent Error Contour')
xlabel('b')
ylabel('a')
legend('MSE','H^')
%Plot the Surface of MSE
figure                  
surf(MESH)              
title('Steepest Descent Error Surface')
xlabel('b')
ylabel('a')
zlabel('MSE')
%Plot the real and predicted output
figure
plot(d,'r'); hold on ;
plot(d_hat,'b');
title('System output') ;
legend('Real Output','Predicted Output')
%Plot parameter b
figure
subplot(2,1,1)
plot(H_vector(1,1:N));
title('Parameter Estimates (b = coefficient of u(k) ') ;
xlabel('Time Steps')
ylabel('b-coefficient of x(k) ')
legend('b='+ string(H(1)))
%Plot parameter a
subplot(2,1,2)
plot(H_vector(2,1:N));
title('Parameter Estimates (a=coefficient of d(k)') ;
xlabel('Time Steps')
ylabel('a-coefficient of d(k)  ')
legend('a='+ string(H(2)))
toc