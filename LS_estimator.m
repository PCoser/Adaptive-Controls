clc
clear all
close all
tic             %measure the processing time
%Generate the original signal
n = 400;                             %number of points
u = 2*randi(2,1,n) - 3;              %generate a random binary signal u(k) = +-1 with n elements
Theta = transpose([1.5 -0.7 1]);     %input of parameters
[Y,PHI] = constructor2(u,Theta,0);   %Build the signal
%Corrupt signal with noise
for i = 1:n
    %Ym(i) = Y(i) + sqrt(4)*randn();     %Uncomment for variance = 4
    Ym(i) = Y(i) + sqrt(9)*randn();    %Uncomment for variance = 9
    Ym = Ym';
end
%Build the corrupted PHI to use in LS estimates
for i = 3:n
    phim = [Ym(i-1) Ym(i-2) u(i-1)];
    for j = 1:length(phim)
        PHIm(j,i) = phim(j);
    end
end

%First estimative of parameters using Least Squares
%Estimate the parameters
Theta_hat = (inv(PHIm*transpose(PHIm))*PHIm)*Ym;         
%Predict Output Y^
Y_hat0 = zeros([n 1]);
for k = 3:n
    Y_hat0(k) = [Y_hat0(k-1) Y_hat0(k-2) u(k-1)]*Theta_hat;
end
%Computing cost function
E = 0;
for k = 1:n
    E = E + (Y(k) - Y_hat0(k))^2;
end
J(1) = (1/(n*2))*E;     %J: vector that hosts the values for cost function in all iterations.

%Using now the LS Instrumental Approach
iter = 10;                   %number of iterations
Thetainst = zeros([3 iter]); %Initialize the matrix that hosts the values of Theta in each iteration
Thetainst(:,1) = Theta_hat;  %On the first iteration Thetainst = Theta from LS estimate
Y_hat = zeros([n iter]);     %Initialize matrix of predict outputs in each iteration
Y_hat(:,1) = Y_hat0;         %First output is the prediction from LS algorithm
%Iteration
for i = 1:iter - 1          % -1 because the predictor predicts the next Theta and Y^
    Thetainst(:,i+1) = inst_estimator(Thetainst(:,i),Ym,PHIm,u); %Instrument approach alg.
    %Predict Output Y^ using the parameter estimated above
    for k = 3:n
        Y_hat(k,i+1) = [Y_hat(k-1,i+1) Y_hat(k-2,i+1) u(k-1)]*Thetainst(:,i+1);
    end
    %Computing cost function of Y^ 
    E = 0;
    for k = 1:n
        E = E + (Y(k) - Y_hat(k,i))^2; %Sum up squared errors
    end
    J(i+1) = (1/(n*2))*E;  %J: vector that hosts the values for cost function in all iterations.
end

%Plot the original signal, corrupted and predicted
figure()
plot([1:n],Y,[1:n],Ym,[1:n],Y_hat(:,1),[1:n],Y_hat(:,2),[1:n],Y_hat(:,3),[1:n],Y_hat(:,4),[1:n],Y_hat(:,5))
title('Output for several iterations') ;
xlabel('Time Steps')
ylabel('Y')
legend('Y','Ym','Yinst 1','Yinst 2','Yinst 3','Yinst 4','Yinst 5')
%Plot each iteration in different graphics
for i = 1:5
    figure()
    plot([1:n],Y,[1:n],Ym,[1:n],Y_hat(:,i));
    title('Output for several iterations') ;
    xlabel('Time Steps')
    ylabel('Y')
    legend('Y','Ym','Yinst '+string(i))
end
%Plot Cost Function for each iteration
figure()
plot(J);
title('Cost Function for each iteration') ;
xlabel('Iteration')
ylabel('Cost Function')
legend('J')
toc         %measure the processing time
