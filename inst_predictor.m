%Instrumental Variable Approach 
function Thetainst = inst_estimator(Theta_prev,Y,PHI,u) %Parameter estimator using Instruments
n = length(Y);   
X = zeros([n 1]);    %Initialize the instrument for 1 iteration
 
for k = 3:n
    zheta = [X(k-1) X(k-2) u(k-1)];  %Zheta    
    for i = 1:length(zheta)         
        %build the Zheta dataset
        ZHETA(i,k) = zheta(i);
    end
    X(k) = zheta*Theta_prev; %Generate the instruments
end

Thetainst = inv(ZHETA*transpose(PHI))*(ZHETA*Y); %Estimate the parameter