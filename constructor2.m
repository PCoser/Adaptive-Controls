function [Y,PHI] = constructor2(u,Theta,db)
    n = length(u);
    Y = zeros([n 1]);
    PHI = zeros([3 n]);
    for k = 3:n
        phi = [Y(k-1) Y(k-2) u(k-1)];
        for i = 1:length(phi)
            %build the PHI dataset
            PHI(i,k) = phi(i);
        end
        %Build the Y dataset
        Y(k) = phi*Theta + sqrt(db)*randn(); %(wgn(1,1,0)*db/4) %adding white gaussian noise(wgn), db is variance
    end
end