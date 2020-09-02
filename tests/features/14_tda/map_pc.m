function [PC] = map_pc(A,t,d,numPts)

% x_i x_{i+tau}
PC = zeros(numPts-t,d);

for j=1:numPts-t
    for k = 1:d
        PC(j,k) = A(j+(k-1)*t);
    end
end


end

