%construct a stacked diagonal matrix from arrays organized as columns
%for arrays X = [x1 x2] and Y = [y1 y2] 
%the result is [x1 0; 0 y1; x2 0; 0 y2]
function X = array_to_stacked(M)
    %extract the number of channels (FFT bins)
    %extract the number of antenna elements
    [N, L] = size(M);
    %pre-allocate X
    X = [];
    
    for i = 1:N
        %append the column for new array to previous array columns
        X = [X sparse((i:L:L*N), ones(N,1), M(:,i), L*N, 1)];
    end
end