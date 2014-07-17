%Convert square matrices to stacked diagonal matrix
%for two square matrices A = [a b; c d] and E = [e f; g h] 
%the result is [a 0 b 0; 0 e 0 f; c 0 d 0; 0 g 0 h]
%multiplying the result by a stacked array matrix 
%for arrays X = [x1 x2] and Y = [y1 y2] ([x1 0; 0 y1; x2 0; 0 y2])
%yields a stacked array matrix equivalent to A*X and E*Y
%[a*x1+b*x2 0; 0 e*y1+f*y2; c*x1+d*x2 0; 0 g*y1+h*y2]
function D = square_to_stacked(M)
    %extract the number of channels (FFT bins)
    %extract the number of antenna elements
    [L, ~, N] = size(M);
    %pre-allocate the stacked diagonal matrix
    D = sparse(L*N,L*N);
    %for each square matrix in M
    for i = 1:N
        %column locations for current square matrix
        col = repmat((i:N:L*N),1,L);
        %row locations for current square matrix
        row = repmat((i:N:L*N),L,1);
        %Combine new stacked matrix with previous
        %since new stacked matrix entries are zero in previous matrix
        %addition can be used
        m = M(:,:,i);
        D = D + sparse(col, row, m(:), L*N, L*N);
    end
end