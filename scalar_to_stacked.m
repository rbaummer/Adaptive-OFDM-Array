%Takes a NxN diagonal matrix and makes a L*NxL*N stacked matrix
function D = scalar_to_stacked(M,L)
    %extract the number of channels (FFT bins)
    %extract the number of antenna elements
    [~, N] = size(M);
    %pre-allocate the stacked diagonal matrix
    D = sparse(L*N,L*N);
    %for each square matrix in M
    for i = 1:N
        %column locations for current square matrix
        col = repmat((i:N:L*N),L,1);
        %row locations for current square matrix
        row = repmat((i:N:L*N),1,L);
        %Combine new stacked matrix with previous
        %since new stacked matrix entries are zero in previous matrix
        %addition can be used
        m = M(i,i)*ones(L^2,1);
        D = D + sparse(row, col(:), m(:), L*N, L*N);
    end
end