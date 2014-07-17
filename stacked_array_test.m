%Test square matrix multiplication for block matrices
function R = stacked_array_test

%Nx1 = NxN * Nx1

%array 1
x1 = (1:3)';
%array 2
x2 = (10:12)';
%array 3
x3 = (20:22)';

%matrix 1
m1 = [1 2 3; 3 4 5; 1 2 5];
%matrix 2
m2 = [10 11 12; 12 13 14; 13 11 12];
%matrix 3
m3 = [5 6 7; 8 9 10; 5 6 17];

M{1} = m1;
M{2} = m2;
M{3} = m3;

%result 1
r1 = m1*x1 %#ok<*NOPRT>
%result 2
r2 = m2*x2
%result 3
r3 = m3*x3

%stacked diagonal matrices
D = square_to_stacked(M);

X = [upsample(x1,3) circshift(upsample(x2,3),1) circshift(upsample(x3,3),2)];
%array length
L = 3;
%# arrays
N = 3;

%another way of constructing X
X = sparse((1:L:L*N),ones(L,1),x1, L*N, 1);
X = [X sparse((2:L:L*N),ones(L,1),x2, L*N, 1)];
X = [X sparse((3:L:L*N),ones(L,1),x3, L*N, 1)];

%result
R = D*X


end

%Convert square matrices to stacked diagonal matrix
function D = square_to_stacked(M)
    D=[upsample(M{1}(1,:),3); circshift(upsample(M{2}(1,:),3),[0 1]); circshift(upsample(M{3}(1,:),3),[0 2]);upsample(M{1}(2,:),3); circshift(upsample(M{2}(2,:),3),[0 1]); circshift(upsample(M{3}(2,:),3),[0 2]); upsample(M{1}(3,:),3); circshift(upsample(M{2}(3,:),3),[0 1]); circshift(upsample(M{3}(3,:),3),[0 2])];
      
    N = length(M);
    L = length(M{1});
    %another way of constructing D
    D = sparse(L*N,L*N);
    for i = 1:length(M)
        col = repmat((i:L:L*N),1,L);
        row = repmat((i:L:L*N),L,1);
        D = D + sparse(col, row, M{i}(:), L*N, L*N);
    end
end