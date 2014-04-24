%Author: Rob Baummer
%
%Form a circulant matrix from the vector r with r as the first column
function R = circ(r)
%make r a column matrix
r = r(:);

%length of r
l = length(r);

%Preallocate R for speed
R = zeros(l,l);
%Set first column of circulant matrix R
R(:,1) = r;
%Set remaining columns of circulant matrix R
for i = 1:l-1
    R(:,i+1) = [r(l-(i-1):end); r(1:l-i)];
end
