%matrization tester

help tensor_toolbox

Y = tenrand([5 4 3]);
A = rand(4,5); B = rand(3,4); C = rand(2,3); U = {A,B,C};

X = ttm(Y,U) %<-- X = Y x_1 A x_2 B x_3 C
Xm1 = kron(B,A)*tenmat(Y,[1 2])*C' %<-- Kronecker product version
Xm2 = tenmat(X,[1 2]) %<-- Matriczed version
norm(Xm1 - Xm2)  % <-- should be zero
Xm1 = B * tenmat(Y,2,[3 1]) * kron(A,C)'; %<-- Kronecker product version
Xm2 = tenmat(X,2,[3 1]); %<-- Matricized version
norm(Xm1 - Xm2) % <-- should be zero
Xm1 = tenmat(Y,[],[1 2 3]) * kron(kron(C,B),A)'; %<-- Vectorized via Kronecker
Xm2 = tenmat(X,[],[1 2 3]); %<-- Vectorized via matricize
norm(Xm1 - Xm2)


for n = 1:ndims(Y)
  X2 = ttm(Y,U,n); %<-- X = Y x_n U{n}
  Xn = U{n} * tenmat(Y,n); %<-- Xn = U{n} * Yn
  norm(tenmat(X2,n) - Xn)  % <-- should be zero
end