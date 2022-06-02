

mu = [9, 8; 6, 7];
sigma_2 = [6, 5; 3, 4];
gauss_mix_weights = [0.6, 0.4];
dimension = 2;
degree = 3;
test_device = DeviceSpecs(mu, sigma_2, gauss_mix_weights, dimension, degree);

B = test_device.xiGenerate(9).phi;


sample_count = 9; %sample count
R = 3; %R count


%initialization
y = tensor(zeros(sample_count, 1));

%Put sample outputs in y vector, right now randomized
y = tensor(rand(sample_count, 1));

U = {};

eta = zeros(R, 1);

q = 0.5;

lambda_0 = 1;
lambda = lambda_0;
Lambda = kron(diag(ones(R, 1)), eye(degree + 1));

epsilon = 1e-10;

for k = 1:dimension
    Phi = tensor(zeros(sample_count, R * (degree + 1)));
    PhiT = tensor(zeros(R * (degree + 1), sample_count));
    for n = 1:sample_count
        O = tenzeros([degree + 1 R]);

        for i = 1 : R
            for j = 1:degree + 1
                O(j, i) = polyval(B(j, :, k).data, y(n)); 
            end
        end
        Phi(n, :) = O(:);
        PhiT(:, n) = O(:);
    end

    
    u_vector = inv(PhiT.data * Phi.data  + lambda * Lambda) * PhiT.data * y.data;
    U{k} = reshape(u_vector, degree + 1, R);
    
end

v = zeros(R, 1);

for r = 1:R
    for k = 1:dimension
        v(r) = v(r) + norm(U{k}(:, r))^2;
    end
    v(r) = sqrt(v(r));
end
v

for r = 1:R
    eta(r) = power(v(r), 2-q) * power(norm(v, q), q) + epsilon;
end
eta

while true
    %Lambda update
    
    Lambda = kron(diag(ones(R, 1) ./ eta), eye(degree + 1));
    
    
    %U(k) update
    
    for k = 1:dimension
        
        for n = 1:N
            
            %Evaluate B tensor at nth sample
            
            B = tensor();
            
            O = mttkrp(
        end
        
        u_vector = inv(PhiT.data * Phi.data  + lambda * Lambda) * PhiT.data * y.data;
        U{k} = reshape(u_vector, degree + 1, R);
    end
    
    %eta update
    
    
    %lambda update
    lambda = lambda_0 * max(eta);
    
    
end


