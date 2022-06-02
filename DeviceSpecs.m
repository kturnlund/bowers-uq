
%%%%%%%%%%%%%%%%%%%
%
% CURRENT ASSUMPTIONS:
%
% - Limited to two parameter dimensions for proof of concept  8-27-21
% 
%
%%%%%%%%%%%%%%%%%%%


classdef DeviceSpecs < handle
    % DeviceSpecs  Object to hold device properties

    
    properties
        dim                 % positive integer indicating dimensionality of system. data type: positive int
        mu                  % vector containing mu values for each parameter.   data type: scalar double
        sigma_2             % vector containing sigma^2 values for each parameter.   data type: scalar double
        gauss_mix_weights   % vector containing Gaussian mixture weight values for each parameter.   data type: scalar double
        p                   % vector containing probability distributions.   data type: symbolic function
        phi
        degree              % Max degree for alpha
        m                   % # of quadrature points.
        C
        xi_ij
        gPCE_weights
        H_weights
        H_xj
    end 
    

    methods (Access = public)
        
        function HG = DeviceSpecs(mu_in, sigma_in, weights_in, dim, degree)
            % HG = DeviceSpecs(mu, sigma, weights) creates an object with 
            % a Gaussian mixture distribution p given mu vector, sigma matrix,
            % and gaussian weights vector
            HG.mu = mu_in;
            HG.sigma_2 = sigma_in;
            HG.gauss_mix_weights = weights_in;
            HG.dim = dim;
            HG.degree = degree;
            HG = MarginPDF(HG);
        end
        
        function HG = xiGenerate(HG, p)
            % Generates the xi values for each parameter that are then used
            % to evaluate the simulation 
            HG.m = p;
            [HG.H_weights, HG.H_xj] = HermiteGaussNodes(HG, HG.m + 1);
            HG.PhiBasisConstruction();
            %test_check_phi = 
             
        end
        
        function HG = cGenerate(HG, array_in) %array_in is results from simulations
            % Calculates coefficients as a tensor product
            
            HG.C = tensor(zeros(HG.degree + 1, HG.degree + 1));

            for i=1:HG.degree
                for j=1:HG.degree
                    
                    p_eval = 1; %placeholder for p(xi_vector)/(all p's multiplied)
                    
                    psi_eval = polyval(HG.phi(i, :, 1).data, HG.xi_ij(1, i)) * polyval(HG.phi(j, :, 2).data, HG.xi_ij(2, j));
                    
                    
                    psi_norm_eval = norm(polyval((HG.phi(i, :, 1).data), HG.xi_ij(1, i)))^2 * ...
                        norm(polyval(flip(HG.phi(j, :, 2).data), HG.xi_ij(2, j)))^2;
                    
                    
                    u_Psi_innerProd = array_in(i, j) * HG.gPCE_weights(1, i) * HG.gPCE_weights(2, j) * psi_eval * p_eval;
                    
                    HG.C(i, j) = u_Psi_innerProd/psi_norm_eval;
                end
            end
            
            %Now we have C and psi's, need to store for u-vals
            
        end
        
        function return_value = u_eval(HG, coordinates)
            %Evaluates u at a given set of points
            
            return_value = zeros(length(coordinates));
            for m = 1 : length(coordinates)
                for k = 1 : length(coordinates)
                    for i = 1 : HG.m+1
                        for j = 1: HG.m+1
                            %evaluate each parameter
                            psi_eval = polyval(HG.phi(i, :, 1).data, coordinates(1, m)) * ...
                                polyval(HG.phi(j, :, 2).data, coordinates(2, k));

                            return_value(m, k) = return_value(m, k) + psi_eval * HG.C(i,j);

                        end
                    end
                end
            end
        end
         
    end
    
    
     methods (Access = private)
        
        function HG = MarginPDF(HG) 
        % HG = MarginPDF(HG) is a private nested function that initializes the 
        % Gaussian probability distributions contained in the vector p
        % corresponding to each parameter.
        
            syms xi;
            
            HG.p = {};
            
            for i = 1 : HG.dim
                HG.p{i, 1} = @(xi) (HG.gauss_mix_weights(1) * exp((xi - HG.mu(i, 1))^2 ...
                     /(-2 * HG.sigma_2(i, 1)))) / (sqrt(2 * pi * HG.sigma_2(i, 1))) ...
                     + (HG.gauss_mix_weights(2) * exp((xi - HG.mu(i, 2))^2 ...
                     /(-2 * HG.sigma_2(i, 2)))) / (sqrt(2 * pi * HG.sigma_2(i, 2)));
            end
            
        end 
       
        function HG = PhiBasisConstruction(HG) % dim is number of phi's to generate, w/ largest phi degree being degree
            %Constructs the phi basis functions for each parameter
            HG.phi = tenzeros(HG.degree + 1, HG.degree + 1, HG.dim);
           
            for i = 1 : HG.dim
                HG.phi(1, HG.degree + 1, i) = 1.0;
                
                temp_poly_list = tenzeros(HG.degree + 1);
                
                temp_poly_list(HG.degree:(HG.degree+1)) = [1, -1 * innerProdXi(HG, 1, i)/innerProd(HG, 1, i)];
                
                HG.phi(2, :, i) = temp_poly_list;

                for j = 3 : (HG.degree + 1)
                    temp_poly_list = tenzeros(HG.degree + 1);
                    
                    a = aGen(HG, j - 1, i);
                    b = bGen(HG, j - 1, i);
                    
                    temp_poly_list(1:(HG.degree)) = HG.phi(j - 1, 2:(HG.degree + 1),  i); %xi * phi_previous
                    
                    temp_poly_list = temp_poly_list + ...
                        HG.phi(j-1, :, i) * (-1) * a + ...
                        HG.phi(j-2, :, i) * (-1) * b;
                    
                    HG.phi(j, :, i) = temp_poly_list;
                end
            end
        end
        
        
        function J = JConstruct(HG, p_index)
            %Creates the Golub-Welsch matrix for determining xi values
            J = zeros(HG.m, HG.m);
            for j = 1 : HG.m
                for k = 1 : HG.m
                    if k == j
                        J(j, k) = aGen(HG, 1, p_index);
                    elseif k == j - 1
                        J(j, k) = sqrt(bGen(HG, j, p_index));
                    elseif k == j + 1
                        J(j, k) = sqrt(bGen(HG, k, p_index));
                    end
                end
            end
        end 

        function [w, xi_out] = cWeights(HG, q)
            % Calculates the xi values and corresponding weights using the
            % Golub-Welsch algorithm
            
            syms xi;
            %Different mu than the Gaussian mixture one, this is mu_ij
            mu = zeros(HG.dim, 1);
            w = zeros(HG.dim, q);
            xi_out = zeros(HG.dim, q);
            
            for i=1:HG.dim

                mu(i) = int(HG.p{i}(xi), -Inf, Inf);
                
                J = JConstruct(HG, i);
                [J_eigenvectors, J_eigenvalues] = eig(J);
                J_eig_norm = J_eigenvectors ./ norm(J_eigenvectors);
                
                for j = 1 : q
                    w(i, j) = mu(i) * J_eig_norm(1, j)^2;
                    xi_out(i, j) = J_eigenvalues(j,j);
                end
            end 

        end 

        function result = innerProdXi(HG, phi_index, i) %i is p index
            % Calculates the inner product of phi and phi * xi using a Hermite-Gauss
            % quadrature
            syms xi;
                sum_terms = zeros(length(HG.gauss_mix_weights), 1);
                
                for k = 1 : length(HG.gauss_mix_weights)
                    for j = 1 : length(HG.H_weights)
                            sum_terms(k) = sum_terms(k) + ...
                                HG.H_weights(j) * (HG.mu(i,k) + HG.sigma_2(i, k) * HG.H_xj(j)) ...
                                * polyval(HG.phi( phi_index, :, i).data, (HG.mu(i,k) + HG.sigma_2(i, k) * HG.H_xj(j)))^2;                       
                    end
                    sum_terms(k) = sum_terms(k) * HG.gauss_mix_weights(k);
                end
                
                result = sum(sum_terms);
        end
        
        function result = innerProd(HG, phi_index, i) %i is p index
            % Calculates the inner product of phi using a Hermite-Gauss
            % quadrature
            syms xi;
            sum_terms = zeros(length(HG.gauss_mix_weights), 1);
                for k = 1 : length(HG.gauss_mix_weights)
                    for j = 1 : length(HG.H_xj)
                            sum_terms(k) = sum_terms(k) + ...
                                HG.H_weights(j) * polyval(HG.phi(phi_index, :, i).data, (HG.mu(i,k) + HG.sigma_2(i, k) * HG.H_xj(j)))^2;
                    end
                    sum_terms(k) = sum_terms(k) * HG.gauss_mix_weights(k);
                end
                
            result = sum(sum_terms);
            
        end

        function a = aGen(HG, phi_index, i)
            % Calculates "a" coefficient for three-term recurrence
                   a = innerProdXi(HG, phi_index, i)/innerProd(HG, phi_index, i);
        end 
            
        function b = bGen(HG, phi_index, i)
            % Calculates "b" coefficient for three-term recurrence
                b = innerProd(HG, phi_index, i)/innerProd(HG, phi_index - 1, i);
        end 
        
        function [weights, x_j] = HermiteGaussNodes(HG, q)
            % Calculates Hermite-Gauss nodes and abcissae used in quadrature
            % calculations
            syms x;
            weights = zeros(q, 1);
            
            H = HermiteGauss(HG, q);
            
            x_j = double(roots(sym2poly(H)));
            
            H = simplify(HermiteGauss(HG, q-1));
            
            for i = 1 : q
                weights(i) = power(2, q - 1) * factorial(q) * sqrt(pi) ...
                    / (q^2 * subs(H,x,x_j(i))^2);
            end
        end
        
        function H = HermiteGauss(HG, n)
            % Calculates nth Hermite-Gauss function
            syms x;
            H = (-1)^n * exp(x^2) * diff(exp(-1 * x^2), n); 
        end
     end
end