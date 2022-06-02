

%%%%%%
% Define parameter count

parameter_count = 10;

monte_carlo_count = 7; %should be smaller than # of parameters intially

%Build a Monte Carlo w/  certain # of initial samples
%Generate a random array of correct size

rng("default");
rand_parameters = rand(monte_carlo_count, parameter_count);

rand_parameters = tensor(rand_parameters);
% **************************************
% multiply parameters by correct scaling
% **************************************

%Save to folder for Lumerical to use

%save("monte_carlo_output_lumerical.mat", rand_parameters);


%Take in Lumerical sampling results

% load("lumerical_results.mat");
% *************************************
%BUILDING PLACEHOLDER: RNG results array in same shape
% *************************************

lum_results = rand(7, 1);



%


rng default;
x = rand([1 10]);
y = rand([1 10]);
[V,C] = voronoin([rand_parameters.data(:, 1), lum_results])

xCenter = cellfun(@(index) mean(V(index,1)),C);
yCenter = cellfun(@(index) mean(V(index,2)),C);


