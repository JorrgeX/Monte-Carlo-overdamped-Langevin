sig = fminbnd(@g_sigma,-5,5);
fprintf('sigmin = %d\n', sig);
fprintf('cmin = %d\n', g_sigma(sig));

% Define the range of sigma values
sigma_values = -2:0.001:2;

% Initialize an array to store the function values for each sigma
result_values = zeros(size(sigma_values));

% Calculate the function values for each sigma
for i = 1:length(sigma_values)
    result_values(i) = g_sigma(sigma_values(i));
end

% Plot the function
figure;
plot(sigma_values, result_values, '-o', 'LineWidth', 2);
title('Plot of c(\sigma)');
xlabel('\sigma');
ylabel('c(\sigma)');
grid on;