function AcceptRejectMethod()
close all;
%% Generate eta ~ N(0,1) using exp(-xi)
% see http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-ARM.pdf
N = 1e6; % the number of samples

% get Z
a = 3;
x = linspace(-a,a,N);
Z = (2*a)*sum(exp(-(x.^4 - 2*x.^2 + 1))) / N;
fprintf('Z = %d\n',Z);

sigma = fminbnd(@g_sigma,-5,5);
fprintf('Sigma = %d\n',sigma);
v = randn(N, 1);
xi = v.*sigma;
% xi = -log(v);

u = rand(N,1);
log_u = log(u);

% generate signs for eta
b = randn(N,1);
ind = b < 0;
s = ones(N,1);
s(ind) = -1;

% log_ratio = -0.5*(xi - 1).^2;
log_ratio = -xi.^4 + 2*(1+1/(4*sigma.^2))*xi.^2 + (1+1/(4*sigma.^2).^2 - 2);

ind = find(log_u <= log_ratio);
Na = length(ind); % the number of accepted RVs
eta = xi(ind).*s(ind);

fprintf('N/Na = %d, C = %d\n',N/Na,g_sigma(sigma)/Z);
% fprintf('N/Na = %d, C = sqrt(2*e/pi) = %d\n',N/Na,sqrt(2*exp(1)/pi));

%% plot a histogram to test the distribution
nbins = 500; % the number of bins
etamax = max(eta);
etamin = min(eta);
nb1 = nbins + 1;
x = linspace(etamin,etamax,nb1);
h = x(2) - x(1); % bin width

% xc = centers of bins
xc = linspace(etamin + 0.5*h,etamax - 0.5*h,nbins);
hh = zeros(nbins,1); % heights of the bins

for i = 1 : nbins
    ind = find(eta >= x(i) & eta < x(i + 1));
    hh(i) = length(ind);
end

hh = hh/(Na*h); % scale the histogram
fprintf('MC: E[|x|] = %d\n',2*sum(hh)/length(hh));
% f = exp(-0.5*x.^2)/sqrt(2*pi);
f = exp(-(x.^4 - 2*x.^2 + 1));
fprintf('True: E[|x|] = %d\n', sum(f.*x)/length(x));
figure;

plot(x,f./Z,'r','Linewidth',2);
hold on;
plot(xc,hh,'b','Linewidth',2);
grid;
set(gca,'Fontsize',20);
xlabel('x','Fontsize',20);
ylabel('f(x)','Fontsize',20);
legend('True f','Generated f');

end