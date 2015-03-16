function p  = Psi(dl, dm, sigma)
%sigma = standard deviation
%dl = overlap patch 1
%dm = overlap patch 2
%p = sum value
p = sum(sum(exp(-(abs(dl-dm))./(2.*(sigma^2))))); % ?? Euclidean Distance ??
end