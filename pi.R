n = 10000000
r = 1
x = runif(n,-1,1)
y = runif(n,-1,1)
ED = sqrt(x^2+y^2)
count = length(ED[ED<=1]) #count/n = pi/4
pi_estimate = count/n*4
pi_estimate