using Plots
function stretching_y_function(x)
  gamma1 = 2.5
  -tanh.(gamma1 .* (x)) ./ tanh.(gamma1)

end



L =0.5
N=100

x = LinRange(-L,L,N+1)
z = zeros(N+1)
S = 1/(2*stretching_y_function(D))
y = stretching_y_function(x)
plot(z,y,marker=(:ro))

