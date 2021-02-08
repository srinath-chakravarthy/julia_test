# Find curvature of top surface
using CSV, DataFrames, Statistics
using Dierckx
using Plots
using LaTeXStrings
## ---- Create spline for each surface
xx = collect(0:0.001:5)
dw = 2.0
dh = 4.0
dInt = 0.1
yy = similar(xx)
# for x in xx
yprev = 0
for (i,x) in enumerate(xx)
    global yprev
    if x <= dw
        yy[i] = -dh/2.0*(1.0 + cos(x/dw))
        yprev = yy[i]
    else
        yy[i] = yprev
    end
end
# create the spline
spl = Spline1D(xx,yy)
# Now compute the spline normal
yp = derivative(spl, xx)
xnew = similar(xx)
ynew = similar(xx)
den = similar(xx)
den = sqrt.(yp.^2 .+ 1.0)
xnew .= xx .- dInt .* yp ./ (sqrt.(1 .+ yp .^ 2))
ynew .= yy .+ dInt ./ (sqrt.(1 .+ yp .^ 2))
spl1 = Spline1D(xnew,ynew)

# xnew .= xx .+ yp .* dInt ./den
# ynew .= yy .+ dInt ./ den
