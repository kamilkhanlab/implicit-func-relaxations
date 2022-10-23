using McCormick, DataFrames, IntervalArithmetic
using Ipopt, JuMP
using Plots, LaTeXStrings
pyplot()

include("compute_ODE_original.jl")
include("compute_ODE_relaxation.jl")
include("compute_ODE_bound.jl")

MCmode = NS
IC = 9.0
N = 20
dT = 1/N
f(x,p) = - x^2 + p
fimp(x1,x2,p) = x2 - x1  - dT*f(x2,p)

xL = 0.1*ones(N+1)
xU = 9.0*ones(N+1)

pL = -1
pU = 1
pI = Interval(pL,pU)
p_range = range(pL,stop=pU,length=10)
#
xL_1, xU_1 = compute_ODE_bound(xL,xU)

function compute_cvcc(xL,xU)
    # create a data frame to store output data
    df = DataFrame(p = Float64[], x = Float64[], cv = Float64[], cc = Float64[])

    for p in p_range
        x = compute_ODE_original(xL,xU,p)
        xcv, xcc = compute_ODE_relaxation(xL,xU,p)
        save_tuple = (p, x[end], xcv[end], xcc[end])
        push!(df, save_tuple)
    end
    return df
end

df = compute_cvcc(xL, xU)
df2 = compute_cvcc(xL_1, xU_1)

plt = plot(legend = (0.65, 0.15), legendfontsize = 16, orientation = :h)
plt = plot!(guidefontsize = 16, tickfontsize = 16)
plt = plot!(df.p, df.x, label=L"z^{20}", ls=:solid, lw=2, lc=:orange)
plt = plot!(df.p, df.cv, label=L"k=0", ls=:dash, lw = 4, lc=:green)
plt = plot!(df.p, df.cc, label="", ls=:dash, lw = 4, lc=:green)
plt = plot!(df.p, df2.cv, label=L"k=1", ls=:dot, lw = 4, lc=:red)
plt = plot!(df.p, df2.cc, label="", ls=:dot, lw = 4, lc=:red)
plt = plot!(xlabel=L"p", ylabel=L"z^{20}")
display(plt)
