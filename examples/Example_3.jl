using McCormick, DataFrames, IntervalArithmetic
using Ipopt, JuMP
MCmode = NS

function optimize(xL,xU,p)
    cv_model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(cv_model, "print_level", 0)
    @variable(cv_model, xL <= x <= xU)
    @NLobjective(cv_model, Min, x)
    JuMP.register(cv_model, :cv, 2, MCcv, MCcv_grad)
    JuMP.register(cv_model, :cc, 2, MCcc, MCcc_grad)
    @NLconstraint(cv_model, con_cv, cv(x,p) <= 0)
    @NLconstraint(cv_model, con_cc, cc(x,p) >= 0)
    optimize!(cv_model)

    cc_model = Model(Ipopt.Optimizer)
    @variable(cc_model, xL <= x <= xU)
    set_optimizer_attribute(cc_model, "print_level", 0)
    @NLobjective(cc_model, Max, x)
    JuMP.register(cc_model, :cv, 2, MCcv, MCcv_grad)
    JuMP.register(cc_model, :cc, 2, MCcc, MCcc_grad)
    @NLconstraint(cc_model, con_cv, cv(x,p) <= 0)
    @NLconstraint(cc_model, con_cc, cc(x,p) >= 0)
    optimize!(cc_model)

    return objective_value(cv_model), objective_value(cc_model)
end

function MCcv(y...)
    x, p = y
    # construct MC objects
    xMC = MC{2,MCmode}(x,Interval(xL,xU),1)
    pMC = MC{2,MCmode}(p,Interval(pL,pU),2)
    fMC = f(xMC,pMC)
    return fMC.cv
end

function MCcc(y...)
    x, p = y
    # construct MC objects
    xMC = MC{2,MCmode}(x,Interval(xL,xU),1)
    pMC = MC{2,MCmode}(p,Interval(pL,pU),2)
    fMC = f(xMC,pMC)
    return fMC.cc
end

function MCcv_grad(grad, y...)
    x, p = y
    # construct MC objects
    xMC = MC{2,MCmode}(x,Interval(xL,xU),1)
    pMC = MC{2,MCmode}(p,Interval(pL,pU),2)
    fMC = f(xMC,pMC)
    grad[1:2] = fMC.cv_grad[1:2]
end

function MCcc_grad(grad, y...)
    x, p = y
    # construct MC objects
    xMC = MC{2,MCmode}(x,Interval(xL,xU),1)
    pMC = MC{2,MCmode}(p,Interval(pL,pU),2)
    fMC = f(xMC,pMC)
    grad[1:2] = fMC.cc_grad[1:2]
end

f(x,p) = x^2 - (sqrt(max(1e-16, p^2 - p)) - 2)^4

pL = -3
pU = 3
pI = Interval(pL,pU)
p_range = range(pL, stop=pU, length=151)

function compute_cvcc(xL,xU)
    # create a data frame to store output data
    df = DataFrame(p = Float64[], x1 = Float64[], x2 = Float64[], cv = Float64[], cc = Float64[])
    X = Interval(xL,xU)

    for p in p_range
        if p^2 - p >=0
            x1 = (sqrt(p^2-p)-2)^2
            x2 = -(sqrt(p^2-p)-2)^2
        else
            x1 = NaN
            x2 = NaN
        end
        xcv, xcc = optimize(xL,xU,p)
        save_tuple = (p, x1, x2, xcv, xcc)
        push!(df, save_tuple)
    end
    return df
end

xL = -10
xU = 10
df = compute_cvcc(xL,xU)

using Plots
pyplot()
plt = plot(legend = (0.65,0.05), legendfontsize = 12, orientation = :h)
plt = plot!(guidefontsize = 16, tickfontsize = 16)
plt = plot!(df.p, df.x1, label="x(p)", ls=:solid, lw = 2, lc=:green)
plt = plot!(df.p, df.x2, label="", ls=:solid, lw = 2, lc=:green)
plt = plot!(df.p, df.cv, label="Relaxations", ls=:dash, lw = 2, lc=:red)
plt = plot!(df.p, df.cc, label="", ls=:dash, lw = 2, lc=:red)
# plt = plot!(title="Convex relaxations of an implicit function")
plt = plot!(xlabel="p", ylabel="x")
plt = plot!(ylims=(-10,10))
display(plt)
