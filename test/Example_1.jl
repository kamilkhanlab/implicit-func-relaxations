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

f(x,p) = x^2 + p*x + 4

pL = 6
pU = 9
pI = Interval(pL,pU)
p_range = range(pL,stop=pU,length=20)

function compute_cvcc(xL,xU)
    df = DataFrame(p = Float64[], x = Float64[], cv = Float64[], cc = Float64[], lb = Float64[], ub = Float64[])

    X = Interval(xL,xU)
    for p in p_range
        x = (-p-sqrt(p^2-16))/2
        xcv, xcc = optimize(xL,xU,p)
        save_tuple = (p, x, xcv, xcc, xL, xU)
        push!(df, save_tuple)
    end
    return df
end

# xL = -0.78
# xU = -0.4
xL = -10
xU = -5
df = compute_cvcc(xL,xU)

## Plot relaxations
using Plots
pyplot()
using LaTeXStrings
data_mc = [df.x, df.cv, df.cc, ]
plt = plot(legend = (0.0,0.02), legendfontsize = 16, orientation = :h)
plt = plot!(guidefontsize = 16, tickfontsize = 16)
plt = plot!(p_range, df.x, label=L"x^‡(p)", ls=:solid, lw = 2)
plt = plot!(p_range, df.cv, label="Relaxations", ls=:dot, lw = 3, lc=:green)
plt = plot!(p_range, df.cc, label="", ls=:dot, lw = 3, lc=:green)
plt = plot!(p_range, df.lb, label="Bounds", ls=:dash, lw = 2, lc=:grey)
plt = plot!(p_range, df.ub, label="", ls=:dash, lw = 2, lc=:grey)
plt = plot!(xlabel="p", ylabel=L"x^‡")
plt = plot!(ylims=(-10.2,-4.8))
display(plt)
