using McCormick, DataFrames, IntervalArithmetic
using Ipopt, JuMP
using NLsolve
MCmode = NS

# constants for carbon dioxide
a = 3.610 # atm dm^6 mol^-2
b = 0.0429 # dm^3 mol^-1
R = 8.20574e-2 # L atm K^−1 mol^−1
Tc = 8*a/(27*R*b)
Vc = 3*b
Pc = a/(27*b^2)
T = 0.98*Tc

# pressure as parameter
pL = 59
pU = 107
pI = Interval(pL,pU)
p_range = range(pL,stop=pU,length=101)

xL = 0.6*Vc
xU = 2.0*Vc
V_range = range(xL,stop=xU,length=101)

f(V,P) = (P+a/V^2)*(V-b) - R*T
f_P(V) = R*T/(V-b) - a/V^2

function solve_original(xL,xU,p)
    # solve for the original system
    function f_vec!(F, x)
        F[1] = (p+a/x[1]^2)*(x[1]-b) - R*T
    end
    initial = [0.5]
    sol = nlsolve(f_vec!, initial, autodiff = :forward)

    return sol.zero[1]
end

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

function compute_original(xL,xU)
    # create a data frame to store output data
    df = DataFrame(p = Float64[], x=Float64[])

    for V in V_range
        x = f_P(V)
        save_tuple = (V, x)
        push!(df, save_tuple)
    end
    return df
end

function compute_cvcc(xL,xU)
    # create a data frame to store output data
    df = DataFrame(p = Float64[], cv = Float64[], cc = Float64[])

    for p in p_range
        xcv, xcc = optimize(xL,xU,p)
        save_tuple = (p, xcv, xcc)
        push!(df, save_tuple)
    end
    return df
end

# implicit function of V(P)
df = compute_cvcc(xL,xU)
df_ori = compute_original(xL,xU)

using Plots
pyplot()
data_mc = [df_ori.x, df.cv, df.cc]
plt = plot(legend = (0.52,0.3), legendfontsize = 16, orientation = :h)
plt = plot!(guidefontsize = 16, tickfontsize = 16)
plt = plot!(data_mc[1], V_range, label="V", ls=:solid, lw = 2)
plt = plot!(p_range, data_mc[2], label="Relaxations", lc=:red, ls=:dash, lw = 2)
plt = plot!(p_range, data_mc[3], label="", ls=:dash, lc=:red, lw = 2)
plt = plot!(xlabel="P (atm)", ylabel="V (L)")
display(plt)
