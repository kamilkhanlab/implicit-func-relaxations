
"""
    relaxation
"""
function compute_ODE_relaxation(xL,xU,p)
    # # original
    # model = Model(Ipopt.Optimizer)
    # set_optimizer_attribute(model, "print_level", 0)
    # @variable(model, lb[i] <= x[i=1:N+1] <= ub[i])
    # @NLobjective(model, Min, 0)
    # @constraint(model, x[1] == IC)
    # JuMP.register(model, :f, 2, f, autodiff=true)
    # @NLconstraint(model, con[i=1:N], x[i+1] == x[i] + dT*f(x[i+1],p))
    # optimize!(model)

    # convex
    cv_model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(cv_model, "print_level", 0)
    @variable(cv_model, xL[i] <= xcv[i=1:N+1] <= xU[i])
    @NLobjective(cv_model, Min, xcv[N+1])
    @constraint(cv_model, xcv[1] == IC)
    JuMP.register(cv_model, :cv, 8, MCcv, MCcv_grad)
    JuMP.register(cv_model, :cc, 8, MCcc, MCcc_grad)
    @NLconstraint(cv_model, con_cv[i=1:N],
        cv(xcv[i],xcv[i+1],xL[i],xL[i+1],xU[i],xU[i+1],p,i) <= 0)
    @NLconstraint(cv_model, con_cc[i=1:N],
        cc(xcv[i],xcv[i+1],xL[i],xL[i+1],xU[i],xU[i+1],p,i) >= 0)
    optimize!(cv_model)

    cc_model = Model(Ipopt.Optimizer)
    @variable(cc_model, xL[i] <= xcc[i=1:N+1] <= xU[i])
    set_optimizer_attribute(cc_model, "print_level", 0)
    @NLobjective(cc_model, Max, xcc[N+1])
    @constraint(cc_model, xcc[1] == IC)
    JuMP.register(cc_model, :cv, 8, MCcv, MCcv_grad)
    JuMP.register(cc_model, :cc, 8, MCcc, MCcc_grad)
    @NLconstraint(cc_model, con_cv[i=1:N],
        cv(xcc[i],xcc[i+1],xL[i],xL[i+1],xU[i],xU[i+1],p,i) <= 0)
    @NLconstraint(cc_model, con_cc[i=1:N],
        cc(xcc[i],xcc[i+1],xL[i],xL[i+1],xU[i],xU[i+1],p,i) >= 0)
    optimize!(cc_model)

    return  JuMP.value.(xcv), JuMP.value.(xcc)
end

function MCcv(y...)
    x1, x2, xL1, xL2, xU1, xU2, p, i = y
    i = Int(i)
    # construct MC objects
    x1MC = MC{3,MCmode}(x1,Interval(xL1,xU1),1)
    x2MC = MC{3,MCmode}(x2,Interval(xL2,xU2),2)
    pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
    fMC = fimp(x1MC, x2MC, pMC)
    return fMC.cv
end

function MCcc(y...)
    x1, x2, xL1, xL2, xU1, xU2, p, i = y
    i = Int(i)
    # construct MC objects
    x1MC = MC{3,MCmode}(x1,Interval(xL1,xU1),1)
    x2MC = MC{3,MCmode}(x2,Interval(xL2,xU2),2)
    pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
    fMC = fimp(x1MC, x2MC, pMC)
    return fMC.cc
end

function MCcv_grad(grad, y...)
    x1, x2, xL1, xL2, xU1, xU2, p, i = y
    i = Int(i)
    # construct MC objects
    x1MC = MC{3,MCmode}(x1,Interval(xL1,xU1),1)
    x2MC = MC{3,MCmode}(x2,Interval(xL2,xU2),2)
    pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
    fMC = fimp(x1MC, x2MC, pMC)
    grad[1:3] = fMC.cv_grad[1:3]
end

function MCcc_grad(grad, y...)
    x1, x2, xL1, xL2, xU1, xU2, p, i = y
    i = Int(i)
    # construct MC objects
    x1MC = MC{3,MCmode}(x1,Interval(xL1,xU1),1)
    x2MC = MC{3,MCmode}(x2,Interval(xL2,xU2),2)
    pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
    fMC = fimp(x1MC, x2MC, pMC)
    grad[1:3] = fMC.cc_grad[1:3]
end
