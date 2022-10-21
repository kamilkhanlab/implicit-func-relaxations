
"""
    Bounds
"""

function compute_ODE_bound(xL, xU)

    function optimize(xL,xU)
        xLTraj = deepcopy(xL)
        xUTraj = deepcopy(xU)
        xLTraj[1] = IC
        xUTraj[1] = IC

        for j in 1:N
            # convex
            cv_model = Model(Ipopt.Optimizer)
            set_optimizer_attribute(cv_model, "print_level", 0)
            @variable(cv_model, xL[i] <= xcv[i=1:N+1] <= xU[i])
            @variable(cv_model, pL <= p <= pU)
            @NLobjective(cv_model, Min, xcv[j+1])
            @constraint(cv_model, xcv[1] == IC)
            JuMP.register(cv_model, :cv, 4, MClb, MClb_grad)
            JuMP.register(cv_model, :cc, 4, MCub, MCub_grad)
            @NLconstraint(cv_model, con_cv[i=1:N], cv(xcv[i],xcv[i+1],p,i) <= 0)
            @NLconstraint(cv_model, con_cc[i=1:N], cc(xcv[i],xcv[i+1],p,i) >= 0)
            optimize!(cv_model)

            cc_model = Model(Ipopt.Optimizer)
            @variable(cc_model, xL[i] <= xcc[i=1:N+1] <= xU[i])
            @variable(cc_model, pL <= p <= pU)
            set_optimizer_attribute(cc_model, "print_level", 0)
            @NLobjective(cc_model, Max, xcc[j+1])
            @constraint(cc_model, xcc[1] == IC)
            JuMP.register(cc_model, :cv, 4, MClb, MClb_grad)
            JuMP.register(cc_model, :cc, 4, MCub, MCub_grad)
            @NLconstraint(cc_model, con_cv[i=1:N], cv(xcc[i],xcc[i+1],p,i) <= 0)
            @NLconstraint(cc_model, con_cc[i=1:N], cc(xcc[i],xcc[i+1],p,i) >= 0)
            optimize!(cc_model)

            xLTraj[j+1] = JuMP.objective_value(cv_model)
            xUTraj[j+1] = JuMP.objective_value(cc_model)
        end

        # return  JuMP.value.(xcv), JuMP.value.(xcc)
        return xLTraj, xUTraj
    end

    function MClb(y...)
        x1, x2, p, i = y
        i = Int(i)
        # construct MC objects
        x1MC = MC{3,MCmode}(x1,Interval(xL[i],xU[i]),1)
        x2MC = MC{3,MCmode}(x2,Interval(xL[i+1],xU[i+1]),2)
        pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
        fMC = fimp(x1MC, x2MC, pMC)
        return fMC.cv
    end

    function MCub(y...)
        x1, x2, p, i = y
        i = Int(i)
        # construct MC objects
        x1MC = MC{3,MCmode}(x1,Interval(xL[i],xU[i]),1)
        x2MC = MC{3,MCmode}(x2,Interval(xL[i+1],xU[i+1]),2)
        pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
        fMC = fimp(x1MC, x2MC, pMC)
        return fMC.cc
    end

    function MClb_grad(grad, y...)
        x1, x2, p, i = y
        i = Int(i)
        # construct MC objects
        x1MC = MC{3,MCmode}(x1,Interval(xL[i],xU[i]),1)
        x2MC = MC{3,MCmode}(x2,Interval(xL[i+1],xU[i+1]),2)
        pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
        fMC = fimp(x1MC, x2MC, pMC)
        grad[1:3] = fMC.cv_grad[1:3]
    end

    function MCub_grad(grad, y...)
        x1, x2, p, i = y
        i = Int(i)
        # construct MC objects
        x1MC = MC{3,MCmode}(x1,Interval(xL[i],xU[i]),1)
        x2MC = MC{3,MCmode}(x2,Interval(xL[i+1],xU[i+1]),2)
        pMC = MC{3,MCmode}(p,Interval(pL,pU),3)
        fMC = fimp(x1MC, x2MC, pMC)
        grad[1:3] = fMC.cc_grad[1:3]
    end

    return optimize(xL,xU)
end
