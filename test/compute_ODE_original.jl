function compute_ODE_original(xL,xU,pv)
    # original
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    @variable(model, xL[i] <= x[i=1:N+1] <= xU[i])
    @NLobjective(model, Min, 0)
    @constraint(model, x[1] == IC)
    JuMP.register(model, :f, 2, f, autodiff=true)
    @NLconstraint(model, con[i=1:N], x[i+1] == x[i] + dT*f(x[i+1],pv))
    optimize!(model)

    return  JuMP.value.(x)
end
