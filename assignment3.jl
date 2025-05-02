using Plots
using Printf
using LinearAlgebra

function question1()
    times_to_save = [5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100, 200, 300, 400, 500]

    Svals = Float64[]    
    tvals = Float64[]    
    Beta_Tvals = Float64[] 
    Ivals = Float64[]
    IncidenceVals = Float64[]

    # constants
    B = 100 # 100 people / day
    mu_S = 0.01 # people / day
    smax = 10
    N = 100
    sval = range(0.0, smax, length=N)
    dt = 1
    ds = 1
    t = 0
    tmax = 200

    # functions and IC
    gamma(s) = max(0.1 * (s-3), 0)
    beta(s) = 0.0001 * ((3.5 <= s && s <= 6.5) ? 1 : 0)
    nI(s) = max(10-s,0)  # initial condition for n(t=0,s)
    S = 10000 # initial condition for S(t=0)
    n = nI.(sval)

    # save initial plot
    p = plot(
        sval, n,
        xlabel = "Infection-age s (days)",
        ylabel = "Density n(t, s)",
        ylim = (0, 400),
        title = @sprintf("Infection-Age Density at t = %.0f days", t),
        lw = 2,
        legend = false
    )
    filename = "plot_output_t$(t)_assignment3_q1.png"
    savefig(p, filename)

    while t < tmax-dt
        t = t + dt

      # could be done in a vector
      for j=1:N
        n[j] = n[j] / (1+dt*(gamma(sval[j])))  # solve first part of n DE.  del_t(n) + gamma(s)*n = 0
      end

      Beta_T = sum(n .* beta.(sval))*ds  # This it total infected people at time at current time
      S = S + dt * (B - mu_S * S - Beta_T * S)
      n0 = Beta_T * S

      # this is the upwind aging out.  Solves the second part of n DE.  del_t(n) + del_s(n) = 0
      nold = n
      n = [n0; nold[1:end-1]]  
      Ival = sum(n) * ds  # total number of infected at time t
      incidence = Beta_T * S

      p = plot(
        sval, n,
        xlabel = "Infection-age s (days)",
        ylabel = "Density n(t, s)",
        ylim = (0, 400),
        title = @sprintf("Infection-Age Density at t = %.0f days", t),
        lw = 2,
        legend = false
    )
      display(p)
      sleep(0.1)

      save_time = round(Int, t)
      if save_time in times_to_save
          filename = "plot_output_t$(save_time)_assignment3_q1.png"
          savefig(p, filename)
      end

      push!(Svals, S)
      push!(tvals, t)
      push!(Beta_Tvals, Beta_T)
      push!(Ivals, Ival)
      push!(IncidenceVals, incidence)

    end

    p = plot(
        tvals, Svals,
        xlabel = "Time (days)",
        ylabel = "Population",
        title = "S(t) and I(t) over Time",
        label = "Susceptibles S(t)",
        lw = 2
    )
    
    plot!(
        tvals, Ivals,
        label = "Infected I(t)",
        lw = 2,
        linestyle = :dash
    )
    
    savefig(p, "assignment3_q1_S_and_I.png")

    p_beta = plot(
    tvals, Beta_Tvals,
    xlabel = "Time (days)",
    ylabel = "Force of Infection (Beta_T)",
    title = "Beta_T(t): Force of Infection",
    ylim = (0, 1),
    lw = 2,
    label = "Beta_T(t)"
    )
    
    savefig(p_beta, "assignment3_q1_Beta_T.png")

    p_inc = plot(
    tvals, IncidenceVals,
    xlabel = "Time (days)",
    ylabel = "New Infections per Day",
    title = "Incidence: Beta_T(t) × S(t)",
    lw = 2,
    label = "Incidence"
    )

    savefig(p_inc, "assignment3_q1_Incidence.png")


end

function question2_integrate_and_fire_neuron()
    Rm = 0.04e9
    Cm = 0.5e-9
    V_rest = -70e-3
    V_threshold = -50e-3
    I_inj = 3e-10
    pulse_length = 100e-3  
    dt = 10e-3         
    tmax = 0.3   
    tau = Rm * Cm   

    # arrays
    Vvals = Float64[]
    Tvals = Float64[]

    # Initial conditions
    t = 0.0
    V = V_rest
    push!(Vvals, V * 1e3)  
    push!(Tvals, t * 1e3)  

    while t < tmax - dt
        t += dt

        I = I_inj
        if t > pulse_length
            I = 0.0
        end

        if V >= V_threshold
            V = V_rest
        else
            V = (V + (dt / tau) * (V_rest + Rm * I)) / (1 + dt / tau)
        end

        push!(Vvals, V * 1e3)
        push!(Tvals, t * 1e3)
    end

    p = plot(
        Tvals, Vvals,
        xlabel = "Time (ms)",
        ylabel = "Voltage (mV) of Membrane",
        title = "Leaky integrate and fire neuron simulation - 0p3 nA",
        lw = 2,
        legend = false
    )
      display(p)
      savefig(p, "assignment3_leaky_neuron_0p3 nA.png")

end

function question3_membrane_dpolarization_cable_model()
    a = -300e-6      # micrometres 
    b = 800e-6       # micrometres
    lambda = 50e-6  # micrometres
    tau_m = 0.02e-3  # millisecond
    M = 11e-6        
    n = 100 
    dx = (b-a) / n
    t = 0
    tmax = 10e-3           # 10 milliseconds total simulation time
    time_points = 1000   # 10,000 time steps
    dt = tmax / time_points  # = 1e-6 seconds (1 microsecond per step)
    D = (lambda^2) / tau_m

    x = range(a, b; length=n)
    V = zeros(n)

    e = ones(n)
    L = diagm(-1 => e[1:end-1], 0 => -2e, 1 => e[1:end-1])
    L[1, end] = 1
    L[end, 1] = 1
    L /= dx^2

    # Initial Conditions
    V[1 + (n ÷ 2)] = M / dx

    frame_num = 0  # Frame counter

    p = plot(x .* 1e6, V .* 1e3,  # µm, mV
    label = "Voltage across membrane in space",
    xlabel = "Position (µm)",
    ylabel = "Membrane Voltage (mV)",
    ylim = (0, 100),
    title = @sprintf("t = %.3f ms", t * 1e3),
    color = :red,
    legend = :topright)

    # display(p)
    filename = "plot_output_INITIAL_assignment3_q3b.png"
    savefig(p, filename)

    center_idx = 1 + n ÷ 2
    target_idx = center_idx + 10

    # simulate solution
    while t < tmax-dt
        t = t + dt
        V_o = V

        # solve PDE for diffusion component only.
        V = (I(n)-dt*D*L)\V_o
        V_o = V

        # solve PDE for reaction component
        # V = V_o / (1 + (dt/tau_m))

        p = plot(x .* 1e6, V .* 1e3,  # µm, mV
        label = "Voltage across membrane in space",
        xlabel = "Position (µm)",
        ylabel = "Membrane Voltage (mV)",
        ylim = (0, 100),
        title = @sprintf("t = %.3f ms", t * 1e3),
        color = :red,
        legend = :topright)

        display(p)
        @printf("t = %.3f ms, V(x[%d]) = %.3f mV\n", t * 1e3, target_idx, V[target_idx] * 1e3)

        # sleep(0.05)

    # Save first 20 frames only
    if frame_num < 20
        filename = @sprintf "plot_output_t%03d_assignment3_q3b.png" frame_num
        savefig(p, filename)
    end

    frame_num += 1
            
    end

    p = plot(x .* 1e6, V .* 1e3,  # µm, mV
    label = "Voltage across membrane in space",
    xlabel = "Position (µm)",
    ylabel = "Membrane Voltage (mV)",
    ylim = (0, 100),
    title = @sprintf("t = %.3f ms", t * 1e3),
    color = :red,
    legend = :topright)

    display(p)
    filename = @sprintf "plot_output_FINAL_t_assignment3_q3b.png"
    savefig(p, filename)

end





function main()
    # question1()
    # question2_integrate_and_fire_neuron()
    question3_membrane_dpolarization_cable_model()
end

main()