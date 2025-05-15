using Plots
using Printf
using LinearAlgebra

function question2()
    k1 = 0.1
    kp = 1
    kq = 0.2
    kappa = 0.05
    CQ = 1
    tmax = 40
    dt = 0.1

    N = Int(tmax / dt)
    times = collect(0:dt:tmax)

    S_vals = range(0.0051, stop=0.0499, length=5)

    plt = plot(
        title = "Phosphorylation Dynamics for Varying S",
        xlabel = "Time t",
        ylabel = "Concentration of Protein P",
        legend = :outertopright,
        size = (900, 600),
        ylim = (0, 1),
        yticks = 0:0.2:1
    )

    # functions for equilibrium analysis
    # comparing at different values of S the equilibrium point of CRP compared to the ODE solver.
    Beta_S(S) = 50S * (kappa - 1) + (1 + kappa)
    C_S(S) = (-Beta_S(S) + sqrt(Beta_S(S)^2 + 4 * (50S - 1) * (50S * kappa))) / (2 * (50S - 1))

    for S in S_vals
        V1 = S * kp / k1
        V2 = kq * CQ

        CRP = zeros(N + 1)
        CRP[1] = 0.0
        
        i = 1

        while i < N + 1
            CRP_t = CRP[i]
            term1 = V1 * (1 - CRP_t) / ((1 - CRP_t) + kappa)
            term2 = V2 * CRP_t / (CRP_t + kappa)
            CRP[i+1] = CRP_t + dt * (term1 - term2)
            i += 1
        end

        plot!(plt, times, CRP, lw=1, label="ODE solved S = $(round(S, digits=5))")

        # analytical solution
        CRP_eq = C_S(S)
        hline!(plt, [CRP_eq], linestyle=:dash, label="analytical sol @ S=$(round(S, digits=5))")
    end

    display(plt)
    savefig(plt, "assignment4_q2.png")
end

using Plots

function question3()
    # Parameters
    k1 = 0.1
    k1_p = 0.5
    kp = 1
    kq = 0.2
    kappa = 0.05
    CQ = 1
    tmax = 100
    dt = 0.1
    N = Int(tmax / dt)
    S_vals = collect(0.001:0.001:0.1)  # going forward in small dS steps up to 0.1
    S_vals_rev = reverse(S_vals)  # going backward in small dS steps back down to 0

    # storing protein P and protein R phosphorylated equilibrium values.
    CP_forward = Float64[]
    CRP_forward = Float64[]
    CP_backward = Float64[]
    CRP_backward = Float64[]

    # adding the unphosphorylated version here as its going to help tell the story about hysterisis
    CR_forward = Float64[]
    CR_backward = Float64[]

    CP0 = S_vals[1] / k1  
    CRP0 = 0.0

    # solving ODE going forward in dS
    for S in S_vals
        V2 = kq * CQ
        CP = CP0
        CRP = CRP0

        for _ in 1:N
            CP_new = CP + dt * (S - (k1 + k1_p * (1 - CRP)) * CP)
            V1 = kp * CP_new
            CRP_new = CRP + dt * (
                V1 * (1 - CRP) / ((1 - CRP) + kappa) -
                V2 * CRP / (CRP + kappa)
            )
            CP, CRP = CP_new, CRP_new
        end

        push!(CP_forward, CP)
        push!(CRP_forward, CRP)
        push!(CR_forward, 1 - CRP)
        CP0, CRP0 = CP, CRP  
    end

    # going backwards
    CP0 = CP_forward[end]
    CRP0 = CRP_forward[end]

    for S in S_vals_rev
        V2 = kq * CQ
        CP = CP0
        CRP = CRP0

        for _ in 1:N
            CP_new = CP + dt * (S - (k1 + k1_p * (1 - CRP)) * CP)
            V1 = kp * CP_new
            CRP_new = CRP + dt * (
                V1 * (1 - CRP) / ((1 - CRP) + kappa) -
                V2 * CRP / (CRP + kappa)
            )
            CP, CRP = CP_new, CRP_new
        end

        push!(CP_backward, CP)
        push!(CRP_backward, CRP)
        push!(CR_backward, 1 - CRP)
        CP0, CRP0 = CP, CRP
    end

    # PLOTS
    layout = @layout [a{0.33h}; b{0.33h}; c{0.33h}]
    plt = plot(layout=layout, size=(900, 800),
        title="Equilibrium Concentrations vs S",
        legend=:outertopright)

    plot!(plt[2], S_vals, CP_forward, label="CP forward", lw=2)
    plot!(plt[2], S_vals, reverse(CP_backward), label="CP backward", lw=2, linestyle=:dash)
    ylabel!(plt[2], "C_P")
    xlabel!(plt[2], "Synthesis rate S")

    plot!(plt[1], S_vals, CRP_forward, label="CRP forward", lw=2)
    plot!(plt[1], S_vals, reverse(CRP_backward), label="CRP backward", lw=2, linestyle=:dash)
    ylabel!(plt[1], "C_RP")
    xlabel!(plt[1], "")

    plot!(plt[3], S_vals, CR_forward, label="C_R forward", lw=2)
    plot!(plt[3], S_vals, reverse(CR_backward), label="C_R backward", lw=2, linestyle=:dash)
    ylabel!(plt[3], "C_R (unphosphorylated R)")
    xlabel!(plt[3], "Synthesis rate S")

    display(plt)
    savefig(plt, "assignment4_q3_hysteresis.png")

end





function main()

    # question2()
    question3()
end

main()