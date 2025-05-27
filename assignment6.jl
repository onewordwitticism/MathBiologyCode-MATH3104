using DifferentialEquations
using Plots
using Statistics
using Printf

# Define parameters as constants
# const S = 0.5  # S = substrate MAX transcription rate per SECOND
const K = 40   # K = concentration to hit half synth rate
const n = 2    #  n = Hill coefficient
const avg_mrna_t = 2  # average length of mrna survival in minutes
const avg_protein_t = 10  # average length of protein survival in minutes
const trans_efficiency = 20 * 60 # 20 proteins per minute per mRNA

const k_mrna = log(2) / avg_mrna_t  # deg rate for mrna's
const k_p = log(2) / avg_protein_t  # deg rate for proteins
const alpha = trans_efficiency / avg_mrna_t  # 10 proteins per min*mrna

function question1()
    # Initial conditions: [A, B, C, TFA, TFB, TFC]
    u0 = [8, 5, 15, 15, 25, 6]
    tspan = (0.0, 1000)

    S_values = 0.1:0.1:3.0 .* 60  # transcripts/min - originally 0.5 was in transcripts/second !!
    TFA_ss = Float64[]; TFB_ss = Float64[]; TFC_ss = Float64[]
    A_ss = Float64[]; B_ss = Float64[]; C_ss = Float64[]

    # Loop over each S value and solve the system to steady state
    for S in S_values
        function odes!(du, u, p, t)
            A, B, C, TFA, TFB, TFC = u
            du[1] = S * (K^n / (K^n + TFC^n)) - k_mrna * A  # TFC represses A
            du[2] = S * (K^n / (K^n + TFA^n)) - k_mrna * B  # TFA represses B
            du[3] = S * (K^n / (K^n + TFB^n)) - k_mrna * C  # TFB represses C
            du[4] = (alpha * A) - k_p * TFA
            du[5] = (alpha * B) - k_p * TFB
            du[6] = (alpha) * C - k_p * TFC
        end

        prob = ODEProblem(odes!, u0, tspan)  ### closest thing in Julia to ode45() in matlab
        sol = solve(prob, Tsit5())  ### closest thing in Julia to ode45() in matlab

        # Plot protein vs time if S == 0.5 * 60 (i.e. original value)
        if S == 0.5 * 60
            t = sol.t
            TFA = sol[4, :]
            TFB = sol[5, :]
            TFC = sol[6, :]

            p_time = plot(t, [TFA, TFB, TFC],
                label = ["TFA" "TFB" "TFC"],
                xlabel = "Time (min)",
                ylabel = "Protein concentration",
                lw = 2,
                title = "Protein Concentrations vs Time (S = 0.5)")

            display(p_time)  
            savefig(p_time, "protein_vs_time_S_0.5_QuestionA.png")
        end

        # Steady-state values
        push!(A_ss, sol[1, end])
        push!(B_ss, sol[2, end])
        push!(C_ss, sol[3, end])
        push!(TFA_ss, sol[4, end])
        push!(TFB_ss, sol[5, end])
        push!(TFC_ss, sol[6, end])
    end

    p1 = plot(S_values, [TFA_ss, TFB_ss, TFC_ss],
        label = ["TFA" "TFB" "TFC"],
        xlabel = "S (max transcription rate)",
        ylabel = "Protein steady-state concentration",
        lw = 2,
        title = "Protein Concentrations at Steady State",
        )

    # Plot steady-state mRNA concentrations
    p2 = plot(S_values, [A_ss, B_ss, C_ss],
        label = ["A (mRNA)" "B (mRNA)" "C (mRNA)"],
        xlabel = "S (max transcription rate)",
        ylabel = "mRNA steady-state concentration",
        lw = 2,
        title = "mRNA Concentrations at Steady State",
        ylims = (0, maximum([A_ss; B_ss; C_ss]) * 1.1))

    vline!(p1, [0.5 * 60], linestyle = :dash, color = :black, label = "")
    vline!(p2, [0.5 * 60], linestyle = :dash, color = :black, label = "")

    plot(p1, p2, layout = (2,1), size = (800, 800))
    full_plot = plot(p1, p2, layout = (2,1), size = (800, 800))
    savefig(full_plot, "protein_mrna_vs_S_fullrange_QuestionA.png")

    mask = S_values .<= 20
    max_zoom = maximum(vcat(TFA_ss[mask], TFB_ss[mask], TFC_ss[mask]))

    p3 = plot(S_values, [TFA_ss, TFB_ss, TFC_ss],
        label = ["TFA" "TFB" "TFC"],
        xlabel = "S (max transcription rate)",
        ylabel = "Protein concentration",
        lw = 2,
        title = "Zoomed Protein Concentrations (S = 0 to 20)",
        xlims = (0, 20),
        ylims = (0, max_zoom * 1.1)
    )

    p4 = plot(S_values, [A_ss, B_ss, C_ss],
        label = ["A (mRNA)" "B (mRNA)" "C (mRNA)"],
        xlabel = "S (max transcription rate)",
        ylabel = "mRNA concentration",
        lw = 2,
        title = "Zoomed mRNA Concentrations (S = 0 to 20)",
        xlims = (0, 20))

    plot(p3, p4, layout = (2,1), size = (800, 800))
    zoomed_plot = plot(p3, p4, layout = (2,1), size = (800, 800))
    savefig(zoomed_plot, "protein_mrna_vs_S_zoomed_S_0_to_20_QuestionA.png")

end

function question2()
    # all parameters are NON DIMENSIONAL
    S0 = 1
    n = 4
    K = 100
    avg_mrna_t = 1
    avg_protein_t = 100
    trans_efficiency = 10
    tmax = 1000

    num_realizations = 10
    all_T = Vector{Vector{Float64}}()
    all_M = Vector{Vector{Float64}}()
    all_P = Vector{Vector{Float64}}()

    for run in 1:num_realizations
        P = Float64[]
        M = Float64[]
        T = Float64[]

        k_mrna = log(2) / avg_mrna_t  # deg rate for mrna's
        k_p = log(2) / avg_protein_t  # deg rate for proteins
        alpha = trans_efficiency / avg_mrna_t  # 10 proteins per min*mrna

        # set IC
        push!(P, 1)
        push!(M, 10)
        push!(T, 0)

        while T[end] < tmax
            v1 =  S0 * (K^n / (K^n + P[end]^n))  # repression by P
            v2 = k_mrna * M[end]
            v3 = alpha * M[end]
            v4 = k_p * P[end]

            total_flux = v1 + v2 + v3 + v4  # total probability/time that one reaction occurs out of our set of reactions
            tau = -log(rand()) / total_flux  # this is the exponential waiting time; we sample and get a time tau
            event = rand()
            scaled_event = event * total_flux  # scale the event to the total flux

            # note this if loop stuff only works because its checking in order.  Otherwise the logic would break.
            if scaled_event < v1  # mrna synthesis
                push!(M, M[end] + 1)
                push!(P, P[end])
            elseif scaled_event < v1 + v2  # mrna degradation
                push!(M, M[end] - 1)
                push!(P, P[end])
            elseif scaled_event < v1 + v2 + v3  # protein synthesis
                push!(M, M[end])
                push!(P, P[end] + 1)
            else  # protein degradation
                push!(M, M[end])
                push!(P, P[end] - 1)
            end

            push!(T, T[end] + tau)
        end

        push!(all_T, T)
        push!(all_M, M)
        push!(all_P, P)
    end

    p1 = plot(all_T[1], all_P[1],
        label = "Protein",
        xlabel = "Time",
        ylabel = "Molecule Count",
        title = "Negative Autoregulation: Protein and mRNA vs Time",
        lw = 1.5)

    plot!(all_T[1], all_M[1], label = "mRNA", lw = 1.5, ls = :dash)

    for i in 2:num_realizations
        plot!(all_T[i], all_P[i], label = "", lw = 1.5)
        plot!(all_T[i], all_M[i], label = "", lw = 1.5, ls = :dash)
    end

    savefig(p1, "negative_autoregulation_protein_mrna_vs_time_QuestionB.png")
end


function noise_given_lifetime(lifetime)
    avg_protein_t = lifetime
    k_p = log(2) / avg_protein_t
    S0 = 1
    n = 4
    K = 100
    avg_mrna_t = 1
    trans_efficiency = 10
    tmax = 100000

    k_mrna = log(2) / avg_mrna_t  
    k_p = log(2) / avg_protein_t  
    alpha = trans_efficiency / avg_mrna_t  # 10 proteins per min*mrna

    P = [1.0]
    M = [10.0]
    T = [0.0]

    while T[end] < tmax
        v1 = S0 * (K^n / (K^n + P[end]^n))
        v2 = k_mrna * M[end]
        v3 = alpha * M[end]
        v4 = k_p * P[end]
        total_flux = v1 + v2 + v3 + v4
        tau = -log(rand()) / total_flux
        r = rand() * total_flux

        if r < v1
            push!(M, M[end] + 1)
            push!(P, P[end])
        elseif r < v1 + v2
            push!(M, M[end] - 1)
            push!(P, P[end])
        elseif r < v1 + v2 + v3
            push!(M, M[end])
            push!(P, P[end] + 1)
        else
            push!(M, M[end])
            push!(P, P[end] - 1)
        end

        push!(T, T[end] + tau)
    end

    burnin = Int(length(T) ÷ 2)
    P_steady = P[burnin:end]
    return std(P_steady) / mean(P_steady)
end


function main()

    # question1()
    question2()

    # # last part of q2
    # lifetimes = [10, 100, 1000]
    # noises = [noise_given_lifetime(L) for L in lifetimes]

    # # Plot the noise vs protein lifetime
    # p2 = plot(lifetimes, noises,
    #     xscale = :log10,
    #     xlabel = "Protein Lifetime",
    #     ylabel = "Noise (σ / μ)",
    #     title = "Noise vs Protein Lifetime",
    #     marker = :o,
    #     lw = 2,
    #     legend = false)
    # savefig(p2, "noise_vs_protein_lifetime_QuestionB.png")
end

main()