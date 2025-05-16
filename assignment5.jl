using CairoMakie

function question1()
    km = 0.02  # k1_act
    kd = 0.01  # k1_deact

    kn = 0.01  # k2_act
    kr = 0.02  # k2_deact

    tol=1e-10
    maxiter=30
    e = 1e-12

    # need epsilon otherwise divides through by 0 breaking stuff
    P1(S1, P) = (S1*P*kd) / (max(km + P - S1*P, e))
    P2(S2, P) = (S2*P*kr) / (max(kn + P - S2*P, e))
    C(P, S1, S2) = P + P1(S1, P) + P2(S2, P) - 1
    C_deriv(P, S1, S2) = 1 + ((S1*kd*(km + P)) / max(km + P - S1*P, e)^2) + ((S2*kr*(kn + P)) / max(kn + P - S2*P, e)^2)

    # need 0 < S1 and S2 < 2.  Probably assume S2 can't go below 0.

    # build the heat-map arrays
    S1_vals = range(0, 2, length=120)   
    S2_vals = range(0, 2, length=120) 

    P1_grid = Matrix{Float64}(undef, length(S1_vals), length(S2_vals))
    P2_grid = similar(P1_grid)

    for (i, s1) in enumerate(S1_vals), (j, s2) in enumerate(S2_vals)
        P = 0.5

        for _ in 1:maxiter
            FP = C(P, s1, s2)
            dFP  = C_deriv(P, s1, s2)
            Pn = clamp(P - FP/dFP, 0.0, 1.0)
            
            if abs(Pn - P) ≤ tol
                P = Pn
                break
            end
            P = Pn
        end

        val1 = P1(s1, P)
        val2 = P2(s2, P)
        
        P1_grid[i, j] = clamp(val1, 0.0, 1.0)
        P2_grid[i, j] = clamp(val2, 0.0, 1.0)
    end

    clims = (0.0, 1.0)         # 0 = blue, 1 = yellow

    fig = Figure(resolution = (1050, 400))

    Label(fig[0, 1:3], "Saturated-constant regime  (Kₘ = 0.02, K_d = 0.01, Kₙ = 0.01, K_r = 0.02)",
      fontsize = 18, padding = (0,0,10,0), halign = :center)

    ax1 = Axis(fig[1,1], title = "P₁ steady state", xlabel = "S₁", ylabel = "S₂")
    ax2 = Axis(fig[1,2], title = "P₂ steady state", xlabel = "S₁", ylabel = "S₂")
    
    h1 = heatmap!(ax1, S1_vals, S2_vals, P1_grid; colormap = :viridis, colorrange = clims)
    h2 = heatmap!(ax2, S1_vals, S2_vals, P2_grid; colormap = :viridis, colorrange = clims)
    
    Colorbar(fig[1,3], h1, label = "fraction", width = 15)   
    save("assignment5_q1_saturated.png", fig)  
    fig

end


function main()

    question1()
end

main()