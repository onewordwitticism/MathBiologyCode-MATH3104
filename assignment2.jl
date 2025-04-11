using Plots
using PlotUtils
using LinearAlgebra
using JSON

# diffusion; just some diffusion term D*dxx*c(x,t) - 1D heat Eq
# diffusion-reaction - diffusion and a creation/depletion term (ie. logistic growth) - Fisher KPP
# diffusion-drift - diffusion term plus a drift term towards something.  IE.  there is a deterministic component - Keller Segel

function read_config(file_path::String)
    return JSON.parsefile(file_path)
end

function q1_KS_uncoupled(a, b, D, chi, dx, dt, no_cells, tmax)
    times_to_save = [2.5, 5, 10, 15, 25, 35, 50, 75, 100, 125, 150, 200]

    # functions of use for solving the DE
    u(x) = (2pi .* (x .- a .- 3)) ./ (b - a)
    S(x) = 2 .* cos.(u(x)) .+ 2
    Sx(x) = (-4pi .* sin.(u(x))) ./ (b - a)
    Sxx(x) = -8 .* cos.(u(x)) .* (pi / (b - a))^2
    # --------------------------------------

    # creating discretization
    n = Int((b - a) / dx)  # number of discrete points on my x grid; indicator is 1,2, ... -> n
    x = range(start=a+dx, stop=b, step=dx)
    t = 0

    # 2nd derivative (dxx) matrix operator as attained by finite difference method; when applied to a vector, it returns the 2nd derivative of the vector
    # typically the vector is a vector of points (e.g. x; and then we get dxx at those points)
    # it incorporates movement along the x grid.  xj - xj-1 etc.
    # each rho derivate at j is being computed at rhoj - rhoj-1 / dx etc.  Iterating through the grid for rho.
    M = diff(diff([zeros(1, n-1) 1; I(n); 1 zeros(1, n-1)], dims=1), dims=1) / dx^2

    # set IC
    rho = zeros(n)
    rho_ic_lb = a + 0.75*(b-a)
    rho_ic_up = b
    rho .= (4 * no_cells) / (b - a)
    rho[x .< rho_ic_lb] .= 0
    rho[x .>= rho_ic_up] .= 0

    # graphiing IC
    p = plot(x, rho, label="ρ(x) - Cell Concentration", xlabel="Space(x)", ylabel="ρ(x)", color=:blue, ylim=(0, 40), legend=:topright)

    plot!(p, x, S.(x), right_y=true,
    ylabel="Concentration ((cells or attractant)/dx)",
        label="S(x) - Chemoattractant Concentration",
        color=:red)

    x_pos = x[end] - (x[end] - x[1]) * 0.1  # Fixed x position near the right edge
    y_pos_time = 3  # Move the Time annotation slightly higher
    y_pos_D = 1     # Fixed y position for the D annotation
    
    annotate!(p, (x_pos, y_pos_time, text("Time = $(round(t, digits=2))", :right, 10)))
    annotate!(p, (x_pos, y_pos_D, text("D = $(D)", :right, 10)))
    display(p)

    filename = "plot_output_t$(round(Int, t))_assignment2_q1_t0.png"
    savefig(p, filename)  

    k = dt*chi/2*dx
    Sx_jm1 = circshift(Sx(x), 1)
    Sx_jp1 = circshift(Sx(x), -1)

    # simulate solution
    while t < tmax-dt
        t = t + dt
        rho_o = rho

        # solve PDE for diffusion component only.
        rho = (I(n)-dt*D*M)\rho_o

        # still at same time step but now treat rho_n after the diffusion update as 'rho_n-1'
        # we are now solving for the ODE drift component using Lax Friedrichs method
        rho_o = rho
        rho_jm1 = circshift(rho_o, 1)
        rho_jp1 = circshift(rho_o, -1)

        # updated rho with drift.
        rho = (0.5 .* (rho_jp1 .+ rho_jm1)) .- (k .* ((rho_jp1.*Sx_jp1) - (rho_jm1.*Sx_jm1)))

        # end computation --------------------------------------------------

        # plots
        p = plot(x, rho,
         label="ρ(x) - Cell Concentration", xlabel="Space(x)", ylabel="ρ(x)", color=:blue, ylim=(0, 40), legend=:topright)

        plot!(p, x, S.(x), right_y=true,
         ylabel="Concentration ((cells or attractant)/dx)",
          label="S(x) - Chemoattractant Concentration",
           color=:red)

        annotate!(p, (x_pos, y_pos_time, text("Time = $(round(t, digits=2))", :right, 10)))
        annotate!(p, (x_pos, y_pos_D, text("D = $(D)", :right, 10)))

        plot!(p, label="Time = $(round(t, digits=2)), D = $(D)")

        display(p)

        if round(Int, t) in times_to_save
            filename = "plot_output_t$(round(Int, t))_assignment2_q1.png"
            savefig(p, filename)  
        end

    end

end

function q2_KS_coupled(a, b, no_cells, D, Ds, chi, alpha, beta, dx, dt, tmax)
    times_to_save = [1,2,3,4,5,6,7,8,9, 10, 15, 25, 35, 50, 75, 100, 125, 150, 200]

    # creating discretization
    n = Int((b - a) / dx)  # number of discrete points on my x grid; indicator is 1,2, ... -> n
    x = range(start=a+dx, stop=b, step=dx)
    t = 0

    # 2nd derivative (dxx) matrix operator as attained by finite difference method; when applied to a vector, it returns the 2nd derivative of the vector
    # typically the vector is a vector of points (e.g. x; and then we get dxx at those points)
    # it incorporates movement along the x grid.  xj - xj-1 etc.
    # each rho derivate at j is being computed at rhoj - rhoj-1 / dx etc.  Iterating through the grid for rho.
    # we are going to use M for the second derivative for rho and S
    M = diff(diff([zeros(1, n-1) 1; I(n); 1 zeros(1, n-1)], dims=1), dims=1) / dx^2

    # set IC
    rho = zeros(n)
    S = zeros(n)
    S .= no_cells/(b-a)  # 100 cells divided equiprobable.  

    rho_ic_lb = a + (b-a)/4
    rho_ic_up = a + 0.75*(b-a)
    rho .= (2 * no_cells) / (b - a)
    rho[x .< rho_ic_lb] .= 0
    rho[x .>= rho_ic_up] .= 0

    p = plot(x, rho, label="ρ(x) - Cell Concentration", xlabel="Space(x)", ylabel="ρ(x)", color=:blue, ylim=(0, 20), legend=:topright)

    plot!(p, x, S, right_y=true,
    ylabel="Concentration ((cells or attractant)/dx)",
        label="S(x) - Chemoattractant Concentration",
        color=:red)

    x_pos = x[end] - (x[end] - x[1]) * 0.1  
    y_pos_time = 3  
    y_pos_D = 1    
    
    annotate!(p, (x_pos, y_pos_time, text("Time = $(round(t, digits=2))", :right, 10)))
    annotate!(p, (x_pos, y_pos_D, text("D = $(D)", :right, 10)))
    display(p)

    filename = "plot_output_t$(t)_assignment2_q2_D=$(D).png"
    savefig(p, filename)

    k = dt*chi/2*dx

    # simulate solution
    while t < tmax-dt
        t = t + dt
        rho_o = rho

        # solve PDE for diffusion component only.
        # step 1. diffusion
        rho = (I(n)-dt*D*M)\rho_o

        # still at same time step but now treat rho_n after the diffusion update as 'rho_n-1'
        # we are now solving for the ODE drift component using Lax Friedrichs method
        rho_o = rho
        rho_jm1 = circshift(rho_o, 1)
        rho_jp1 = circshift(rho_o, -1)
        
        # solve linear system for S(x)
        S = (M - beta * I(n))\(-alpha.*rho)

        # Get derivative of S wrt x; for use in Lax-Friedrichs below
        Sx = (circshift(S, -1) .- circshift(S, 1)) ./ (2dx)

        Sx_jm1 = circshift(Sx, 1)
        Sx_jp1 = circshift(Sx, -1)

        # # updated rho with chemotactic drift.  (Lax-Friedrich eq)
        rho = (0.5 .* (rho_jp1 .+ rho_jm1)) .- (k .* ((rho_jp1.*Sx_jp1) - (rho_jm1.*Sx_jm1)))
        # end computation --------------------------------------------------

        # plots
        p = plot(x, rho,
         label="ρ(x) - Cell Concentration", xlabel="Space(x)", ylabel="ρ(x)", color=:blue, ylim=(0, 20), legend=:topright)

        plot!(p, x, S, right_y=true,
         ylabel="Concentration ((cells or attractant)/dx)",
          label="S(x) - Chemoattractant Concentration",
           color=:red)

        annotate!(p, (x_pos, y_pos_time, text("Time = $(round(t, digits=2))", :right, 10)))
        annotate!(p, (x_pos, y_pos_D, text("D = $(D)", :right, 10)))
        
        display(p)

        save_time = round(Int, t)
        if save_time in times_to_save
            filename = "plot_output_t$(save_time)_assignment2_q2_D=$(D).png"
            savefig(p, filename)
        end
    end
end

function main()
    config = read_config("config_assignment2.json")
    
    q1_params = config["q1_parameters"]
    q2_params = config["q2_parameters"]

    q1_KS_uncoupled(
        q1_params["x_lb"],
        q1_params["x_up"],
        q1_params["D_diffusion_coefficient"],
        q1_params["chi_chemo_sens"],
        q1_params["dx"],
        q1_params["dt"],
        q1_params["no_cells"],
        q1_params["tmax"]
        )

    q2_KS_coupled(
        q2_params["x_lb"],
        q2_params["x_up"],
        q2_params["no_cells"],
        0.5,
        q2_params["Ds"],
        q2_params["chi"],
        q2_params["alpha"],
        q2_params["beta"],
        q2_params["dx"],
        q2_params["dt"],
        q2_params["tmax"]
    )

    q2_KS_coupled(
        q2_params["x_lb"],
        q2_params["x_up"],
        q2_params["no_cells"],
        3,
        q2_params["Ds"],
        q2_params["chi"],
        q2_params["alpha"],
        q2_params["beta"],
        q2_params["dx"],
        q2_params["dt"],
        q2_params["tmax"]
    )

end

main()
