using Plots
using LinearAlgebra

function q2_diffusion(x_lb, x_rb, dx, dt, ic_p1, ic_p2, tmax)
    times_to_save = [25,50,75,100,125,150]
    n = Int((x_rb - x_lb) / dx)
    x = range(start=x_lb+dx, stop=x_rb, step=dx)
    D = dx^2 / (2*dt)
    t = 0

    # 2nd derivative (dxx) matrix operator as attained by finite difference method; when applied to a vector, it returns the 2nd derivative of the vector
    # typically the vector is a vector of points (e.g. x; and then we get dxx at those points)
    M = diff(diff([zeros(1, n-1) 1; I(n); 1 zeros(1, n-1)], dims=1), dims=1) / dx^2

    # set IC
    # note that dx * sum(rho) = 1 if pdf; or total_cells if there is a count population (just renormalizing the pdf)
    rho = zeros(n)
    rho[floor(Int, ((ic_p1[2] + 20) / 40)*n)] =  ic_p1[1] / dx
    rho[floor(Int, ((ic_p2[2] + 20) / 40)*n)] =  ic_p2[1] / dx

    # save initial 50 figure
    p = plot(x, rho, ylim=(0.0, 5000), xlabel="Space", ylabel="Cell concentration - Number of Cells / dx",
    label="Time = $(round(t, digits=1))", legend=:topright, show=true)
    filename = "plot_output_t$(round(Int, t))_q2.png"
    savefig(p, filename)  

    println("Total sum of the population at the START (sum(rho*dx)): $(sum(rho*dx))")

    while t < tmax-dt
        t = t + dt
        rho_old = rho
        rho = (I(n)-dt*D*M)\rho_old
        p = plot(x, rho, xlabel="Space", ylabel="Cell concentration - Number of Cells / dx",
        label="Time = $(round(t, digits=1))", legend=:topright, show=true)

        display(p)

        if round(Int, t) in times_to_save
            filename = "plot_output_t$(round(Int, t))_q2.png"
            savefig(p, filename)  
        end
    end
    println("Total sum of the population at the END (sum(rho*dx)): $(sum(rho*dx))")

end

function q3_fisher_kpp(x_lb, x_ub, dx, dt, D, alpha, K, tmax)
    # transformation
    t_c = 1 / alpha
    rho_c = K
    x_c = sqrt(D/alpha)
    ttilda=0
    t=0

    xtilda_lb = x_lb / x_c
    xtilda_ub = x_ub / x_c
    dxtilda = dx / x_c
    dttilda = dt / t_c
    tmaxtilda = tmax / t_c

    xtilda = range(start=xtilda_lb+dxtilda, stop=xtilda_ub, step=dxtilda)
    n = length(xtilda)
    ttilda = 0
    rho = zeros(n)
    total_space_length = xtilda_ub - xtilda_lb
    
    # set ic of rho
    index1 = floor(Int, (xtilda_ub - 0.5/x_c)*n/total_space_length)
    index2 = floor(Int, (xtilda_ub + 0.5/x_c)*n/total_space_length)
    rho[index1:index2] .= 1/4

    M = diff(diff([zeros(1, n-1) 1; I(n); 1 zeros(1, n-1)], dims=1), dims=1) / dxtilda^2

    p = plot(xtilda, rho, ylim=(0.0, 1),
    xlabel="x̃ (dimensionless distance)",
    ylabel="ρ (tumor density, population / x̃)\nK=1, acts like PDF", legend=false, label="Time = $(round(t, digits=1))", show=false)
    filename = "plot_output_ttilda$(round(Int, t))_q3.png"
    savefig(p, filename)  

    # Decide how often to plot
    num_plots = 10  
    update_interval = 150  
    step = 0

    # this is specifically for question 3c.  Tracking the front of the tumour growth.
    # Nothing to do with the iteration method.
    front_positions_right = Float64[]
    threshold = 1e-4

    while ttilda < tmaxtilda-dttilda
        ttilda = ttilda + dttilda
        step += 1
        
        # calculate diffusion component of rho_n (rho); using rho_n-1 (rho_old)
        rho_old = rho
        rho = (I(n)-dttilda*D*M)\rho_old

        # still at same time step but now treat rho_n after the diffusion update as 'rho_n-1'
        rho_old = rho
        rho = (rho_old .* (1 .+ dttilda)) ./ (1 .+ dttilda .* rho_old)

        # specifically for q3c; nothing to do with iteration method which is just the lines above.
        right_index = findlast(rho .> threshold)

        if right_index !== nothing
            x_right = xtilda[right_index]
            push!(front_positions_right, x_right)
        end

        # Only plot every 'update_interval' steps
        if step % update_interval == 0
            plot!(xtilda, rho, xlabel="x̃ (dimensionless distance)",
            ylabel="ρ (tumor density, population / x̃)\nK=1, acts like PDF", ylim=(0.0, 1), legend=false, label="Time = $(round(ttilda, digits=1))", show=false)
            display(p)  # Ensure the plot updates

            filename = "plot_output_ttilda=$(round(Int, ttilda))_q3.png"
            savefig(p, filename)  
        end

    end

    # Compute discrete growth rate (cm/day)
    # mine appears closer to 0.1 cm/day as opposed to the theoretical 0.08944 cm/day = 2 * sqrt(D/alpha)
    front_growth_rate = (diff(front_positions_right) .* x_c) ./ (dttilda * t_c)

    # Plot and save
    p_growth = plot(front_growth_rate, xlabel="Step", ylabel="Growth rate (cm/day)",
    title="Tumor front growth rate (right side)", legend=false)
    savefig(p_growth, "tumor_front_growth_rate_cm_per_day.png")
    display(p_growth)

end

function main()
    dx_q2 = 0.2
    dt_q2 = 0.1
    x_lb_q2 = -20
    x_rb_q2 = 20
    ic_p1_q2 = [1000, 10]
    ic_p2_q2 = [500, -10]
    tmax_q2 = 150

    x_lb_q3 = -5
    x_ub_q3 = 5
    dx_q3 = 0.1
    dt_q3 = 0.1
    D = 0.02
    alpha = 0.1
    K = 1
    tmax_q3 = 300

    # just run one at once; comment out whichever one you don't want to show.
    # q2_diffusion(x_lb_q2, x_rb_q2, dx_q2, dt_q2, ic_p1_q2, ic_p2_q2, tmax_q2)
    q3_fisher_kpp(x_lb_q3, x_ub_q3, dx_q3, dt_q3, D, alpha, K, tmax_q3)

end

main()
