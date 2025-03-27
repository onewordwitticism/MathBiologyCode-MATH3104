using Plots
using LinearAlgebra

function forward_diff_firstderiv(f, x, dx)
    return (f(x + dx) - f(x)) / dx
end

function backward_diff_firstderiv(f, x, dx)
    return (f(x) - f(x - dx)) / dx
end

function central_diff_firstderiv(f, x, dx)
    return (f(x + dx) - f(x - dx)) / (2 * dx)
end

function secondderiv_approx(f, x, dx)
    # approx 1 and 2 are equivalent - just writing down both for learning purposes
    approx_1 = forward_diff_firstderiv(f, x, dx) - backward_diff_firstderiv(f, x, dx) / dx
    approx_2 = (f(x + dx) - 2 * f(x) + f(x - dx)) / dx^2
    return approx_2
end


function solve_basic_ODE_explicit(x0, xn, n, y0)
    # solving dy/dx = x^2 - y
    y = zeros(n)
    dx = (xn - x0) / n
    number_of_steps = n

    x = range(start=x0 + dx, stop=xn, step=dx)
    y[1] = y0

    # note i is the counter; we have an (xi, yi) pair at each step
    # one could also 'calculate' the xi at each particular location
    # xi = x0 + i * dx
    for i in 1:number_of_steps-1
        y[i + 1] = y[i] + dx * (x[i]^2 - y[i])
    end

    display(plot(x, y, label="explicit", xlabel="x", ylabel="y", title="Basic ODE solution", legend=:topleft))

end

function solve_basic_ODE_implicit(x0, xn, n, y0)
    # solving dy/dx = x^2 - y
    y = zeros(n)
    dx = (xn - x0) / n
    number_of_steps = n

    x = range(start=x0 + dx, stop=xn, step=dx)
    y[1] = y0

    for i in 1:number_of_steps-1
        y[i + 1] = (y[i] + dx * (x[i + 1]^2)) / (1 + dx)
    end

    display(plot(x, y, label="implicit", xlabel="x", ylabel="y", title="Basic ODE solution", legend=:topleft))

end

function solve_basic_ODE_explicit_q2(t0, tn, n, y0, tau, ybar)
    # solving dy/dt = (1/tau) * (ybar - y)
    y = zeros(n)
    dt = (tn - t0) / n
    number_of_steps = n

    t = range(start=t0 + dt, stop=tn, step=dt)
    y[1] = y0

    for i in 1:number_of_steps-1
        y[i + 1] = y[i] + (dt/tau) * (ybar - y[i])
    end

    display(plot(t, y, label="explicit", xlabel="t", ylabel="y", title="Basic ODE solution", legend=:topleft))

end

function solve_basic_ODE_implicit_q2(t0, tn, n, y0, tau, ybar)
    # solving dy/dt = (1/tau) * (ybar - y)
    y = zeros(n)
    dt = (tn - t0) / n
    number_of_steps = n

    t = range(start=t0 + dt, stop=tn, step=dt)
    y[1] = y0

    for i in 1:number_of_steps-1
        y[i + 1] = (y[i] + (dt*ybar/tau)) / (1 + dt/tau)
    end

    display(plot(t, y, label="implicit", xlabel="t", ylabel="y", title="Basic ODE solution", legend=:topleft))

end



function solve_basic_ODE_chatgpt(t0, tn, n, y0, tau, ybar)
    # Solving dy/dt = (1/tau) * (ybar - y) using both Explicit and Implicit Euler

    dt = (tn - t0) / n
    number_of_steps = n

    # Time grid
    t = range(start=t0 + dt, stop=tn, step=dt)

    # Initialize solution arrays
    y_explicit = zeros(n)
    y_implicit = zeros(n)

    # Initial conditions
    y_explicit[1] = y0
    y_implicit[1] = y0

    # Compute both Explicit and Implicit solutions
    for i in 1:number_of_steps-1
        y_explicit[i + 1] = y_explicit[i] + (dt/tau) * (ybar - y_explicit[i])
        y_implicit[i + 1] = (y_implicit[i] + (dt * ybar / tau)) / (1 + dt/tau)
    end

    # Compute the analytical solution
    y_analytical = (y0 - ybar) .* exp.(-t ./ tau) .+ ybar

    # Plot all three solutions
    plot(t, y_explicit, label="Explicit Euler", xlabel="t", ylabel="y", title="Explicit vs. Implicit vs. Analytical", legend=:bottomright, size=(900, 600))
    plot!(t, y_implicit, label="Implicit Euler", linestyle=:dash)
    plot!(t, y_analytical, label="Analytical Solution", linestyle=:dot, lw=2, color=:black)  # Dashed black line for exact solution

    display(plot!())  # Show the plot
end



function main()
    f(x) = sin(x)  # given function
    x = π / 4     # point of differentiation
    dx = 0.01      # Step-size

    println("First derivative of f at x = $x using forward difference: $(forward_diff_firstderiv(f, x, dx))")
    println("First derivative of f at x = $x using backward difference: $(backward_diff_firstderiv(f, x, dx))")
    println("First derivative of f at x = $x using central difference: $(central_diff_firstderiv(f, x, dx))")
    println("Second derivative of f at x = $x using central difference: $(secondderiv_approx(f, x, dx))")

    x0 = 0.0
    xn = 10.0
    n = 100
    y0 = 0

    # solve_basic_ODE_explicit(x0, xn, n, y0)
    # solve_basic_ODE_implicit(x0, xn, n, y0)

    t0 = 0
    tn = 10
    tau = 0.1
    ybar = 10

    # solve_basic_ODE_explicit_q2(t0, tn, n, y0, tau, ybar)
    # solve_basic_ODE_implicit_q2(t0, tn, n, y0, tau, ybar)

    solve_basic_ODE_chatgpt(t0, tn, n, y0, tau, ybar)

end

main()