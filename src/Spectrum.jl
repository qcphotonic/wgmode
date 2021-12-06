#= Calculate the Wispering Gallery eigenmodes (n, l, m). Generally speaking, WGM refers to the mode where n=1. We can compute modes of 
any order n utilizing spectrum function. 
    
pre_loop(l, lambda, mode, n, R)
    l::Int 
        l mode number
    lambda::Array{float64}
        wavelenght searching range
    mode::String ("TE" or "TM")
        mode type
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    
    return plot

    Given a l mode number, plot all modes in the wavelenght range.

spectrum(lambda, mode, n_num, n, R; Q_factor=18, option="n_num depend")
    lambda::Array{float64}
        wavelenght searching range
    mode::String ("TE" or "TM")
        mode type
    n_num::Int
        n mode number
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    Q_factor{::Int, "open"}
        quality factor option, calculate quality factor with "open" or set an integer.
    option::String ("n_num depend" or "all")
        calculate radial number up to n_num or all possible modes

    return DataFrame

    Given electric and geometric parameters, focus on the interest wavelength region and mode, calcualte modes to order n.

view_spectrum(lambda, data, order; view_mode = "")
    lambda::Array{float64}
        wavelenght searching range
    data::DataFrame
        spectrum info
    order::Int
        n mode select
    view_mode::String ("" or "details")
        view mode option, "details" is all information display mode
    
    return plot

    Given wavelength range and spectrum info, plot the spectrum
=#

using SpecialFunctions
using Gadfly
using Compose
using Colors
using Roots
using Zygote
using Optim

using Distributed
using ProgressMeter

using DataFrames
using CSV

using Printf: @sprintf


# bessely(l-1/2, x0)/bessely(l+1/2, x0) overflows when l is much larger than fixed x0. Possible solution but not supported!
#= function besselho(l, x)
        try
            p = besselh(l-1/2, 2, x)/besselh(l+1/2, 2, x)
        catch e
            besselhr(l, x) = l < 1 ? besselh(l+1, 2, x)/besselh(l, 2, x) : 2*l/x - 1/besselhr(l - 1, x)
            p = 1/besselhr(l-1/2, x)
        end
    end
=#
# abs^2 of the complex critical function
function cf_modulus(x0, kappa, l, mode, n)
    nk = n+im*10.0^kappa
    if mode == "TE"
        P = nk
    else
        P = 1/nk
    end
    x = nk*x0
    val = besselh(l-1/2, 2, x0)/besselh(l+1/2, 2, x0)-P*besselj(l-1/2, x)/besselj(l+1/2, x)-l*(1/x0-P/x)
    return abs(val)^2
end

#= function besselyo(l, x)
        try
            p = bessely(l-1/2, x)/bessely(l+1/2, x)
        catch e
            besselyr(l, x) = l < 1 ? bessely(l+1, x)/bessely(l, x) : 2*l/x - 1/besselyr(l - 1, x)
            p = 1/besselyr(l-1/2, x)
        end
    end
=#
# the real critical function
function cf_real(x0, l, mode, n)
    nk = n
    if mode == "TE"
        P = nk
    else
        P = 1/nk
    end
    x = nk*x0
    val = bessely(l-1/2, x0)/bessely(l+1/2, x0)-P*besselj(l-1/2, x)/besselj(l+1/2, x)-l*(1/x0-P/x)
    return val
end

# find roots function utilizing bisection method
function zeros_bsolver(f, range; accuracy=1e-3)
    xmin, xmax = range
    x_range = xmin:accuracy:xmax
    xdata = zeros(0)
    u = 0
    for i in x_range
        d = u
        u = f(i)
        if u*d < 0
            root = find_zero(f, (i-accuracy, i))
            append!(xdata, root)
        end
    end
    return xdata
end

function first_zeros_bsolver(f, range; accuracy=1e-3)
    xmin, xmax = range
    x_range = xmin:accuracy:xmax
    u = 0
    for i in x_range
        d = u
        u = isnan(f(i)) ? 0 : f(i)
        if u*d < 0
            root = find_zero(f, (i-accuracy, i))
            return root
        end
    end
    return nothing
end

# x derivative of the cf_modulus function
cf_modulus_dx(x, kappa, l, mode, rn) = gradient(x0->cf_modulus(x0, kappa, l, mode, rn), x)[1]

# material property: absorption curve 
function absorption_curve(lambda)
    Nx = 5.3e14
    NY1 = 1.5e21
    NY2 = 10e21
    alpha_sp = 1.1e-26*lambda*exp(-10.7*(1-780/lambda)^2)*Nx
    alpha_bp = 5.5e-27*lambda*exp(-9.6*(1-500/lambda)^2)*NY1
    alpha_fe = 1.6e-23*(exp(-32*(1-480/lambda)^2)+exp(-15.4*(1-1100/lambda)^2)/16)*NY2
    alpha_Bc = exp(19*(1-lambda/379))
    alpha_Rc = exp(16.7*(1-6800/lambda))
    return alpha_sp+alpha_bp+alpha_fe+alpha_Bc+alpha_Rc
end

# absorption Quality factor calculation
function Q_ab(lambda)
    return 1e7*2*pi/(lambda*absorption_curve(lambda))
end

# tatal Quality factor calculation
function Q_total(n, Qrad, lambda)
    Q_ab = 1e7*2*pi*n/(lambda*absorption_curve(lambda))
    Qrad = 10^Qrad
    Q = 1/(1/Q_ab+1/Qrad)
    return log10(Q)
end

# test before computing, make sure the appropriate match between l and wavelenght range
function pre_loop(l, lambda, mode, n, R)
    b, a = [2*pi*R*1e3/i for i in lambda]
    kappa = -10
    Gadfly.plot(y=[x0->cf_modulus(x0, kappa, l, mode, n)], xmin=[a], xmax=[b], Coord.cartesian(ymin=0, ymax=10), Stat.func(1000), Geom.line)
end

# return spectrum, Q_radiation in the given wavelength window and given l number
function spectrum_l(l, lambda, mode, n, R; Q_factor="open")
    b, a = [2*pi*R*1e3/i for i in lambda]
    kappa = -10
    bar = minimum([cf_modulus(i, kappa, 0, mode, n) for i in a:1e-2:b])
    x_zeros = zeros_bsolver(x -> cf_modulus_dx(x, kappa, l, mode, n), [a, b])
    xdata = zeros(0)
    for x in x_zeros
        if cf_modulus(x, kappa, l, mode, n) < bar + 10
            append!(xdata, x)
        end
    end
    # ...
    if Q_factor == "open"
        kappa_data = zeros(0)
        for x0 in xdata
            opt = optimize(kappa->cf_modulus(x0, kappa, l, mode, n), -20, 0)
            append!(kappa_data, Optim.minimizer(opt))
        end
        Qrad = [log10(n)-i for i in kappa_data]
    else
        Qrad = fill(Q_factor, (length(xdata),))
    end
    spectrum = [2*pi*R*1e3/i for i in xdata]
    return spectrum, Qrad
end

# spectrum calculation function
function spectrum(lambda, mode, n_num, n, R; Q_factor=18, option="n_num depend")

    rn = n
    # set sweep parameters to find l_max
    sweep_step = 50
    start = 0

    # set progressbar
    p = ProgressMeter.Progress(101, 0.1, "Initializing... ")
    msg = @sprintf "%s" p.desc 
    ProgressMeter.printover(p.output, msg, p.color)
    ProgressMeter.printvalues!(p, (); color = p.color, truncate = false)
    sleep(0.1)

    # sweep to find l_max=start!
    while true
        l_zeros = zeros_bsolver(l->cf_real(2*pi*R*1e3/lambda[1], l, mode, rn), [start, start+2*sweep_step])
        if length(l_zeros) > 4 && start+2*sweep_step-l_zeros[end] > l_zeros[end]-l_zeros[end-4]
            start = floor(l_zeros[end])
            break
        elseif length(l_zeros) == 0
            l_zeros = zeros_bsolver(l->cf_real(2*pi*R*1e3/lambda[1], l, mode, rn), [start-sweep_step, start])
            start = floor(l_zeros[end])
            break
        end
        start += sweep_step
    end
    p.desc = "Computing...    "
    
    # decrease l from l_max to the given n order precision, calcualte all modes
    var =  0
    lambda_end = lambda[2]
    N = zeros(0)
    L = zeros(0)
    W = zeros(0)
    Q = zeros(0)
    var1 = 1
    l = start
    dlambda = lambda[2]-lambda[1]
    while var == 0
        lambda_sweep = [lambda[1], lambda_end]
        spectrum, Qrad = spectrum_l(l, lambda_sweep, mode, rn, R, Q_factor=Q_factor)
        Qtt = zeros(0)
        for (i, q) in enumerate(Qrad)
            append!(Qtt, Q_total(n, q, spectrum[i]))
        end
        n = length(spectrum)
        if n != 0
            if n > n_num && option=="n_num depend"
                if spectrum[n_num]>lambda[2]
                    var = 1
                end
            end
            if n > 2 && var == 0
                if lambda_end-spectrum[1] < spectrum[1]-spectrum[3]
                    lambda_end += 2*(spectrum[1]-spectrum[2])
                end
            elseif lambda_end - spectrum[1] < spectrum[1] - lambda[1] && var == 0
                lambda_end += dlambda
            end
            for i in 1:n
                if spectrum[i] <= lambda[2] && (i<=n_num || option=="all") && var == 0
                    append!(N, floor(Int, i))
                    append!(L, floor(Int, l))
                    append!(W, spectrum[i])
                    append!(Q, Qtt[i])
                end
            end
            if n > 1 && option == "n_num depend"
                dsp = spectrum[1]-spectrum[end]
                percentage = n*dsp/(dlambda*(n-1)+n_num*dsp)
                if percentage > 1
                    percentage = 1
                end
                ProgressMeter.update!(p, floor(Int, 100*percentage))
            elseif option == "n_num depend"
                percentage = 1/(n_num+dlambda/(spectrum[1]-lambda[1]))
                ProgressMeter.update!(p, floor(Int, 100*percentage))
            end
        else
            var1 += 1
        end
        if var1 == 2
            println("Error: l_max found not accurate! Suggestion: pre_loop function")
            break
        end
        l -= 1
        if option == "all"
            var = 0
            percentage=(1-l/start)^2
            ProgressMeter.update!(p, floor(Int, 100*percentage)>99 ? 100 : floor(Int, 100*percentage))
            if l < 0
                break
            end
        end
    end
    p.desc = "Finished ✓      "
    ProgressMeter.update!(p, 101)

    df = DataFrame(n = Int.(N), l = Int.(L), wavelength = W, Q = Q)
    return sort!(df)
end

# plot settings for specific order
function select_order(n, N)
    C = zeros(N)
    C[n] = 1
    return C
end

# plot the spectrum
function view_spectrum(lambda, data, order; view_mode = "")
    n = data.n 
    l = data.l 
    w = data.wavelength
    q = data.Q

    xp = zeros(0)
    yp = zeros(0)
    index = zeros(0)
    np = zeros(0)
    lp = zeros(0)
    wp = zeros(0)
    qp = zeros(0)
    for i in 1:length(n)
        append!(xp, w[i], w[i])
        append!(yp, 0, q[i])
        append!(index, i, i)
        append!(np, n[i], n[i])
        append!(lp, l[i], l[i])
        append!(wp, w[i], w[i])
        append!(qp, q[i], q[i])
    end
    df = DataFrame(index=index, x_col=xp, y_col=yp, n=np, l=lp, wavelength=wp, Qtt=qp);
    gd = groupby(df, :n)
    N = length(gd)
    if order > N
        println("Error, order selected is out of range of calculation")
    else
        datag = groupby(data, :n)
        C_default = [1/i for i=1:N]
        C = (order == 0) ? C_default : select_order(order, N)
        if order != 0 && view_mode == "details"
            fn = datag[order].wavelength[1]/140
            L=layer(gd[order], x=:x_col, y=:y_col, group=index, Geom.line, Theme(default_color=HSV(210, 1, C[order]), line_width=0.1mm))
            pl = Gadfly.plot(L, Coord.cartesian(xmin=lambda[1], xmax=lambda[2], ymin=0), Guide.xlabel("Wavelength/nm"), Guide.ylabel("Q/-log10"), 
                 Guide.annotation(Compose.compose(context(), Compose.text([i for i=datag[order].wavelength], [i-1 for i=datag[order].Q], 
                 [string(i) for i=datag[order].l]), stroke("gray50"), fontsize(fn*1pt))), Guide.xticks(ticks=datag[order].wavelength, orientation=:vertical), 
                 Theme(background_color=HSV(0, 1, 0), grid_line_width=0mm))
        else
            L=[layer(gd[i], x=:x_col, y=:y_col, group=index, Geom.line, Theme(default_color=HSV(210, 1, C[i]), line_width=0.1mm)) for i in 1:N]
            pl = Gadfly.plot(L..., Coord.cartesian(xmin=lambda[1], xmax=lambda[2], ymin=0), Guide.xlabel("Wavelength/nm"), Guide.ylabel("Q/-log10"), 
                 Theme(background_color=HSV(0, 1, 0), grid_line_width=0mm))
        end
        set_default_plot_size(28cm, 13cm)
        display(pl)
    end
end
