#= Calculate the nonlinearity coupling strength gnl.
    
overlap_nonlinearity(field_parameters, n, R, d; digit=3)
    field_parameters::Vector{Vector{Any}}
        selected mode number, mode type for two fields
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    d::Vector{float64}
        nonlinearity coefficient
    digit::Int
        effective number of the result
    
    return float64

    Given the info of selected two fields and structural parameters, calculate the nonlinearity coupling strenght.
=#

using DataFrames
using ProgressMeter

function detuning(data_f, data_shg; threshold = "")
    n_f = data_f.n
    l_f = data_f.l
    Q_f = data_f.Q
    wav_f = data_f.wavelength
    mode_f = data_f.mode
    n_shg = data_shg.n
    l_shg = data_shg.l
    Q_shg = data_shg.Q
    wav_shg = data_shg.wavelength
    mode_shg = data_shg.mode

    wav_eff_f = zeros(0)
    wav_eff_shg = zeros(0)

    n_f_new = zeros(0)
    l_f_new = zeros(0)
    wav_f_new = zeros(0)
    Q_f_new = zeros(0)
    mode_f_new = String[]
    n_shg_new = zeros(0)
    l_shg_new = zeros(0)
    wav_shg_new = zeros(0)
    Q_shg_new = zeros(0)
    Ratio = zeros(0)
    mode_shg_new = String[]

    for (i, w_shg) in enumerate(wav_shg)
        for (j, w_f) in enumerate(wav_f)
            G = (1-2*w_shg/w_f)*10^Q_shg[i]
            ratio = 1/(1+G^2)
            if typeof(threshold) == Float64
                if ratio > threshold 
                    append!(wav_eff_shg, w_shg)
                    append!(wav_eff_f, w_f)
                end
            else
                append!(Ratio, ratio)
                append!(n_f_new, n_f[j])
                append!(l_f_new, l_f[j])
                append!(wav_f_new, wav_f[j])
                append!(Q_f_new, Q_f[j])
                append!(mode_f_new, [mode_f[j]])
                append!(n_shg_new, n_shg[i])
                append!(l_shg_new, l_shg[i])
                append!(wav_shg_new, wav_shg[i])
                append!(Q_shg_new, Q_shg[i])
                append!(mode_shg_new, [mode_shg[i]])
            end
        end
    end

    if typeof(threshold) == Float64
        return wav_eff_f, wav_eff_shg
    else
        df = DataFrame(ratio_g = Ratio, n_f = Int.(n_f_new), l_f = Int.(l_f_new), mode_f = mode_f_new, wavelength_f = wav_f_new, Q_f = Q_f_new, 
            n_shg = Int.(n_shg_new), l_shg = Int.(l_shg_new), mode_shg = mode_shg_new, wavelength_shg = wav_shg_new, Q_shg = Q_shg_new)
        return sort!(df, rev=true)
    end
end

function overlap_array(df, n, R, d, limitation, cutoff; pgbar = "on")
    L = length(detune.ratio_g)
    if cutoff > L
        println("Error! cutoff=$cutoff should be smaller than mode pairs=$L")
        return 0
    end
    delta_f, delta_shg = limitation

    ratio_new = zeros(0)

    n_f_new = zeros(0)
    l_f_new = zeros(0)
    m_f_new = zeros(0)
    wav_f_new = zeros(0)
    Q_f_new = zeros(0)
    mode_f_new = String[]

    n_shg_new = zeros(0)
    l_shg_new = zeros(0)
    m_shg_new = zeros(0)
    wav_shg_new = zeros(0)
    Q_shg_new = zeros(0)
    mode_shg_new = String[]

    g_new = zeros(0)
    gtt_new = complex(zeros(0))
    if pgbar == "on"
        r = Progress(cutoff+1, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
        r.desc = "Computing...    "
    end
    for i = 1: cutoff
        slice = df[i, :]
        
        lambda1 = slice.wavelength_f
        n_num1 = slice.n_f
        l_num1 = slice.l_f
        Q_f = slice.Q_f
        mode1 = slice.mode_f

        lambda2 = slice.wavelength_shg
        n_num2 = slice.n_shg
        l_num2 = slice.l_shg
        Q_shg = slice.Q_shg
        mode2 = slice.mode_shg

        ratio = slice.ratio_g

        for m_num1 = l_num1:-1:l_num1-delta_f+1
            for m_num2 = l_num2:-1:l_num2-delta_shg+1
                if m_num2-3 <= 2*m_num1 <= m_num2+3 && 1 <= delta_f <= 2*l_num1 && abs(l_num2-delta_shg+1) <= l_num2
                    field_parameters_nonlinearity = [[lambda1, l_num1, m_num1, mode1], [lambda2, l_num2, m_num2, mode2]]
                    g, contribution = overlap_nonlinearity(field_parameters_nonlinearity, n, R, d; digit=1, category="sweep")
                    
                    append!(ratio_new, ratio)
                    append!(g_new, abs(g))
                    append!(gtt_new, g*ratio)

                    append!(n_f_new, n_num1)
                    append!(l_f_new, l_num1)
                    append!(m_f_new, m_num1)
                    append!(wav_f_new, lambda1)
                    append!(Q_f_new, Q_f)
                    append!(mode_f_new, [mode1])

                    append!(n_shg_new, n_num2)
                    append!(l_shg_new, l_num2)
                    append!(m_shg_new, m_num2)
                    append!(wav_shg_new, lambda2)
                    append!(Q_shg_new, Q_shg)
                    append!(mode_shg_new, [mode2])
                end
            end
        end
        if pgbar == "on"
            sleep(0.02)
            ProgressMeter.next!(r)
        end
    end
    if pgbar == "on"
        r.desc = "Finished ???      "
        ProgressMeter.next!(r)
    end 
    df = DataFrame(g = g_new, g_total = gtt_new, ratio_g = ratio_new, n_f = Int.(n_f_new), l_f = Int.(l_f_new), m_f = Int.(m_f_new), mode_f = mode_f_new, wavelength_f = wav_f_new, 
         Q_f = Q_f_new, n_shg = Int.(n_shg_new), l_shg = Int.(l_shg_new), m_shg = Int.(m_shg_new), mode_shg = mode_shg_new, wavelength_shg = wav_shg_new, Q_shg = Q_shg_new)
    return sort!(df, [:g], rev=true)
end

function overlap_nonlinearity(field_parameters, n, R, d; digit=3, category = "")
    error = 10.0^(-floor(Int, digit/2)-1)
    rtol = 10.0^(-digit)
    contribution = complex(zeros(0))
    lambda1, l_num1, m_num1, mode1 = field_parameters[1]
    lambda2, l_num2, m_num2, mode2 = field_parameters[2]
    region1_r = integral_region(lambda1, l_num1, m_num1, n, R, mode1, "E_r", error=error)
    region1_theta = integral_region(lambda1, l_num1, m_num1, n, R, mode1, "E_theta", error=error)
    region1_phi = integral_region(lambda1, l_num1, m_num1, n, R, mode1, "E_phi", error=error)
    region2_r = integral_region(lambda2, l_num2, m_num2, n, R, mode2, "E_r", error=error)
    region2_theta = integral_region(lambda2, l_num2, m_num2, n, R, mode2, "E_theta", error=error)
    region2_phi = integral_region(lambda2, l_num2, m_num2, n, R, mode2, "E_phi", error=error)
    region = [region1_r, region1_theta, region1_phi, region2_r, region2_theta, region2_phi]
    coupling = ["rrr", "rrt", "rrp", "rtr", "rtt", "rtp", "rpr", "rpt", "rpp", "ttr", "ttt", "ttp", "tpr", "tpt", "tpp", "ppr", "ppt", "ppp"]
    field_type = ["E_r", "E_theta", "E_phi"]
    if category != "sweep"
        p = Progress(20, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        p.desc = "Computing...    "
    end 
    for (i, coupling_types) in enumerate(coupling)
        k1 = floor(Int, floor((i-1)/9)+1+floor((i-1)/15))
        k2 = floor(Int, floor((i-1)/3)%3+1+floor((i-1)/9)-floor((i-1)/15))
        k3 = floor(Int, (i-1)%3+4)
        k4 = floor(Int, (i-1)%3+1)
        integral = overlap_nonlinearity_calculation(field_parameters, n, R, d, region[k1], region[k2], region[k3], field_type[k1], field_type[k2], field_type[k4], coupling_types, rtol)
        append!(contribution, integral)
        if category != "sweep"
            ProgressMeter.next!(p)
        end
    end
    coefficient = 12*pi^2*c*sqrt(6.62607004*1e-34*c/(2*8.8541878*1e-12))*1e9*sqrt(1e9)/(lambda1*sqrt(lambda2)*n^3)
    G1 = overlap_field([vcat(field_parameters[1], ["E"]), vcat(field_parameters[1], ["E"])], n, R)
    if category != "sweep"
        ProgressMeter.next!(p)
    end
    G2 = overlap_field([vcat(field_parameters[2], ["E"]), vcat(field_parameters[2], ["E"])], n, R)
    if category != "sweep"
        p.desc = "Finished ???      "
        ProgressMeter.next!(p)
    end
    gnl_complex = 1e-12*sum(contribution)*coefficient/(G1*sqrt(abs(G2)))/(2*pi) # 2pi for unit Hz
    # println((G1*sqrt(abs(G2))))
    if real(gnl_complex)==0 || imag(gnl_complex)==0
        gnl = abs(gnl_complex)
        if category != "sweep"
            println("gnl = $gnl")
        end
    else
        gnl = gnl_complex
        if category != "sweep"  
            println("Waring! gnl_complex = $gnl_complex")
        end
    end
    return gnl, contribution
end

function overlap_nonlinearity_calculation(field_parameters, n, R, d, region1, region2, region3, field_type1, field_type2, field_type3, coupling_types, rtol)
    lambda1, l_num1, m_num1, mode1 = field_parameters[1]
    lambda2, l_num2, m_num2, mode2 = field_parameters[2]
    if region1[1] !== nothing && region1[2] !== nothing && region2[1] !== nothing && region2[2] !== nothing && region3[1] !== nothing && region3[2] !== nothing
        mregion = mutual_region(region3, mutual_region(region1, region2))
        region = [mregion[1], mregion[2], R, pi/2]
        f1(p) = field(p[1], p[2], lambda1, l_num1, m_num1, n, R, mode1, field_type1)
        f2(p) = field(p[1], p[2], lambda1, l_num1, m_num1, n, R, mode1, field_type2)
        f3(p) = conj(field(p[1], p[2], lambda2, l_num2, m_num2, n, R, mode2, field_type3))
        f4(p) = p[1]^2*f1(p)*f2(p)*f3(p)*theta_depends(p[2], m_num1, m_num2, d, coupling_types)
        integration = abs(2*m_num1-m_num2) <= 3 ? hcubature(f4, region[1:2], region[3:4], rtol=rtol)[1]*4 : 0
        return integration*1e-18
    else
        return 0
    end
end

function d_rotation(d, thetax, thetay, thetaz)
    r11 = cos(thetay)^2*cos(thetaz)^2
    r12 = (sin(thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(thetax))^2
    r13 = (-sin(thetax)*sin(thetaz) + sin(thetay)*cos(thetax)*cos(thetaz))^2
    r14 = sin(thetax)^2*sin(thetay)*sin(2*thetaz)/2 - sin(thetax)*sin(thetay)^2*cos(thetax)*cos(thetaz)^2 + sin(2*thetax)*sin(thetaz)^2/2 - sin(thetay)*sin(thetaz)*cos(thetax)^2*cos(thetaz)
    r15 = (sin(thetax)*sin(thetaz) - sin(thetay)*cos(thetax)*cos(thetaz))*cos(thetay)*cos(thetaz)
    r16 = (sin(thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(thetax))*cos(thetay)*cos(thetaz)
    r21 = sin(thetaz)^2*cos(thetay)^2
    r22 = (-sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(thetaz))^2
    r23 = (sin(thetax)*cos(thetaz) + sin(thetay)*sin(thetaz)*cos(thetax))^2
    r24 = (-sin(thetay)^2*sin(thetaz)^2 + cos(thetaz)^2)*sin(thetax)*cos(thetax) - sin(thetax)^2*sin(thetay)*sin(thetaz)*cos(thetaz) + sin(thetay)*sin(thetaz)*cos(thetax)^2*cos(thetaz)
    r25 = -(sin(thetax)*cos(thetaz) + sin(thetay)*sin(thetaz)*cos(thetax))*sin(thetaz)*cos(thetay)
    r26 = (sin(thetax)*sin(thetay)*sin(thetaz) - cos(thetax)*cos(thetaz))*sin(thetaz)*cos(thetay)
    r31 = sin(thetay)^2
    r32 = sin(thetax)^2*cos(thetay)^2
    r33 = cos(thetax)^2*cos(thetay)^2
    r34 = -sin(thetax)*cos(thetax)*cos(thetay)^2
    r35 = sin(2*thetay)*cos(thetax)/2
    r36 = -sin(thetax)*sin(thetay)*cos(thetay)
    r41 = -2*sin(thetay)*sin(thetaz)*cos(thetay)
    r42 = 2*(sin(thetax)*sin(thetay)*sin(thetaz) - cos(thetax)*cos(thetaz))*sin(thetax)*cos(thetay)
    r43 = 2*(sin(thetax)*cos(thetaz) + sin(thetay)*sin(thetaz)*cos(thetax))*cos(thetax)*cos(thetay)
    r44 = (-sin(thetax)^2*cos(thetaz) - 2*sin(thetax)*sin(thetay)*sin(thetaz)*cos(thetax) + cos(thetax)^2*cos(thetaz))*cos(thetay)
    r45 = sin(thetax)*sin(thetay)*cos(thetaz) - sin(thetaz)*cos(thetax)*cos(2*thetay)
    r46 = sin(thetax)*sin(thetaz)*cos(2*thetay) + sin(thetay)*cos(thetax)*cos(thetaz)
    r51 = sin(2*thetay)*cos(thetaz)
    r52 = -2*(sin(thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(thetax))*sin(thetax)*cos(thetay)
    r53 = 2*(sin(thetax)*sin(thetaz) - sin(thetay)*cos(thetax)*cos(thetaz))*cos(thetax)*cos(thetay)
    r54 = (sin(2*thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(2*thetax))*cos(thetay)
    r55 = sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(2*thetay)*cos(thetaz)
    r56 = (sin(thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(thetax))*sin(thetay) - sin(thetax)*cos(thetay)^2*cos(thetaz)
    r61 = -2*sin(thetaz)*cos(thetay)^2*cos(thetaz)
    r62 = -2*sin(thetax)^2*sin(thetay)^2*sin(thetaz)*cos(thetaz) + 2*sin(thetax)*sin(thetay)*cos(thetax)*cos(2*thetaz) + 2*sin(thetaz)*cos(thetax)^2*cos(thetaz)
    r63 = sin(thetax)^2*sin(2*thetaz) - 2*sin(thetax)*sin(thetay)*cos(thetax)*cos(thetaz)^2 + sin(2*thetax)*sin(thetay)*sin(thetaz)^2 - 2*sin(thetay)^2*sin(thetaz)*cos(thetax)^2*cos(thetaz)
    r64 = -(cos(2*thetay)/2 - 3/2)*sin(2*thetax)*sin(2*thetaz)/2 + sin(thetax)^2*sin(thetay)*cos(2*thetaz) - sin(thetay)*cos(thetax)^2*cos(2*thetaz)
    r65 = (-sin(thetax)*sin(thetaz)^2 + sin(thetax)*cos(thetaz)^2 + sin(thetay)*sin(2*thetaz)*cos(thetax))*cos(thetay)
    r66 = (-2*sin(thetax)*sin(thetay)*sin(thetaz)*cos(thetaz) + cos(thetax)*cos(2*thetaz))*cos(thetay)
    Rb = [[r11, r12, r13, r14, r15, r16], [r21, r22, r23, r24, r25, r26], [r31, r32, r33, r34, r35, r36], [r41, r42, r43, r44, r45, r46], [r51, r52, r53, r54, r55, r56], [r61, r62, r63, r64, r65, r66]]
    Ra = [[cos(thetay)*cos(thetaz), -sin(thetaz)*cos(thetay), sin(thetay)], 
         [sin(thetax)*sin(thetay)*cos(thetaz) + sin(thetaz)*cos(thetax), -sin(thetax)*sin(thetay)*sin(thetaz) + cos(thetax)*cos(thetaz), -sin(thetax)*cos(thetay)], 
         [sin(thetax)*sin(thetaz) - sin(thetay)*cos(thetax)*cos(thetaz), sin(thetax)*cos(thetaz) + sin(thetay)*sin(thetaz)*cos(thetax), cos(thetax)*cos(thetay)]]
    Ra = hcat(Ra...)'
    Rb = hcat(Rb...)'
    d = hcat(d...)'
    d = Ra*d*Rb
    d = [[d[1, 1], d[1, 2], d[1, 3], d[1, 4], d[1, 5], d[1, 6]], [d[2, 1], d[2, 2], d[2, 3], d[2, 4], d[2, 5], d[2, 6]], [d[3, 1], d[3, 2], d[3, 3], d[3, 4], d[3, 5], d[3, 6]]]
    return d
end

function shade(f, x0; epsilon=1e-3, steps = 30)
    k = pi/10
    delta = 2
    maxi = 0
    while delta>epsilon && steps > 0
        a1 = f(x0[1], x0[2], x0[3])
        a = [f(x0[1]+k, x0[2], x0[3]), f(x0[1]-k, x0[2], x0[3]), f(x0[1], x0[2]+k, x0[3]), f(x0[1], x0[2]-k, x0[3]), f(x0[1], x0[2], x0[3]+k), f(x0[1], x0[2], x0[3]-k)]
        m = maximum(a)
        i = argmax(a)
        if m>a1
            if i==1
                x0[1]+=k
            elseif i==2
                x0[1]-=k
            elseif i==3
                x0[2]+=k
            elseif i==4
                x0[2]-=k
            elseif i==5
                x0[3]+=k
            else
                x0[3]-=k
            end
            delta = abs(m-a1)/a1
            steps-=1
            println("maximum found is $a1, at position $x0, delta=$delta, steps=$steps")
        else
            k = 0.2*k
            println(k)
        end
        maxi = a1
    end
    return maxi, x0, delta
end