using DataFrames

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

function overlap_array(df, n, R, d, limitation, cutoff)
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
    Ratio = zeros(0)
    mode_shg_new = String[]
    g_new = complex(zeros(0))
    gtt_new = complex(zeros(0))

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
                if m_num2-3 <= 2*m_num1 <= m_num2+3
                    field_parameters_nonlinearity = [[lambda1, l_num1, m_num1, mode1], [lambda2, l_num2, m_num2, mode2]]
                    g, contribution = overlap_nonlinearity(field_parameters_nonlinearity, n, R, d; digit=1)
                    
                    append!(ratio_new, ratio)
                    append!(g_new, g)
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
    end

    df = DataFrame(g = g_new, g_total = gtt_new, ratio_g = ratio_new, n_f = Int.(n_f_new), l_f = Int.(l_f_new), m_f = Int.(m_f_new), mode_f = mode_f_new, wavelength_f = wav_f_new, 
         Q_f = Q_f_new, n_shg = Int.(n_shg_new), l_shg = Int.(l_shg_new), m_shg = Int.(m_shg_new), mode_shg = mode_shg_new, wavelength_shg = wav_shg_new, Q_shg = Q_shg_new)
    return df
end

# sort!(df, rev=true)    
