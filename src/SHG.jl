using DataFrames

function detuning(data_f, data_shg; threshold = "")
    n_f = data_f.n
    l_f = data_f.l
    Q_f = data_f.Q
    wav_f = data_f.wavelength
    n_shg = data_shg.n
    l_shg = data_shg.l
    Q_shg = data_shg.Q
    wav_shg = data_shg.wavelength

    wav_eff_f = zeros(0)
    wav_eff_shg = zeros(0)

    n_f_new = zeros(0)
    l_f_new = zeros(0)
    wav_f_new = zeros(0)
    Q_f_new = zeros(0)
    n_shg_new = zeros(0)
    l_shg_new = zeros(0)
    wav_shg_new = zeros(0)
    Q_shg_new = zeros(0)
    Ratio = zeros(0)

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
                append!(n_shg_new, n_shg[j])
                append!(l_shg_new, l_shg[j])
                append!(wav_shg_new, wav_shg[j])
                append!(Q_shg_new, Q_shg[j])
            end
        end
    end

    if typeof(threshold) == Float64
        return wav_eff_f, wav_eff_shg
    else
        df = DataFrame(ratio_g = Ratio, n_f = Int.(n_f_new), l_f = Int.(l_f_new), wavelength_f = wav_f_new, Q_f = Q_f_new, 
            n_shg = Int.(n_shg_new), l_shg = Int.(l_shg_new), wavelength_shg = wav_shg_new, Q_shg = Q_shg_new)
        return sort!(df, rev=true)
    end
end

