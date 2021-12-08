using DataFramesMeta
using DataFrames
using CSV

function query(n_num, l_num; option)
    error = 0
    if typeof(option)==String
        DF = CSV.File(String)
    elseif typeof(option)==DataFrame
        DF = option
    else
        println("Error! Option should be path or dataframe.")
        error = 1
    end
    if error==0
        lambda_df = @chain DF begin
            @rsubset :n == n_num && :l == l_num
            @select(:wav = :wavelength)
        end
        lambda_select = lambda_df.wav[1]
    end
    return lambda_select
end
