#= Calculate the field distribution for each mode selected. View the field in the r-theta plane.
    
field(r, theta, lambda, l_num, m_num, n, R, mode, field; phi=0)
    r::float64
        radial parameter
    theta::float64
        polar parameter
    phi::float64
        azimuthal parameter
    lambda::float64
        wavelength
    l_num::Int 
        l mode number
    m_num::Int 
        m mode number
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    mode::String ("TE" or "TM")
        mode type
    field::String ("E_r", "E_theta", "E_phi", "E", "E^2", "H_r", "H_theta", "H_phi", "H", "H^2",
                    "nE_r", "nE_theta", "nE_phi", "nH_r", "nH_theta", "nH_phi")
        field type
    
    return float64

    Given all information above, return unnormalized field value at (r, theta, phi) for the chosen field type.


view_field(data, n_num, l_num, m_num, n, R, mode, field_tp; quality="coarse", scale="normal", R_region=R, half_angle="no")
    data{::DataFrames, ::float64, 1}
        spectrum data or radial filled ratio (0<data<=1)
    n_num::Int 
        n mode number
    l_num::Int 
        l mode number
    m_num::Int 
        l mode number
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    mode::String ("TE" or "TM")
        mode type
    field_tp::String ("E_r", "E_theta", "E_phi", "E", "E^2", "H_r", "H_theta", "H_phi", "H", "H^2",
                    "nE_r", "nE_theta", "nE_phi", "nH_r", "nH_theta", "nH_phi")
        field type
    qualtiy::String ("coarse" or "fine")
        image quality
    scale::String ("normal" or "log")
        colarbar scaling
    R_region::float64
        radius range of viewing, R_region could be less than R, equal to R, or larger than R.
    half_angle::String ("yes" or "no")
        determine whether return the field distribution in the theta range or not
    
    return {plot, ?::float64}
    
    Given all information above, viewing the field distribution in r-theta plane.


overlap_field(field_parameters, n, R; R_region=R, digit=3)
    field_parameters::Vector{Vector{Any}}
        selected mode number, mode type and field type for two fields
    n::float64
        refractive index
    R::float64
        radius of the dielectric sphere
    R_region::float64
        radius range of viewing, R_region could be less than R, equal to R, or larger than R.
    digit::Int
        effective number of the result

    Given the info of selected two fields and structural parameters, calculate the field overlapping.
=#

using VectorSphericalHarmonics
using SpecialFunctions
using LinearAlgebra
using GR
using HCubature

using ProgressMeter
using DataFramesMeta

c = 299792458


function sphericalbesselj_d(l, x)
    return sqrt(pi*x/2)*(besselj(l+1/2, x)/(2*x)+besselj(l-1/2, x)/2-besselj(l+3/2, x)/2)
end

function sphericalbesselh_d(l, x)
    return sqrt(pi*x/2)*(besselh(l+1/2, 2, x)/(2*x)+besselh(l-1/2, 2, x)/2-besselh(l+3/2, 2, x)/2)
end

function E_TE(r, theta, phi, lambda, l, m, n, R)
    xr0 = 2*pi*r*1e3/lambda
    xr = xr0*n
    xR0 = 2*pi*R*1e3/lambda
    xR = xR0*n
    if r<=R
        fr = sqrt(pi/(2*xr))*besselj(l+1/2, xr)
    else
        fr = sqrt(pi/(2*xr))*besselj(l+1/2, xR)*besselh(l+1/2, 2, xr0)/besselh(l+1/2, 2, xR0)
    end
    fo = -1im*r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 0, theta, phi)
    return fr*fo
end

function H_TE(r, theta, phi, lambda, l, m, n, R)
    xr0 = 2*pi*r*1e3/lambda
    xr = xr0*n
    xR0 = 2*pi*R*1e3/lambda
    xR = xR0*n
    if r<=R
        f1 = l*(l+1)*sqrt(pi*xr/2)*besselj(l+1/2, xr)*vshbasis(Hansen(), Polar(), l, m, -1, theta, phi)/xr^2
        f2 = sphericalbesselj_d(l, xr)*r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 1, theta, phi)/xr
        f_tol = (1im*n/c)*(f1+f2)
    else
        f1 = l*(l+1)*sqrt(pi*xr0/2)*besselh(l+1/2, 2, xr0)*vshbasis(Hansen(), Polar(), l, m, -1, theta, phi)/xr0^2
        f2 = sphericalbesselh_d(l, xr0)*r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 1, theta, phi)/xr0
        f_tol = (1im/(n*c))*sqrt(n)*(besselj(l+1/2, xR)/besselh(l+1/2, 2, xR0))*(f1+f2)
    end
    return f_tol
end

function E_TM(r, theta, phi, lambda, l, m, n, R)
    xr0 = 2*pi*r*1e3/lambda
    xr = xr0*n
    xR0 = 2*pi*R*1e3/lambda
    xR = xR0*n
    if r<=R
        f1 = l*(l+1)*sqrt(pi*xr/2)*besselj(l+1/2, xr)*vshbasis(Hansen(), Polar(), l, m, -1, theta, phi)/xr^2
        f2 = sphericalbesselj_d(l, xr)*r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 1, theta, phi)/xr
        f_tol = f1+f2
    else
        f1 = l*(l+1)*sqrt(pi*xr0/2)*besselh(l+1/2, 2, xr0)*vshbasis(Hansen(), Polar(), l, m, -1, theta, phi)/xr0^2
        f2 = sphericalbesselh_d(l, xr0)*r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 1, theta, phi)/xr0
        f_tol = sqrt(n)*(besselj(l+1/2, xR)/besselh(l+1/2, 2, xR0))*(f1+f2)
    end
    return f_tol
end

function H_TM(r, theta, phi, lambda, l, m, n, R)
    xr0 = 2*pi*r*1e3/lambda
    xr = xr0*n
    xR0 = 2*pi*R*1e3/lambda
    xR = xR0*n
    if r<=R
        fr = sqrt(pi/(2*xr))*besselj(l+1/2, xr)
    else
        fr = sqrt(pi/(2*xr))*besselj(l+1/2, xR)*besselh(l+1/2, 2, xr0)/besselh(l+1/2, 2, xR0)
    end
    fo = r*sqrt(l*(l+1))*vshbasis(Hansen(), Polar(), l, m, 0, theta, phi)
    return n*fr*fo/c
end

# field selector
function field(r, theta, lambda, l_num, m_num, n, R, mode, field; phi=0)
    field_type = ["E_r", "E_theta", "E_phi", "E", "E^2", "H_r", "H_theta", "H_phi", "H", "H^2",
                    "nE_r", "nE_theta", "nE_phi", "nH_r", "nH_theta", "nH_phi"]
    type = findfirst(x -> x==field, field_type)
    if mode == "TE"
        if 1 <= type <= 3
            return E_TE(r, theta, phi, lambda, l_num, m_num, n, R)[type]
        elseif type == 4
            return norm(E_TE(r, theta, phi, lambda, l_num, m_num, n, R))
        elseif type == 5
            return norm(E_TE(r, theta, phi, lambda, l_num, m_num, n, R))^2
        elseif 6 <= type <= 8
            return H_TE(r, theta, phi, lambda, l_num, m_num, n, R)[type-5]
        elseif type == 9
            return norm(H_TE(r, theta, phi, lambda, l_num, m_num, n, R))
        elseif type == 10
            return norm(H_TE(r, theta, phi, lambda, l_num, m_num, n, R))^2
        elseif 11 <= type <= 13
            return norm(E_TE(r, theta, phi, lambda, l_num, m_num, n, R)[type-10])
        elseif 14 <= type <= 16
            return norm(H_TE(r, theta, phi, lambda, l_num, m_num, n, R)[type-13])
        else
            println("Error: field_type of $filed doesn't exist")
            return nothing
        end
    elseif mode == "TM"
        if 1 <= type <= 3
            return E_TM(r, theta, phi, lambda, l_num, m_num, n, R)[type]
        elseif type == 4
            return norm(E_TM(r, theta, phi, lambda, l_num, m_num, n, R))
        elseif type == 5
            return norm(E_TM(r, theta, phi, lambda, l_num, m_num, n, R))^2
        elseif 6 <= type <= 8
            return H_TM(r, theta, phi, lambda, l_num, m_num, n, R)[type-5]
        elseif type == 9
            return norm(H_TM(r, theta, phi, lambda, l_num, m_num, n, R))
        elseif type == 10
            return norm(H_TM(r, theta, phi, lambda, l_num, m_num, n, R))^2
        elseif 11 <= type <= 13
            return norm(E_TM(r, theta, phi, lambda, l_num, m_num, n, R)[type-10])
        elseif 14 <= type <= 16
            return norm(H_TM(r, theta, phi, lambda, l_num, m_num, n, R)[type-13])
        else
            println("Error: field_type of $filed doesn't exist")
            return nothing
        end
    else
        println("Error: mode of $mode doesn't exist")
        return nothing
    end
end

# fill the whole plane with objective field
function whole(θ, polardata)
    nq = length(θ)
    θ = collect(θ)
    for (i, v) in enumerate(reverse(θ)[2:nq])
        append!(θ, pi-v)
        polardata = hcat(polardata, polardata[:, nq-i])
    end
    for (i, v) in enumerate(reverse(θ)[2:2*nq-1])
        append!(θ, 2*pi-v)
        polardata = hcat(polardata, polardata[:, 2*nq-1-i])
    end
    return θ, polardata
end

# sampling array
function radius_array(len, R, R_region, n_num, n_max)
    if R_region <= R
        r = sort([density(x, R, R_region, n_num, n_max) for x = range(0, stop=R_region, length=len)])
    else
        l1 = floor(Int, len*R/R_region)
        range_r1 = range(0, stop=R, length=l1)
        p = (R_region-R)/(len-l1)
        range_r2 = R+p:p:R_region
        range_r = vcat(range_r1, range_r2)
        r = sort([density(x, R, R_region, n_num, n_max) for x = range_r])
    end
    return r
end

# sampling function
function density(x, R, R_region, n_num, n_max)
    ratio = 2-n_num/n_max
    if R_region == R
        return (R^ratio-x^ratio)^(1/ratio)
    elseif R_region > R
        return (x <= R) ? (R^ratio-x^ratio)^(1/ratio) : R+(x-R)^2/(R_region-R)
    else
        return (R_region^ratio-x^ratio)^(1/ratio)
    end
end

# log view mode assistance
function log_array(pd)
    s = size(pd)
    pd_log = zeros(s)
    for i=1:s[1]
        for j=1:s[2]
            pd_log[i ,j]=log(pd[i ,j]+1)
        end
    end
    return pd_log
end
    
# display the field distribution in r-theta plane
function view_field(data, n_num, l_num, m_num, n, R, mode, field_tp; quality="coarse", scale="normal", R_region=R, half_angle="no")
    
    if typeof(data) == DataFrame
        lambda_df = @chain data begin
                         @rsubset :n == n_num && :l == l_num
                         @select(:wav = :wavelength)
        end
        lambda = lambda_df.wav[1]
        
        l_max = maximum(data.l)

        n_max = l_max*n_num/(l_max-l_num)
    elseif (typeof(data) == Float64 && 0 < data < 1) || data == 1
        n_max = floor(Int, n_num/data)
    else
        println("Error: $data is not a DataFrame nor float less than 1")
    end

    dsp = [floor(Int, 40*(1+n_num/n_max)), 100] # precision of coarse quality
    θ = range(0, stop=pi/2, length=dsp[2])
    r = radius_array(dsp[1], R, R_region, n_num, n_max)
    
    polardata = zeros((dsp[1], dsp[2]))
    for i=1:length(r)
        for j=1:length(θ)
            pd = field(r[i], θ[j], lambda, l_num, m_num, n, R, mode, field_tp)
            polardata[i, j] = isnan(pd) ? 0 : pd
        end
    end

    if quality == "fine"
        m = findfirst(x -> x > maximum(polardata)/1e4, polardata)[2]
        dspm = [floor(Int, 400*(1+n_num/n_max)), m] # precision of fine quality
        θ = range(0, stop=pi*(m-1)/(2*(dsp[2]-1)), length=m)
        r_fine = radius_array(dspm[1], R, R_region, n_num, n_max)
        polardata1 = zeros((dspm[1], m))

        p = ProgressMeter.Progress(dspm[1]+1, 0.01, "Initializing... ")
        for i=1:length(r_fine)
            for j=1:length(θ)
                pd = field(r_fine[i], θ[j], lambda, l_num, m_num, n, R, mode, field_tp)
                polardata1[i, j] = isnan(pd) ? 0 : pd
            end
            sleep(0.001)
            ProgressMeter.next!(p)
        end
        p.desc = "Computing...    "
        p.counter=p.prev_update_count=0
        s = pi*(m-1)/(2*(dsp[2]-1)):0.001:pi/2
        polardata2= zeros((dspm[1], length(s)))
        for i=1:length(r_fine)
            for j=1:length(s)
                pd = field(r_fine[i], s[j], lambda, l_num, m_num, n, R, mode, field_tp)
                polardata2[i, j] = isnan(pd) ? 0 : pd
            end
            ProgressMeter.next!(p)
        end
        for i=1:length(r_fine)
            c = vcat(polardata1[i,:][1:m-1], polardata2[i,:])
            if i == 1
                polardata_fine = c
            else
                polardata_fine = hcat(polardata_fine, c)
            end
        end
        polardata_fine = transpose(polardata_fine)
        θ_fine = vcat(θ[1:m-1], s)
        θ_fine, polardata_fine = whole(θ_fine, polardata_fine)
        p.desc = "Finished ✓      "
        ProgressMeter.next!(p)
        if scale == "normal"
            pic = nonuniformpolarheatmap(θ_fine,r_fine,polardata_fine)
        elseif scale=="log"
            pic = nonuniformpolarheatmap(θ_fine,r_fine,log_array(polardata_fine))
        else
            println("Error: scale no option of $scale")
        end
        display(pic)
    else
        θ, polardata = whole(θ, polardata)
        if scale == "normal"
            pic = nonuniformpolarheatmap(θ,r,polardata)
        elseif scale=="log"
            pic = nonuniformpolarheatmap(θ,r,log_array(polardata))
        else
            println("Error: scale no option of $scale")
        end
        display(pic)
    end
    
    # ? 1e4 controversial measure
    if half_angle == "yes"
        m = findfirst(x -> x > maximum(polardata)/1e4, polardata)[2]
        return (dsp[2]-m)*pi/(2*(dsp[2]-1))
    else
        return nothing
    end
end

# find the effective integration region
function integral_region(lambda,l_num, m_num, n, R, mode, field_type; error=1e-2)
    r_border = first_zeros_bsolver(r->norm(field(r, pi/2, lambda, l_num, m_num, n, R, mode, field_type))-error, [0, R], accuracy=5)
    if m_num>2
        theta_border = first_zeros_bsolver(theta->norm(field(50, theta, lambda, l_num, m_num, n, R, mode, field_type))-error, [0, pi/2], accuracy=0.5)
    else
        theta_border = 0
    end
    if r_border === nothing && theta_border !== nothing
        r_border = first_zeros_bsolver(r->norm(field(r, pi/2, lambda, l_num, m_num, n, R, mode, "E"))-error, [0, R], accuracy=5)
    end
    return r_border, theta_border
end

function mutual_region(b1, b2)
    return maximum([b1[1], b2[1]]), maximum([b1[2], b2[2]])
end

# Calculate two fields overlaping
function overlap_field(field_parameters, n, R; R_region=R, digit=3)
    error = 10.0^(-floor(Int, digit/2)-1)
    rtol = 10.0^(-digit)
    lambda1, l_num1, m_num1, mode1, field_type1 = field_parameters[1]
    lambda2, l_num2, m_num2, mode2, field_type2 = field_parameters[2]
    region1 = integral_region(lambda1,l_num1, m_num1, n, R, mode1, field_type1, error=error)
    region2 = integral_region(lambda2,l_num2, m_num2, n, R, mode2, field_type2, error=error)
    if region1[1] !== nothing && region1[2] !== nothing && region2[1] !== nothing && region2[2] !== nothing
        mregion = mutual_region(region1, region2)
        region = [mregion[1], mregion[2], R_region, pi/2]
        f1(p) = conj(field(p[1], p[2], lambda1, l_num1, m_num1, n, R, mode1, field_type1))
        f2(p) = field(p[1], p[2], lambda2, l_num2, m_num2, n, R, mode2, field_type2)
        f3(p) = p[1]^2*sin(p[2])*f1(p)*f2(p)
        integration = (m_num1 == m_num2 && l_num1 == l_num2) ? hcubature(f3, region[1:2], region[3:4], rtol=rtol)[1]*2*pi : 0
        return round(4*integration*1e-18, sigdigits=digit)
    else
        return 0
    end
end
