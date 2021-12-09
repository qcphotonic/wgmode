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

include("modesolutions.jl")
using ProgressMeter

c = 299792458

function Irrr(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return sin(theta)*(d33*cos(theta)^3+3*d31*cos(theta)*sin(theta)^2)
    elseif 2*m_a-m_b == -3
        return 1im*d22*sin(theta)^4/2
    elseif 2*m_a-m_b == 3
        return -1im*d22*sin(theta)^4/2
    else
        return 0
    end
end

function Irrt(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return sin(theta)^2*(d31-d33+(3*d31-d33)*cos(2*theta))/2
    elseif 2*m_a-m_b == -3
        return 1im*d22*sin(theta)^3*cos(theta)/2
    elseif 2*m_a-m_b == 3
        return -1im*d22*sin(theta)^3*cos(theta)/2
    else 
        return 0
    end
end

function Irrp(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -(d22+4*d31*cos(theta))*sin(theta)^3/4
    elseif 2*m_a-m_b == -2
        return -1im*(d31*cos(theta)+d22*sin(theta)^2)*sin(theta)^2/2
    elseif 2*m_a-m_b == 2
        return 1im*(d31*cos(theta)+d22*sin(theta)^2)*sin(theta)^2/2
    elseif abs(2*m_a-m_b) == 3
        return -d22*sin(theta)^3/4
    else
        return 0
    end
end

function Irtr(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return sin(theta)^2*(d31-d33+(3*d31-d33)*cos(2*theta))
    elseif 2*m_a-m_b == -3
        return 1im*d22*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == 3
        return -1im*d22*sin(theta)^3*cos(theta)
    else
        return 0
    end
end

function Irtt(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return sin(2*theta)*(d33-d31+(3*d31-d33)*cos(2*theta))/2
    elseif 2*m_a-m_b == -3
        return 1im*d22*sin(2*theta)^2/4
    elseif 2*m_a-m_b == 3
        return -1im*d22*sin(2*theta)^2/4
    else
        return 0
    end
end

function Irtp(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -(d22*cos(theta)+2*d31*cos(2*theta))*sin(theta)^2/2
    elseif 2*m_a-m_b == -2
        return -1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif 2*m_a-m_b == 2
        return 1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif abs(2*m_a-m_b) == 3
        return -d22*cos(theta)*sin(theta)^2/2
    else 
        return 0
    end
end

function Irpr(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -sin(theta)^3*(d22+4*d31*cos(theta))/2
    elseif 2*m_a-m_b == -2
        return -1im*sin(theta)^2*(d31*cos(theta)+d22*sin(theta)^2)
    elseif 2*m_a-m_b == 2
        return 1im*sin(theta)^2*(d31*cos(theta)+d22*sin(theta)^2)
    elseif abs(2*m_a-m_b) == 3
        return -d22*sin(theta)^3/2
    else
        return 0
    end
end

function Irpt(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -sin(theta)^2*(d22*cos(theta)+2*d31*cos(2*theta))/2
    elseif 2*m_a-m_b == -2
        return -1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif 2*m_a-m_b == 2
        return 1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif abs(2*m_a-m_b) == 3
        return -d22*sin(theta)*sin(2*theta)/4
    else
        return 0
    end
end

function Irpp(theta, m_a, m_b, d)
    d22, d31, = d
    if 2*m_a-m_b == 0
        return sin(theta)*(d31*(3*cos(theta)-cos(3*theta))+4*d22*sin(theta)^2)/2
    elseif 2*m_a-m_b == -1
        return -1im*d22*(2*cos(2*theta)-1)*sin(theta)^2/4
    elseif 2*m_a-m_b == 1
        return 1im*d22*(2*cos(2*theta)-1)*sin(theta)^2/4
    elseif abs(2*m_a-m_b) == 2
        return sin(theta)*(d31*cos(theta)+2*d22*sin(theta)^2)/2
    elseif 2*m_a-m_b == -3
        return -1im*d22*sin(theta)^2/4
    elseif 2*m_a-m_b == 3
        return 1im*d22*sin(theta)^2/4
    else
        return 0
    end
end

function Ittr(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return (d33-d31+(3*d31-d33)*cos(2*theta))*sin(2*theta)/4
    elseif 2*m_a-m_b == -3
        return 1im*d22*sin(2*theta)^2/8
    elseif 2*m_a-m_b == 3
        return -1im*d22*sin(2*theta)^2/8
    else
        return 0
    end
end

function Ittt(theta, m_a, m_b, d)
    d22, d31, d33 = d
    if 2*m_a-m_b == 0
        return -(3*d31+d33+(3*d31-d33)*cos(2*theta))*sin(theta)^2/2
    elseif 2*m_a-m_b == -3
        return 1im*d22*cos(theta)^3*sin(theta)/2
    elseif 2*m_a-m_b == 3
        return -1im*d22*cos(theta)^3*sin(theta)/2
    else
        return 0
    end
end

function Ittp(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -(d22*cos(theta)-4*d31*sin(theta)^2)*sin(2*theta)/8
    elseif 2*m_a-m_b == -2
        return 1im*cos(theta)*(d31-d22*cos(theta))*sin(theta)^2/2
    elseif 2*m_a-m_b == 2
        return -1im*cos(theta)*(d31-d22*cos(theta))*sin(theta)^2/2
    elseif abs(2*m_a-m_b) == 3
        return -d22*cos(theta)^2*sin(theta)/4
    else
        return 0
    end
end

function Itpr(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -(d22*cos(theta)+2*d31*cos(2*theta))*sin(theta)^2/2
    elseif 2*m_a-m_b == -2
        return -1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif 2*m_a-m_b == 2
        return 1im*sin(theta)*(d31*cos(2*theta)+d22*sin(theta)*sin(2*theta))/2
    elseif abs(2*m_a-m_b) == 3
        return -d22*cos(theta)*sin(theta)^2/2
    else
        return 0
    end
end

function Itpt(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 1
        return -sin(2*theta)*(d22*cos(theta)-4*d31*sin(theta)^2)/4
    elseif 2*m_a-m_b == -2
        return 1im*cos(theta)*(d31-d22*cos(theta))*sin(theta)^2
    elseif 2*m_a-m_b == 2
        return -1im*cos(theta)*(d31-d22*cos(theta))*sin(theta)^2
    elseif abs(2*m_a-m_b) == 3
        return -d22*cos(theta)*sin(2*theta)/4
    else
        return 0
    end
end

function Itpp(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 0
        return sin(theta)^2*(2*d22*cos(theta)+d31*(cos(2*theta)-2))
    elseif 2*m_a-m_b == -1
        return -1im*d22*cos(3*theta)*sin(theta)/4
    elseif 2*m_a-m_b == 1
        return 1im*d22*cos(3*theta)*sin(theta)/4
    elseif abs(2*m_a-m_b) == 2
        return -sin(theta)^2*(d31-2*d22*cos(theta))/2
    elseif 2*m_a-m_b == -3
        return -1im*d22*sin(theta)*cos(theta)/4
    elseif 2*m_a-m_b == 3
        return 1im*d22*sin(theta)*cos(theta)/4
    else
        return 0
    end
end

function Ippr(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 0
        return sin(theta)*(d31*(3*cos(theta)-cos(3*theta))+4*d22*sin(theta)^2)/4
    elseif 2*m_a-m_b == -1
        return -1im*d22*sin(theta)^2*(2*cos(2*theta)-1)/8
    elseif 2*m_a-m_b == 1
        return 1im*d22*sin(theta)^2*(2*cos(2*theta)-1)/8
    elseif abs(2*m_a-m_b) == 2
        return sin(theta)*(d31*cos(theta)+2*d22*sin(theta)^2)/4
    elseif 2*m_a-m_b == -3
        return -1im*d22*sin(theta)^2/8
    elseif 2*m_a-m_b == 3
        return 1im*d22*sin(theta)^2/8
    else
        return 0
    end
end

function Ippt(theta, m_a, m_b, d)
    d22, d31, = d
    if abs(2*m_a-m_b) == 0
        return sin(theta)^2*(2*d22*cos(theta)+d31*(cos(2*theta)-2))/2
    elseif 2*m_a-m_b == -1
        return 1im*d22*(sin(2*theta)-sin(4*theta))/16
    elseif 2*m_a-m_b == 1
        return -1im*d22*(sin(2*theta)-sin(4*theta))/16
    elseif abs(2*m_a-m_b) == 2
        return -sin(theta)^2*(d31-2*d22*cos(theta))/4
    elseif 2*m_a-m_b == -3
        return -1im*d22*sin(theta)*cos(theta)/8
    elseif 2*m_a-m_b == 3
        return 1im*d22*sin(theta)*cos(theta)/8
    else
        return 0
    end
end

function Ippp(theta, m_a, m_b, d)
    d22, = d
    if abs(2*m_a-m_b) == 1
        return 3*d22*(sin(3*theta)-2*sin(theta))/8
    elseif abs(2*m_a-m_b) == 3
        return d22*sin(theta)/8
    else
        return 0
    end
end

function theta_depends(theta, m_a, m_b, d, coupling_types)
    if coupling_types == "rrr"
        return Irrr(theta, m_a, m_b, d)
    elseif coupling_types == "rrt"
        return Irrt(theta, m_a, m_b, d)
    elseif coupling_types == "rrp"
        return Irrp(theta, m_a, m_b, d)
    elseif coupling_types == "rtr"
        return Irtr(theta, m_a, m_b, d)
    elseif coupling_types == "rtt"
        return Irtt(theta, m_a, m_b, d)
    elseif coupling_types == "rtp"
        return Irtp(theta, m_a, m_b, d)
    elseif coupling_types == "rpr"
        return Irpr(theta, m_a, m_b, d)
    elseif coupling_types == "rpt"
        return Irpt(theta, m_a, m_b, d)
    elseif coupling_types == "rpp"
        return Irpp(theta, m_a, m_b, d)
    elseif coupling_types == "ttr"
        return Ittr(theta, m_a, m_b, d)
    elseif coupling_types == "ttt"
        return Ittt(theta, m_a, m_b, d)
    elseif coupling_types == "ttp"
        return Ittp(theta, m_a, m_b, d)
    elseif coupling_types == "tpr"
        return Itpr(theta, m_a, m_b, d)
    elseif coupling_types == "tpt"
        return Itpt(theta, m_a, m_b, d)
    elseif coupling_types == "tpp"
        return Itpp(theta, m_a, m_b, d)
    elseif coupling_types == "ppr"
        return Ippr(theta, m_a, m_b, d)
    elseif coupling_types == "ppt"
        return Ippt(theta, m_a, m_b, d)
    elseif coupling_types == "ppp"
        return Ippp(theta, m_a, m_b, d)
    else
        println("Error: No $coupling_types found.")
        return 0
    end
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
        p.desc = "Finished ✓      "
        ProgressMeter.next!(p)
    end
    gnl_complex = 1e-12*sum(contribution)*coefficient/(G1*sqrt(abs(G2)))
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


    

