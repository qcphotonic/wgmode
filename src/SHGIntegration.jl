include("modesolutions.jl")
using ProgressMeter

c = 299792458

function Irrr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (2*d33*cos(theta)^3 + (d15 + d24)*sin(theta)*sin(2*theta) + (d31 + d32)*sin(theta)^2*cos(theta))*sin(theta)/2
    elseif 2*m_a-m_b == -1
        return (d13 - 1im*(d23 + 2*d34 + 2*1im*d35))*sin(2*theta)^2/8 + (3*d11 + d12 - 1im*(2*d16 + d21 + 3*d22 + 2*1im*d26))*sin(theta)^4/8
    elseif 2*m_a-m_b == 1
        return (d13 + 1im*d23 + 2*1im*d34 + 2*d35)*sin(2*theta)^2/8 + (3*d11 + d12 + 2*1im*d16 + 1im*d21 + 3*1im*d22 + 2*d26)*sin(theta)^4/8
    elseif 2*m_a-m_b == -2
        return (-1im*d14/2 + d15/2 - d24/2 - 1im*d25/2 + d31/4 - d32/4 - 1im*d36/2)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == 2
        return (1im*d14/2 + d15/2 - d24/2 + 1im*d25/2 + d31/4 - d32/4 + 1im*d36/2)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -3
        return (d11/8 - d12/8 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/8)*sin(theta)^4
    elseif 2*m_a-m_b == 3
        return (d11/8 - d12/8 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/8)*sin(theta)^4
    else
        return 0
    end
end

function Irrt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(d31/2 + d32/2)*sin(theta)^4 + (d15 + d24 - d33)*sin(theta)^2*cos(theta)^2
    elseif 2*m_a-m_b == -1
        return (4*d13 - 4*1im*d23)*sin(theta)*cos(theta)^3/8 + (3*d11 + d12 - 1im*(2*d16 + d21 + 3*d22 + 2*1im*d26 - 8*d34 - 8*1im*d35))*sin(theta)^3*cos(theta)/8
    elseif 2*m_a-m_b == 1
        return (4*d13 + 4*1im*d23)*sin(theta)*cos(theta)^3/8 + (3*d11 + d12 + 1im*(2*d16 + d21 + 3*d22 - 2*1im*d26 - 8*d34 + 8*1im*d35))*sin(theta)^3*cos(theta)/8
    elseif 2*m_a-m_b == -2
        return -((d31 - d32 - 2*1im*d36)*sin(theta)^3 + 2*1im*(d14 + 1im*d15 - 1im*d24 + d25)*sin(theta)*cos(theta)^2)*sin(theta)/4
    elseif 2*m_a-m_b == 2
        return (-d31 + d32 - 2*1im*d36)*sin(theta)^4/4 + (2*1im*d14 + 2*d15 - 2*d24 + 2*1im*d25)*sin(theta)^2*cos(theta)^2/4
    elseif 2*m_a-m_b == -3
        return (d11/8 - d12/8 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/8)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == 3
        return (d11/8 - d12/8 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/8)*sin(theta)^3*cos(theta)
    else
        return 0
    end
end

function Irrp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d13*sin(theta)*cos(theta)^2 - d25*sin(2*theta) + (d11 + d12)*sin(theta)^3)*sin(theta)/2
    elseif 2*m_a-m_b == -1
        return (4*d23*cos(theta)^2 + 8*1im*(d14 + 1im*d15)*sin(theta)^2*cos(theta) + (3*d21 + d22 - 2*1im*d26)*sin(theta)^2)*sin(theta)/8
    elseif 2*m_a-m_b == 1
        return (4*d23*cos(theta)^2 - 8*1im*(d14 - 1im*d15)*sin(theta)^2*cos(theta) + (3*d21 + d22 + 2*1im*d26)*sin(theta)^2)*sin(theta)/8
    elseif 2*m_a-m_b == -2
        return (-2*1im*d24*sin(theta)*cos(theta) + d25*sin(2*theta) + (-d11 + d12 + 2*1im*d16)*sin(theta)^3)*sin(theta)/4
    elseif 2*m_a-m_b == 2
        return (2*1im*d24*sin(theta)*cos(theta) + d25*sin(2*theta) + (-d11 + d12 - 2*1im*d16)*sin(theta)^3)*sin(theta)/4
    elseif 2*m_a-m_b == -3
        return (d21/8 - d22/8 - 1im*d26/4)*sin(theta)^3
    elseif 2*m_a-m_b == 3
        return (d21/8 - d22/8 + 1im*d26/4)*sin(theta)^3
    else
        return 0
    end
end

function Irtr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (d31/2 + d32/2 - d33 + (2*d15 + 2*d24 + d31 + d32 - 2*d33)*cos(2*theta)/2)*sin(theta)^2
    elseif 2*m_a-m_b == -1
        return (-1im*d34 + d35)*sin(theta)*cos(theta)^3 + (3*d11/4 + d12/4 - d13 - 1im*d16/2 - 1im*d21/4 - 3*1im*d22/4 + 1im*d23 + d26/2 + 1im*d34 - d35)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == 1
        return (1im*d34 + d35)*sin(theta)*cos(theta)^3 + (3*d11/4 + d12/4 - d13 + 1im*d16/2 + 1im*d21/4 + 3*1im*d22/4 - 1im*d23 + d26/2 - 1im*d34 - d35)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -2
        return (d31/4 - d32/4 - 1im*d36/2 + (-2*1im*d14 + 2*d15 - 2*d24 - 2*1im*d25 + d31 - d32 - 2*1im*d36)*cos(2*theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == 2
        return (d31/4 - d32/4 + 1im*d36/2 + (2*1im*d14 + 2*d15 - 2*d24 + 2*1im*d25 + d31 - d32 + 2*1im*d36)*cos(2*theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == -3
        return (d11/4 - d12/4 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/4)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == 3
        return (d11/4 - d12/4 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/4)*sin(theta)^3*cos(theta)
    else
        return 0
    end
end

function Irtt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (d15 + d24)*sin(theta)*cos(theta)^3 - (d15 + d24 + d31 + d32 - 2*d33)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -1
        return (3*d11/8 + d12/8 - d13/2 - 1im*d16/4 - 1im*d21/8 - 3*1im*d22/8 + 1im*d23/2 + d26/4 + (3*d11 + d12 - 4*d13 - 2*1im*d16 - 1im*d21 - 3*1im*d22 + 4*1im*d23 + 2*d26 + 8*1im*d34 - 8*d35)*cos(2*theta)/8)*sin(theta)^2
    elseif 2*m_a-m_b == 1
        return (3*d11/8 + d12/8 - d13/2 + 1im*d16/4 + 1im*d21/8 + 3*1im*d22/8 - 1im*d23/2 + d26/4 + (3*d11 + d12 - 4*d13 + 2*1im*d16 + 1im*d21 + 3*1im*d22 - 4*1im*d23 + 2*d26 - 8*1im*d34 - 8*d35)*cos(2*theta)/8)*sin(theta)^2
    elseif 2*m_a-m_b == -2
        return (-d31/8 + d32/8 + 1im*d36/4 + (-2*1im*d14 + 2*d15 - 2*d24 - 2*1im*d25 + d31 - d32 - 2*1im*d36)*cos(2*theta)/8)*sin(2*theta)
    elseif 2*m_a-m_b == 2
        return (-d31/8 + d32/8 - 1im*d36/4 + (2*1im*d14 + 2*d15 - 2*d24 + 2*1im*d25 + d31 - d32 + 2*1im*d36)*cos(2*theta)/8)*sin(2*theta)
    elseif 2*m_a-m_b == -3
        return (d11/16 - d12/16 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/16)*sin(2*theta)^2
    elseif 2*m_a-m_b == 3
        return (d11/16 - d12/16 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/16)*sin(2*theta)^2
    else
        return 0
    end
end

function Irtp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -d25*sin(theta)^3 + d25*sin(theta)*cos(theta)^2 - (d11 + d12 - 2*d13)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -1
        return (1im*(d14 + 1im*d15)*cos(2*theta) + (3*d21 + d22 - 4*d23 - 2*1im*d26)*cos(theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == 1
        return (-1im*(d14 - 1im*d15)*cos(2*theta) + (3*d21 + d22 - 4*d23 + 2*1im*d26)*cos(theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == -2
        return ((-4*1im*d24 + 4*d25)*cos(2*theta)/8 + (-d11 + d12 + 2*1im*d16)*cos(theta)/8 + (d11 - d12 - 2*1im*d16)*cos(3*theta)/8)*sin(theta)
    elseif 2*m_a-m_b == 2
        return ((4*1im*d24 + 4*d25)*cos(2*theta)/8 + (-d11 + d12 - 2*1im*d16)*cos(theta)/8 + (d11 - d12 + 2*1im*d16)*cos(3*theta)/8)*sin(theta)
    elseif 2*m_a-m_b == -3
        return (d21/4 - d22/4 - 1im*d26/2)*sin(theta)^2*cos(theta)
    elseif 2*m_a-m_b == 3
        return (d21/4 - d22/4 + 1im*d26/2)*sin(theta)^2*cos(theta)
    else
        return 0
    end
end

function Irpr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (-2*d35*cos(theta)^2 - (d11 + d26)*sin(theta)^2 + (d14 + d36)*cos(theta))*sin(theta)^2
    elseif 2*m_a-m_b == -1
        return d34*sin(theta)*cos(theta)^2 - (d15 - 1im*(d25 + 1im*d31 + d36))*sin(theta)^3*cos(theta) - 1im*(d12 + 3*1im*d16 + 1im*d22 + d26)*sin(theta)^3/4
    elseif 2*m_a-m_b == 1
        return d34*sin(theta)*cos(theta)^2 - (d15 + 1im*d25 + d31 + 1im*d36)*sin(theta)^3*cos(theta) + (1im*d12/4 + 3*d16/4 + d22/4 + 1im*d26/4)*sin(theta)^3
    elseif 2*m_a-m_b == -2
        return (d14 - 1im*(d24 + d32 + 1im*d36))*sin(theta)^2*cos(theta)/2 + (-d11 + 1im*d16 + 1im*d21 + d26)*sin(theta)^4/2
    elseif 2*m_a-m_b == 2
        return -(d11 + 1im*(d16 + d21 + 1im*d26))*sin(theta)^4/2 + (d14 + 1im*d24 + 1im*d32 + d36)*sin(theta)^2*cos(theta)/2
    elseif 2*m_a-m_b == -3
        return -1im*(d12 + 1im*d16 - 1im*d22 + d26)*sin(theta)^3/4
    elseif 2*m_a-m_b == 3
        return 1im*(d12 - 1im*d16 + 1im*d22 + d26)*sin(theta)^3/4
    else
        return 0
    end
end

function Irpt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return d14*sin(theta)*cos(theta)^2 - d36*sin(theta)^3 - (d11 + d26 - 2*d35)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -1
        return (-(4*d15 - 4*1im*d25)*cos(theta)^2 + (4*d31 - 4*1im*d36)*sin(theta)^2 + (-1im*d12 + 3*d16 + d22 - 1im*d26 - 4*d34)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == 1
        return (-(4*d15 + 4*1im*d25)*cos(theta)^2 + (4*d31 + 4*1im*d36)*sin(theta)^2 + (1im*d12 + 3*d16 + d22 + 1im*d26 - 4*d34)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == -2
        return -(d11 - 1im*(d16 + d21 - 1im*d26))*sin(theta)^3*cos(theta)/2 + (d14 - 1im*d24)*sin(theta)*cos(theta)^2/2 + 1im*(d32 + 1im*d36)*sin(theta)^3/2
    elseif 2*m_a-m_b == 2
        return -(d11 + 1im*(d16 + d21 + 1im*d26))*sin(theta)^3*cos(theta)/2 + (d14 + 1im*d24)*sin(theta)*cos(theta)^2/2 - (1im*d32 + d36)*sin(theta)^3/2
    elseif 2*m_a-m_b == -3
        return -1im*(d12 + 1im*d16 - 1im*d22 + d26)*sin(theta)^2*cos(theta)/4
    elseif 2*m_a-m_b == 3
        return (1im*d12/4 + d16/4 - d22/4 + 1im*d26/4)*sin(theta)^2*cos(theta)
    else
        return 0
    end
end

function Irpp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (-d16 - d21)*sin(theta)^3 + (2*d15*sin(theta)^2 + d24)*sin(theta)*cos(theta)
    elseif 2*m_a-m_b == -1
        return (-1im*d22 + 3*d26 + (4*d11 - 4*1im*d16)*sin(theta)^2 - (4*d14 + 4*d25)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == 1
        return (1im*d22 + 3*d26 + (4*d11 + 4*1im*d16)*sin(theta)^2 - (4*d14 + 4*d25)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == -2
        return d24*sin(theta)*cos(theta)/2 + 1im*(d12 + 1im*d16 + 1im*d21 + d26)*sin(theta)^3/2
    elseif 2*m_a-m_b == 2
        return d24*sin(theta)*cos(theta)/2 - 1im*(d12 - 1im*(d16 + d21 + 1im*d26))*sin(theta)^3/2
    elseif 2*m_a-m_b == -3
        return (-1im*d22/4 + d26/4)*sin(theta)^2
    elseif 2*m_a-m_b == 3
        return (1im*d22/4 + d26/4)*sin(theta)^2
    else
        return 0
    end
end

function Ittr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return (d31/2 + d32/2)*sin(theta)*cos(theta)^3 - (d15 + d24 - d33)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -1
        return (4*d13 - 4*1im*d23)*sin(theta)^4/8 + (3*d11 + d12 - 1im*(2*d16 + d21 + 3*d22 + 2*1im*d26 - 8*d34 - 8*1im*d35))*sin(theta)^2*cos(theta)^2/8
    elseif 2*m_a-m_b == 1
        return (4*d13 + 4*1im*d23)*sin(theta)^4/8 + (3*d11 + d12 + 1im*(2*d16 + d21 + 3*d22 - 2*1im*d26 - 8*d34 + 8*1im*d35))*sin(theta)^2*cos(theta)^2/8
    elseif 2*m_a-m_b == -2
        return ((d31 - d32 - 2*1im*d36)*cos(theta)^3 + (1im*d14 - d15 + d24 + 1im*d25)*sin(theta)*sin(2*theta))*sin(theta)/4
    elseif 2*m_a-m_b == 2
        return ((d31 - d32 + 2*1im*d36)*cos(theta)^3 + (-1im*d14 - d15 + d24 - 1im*d25)*sin(theta)*sin(2*theta))*sin(theta)/4
    elseif 2*m_a-m_b == -3
        return (d11/32 - d12/32 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/32)*sin(2*theta)^2
    elseif 2*m_a-m_b == 3
        return (d11/32 - d12/32 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/32)*sin(2*theta)^2
    else
        return 0
    end
end

function Ittt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d33*sin(theta)^3 + (2*d15 + 2*d24 + d31 + d32)*sin(theta)*cos(theta)^2)*sin(theta)/2
    elseif 2*m_a-m_b == -1
        return ((4*d13 - 4*1im*(d23 + 2*d34 + 2*1im*d35))*sin(theta)^2*cos(theta) + (3*d11 + d12 - 1im*(2*d16 + d21 + 3*d22 + 2*1im*d26))*cos(theta)^3)*sin(theta)/8
    elseif 2*m_a-m_b == 1
        return ((4*d13 + 4*1im*d23 + 8*1im*d34 + 8*d35)*sin(theta)^2*cos(theta) + (3*d11 + d12 + 2*1im*d16 + 1im*d21 + 3*1im*d22 + 2*d26)*cos(theta)^3)*sin(theta)/8
    elseif 2*m_a-m_b == -2
        return (1im*d14/8 - d15/8 + d24/8 + 1im*d25/8 - d31/16 + d32/16 + 1im*d36/8)*sin(2*theta)^2
    elseif 2*m_a-m_b == 2
        return (-1im*d14/8 - d15/8 + d24/8 - 1im*d25/8 - d31/16 + d32/16 - 1im*d36/8)*sin(2*theta)^2
    elseif 2*m_a-m_b == -3
        return (d11/8 - d12/8 - 1im*(2*d16 + d21 - d22 - 2*1im*d26)/8)*sin(theta)*cos(theta)^3
    elseif 2*m_a-m_b == 3
        return (d11/8 - d12/8 + 1im*(2*d16 + d21 - d22 + 2*1im*d26)/8)*sin(theta)*cos(theta)^3
    else
        return 0
    end
end

function Ittp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d13*sin(theta)^3 + d25*sin(2*theta) + (d11 + d12)*sin(theta)*cos(theta)^2)*sin(theta)/2
    elseif 2*m_a-m_b == -1
        return ((4*d23 + 4*(-2*1im*d14 + 2*d15)*cos(theta))*sin(theta)^2 + (3*d21 + d22 - 2*1im*d26)*cos(theta)^2)*sin(theta)/8
    elseif 2*m_a-m_b == 1
        return ((4*d23 + 4*(2*1im*d14 + 2*d15)*cos(theta))*sin(theta)^2 + (3*d21 + d22 + 2*1im*d26)*cos(theta)^2)*sin(theta)/8
    elseif 2*m_a-m_b == -2
        return -(-2*1im*d24 + 2*d25 + (d11 - d12 - 2*1im*d16)*cos(theta))*sin(theta)^2*cos(theta)/4
    elseif 2*m_a-m_b == 2
        return -(2*1im*d24 + 2*d25 + (d11 - d12 + 2*1im*d16)*cos(theta))*sin(theta)^2*cos(theta)/4
    elseif 2*m_a-m_b == -3
        return (d21/8 - d22/8 - 1im*d26/4)*sin(theta)*cos(theta)^2
    elseif 2*m_a-m_b == 3
        return (d21/8 - d22/8 + 1im*d26/4)*sin(theta)*cos(theta)^2
    else
        return 0
    end
end

function Itpr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -d14*sin(theta)^3 + d36*sin(theta)*cos(theta)^2 - (d11 + d26 - 2*d35)*sin(theta)^3*cos(theta)
    elseif 2*m_a-m_b == -1
        return ((4*d15 - 4*1im*d25)*sin(theta)^2 - (4*d31 - 4*1im*d36)*cos(theta)^2 + (-1im*d12 + 3*d16 + d22 - 1im*d26 - 4*d34)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == 1
        return ((4*d15 + 4*1im*d25)*sin(theta)^2 - (4*d31 + 4*1im*d36)*cos(theta)^2 + (1im*d12 + 3*d16 + d22 + 1im*d26 - 4*d34)*cos(theta))*sin(theta)^2/4
    elseif 2*m_a-m_b == -2
        return -(d11 - 1im*(d16 + d21 - 1im*d26))*sin(theta)^3*cos(theta)/2 - (d14 - 1im*d24)*sin(theta)^3/2 + (-1im*d32 + d36)*sin(theta)*cos(theta)^2/2
    elseif 2*m_a-m_b == 2
        return -((d11 + 1im*(d16 + d21 + 1im*d26))*sin(theta)^2*cos(theta) + (d14 + 1im*d24)*sin(theta)^2 + (-1im*d32 - d36)*cos(theta)^2)*sin(theta)/2
    elseif 2*m_a-m_b == -3
        return -1im*(d12 + 1im*d16 - 1im*d22 + d26)*sin(theta)^2*cos(theta)/4
    elseif 2*m_a-m_b == 3
        return (1im*d12/4 + d16/4 - d22/4 + 1im*d26/4)*sin(theta)^2*cos(theta)
    else
        return 0
    end
end

function Itpt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d35*sin(theta)^2 + (d11 + d26)*cos(theta)^2 + (d14 + d36)*cos(theta))*sin(theta)^2
    elseif 2*m_a-m_b == -1
        return (d34 + (d15 - 1im*(d25 + 1im*d31 + d36))*cos(theta))*sin(theta)^3 - 1im*(d12 + 3*1im*d16 + 1im*d22 + d26)*sin(theta)*cos(theta)^2/4
    elseif 2*m_a-m_b == 1
        return d34*sin(theta)^3 + (d15 + 1im*d25 + d31 + 1im*d36)*sin(theta)^3*cos(theta) + (1im*d12/4 + 3*d16/4 + d22/4 + 1im*d26/4)*sin(theta)*cos(theta)^2
    elseif 2*m_a-m_b == -2
        return -(d14 + (d11 - 1im*(d16 + d21 - 1im*d26))*cos(theta) - 1im*(d24 + d32 + 1im*d36))*sin(theta)^2*cos(theta)/2
    elseif 2*m_a-m_b == 2
        return -(d14 + 1im*d24 + 1im*d32 + d36 + (d11 + 1im*(d16 + d21 + 1im*d26))*cos(theta))*sin(theta)^2*cos(theta)/2
    elseif 2*m_a-m_b == -3
        return -1im*(d12 + 1im*d16 - 1im*d22 + d26)*sin(theta)*cos(theta)^2/4
    elseif 2*m_a-m_b == 3
        return 1im*(d12 - 1im*d16 + 1im*d22 + d26)*sin(theta)*cos(theta)^2/4
    else
        return 0
    end
end

function Itpp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d15*sin(theta)^2 + d24 + (d16 + d21)*cos(theta))*sin(theta)^2
    elseif 2*m_a-m_b == -1
        return -1im*(d22 + 3*1im*d26)*sin(2*theta)/8 + (d14 + d25 + (d11 - 1im*d16)*cos(theta))*sin(theta)^3
    elseif 2*m_a-m_b == 1
        return (d14 + d25)*sin(theta)^3 + (1im*d22 + 3*d26 + (4*d11 + 4*1im*d16)*sin(theta)^2)*sin(theta)*cos(theta)/4
    elseif 2*m_a-m_b == -2
        return (-d24/2 - (-1im*d12 + d16 + d21 - 1im*d26)*cos(theta)/2)*sin(theta)^2
    elseif 2*m_a-m_b == 2
        return (-d24/2 - (1im*d12 + d16 + d21 + 1im*d26)*cos(theta)/2)*sin(theta)^2
    elseif 2*m_a-m_b == -3
        return (-1im*d22/4 + d26/4)*sin(theta)*cos(theta)
    elseif 2*m_a-m_b == 3
        return (1im*d22/4 + d26/4)*sin(theta)*cos(theta)
    else
        return 0
    end
end 

function Ippr(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -d16*sin(theta)^3 + (d31*sin(theta)^3 + d32*sin(theta)/2)*cos(theta)
    elseif 2*m_a-m_b == -1
        return (3*d12 - 1im*d22 - 8*d36*cos(theta) + (4*d11 - 4*1im*d21)*sin(theta)^2)*sin(theta)^2/8
    elseif 2*m_a-m_b == 1
        return (3*d12 + 1im*d22 - 8*d36*cos(theta) + (4*d11 + 4*1im*d21)*sin(theta)^2)*sin(theta)^2/8
    elseif 2*m_a-m_b == -2
        return d32*sin(theta)*cos(theta)/4 - (2*d16 - 2*1im*d26)*sin(theta)^3/4
    elseif 2*m_a-m_b == 2
        return d32*sin(theta)*cos(theta)/4 - (2*d16 + 2*1im*d26)*sin(theta)^3/4
    elseif 2*m_a-m_b == -3
        return (d12/8 - 1im*d22/8)*sin(theta)^2
    elseif 2*m_a-m_b == 3
        return (d12/8 + 1im*d22/8)*sin(theta)^2
    else
        return 0
    end
end

function Ippt(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(d16*sin(2*theta) + 2*d31*sin(theta)^3 + d32*sin(theta))*sin(theta)/2
    elseif 2*m_a-m_b == -1
        return d36*sin(theta)^3 + (3*d12 - 1im*d22 + (4*d11 - 4*1im*d21)*sin(theta)^2)*sin(theta)*cos(theta)/8
    elseif 2*m_a-m_b == 1
        return d36*sin(theta)^3 + (3*d12 + 1im*d22 + (4*d11 + 4*1im*d21)*sin(theta)^2)*sin(theta)*cos(theta)/8
    elseif 2*m_a-m_b == -2
        return (-d32/4 - (2*d16 - 2*1im*d26)*cos(theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == 2
        return (-d32/4 - (2*d16 + 2*1im*d26)*cos(theta)/4)*sin(theta)^2
    elseif 2*m_a-m_b == -3
        return (d12/8 - 1im*d22/8)*sin(theta)*cos(theta)
    elseif 2*m_a-m_b == 3
        return (d12/8 + 1im*d22/8)*sin(theta)*cos(theta)
    else
        return 0
    end
end

function Ippp(theta, m_a, m_b, d)
    d11, d12, d13, d14, d15, d16 = d[1]
    d21, d22, d23, d24, d25, d26 = d[2]
    d31, d32, d33, d34, d35, d36 = d[3]
    if 2*m_a-m_b == 0
        return -(2*d11*sin(theta)^2 + d12 + 2*d26)*sin(theta)^2/2
    elseif abs(2*m_a-m_b) == 1
        return 3*d22*sin(theta)/8 + (8*d16 + 4*d21)*sin(theta)^3/8
    elseif abs(2*m_a-m_b) == 2
        return (-d12/4 - d26/2)*sin(theta)^2
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
