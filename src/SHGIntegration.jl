include("ModeSolutions.jl")

d22 = 1
d31 = 1
d33 = 1

function Irrr(theta, m_a, m_b)
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

function Irrt(theta, m_a, m_b)
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

function Irrp(theta, m_a, m_b)
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

function Irtr(theta, m_a, m_b)
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

function Irtt(theta, m_a, m_b)
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

function Irtp(theta, m_a, m_b)
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

function Irpr(theta, m_a, m_b)
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

function Irpt(theta, m_a, m_b)
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

function Irpp(theta, m_a, m_b)
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

function Ittr(theta, m_a, m_b)
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

function Ittt(theta, m_a, m_b)
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

function Ittp(theta, m_a, m_b)
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

function Itpr(theta, m_a, m_b)
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

function Itpt(theta, m_a, m_b)
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

function Itpp(theta, m_a, m_b)
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

function Ippr(theta, m_a, m_b)
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

function Ippt(theta, m_a, m_b)
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

function Ippp(theta, m_a, m_b)
    if abs(2*m_a-m_b) == 1
        return 3*d22*(sin(3*theta)-2*sin(theta))/8
    elseif abs(2*m_a-m_b) == 3
        return d22*sin(theta)/8
    else
        return 0
    end
end