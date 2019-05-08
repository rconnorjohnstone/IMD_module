#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# INTEGRATION MODULE
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

module Integrate

    using ..imd
    using LinearAlgebra

    export RK4,DP5,Solver,propogate

    const a_dp5 = [0  0            0             0            0          0             0;
                   0  1//5         0             0            0          0             0;
                   0  3//40        9//40         0            0          0             0;
                   0  44//45       -56//15       32//9        0          0             0;
                   0  19372//6561  -25360//2187  64448//6561  -212//729  0             0;
                   0  9017//3168   -355//33      46732//5247  49//176    -5103//18656  0;
                   0  35//384      0             500//1113    125//192   -2187//6784   11//84]
    const b_dp5 = [35//384, 0, 500//1113, 125//192, -2187//6784, 11//84, 0]
    const b2_dp5 = [5179//57600, 0, 7571//16695, 393//640, -92097//339200, 187//2100, 1//40]
    const c_dp5 = [0, 1//5, 3//10, 4//5, 8//9, 1, 1]

    function RK4(x::Vector,t::Float64,h::Float64,f::Function,args...)
        a = [0  0     0     0;
             0  1//2  0     0;
             0  0     1//2  0;
             0  0     0     1]
        b = [1//6, 1//3, 1//3, 1//6]
        c = [0, 1//2, 1//2, 1]
        k1, k2, k3, k4 = zero(x), zero(x), zero(x), zero(x)
        k1 = f(x + h*[k1 k2 k3 k4]*a[1,:], t + h*c[1], args...)
        k2 = f(x + h*[k1 k2 k3 k4]*a[2,:], t + h*c[2], args...)
        k3 = f(x + h*[k1 k2 k3 k4]*a[3,:], t + h*c[3], args...)
        k4 = f(x + h*[k1 k2 k3 k4]*a[4,:], t + h*c[4], args...)
        return t + h, x + h*(b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4)
    end

    function RK4(x::Float64,t::Float64,h::Float64,f::Function,args...)
        a = [0  0     0     0;
             0  1//2  0     0;
             0  0     1//2  0;
             0  0     0     1]
        b = [1//6, 1//3, 1//3, 1//6]
        c = [0, 1//2, 1//2, 1]
        k1, k2, k3, k4 = zero(x), zero(x), zero(x), zero(x)
        k1 = f(x + (h*[k1 k2 k3 k4]*a[1,:])[1], t + h*c[1], args...)
        k2 = f(x + (h*[k1 k2 k3 k4]*a[2,:])[1], t + h*c[2], args...)
        k3 = f(x + (h*[k1 k2 k3 k4]*a[3,:])[1], t + h*c[3], args...)
        k4 = f(x + (h*[k1 k2 k3 k4]*a[4,:])[1], t + h*c[4], args...)
        return t + h, x + h*(b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4)
    end

    function DP5(x::Array{Float64,1},t::Float64,h::Float64,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1,
        a::Array{Rational{Int64},2}=a_dp5,b::Array{Rational{Int64},1}=b_dp5,b2::Array{Rational{Int64},1}=b2_dp5,
        c::Array{Rational{Int64},1}=c_dp5)

        k1, k2, k3, k4, k5, k6, k7 = zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), zero(x)
        k1 = x
        k3 = f(x + h*[k1 k2 k3 k4 k5 k6 k7]*a[3,:], t + h*c[3], args...)
        k4 = f(x + h*[k1 k2 k3 k4 k5 k6 k7]*a[4,:], t + h*c[4], args...)
        k5 = f(x + h*[k1 k2 k3 k4 k5 k6 k7]*a[5,:], t + h*c[5], args...)
        k6 = f(x + h*[k1 k2 k3 k4 k5 k6 k7]*a[6,:], t + h*c[6], args...)
        k7 = f(x + h*[k1 k2 k3 k4 k5 k6 k7]*a[7,:], t + h*c[7], args...)
        soln = x + h*(b[1]*k1 + b[3]*k3 + b[4]*k4 + b[5]*k5 + b[6]*k6 + b[7]*k7)
        soln2 = x + h*(b2[1]*k1 + b2[3]*k3 + b2[4]*k4 + b2[5]*k5 + b2[6]*k6 + b2[7]*k7)
        err = soln - soln2
        # println(norm(err), " ", t)
        tol = min(reltol*norm(x),abstol)
        h_plus = h * min(max((tol/norm(err))^0.99,0.3),100)
        return t + h, soln, h_plus
    end

    function DP5(x::Float64,t::Float64,h::Float64,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1)
        a = [0  0            0             0            0          0             0;
             0  1//5         0             0            0          0             0;
             0  3//40        9//40         0            0          0             0;
             0  44//45       -56//15       32//9        0          0             0;
             0  19372//6561  -25360//2187  64448//6561  -212//729  0             0;
             0  9017//3168   -355//33      46732//5247  49//176    -5103//18656  0;
             0  35//384      0             500//1113    125//192   -2187//6784   11//84]
        b = [35//384, 0, 500//1113, 125//192, -2187//6784, 11//84, 0]
        b2 = [5179//57600, 0, 7571//16695, 393//640, -92097//339200, 187//2100, 1//40]
        c = [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
        k1, k2, k3, k4, k5, k6, k7 = zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), zero(x)
        k1 = x
        # k2 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[2,:])[1], t + h*c[2], args...)
        k3 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[3,:])[1], t + h*c[3], args...)
        k4 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[4,:])[1], t + h*c[4], args...)
        k5 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[5,:])[1], t + h*c[5], args...)
        k6 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[6,:])[1], t + h*c[6], args...)
        k7 = f(x + (h*[k1 k2 k3 k4 k5 k6 k7]*a[7,:])[1], t + h*c[7], args...)
        soln = x + h*(b[1]*k1 + b[3]*k3 + b[4]*k4 + b[5]*k5 + b[6]*k6 + b[7]*k7)
        soln2 = x + h*(b2[1]*k1 + b2[3]*k3 + b2[4]*k4 + b2[5]*k5 + b2[6]*k6 + b2[7]*k7)
        err = soln - soln2
        tol = min(reltol*abs(x),abstol)
        h_plus = 0.9 * h * min(max(tol/norm(err),0.3),2)
        return t + h, soln, h_plus
    end

    function Solver(m::Function,x0::Vector,t0::Float64,tf::Float64,N::Int64,f::Function,args...)
        h = (tf-t0)/N
        d = length(x0)
        ts = zeros(Float64,N+1)
        xs = zeros(Float64,d,N+1)

        ts[1] = t0
        xs[:,1] = x0

        for i in 2:(N+1)
            ts[i],xs[:,i] = m(xs[:,i-1],ts[i-1],h,f,args...)
        end

        return Solution(xs',ts)
    end

    function Solver(m::Function,x0::Float64,t0::Float64,tf::Float64,N::Int64,f::Function,args...)
        h = (tf-t0)/N
        d = length(x0)
        ts = zeros(Float64,N+1)
        xs = zeros(Float64,N+1)

        ts[1] = t0
        xs[1] = x0

        for i in 2:(N+1)
            ts[i],xs[i] = m(xs[i-1],ts[i-1],h,f,args...)
        end

        return Solution(xs,ts)
    end

    function Solver(m::Function,x0::Vector,t0::Real,tf::Real,N::Real,f::Function,args...)
        t0,tf,N = Float64(t0),Float64(tf),Int64(N)
        return Solver(m,x0,t0,tf,N,f,args...)
    end

    function Solver(m::Function,x0::Real,t0::Real,tf::Real,N::Real,f::Function,args...)
        x0,t0,tf,N = Float64(x0),Float64(t0),Float64(tf),Int64(N)
        return Solver(m,x0,t0,tf,N,f,args...)
    end

    function Solver(m::Function,x0::Vector,t0::Float64,tf::Float64,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1)
        h = min(reltol*norm(x0), abstol)
        d = length(x0)

        ts = [t0]
        xs = reshape(x0,(length(x0),1))

        while ts[end] < tf
            t,x,h = m(xs[:,end],ts[end],h,f,args...;reltol=reltol,abstol=abstol)
            push!(ts,t)
            xs = [xs x]
        end

        return Solution(xs',ts)
    end

    function Solver(m::Function,x0::Float64,t0::Float64,tf::Float64,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1)
        h = min(reltol*norm(x0), abstol)
        d = length(x0)

        ts = [t0]
        xs = [x0]

        while ts[end] < tf
            t,x,h = m(xs[end],ts[end],h,f,args...;reltol=reltol,abstol=abstol)
            push!(ts,t)
            push!(xs,x)
        end

        return Solution(xs',ts)
    end

    function Solver(m::Function,x0::Vector,t0::Real,tf::Real,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1)
        t0,tf = Float64(t0),Float64(tf)
        return Solver(m,x0r,t0,tf,f,args...,reltol=reltol,abstol=abstol)
    end

    function Solver(m::Function,x0::Real,t0::Real,tf::Real,f::Function,args...;
        reltol::Float64=1e-1,abstol::Float64=1e-1)
        x0,t0,tf = Float64(x0),Float64(t0),Float64(tf)
        return Solver(m,x0r,t0,tf,f,args...,reltol=reltol,abstol=abstol)
    end

    function propogate(state::imd.State,t0::Real,tf::Real,N::Int64=100000,
        m::Function=imd.RK4,f::Function=imd.ForceModels.EOM_2BP)
        return Solver(m,state.cart_vec,t0,tf,N,f,state.central_body.μ)
    end

end

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# FORCE MODELS MODULE
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

module ForceModels
    using LinearAlgebra, ForwardDiff
    using ..imd

    export EOM_2BP, norm_EOM_3BP, norm_3BP_STM

    function EOM_2BP(x::Vector,t::Real,μ::Real)
        grav = -μ/norm(x[1:3])^3
        return [x[4],x[5],x[6],grav*x[1],grav*x[2],grav*x[3]]
    end

    function norm_EOM_3BP(state::Array{Float64,1},μ::Float64,t::Float64)
        x,y,z,xd,yd,zd = state
        r1 = √((x+μ)^2 + y^2 + z^2)
        r2 = √((x-1+μ)^2 + y^2 + z^2)
        xdd = 2yd + x - (1-μ)*(x+μ)/r1^3 - μ*(x-1+μ)/r2^3
        ydd = -2xd + y - (1-μ)*y/r1^3 - μ*y/r2^3
        zdd = -(1-μ)*z/r1^3 - μ*z/r2^3
        return [xd,yd,zd,xdd,ydd,zdd]
    end

    function norm_3BP_STM(state::Array{Float64,1},μ::Float64,t::Float64)
        Φ = reshape(state[7:end],6,6)
        function differ(state::AbstractArray)
            x,y,z,xd,yd,zd = state
            r1 = √((x+μ)^2 + y^2 + z^2)
            r2 = √((x-1+μ)^2 + y^2 + z^2)
            xdd = 2yd + x - (1-μ)*(x+μ)/r1^3 - μ*(x-1+μ)/r2^3
            ydd = -2xd + y - (1-μ)*y/r1^3 - μ*y/r2^3
            zdd = -(1-μ)*z/r1^3 - μ*z/r2^3
            return [xd,yd,zd,xdd,ydd,zdd]
        end
        xdot = differ(state[1:6])
        A = ForwardDiff.jacobian(differ,state[1:6])
        Φdot = A * Φ
        result = vcat(xdot,vec(Φdot))
        return result
    end

end

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# LAMBERTS MODULE
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

module Lamberts
    using LinearAlgebra
    using ..imd

    export lamberts, costed_lamberts, solve_lamberts

    function solve_lamberts(state1::State,state2::State,tof_req::Real)
        r1 = state1.cart_vec[1:3] ; r1mag = norm(r1)
        r2 = state2.cart_vec[1:3] ; r2mag = norm(r2)
        μ = state1.central_body == state2.central_body ? state1.central_body.μ : error("Not the same planet")

        cos_dθ = dot(r1,r2)/(r1mag*r2mag)
        dθ = atan(r2[2],r2[1]) - atan(r1[2],r1[1])
        dθ = dθ > 2π ? dθ-2π : dθ
        dθ = dθ < 0.0 ? dθ+2π : dθ
        DM = abs(dθ) > π ? -1 : 1
        A = DM * √(r1mag * r2mag * (1 + cos_dθ))
        dθ == 0 || A == 0 && error("Can't solve Lambert's Problem")

        ψ, c2, c3 = 0, 1//2, 1//6
        ψ_down = -4π ; ψ_up = 4π^2
        y = r1mag + r2mag + (A*(ψ*c3 - 1)) / √(c2) ; χ = √(y/c2)
        tof = ( χ^3*c3 + A*√(y) ) / √(μ)

        i = 0
        while abs(tof-tof_req) > 1e-2
            y = r1mag + r2mag + (A*(ψ*c3 - 1)) / √(c2)
            while y/c2 <= 0
                # println("You finally hit that weird issue... ")
                ψ += 0.1
                if ψ > 1e-6
                    c2 = (1 - cos(√(ψ))) / ψ ; c3 = (√(ψ) - sin(√(ψ))) / √(ψ^3)
                elseif ψ < -1e-6
                    c2 = (1 - cosh(√(-ψ))) / ψ ; c3 = (-√(-ψ) + sinh(√(-ψ))) / √((-ψ)^3)
                else
                    c2 = 1//2 ; c3 = 1//6
                end
                y = r1mag + r2mag + (A*(ψ*c3 - 1)) / √(c2)
            end
            χ = √(y/c2)

            tof = ( c3*χ^3 + A*√(y) ) / √(μ)
            tof < tof_req ? ψ_down = ψ : ψ_up = ψ
            ψ = (ψ_up + ψ_down) / 2

            if ψ > 1e-6
                c2 = (1 - cos(√(ψ))) / ψ ; c3 = (√(ψ) - sin(√(ψ))) / √(ψ^3)
            elseif ψ < -1e-6
                c2 = (1 - cosh(√(-ψ))) / ψ ; c3 = (-√(-ψ) + sinh(√(-ψ))) / √((-ψ)^3)
            else
                c2 = 1//2 ; c3 = 1//6
            end

            i += 1
            i > 500 && return [NaN,NaN,NaN],[NaN,NaN,NaN]
        end
        f = 1 - y/r1mag ; g_dot = 1 - y/r2mag ; g = A * √(y/μ)
        v0t = (r2 - f*r1)/g ; vft = (g_dot*r2 - r1)/g
        return v0t, vft
    end

    function solve_lamberts(p1_name::String,p2_name::String,date1::Real,date2::Real)
        state1 = imd.ephemeris(p1_name,date1)
        state2 = imd.ephemeris(p2_name,date2)
        tof_req = (date2-date1)*86400
        return solve_lamberts(state1,state2,tof_req)
    end

    function lamberts(state1::State,state2::State,tof_req::Real;
        outs::Array{String,1}=["v_inf","C3"],type::String="values")

        v_p1 = state1.cart_vec[4:6] ; v_p2 = state2.cart_vec[4:6]
        v_leave, v_arrive = solve_lamberts(state1,state2,tof_req)
        v_inf_out = v_leave-v_p1 ; v_inf_in = v_arrive-v_p2

        if type == "strings"
            outputs = Array{String,1}(undef,0)
            for out in outs
                if "v_inf_in" == out
                    push!(outputs,"v_inf_in = $v_inf_in")
                elseif "v_inf_out" == out
                    push!(outputs,"v_inf_out = $v_inf_out")
                elseif "v_inf_out_norm" == out
                    push!(outputs,"v_inf_out_norm = $(norm(v_inf_out))")
                elseif "v_inf" == out
                    push!(outputs,"v_inf = $(norm(v_inf_in))")
                elseif "C3" == out
                    push!(outputs,"C3 = $(norm(v_inf_out)^2)")
                elseif "RLA" == out
                    push!(outputs,"RLA = $(atand(v_inf_out[2]/v_inf_out[1]))")
                elseif "DLA" == out
                    push!(outputs,"DLA = $(asind(v_inf_out[3]/norm(v_inf_out)))")
                elseif "TOF" == out
                    push!(outputs,"TOF = $(tof_req/86400)")
                end
            end
        elseif type == "values"
            outputs = Array{Any,1}(undef,0)
            for out in outs
                if "v_inf_in" == out
                    push!(outputs,v_inf_in...)
                elseif "v_inf_out" == out
                    push!(outputs,v_inf_out...)
                elseif "v_inf_out_norm" == out
                    push!(outputs,norm(v_inf_out))
                elseif "v_inf" == out
                    push!(outputs,norm(v_inf_in))
                elseif "C3" == out
                    push!(outputs,norm(v_inf_out)^2)
                elseif "RLA" == out
                    push!(outputs,atand(v_inf_out[2]/v_inf_out[1]))
                elseif "DLA" == out
                    push!(outputs,asind(v_inf_out[3]/norm(v_inf_out)))
                elseif "TOF" == out
                    push!(outputs,tof_req/86400)
                end
            end
        else
            error("Invalid output-type declaration")
        end
        return outputs
    end

    function lamberts(p1_name::String,p2_name::String,date1::Real,date2::Real;
        outs::Array{String,1}=["v_inf","C3"],type::String="values")
        state1 = imd.ephemeris(p1_name,date1)
        state2 = imd.ephemeris(p2_name,date2)
        tof_req = (date2-date1)*86400
        return lamberts(state1,state2,tof_req,outs=outs,type=type)
    end

    function lamberts(p1::Planet,p2::Planet,date1::Real,date2::Real;
        outs::Array{String,1}=["v_inf","C3"],type::String="values")

        state1 = imd.ephemeris(p1,date1)
        state2 = imd.ephemeris(p2,date2)
        tof_req = (date2-date1)*86400
        return lamberts(state1,state2,tof_req,outs=outs,type=type)
    end

end
