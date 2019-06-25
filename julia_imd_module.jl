#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# MAIN IMD MODULE
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

module imd
using LinearAlgebra, PlotlyJS

export μs,rs,as,j2s,p_colors,AU,State,Body,Solution,RK4,Planet,Star,Moon,
KepState,CartesianState,solve_lamberts,lamberts,costed_lamberts

# -----------------------------------------------------------------------------
# DEFINING CONSTANTS
# -----------------------------------------------------------------------------

# Gravitational Constants
  const μs = Dict(
  "Sun" => 1.32712440018e11,
  "Mercury" => 2.2032e4,
  "Venus" => 3.257e5,
  "Earth" => 3.986004415e5,
  "Moon" => 4.902799e3,
  "Mars" => 4.305e4,
  "Jupiter" => 1.266865361e8,
  "Saturn" => 3.794e7,
  "Uranus" => 5.794e6,
  "Neptune" => 6.809e6,
  "Pluto" => 9e2)

# Radii
  const rs = Dict(
  "Sun" => 696000.,
  "Mercury" => 2439.,
  "Venus" => 6052.,
  "Earth" => 6378.1363,
  "Moon" => 1738.,
  "Mars" => 3397.2,
  "Jupiter" => 71492.,
  "Saturn" => 60268.,
  "Uranus" => 25559.,
  "Neptune" => 24764.,
  "Pluto" => 1151.)

# Semi Major Axes
  const as = Dict(
  "Mercury" => 57909083.,
  "Venus" => 108208601.,
  "Earth" => 149598023.,
  "Moon" => 384400.,
  "Mars" => 227939186.,
  "Jupiter" => 778298361.,
  "Saturn" => 1429394133.,
  "Uranus" => 2875038615.,
  "Neptune" => 4504449769.,
  "Pluto" => 5915799000.)

# J2 for basic oblateness
  const j2s = Dict(
  "Mercury" => 0.00006,
  "Venus" => 0.000027,
  "Earth" => 0.0010826269,
  "Moon" => 0.0002027,
  "Mars" => 0.001964,
  "Jupiter" => 0.01475,
  "Saturn" => 0.01645,
  "Uranus" => 0.012,
  "Neptune" => 0.004,
  "Pluto" => 0.)

# These are just the colors for plots
  const p_colors = Dict(
  "Sun" => :Electric,
  "Mercury" => :heat,
  "Venus" => :turbid,
  "Earth" => :Blues,
  "Moon" => :Greys,
  "Mars" => :Reds,
  "Jupiter" => :solar,
  "Saturn" => :turbid,
  "Uranus" => :haline,
  "Neptune" => :ice,
  "Pluto" => :matter)

# This part is just for defining what's what later
  stars = ("Sun",)

  planets =
  ("Mercury","Venus","Earth","Mars","Jupiter","Uranus","Neptune","Pluto")

  moons = ("Moon",)

  const AU = 149597870.691 #km

# -----------------------------------------------------------------------------
# DEFINING PLANETARY OBJECTS
# -----------------------------------------------------------------------------

"""
The Body class is the basic class that covers all sorts of planetary bodies.
Each type of Body contains different fields, but they all contain:
 - name
 - gravitational parameter
 - radius
 - color (for plotting)
 - diameter (for plotting purposes, may be scaled)
"""
  abstract type Body end

"""
Planets are pretty self explanatory. Each planetary body also contains a j2 
coefficient.
"""
  struct Planet <: Body

    name::String
    μ::Real
    r::Real
    a::Real
    j2::Real
    col::Symbol
    d_size::Real

  end

"""
Currently, the Sun is the only Star
"""
  struct Star <: Body

    name::String
    μ::Real
    r::Real
    col::Symbol
    d_size::Real

  end

"""
Currently only one Moon (Luna) is supported
"""
  struct Moon <: Body

    name::String
    μ::Real
    r::Real
    a_planet::Real
    a_moon::Real
    j2::Real
    col::Symbol
    d_size::Real

  end

"""
The constructor for the Body class will just call the appropriate smaller
constructor (in other words, Body is abstract only)
"""
  function Body(name::String)

    if name in stars
      func = Star
    elseif name in planets
      func = Planet
    elseif name in moons
      func = Moon
    else
      error("Not a recognized celestial body")
    end
    return func(name)

  end

"""Star constructor"""
  Star(name::String) = Star(
  name,
  μs[name],
  rs[name],
  p_colors[name],
  15*rs[name])

"""Planet constructor"""
  Planet(name::String) = Planet(
  name,
  μs[name],
  rs[name],
  as[name],
  j2s[name],
  p_colors[name],
  rs[name])

"""Moon constructor"""
  Moon(name::String) =Moon(
  name,
  μs[name],
  rs[name],
  as["Earth"],
  as[name],
  j2s[name],
  p_colors[name],
  rs[name])

# ------------------------------------------------------------------------------
# DEFINING STATE OBJECTS
# ------------------------------------------------------------------------------

"""
States are an important part of the imd module. If I did my job correctly, then
the vast majority of the functions in this module should be written to operate
on states. I'm not a big fan of OOP just for the sake of it. The idea behind
enforcing a state class is to ensure that the vector which describes that state
is of the correct type (pos/vel, keplerian, etc.) so that errors don't occur
"""
  abstract type State end

"""
Cartesian states are the most basic kinds of states. Simply define a position
and velocity. Frame choice isn't currently strictly implemented.

TODO: Determine the frame choice and possibly allow for representation in 
multiple frames.
"""
  struct CartesianState <: State

    cart_vec::Vector
    central_body::Body
    name::String

  end

"""
Keplerian states also have a cartesian state vector (for use in the majority
of the functions). However, they are also described by a 6 element Vector of
semi-major axis, eccentricity, inclination, RAAN, argument of periapsis, and
true anomaly. All angular quantities are in radians.
"""
  struct KepState <: State

    cart_vec::Vector
    oe::Array
    central_body::Body
    name::String

  end

"""
Also called a  RIC state in some contexts, this is a basic Hill-frame defined 
orbit. The components of the RΘH frame must fit normal conventions (Θ and H 
components of position and H component of velocity must be zero). Also contains
within it a Cartesian state vector
"""
  struct RθHState <: State

    cart_vec::Vector
    rθh_vec::Vector
    central_body::Body
    name::String

  end

"""
Constructor for a Cartesian State (shortcutted with 'State')
"""
  function State(vec::Vector,central_body::Body,name::String="Spaceraft")

    return CartesianState(vec,central_body,name)

  end

"""
Constructor for a KepState
"""
  function KepState(oe::Vector,central_body::Body,name::String="Spaceraft")

    cart_vec = oe_to_xyz(oe,central_body.μ)
    return KepState(cart_vec,oe,central_body,name)

  end

"""
Constructor for a KepState, starting with a Cartesian State

This provides for very simple conversion
"""
  function KepState(cart_state::CartesianState)

    oe = xyz_to_oe(cart_state.cart_vec,cart_state.central_body.μ)
    return KepState(
    cart_state.cart_vec,
    oe,
    cart_state.central_body,
    cart_state.name)

    end

# ------------------------------------------------------------------------------
# OTHER TYPE DEFINITIONS
# ------------------------------------------------------------------------------

"""
MOSTLY DEFUNCT
Used for numerical integration. I recommend instead using ForwardDiff
"""
  struct Solution

    x::Array
    t::Array{Real,1}

  end

# ------------------------------------------------------------------------------
# FUNCTIONS THAT OTHER MODULES CAN RELY ON
# ------------------------------------------------------------------------------

"""
Provides a really basic, but fast, ephemeris for most planets.
TODO: Update with more accurate ephemeris
"""
  function ephemeris(planet::Planet,date::Real)

    century = (date-2451545)/36525
    if planet.name == "Mercury"
      table = [ 252.250906   149472.6746358  0.00000535     0.000000002 ;
                0.387098310  0               0              0 ;
                0.20563175   0.000020406     0.0000000284   0.00000000017 ;
                7.004986     -0.0059516      0.00000081     0.000000041 ;
                48.330893    -0.1254229      -0.00008833     -0.000000196 ;
                77.456119    0.1588643       -0.00001343     0.000000039]
    elseif planet.name == "Venus"
      table = [ 181.979801   58517.8156760   0.00000165      -0.000000002 ;
                0.72332982   0               0               0 ;
                0.00677188   -0.000047766    0.0000000975    0.00000000044 ;
                3.394662     -0.0008568      -0.00003244     0.000000010 ;
                76.679920    -0.2780080      -0.00014256     -0.000000198 ;
                131.563707   0.0048646       -0.00138232     -0.000005332]
    elseif planet.name == "Earth"
      table = [ 100.466449   35999.3728519   -0.00000568     0 ;
                1.000001018  0               0               0 ;
                0.01670862   -0.000042037    -0.0000001236   0.00000000004 ;
                0            0.0130546       -0.00000931     -0.000000034 ;
                174.873174   -0.2410908      0.00004067      -0.000001327 ;
                102.937348   0.3225557       0.00015026      0.000000478]
    elseif planet.name == "Mars"
      table = [ 355.433275   19140.2993313   0.00000261      -0.000000003 ;
                1.523679342  0               0               0 ;
                0.09340062   0.000090483     -0.0000000806   -0.00000000035 ;
                1.849726     -0.0081479      -0.00002255     -0.000000027 ;
                49.558093    -0.2949846      -0.00063993     -0.000002143 ;
                336.060234   0.4438898       -0.00017321     0.000000300]
    elseif planet.name == "Jupiter"
      table = [ 34.351484    3034.9056746    -0.00008501     0.000000004 ;
                5.202603191  0.0000001913    0               0 ;
                0.04849485   0.000163244     -0.0000004719   -0.00000000197 ;
                1.303270     -0.0019872      0.00003318      0.000000092 ;
                100.464441   0.1766828       0.00090387      -0.000007032 ;
                14.331309    0.2155525       0.00072252      -0.000004590]
    elseif planet.name == "Saturn"
      table = [ 50.077471    1222.1137943    0.00021004      -0.000000019 ;
                9.554909596  -0.0000021389   0               0 ;
                0.05550862   -0.000346818    -0.0000006456   0.00000000338 ;
                2.488878     0.0025515       -0.00004903     0.000000018 ;
                113.665524   -0.2566649      -0.00018345     0.000000357 ;
                93.056787    0.5665496       0.00052809      0.000004882]
    elseif planet.name == "Uranus"
      table = [ 314.055005   429.8640561     0.00030434      0.000000026 ;
                19.218446062 -0.0000000372   0.00000000098   0 ;
                0.04629590   -0.000027337    0.0000000790    0.00000000025 ;
                0.773196     0.0007744       0.00003749      -0.000000092 ;
                74.005947    0.5211258       0.00133982      0.000018516 ;
                173.005159   1.4863784       0.0021450       0.000000433]
    elseif planet.name == "Neptune"
      table = [ 304.348665   219.8833092     0.00030926      0.000000018 ;
                30.110386869 -0.0000001663   0.00000000069   0 ;
                0.00898809   0.000006408     -0.0000000008   -0.00000000005 ;
                1.769952     -0.0093082      -0.00000708     0.000000028 ;
                131.784057   1.1022057       0.00026006      -0.000000636 ;
                48.123691    1.4262677       0.00037918      -0.000000003]
    elseif planet.name == "Pluto"
      table = [ 238.92903833 145.20780515    0               0 ;
                39.48211675  -0.00031596     0               0 ;
                0.24882730   0.00005170      0               0 ;
                17.14001206  0.00004818      0               0 ;
                110.30393684 -0.01183482     0               0 ;
                224.06891629 -0.04062942     0               0]
    else
      return "No such planet"
    end
    poly = [1,century,century^2,century^3]
    L,a,e,i,Ω,P = table * poly
    L,i,Ω,P = [L,i,Ω,P] * (π/180)
    ω = P-Ω
    M = L-P
    θ = M + (2e - 0.25e^3 + (5/96)*e^5)*sin(M) +
    (1.25e^2 - (11/24)*e^4)*sin(2M) + ((13/12)*e^3 - (43/64)*e^5)*sin(3M) +
    (103/96)*e^4*sin(4M) + (1097/960)*e^5*sin(5M)
    oe = [a*AU,e,i,Ω,ω,θ]
    return KepState(oe,Body("Sun"),planet.name)

  end

"""
Overloaded ephemeris function to allow calling with just the name of the planet
"""
  function ephemeris(planet_name::String,date::Real)

    planet = Planet(planet_name)
    return ephemeris(planet,date)

  end

  function norm3BP(
    state::Array{Real,1},t::Real,L::Real,GM1::Real,GM2::Real)

    x = state[1:3]
    v = state[4:6]
    x /= L
    n = √((GM1 + GM2)/L^3)
    t *= n
    v /= n*L
    state = [x; v]
    return state, t

  end

  function norm3BP(
    state::Array{Real,1},L::Real,GM1::Real,GM2::Real)

    x = state[1:3]
    v = state[4:6]
    x /= L
    n = √((GM1 + GM2)/L^3)
    v /= n*L
    state = [x; v]
    return state

  end

  function norm3BP(t::Real,L::Real,GM1::Real,GM2::Real)

    n = √((GM1 + GM2)/L^3)
    t *= n
    return t

  end

  function norm3BP(state::Array{Real,1},t::Real,L::Real,GM1::Real,GM2::Real)

    t = typeof(t) != Real ? Real(t) : t
    L = typeof(L) != Real ? Real(L) : L
    GM1 = typeof(GM1) != Real ? GM1 = Real(GM1) : GM1
    GM2 = typeof(GM2) != Real ? GM2 = Real(GM2) : GM2
    return norm3BP(state,t,L,GM1,GM2)

  end

  function norm3BP(t::Real,L::Real,GM1::Real,GM2::Real)

    t = typeof(t) != Real ? Real(t) : t
    L = typeof(L) != Real ? Real(L) : L
    GM1 = typeof(GM1) != Real ? GM1 = Real(GM1) : GM1
    GM2 = typeof(GM2) != Real ? GM2 = Real(GM2) : GM2
    return norm3BP(state,t,L,GM1,GM2)

  end

  function norm3BP(state::Array{Real,1},L::Real,GM1::Real,GM2::Real)

    L = typeof(L) != Real ? Real(L) : L
    GM1 = typeof(GM1) != Real ? GM1 = Real(GM1) : GM1
    GM2 = typeof(GM2) != Real ? GM2 = Real(GM2) : GM2
    return norm3BP(state,L,GM1,GM2)

  end

  function norm3BP(state::State,t::Real,date::Real,body2::Body)

    body1 = state.central_body
    GM1 = body1.μ
    GM2 = body2.μ
    r1 = body1.name=="Sun" ? [0,0,0] : ephemeris(body1,date).cart_vec[1:3]
    r2 = body2.name=="Sun" ? [0,0,0] : ephemeris(body2,date).cart_vec[1:3]
    r_diff = r2-r1
    L = norm(r_diff)
    return norm3BP(state.cart_vec,t,L,GM1,GM2)

  end

  function norm3BP(state::State,date::Real,body2::Body)

    body1 = state.central_body
    GM1 = body1.μ
    GM2 = body2.μ
    r1 = body1.name=="Sun" ? [0,0,0] : ephemeris(body1,date).cart_vec[1:3]
    r2 = body2.name=="Sun" ? [0,0,0] : ephemeris(body2,date).cart_vec[1:3]
    r_diff = r2-r1
    L = norm(r_diff)
    return norm3BP(state.cart_vec,L,GM1,GM2)

  end

  function norm3BP(state::State,t::Real,date::Real,body2::Body)

    t = typeof(t) != Real ? Real(t) : t
    date = typeof(date) != Real ? Real(date) : date
    return norm3BP(state,t,date,body2)

  end

  function norm3BP(state::State,date::Real,body2::Body)

    date = typeof(date) != Real ? Real(date) : date
    return norm3BP(state,date,body2)

  end

  function JD_to_str(J::Real)

    J = J+0.5
    y = 4716; j = 1401; m = 2; n = 12; r = 4; p = 1461
    v = 3; u = 5; s = 153; w = 2; B = 274277; C = -38
    f = J + j + div(div(4J + B,146097)*3,4)+C
    e = r*f + v
    g = div(e%p,r)
    h = u*g + w
    D = UInt64(div(h%s,u)+1)
    M = UInt64((div(h,s)+m)%n + 1)
    Y = UInt64(div(e,p) - y + div(n + m - M,n))
    return Base.dec(Y,4,false) * "-" * Base.dec(M,2,false) *
    "-" * Base.dec(D,2,false)

  end


# ------------------------------------------------------------------------------
# IMPORT OTHER MODULES
# ------------------------------------------------------------------------------

include("julia_imd_integrate.jl")
using .Integrate,.ForceModels,.Lamberts

# ------------------------------------------------------------------------------
# GENERAL FUNCTIONS
# ------------------------------------------------------------------------------

  function oe_to_rθh(oe::Vector,μ::Real) :: Vector

    a,e,i,Ω,ω,ν = oe
    return [a*(1-e^2)/(1+e*cos(ν)),
    0,
    0,
    (μ/sqrt(μ*a*(1-e^2)))*e*sin(ν),
    (μ/sqrt(μ*a*(1-e^2)))*(1+e*cos(ν)),
    0]

  end

  function rθh_to_xyz(rθh_vec::Vector,oe::Vector)

    a,e,i,Ω,ω,ν = oe
    θ = ω+ν
    cΩ,sΩ,ci,si,cθ,sθ = cos(Ω),sin(Ω),cos(i),sin(i),cos(θ),sin(θ)
    DCM = [cΩ*cθ-sΩ*ci*sθ -cΩ*sθ-sΩ*ci*cθ sΩ*si;
    sΩ*cθ+cΩ*ci*sθ -sΩ*sθ+cΩ*ci*cθ -cΩ*si;
    si*sθ si*cθ ci]
    DCM = kron(Matrix(I,2,2),DCM)
    return DCM*rθh_vec

  end

  function oe_to_xyz(oe::Vector,μ::Real)

    return rθh_to_xyz(oe_to_rθh(oe,μ),oe)

  end

  function xyz_to_oe(cart_vec::Vector,μ::Real)

    r_xyz, v_xyz = cart_vec[1:3],cart_vec[4:6]
    r = norm(r_xyz)
    h_xyz = cross(r_xyz,v_xyz) #km^2/s
    h = norm(h_xyz) #km^2/s
    ξ = dot(v_xyz,v_xyz)/2 - μ/r #km^2 s^-2
    a = -μ/(2ξ) #km
    e = sqrt(1 + (2h^2*ξ)/μ^2)
    e_xyz = cross(v_xyz,h_xyz)/μ - r_xyz/r
    i = acos(h_xyz[3]/h) #rad
    n_xyz = cross([0,0,1],h_xyz)
    Ω = acos(dot(n_xyz,[1,0,0])/norm(n_xyz))
    ω = acos((dot(n_xyz,e_xyz)/(norm(n_xyz)*e)))
    ν = acos((dot(r_xyz,e_xyz))/(r*norm(e_xyz)))
    Ω = dot(n_xyz,[0,1,0]) > 0. ? Ω : -Ω
    ω = dot(e_xyz,[0,0,1]) > 0. ? ω : -ω
    ν = dot(r_xyz,v_xyz) > 0. ? ν : -ν
    return [a,e,i,Ω,ω,ν]

  end

  function rθh_to_vnb(x::Vector)

    fpa = atan(x[4]/x[5])
    dcm = [cos(fpa) -sin(fpa) 0; sin(fpa) cos(fpa) 0; 0 0 1]
    dcm = kron(Matrix(I,2,2),dcm)
    return dcm*x

  end

  function ψ_to_rp(ψ::Real,v_inf::Real,μ::Real)

    ρ = (π - ψ)/2
    return (((1/cos(ρ))-1)*μ)/norm(v_inf)^2

  end

  function bplane(
    v_inf_in::Array{Real,1},v_inf_out::Array{Real,1},μ::Real)

    abs(norm(v_inf_in) - norm(v_inf_out)) > 0.1 &&
    error("Velocities are too different")
    v_inf = norm(v_inf_in)
    ψ = acos(dot(v_inf_in,v_inf_out)/(norm(v_inf_in)*norm(v_inf_out)))
    r_p = ψ_to_rp(ψ,v_inf,μ)

    s_hat = v_inf_in/v_inf
    h_hat = cross(v_inf_in,v_inf_out)/norm(cross(v_inf_in,v_inf_out))
    b_hat = cross(s_hat,h_hat)
    b_mag = (μ/v_inf^2)*√((1 + v_inf^2 * (r_p/μ))^2 - 1)
    b = b_mag * b_hat

    t_hat = cross(s_hat,[0,0,1])/norm(cross(s_hat,[0,0,1]))
    r_hat = cross(s_hat,t_hat)

    return dot(b,t_hat), dot(b,r_hat), ψ, r_p

  end

  function bplane(v_inf_in::Array,v_inf_out::Array,μ::Real)

    abs(norm(v_inf_in) - norm(v_inf_out)) > 0.1 &&
    error("Velocities are too different")
    v_inf = norm(v_inf_in)
    ψ = acos(dot(v_inf_in,v_inf_out)/(norm(v_inf_in)*norm(v_inf_out)))
    r_p = ψ_to_rp(ψ,v_inf,μ)

    s_hat = v_inf_in/v_inf
    h_hat = cross(v_inf_in,v_inf_out)/norm(cross(v_inf_in,v_inf_out))
    b_hat = cross(s_hat,h_hat)
    b_mag = (μ/v_inf^2)*√((1 + v_inf^2 * (r_p/μ))^2 - 1)
    b = b_mag * b_hat

    t_hat = cross(s_hat,[0,0,1])/norm(cross(s_hat,[0,0,1]))
    r_hat = cross(s_hat,t_hat)

    return dot(b,t_hat), dot(b,r_hat), ψ, r_p

  end

  function bplane(v_inf_in::Array{Any,1},v_inf_out::Array{Any,1},μ::Real)

    v_inf = norm(v_inf_in)
    ψ = acos(dot(v_inf_in,v_inf_out)/(norm(v_inf_in)*norm(v_inf_out)))
    r_p = ψ_to_rp(ψ,v_inf,μ)

    s_hat = v_inf_in/v_inf
    h_hat = cross(v_inf_in,v_inf_out)/norm(cross(v_inf_in,v_inf_out))
    b_hat = cross(s_hat,h_hat)
    b_mag = (μ/v_inf^2)*√((1 + v_inf^2 * r_p/μ)^2 - 1)
    b = b_mag * b_hat

    t_hat = cross(s_hat,[0,0,1])/norm(cross(s_hat,[0,0,1]))
    r_hat = cross(s_hat,t_hat)

    return dot(b,t_hat), dot(b,r_hat), ψ, r_p

  end

  function bplane(x::State)

    r = x[1:3] ; v = x[4:6]
    μ = x.central_body.μ

    h_hat = cross(r,v)/norm(cross(r,v))
    e = (1/μ) * ( ( dot(v,v) - μ/norm(r) ) * r - dot(r,v) * v )
    ρ = acos(1/norm(e))
    ψ = 2*asin(1/norm(e))

    s_hat = cos(ρ) * e/norm(e) + sin( cross(h_hat,e) / norm(cross(h_hat,e)) )
    t_hat = cross(s_hat,[0,0,1]) / norm(cross(s_hat,[0,0,1]))
    r_hat = cross(s_hat,t_hat)

    ξ = 0.5*dot(v,v) - μ/norm(r)
    ξ > 0 || error("Must be hyperbolic")
    a = -μ/(2ξ)

    b_mag = a * √(norm(e)^2 - 1)
    b_hat = cross(s_hat,h_hat) / norm(cross(s_hat,h_hat))
    b = b_mag * b_hat

    return dot(b,t_hat), dot(b,r_hat), ψ, r_p

  end

  function patched_conics(v_inf::Real,μ::Real,final_ξ::Real,r::Real)

    flyby_ξ = v_inf^2/2
    flyby_vel = √(2*(flyby_ξ + μ/r))
    orbit_vel = √(2*(final_ξ + μ/r))
    ΔV = orbit_vel - flyby_vel
    return ΔV

  end

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

  function plot_integrated(
    path::Array{Real,2},
    central_body::Body,
    plot_theme::Symbol=:juno,
    title::String="Spacecraft Position")

    N = 32
    θ = collect(range(0,length=N,stop=2π))
    ϕ = collect(range(0,length=N,stop=π))
    x = cos.(θ) * sin.(ϕ)'
    y = sin.(θ) * sin.(ϕ)'
    z = repeat(cos.(ϕ)',outer=[N, 1])
    x_p,y_p,z_p = central_body.d_size .* (x,y,z)

    limit = maximum([x_p y_p z_p])
    limit = 1.3 * max(limit,maximum(path))

    t1 = scatter3d(;x=path[:,1],y=path[:,2],z=path[:,3],
    mode="lines",name="orbit",
    line_color="#db2100", line_width=3)
    t2 = surface(
    ;x=x_p,
    y=y_p,
    z=z_p,
    showscale=false,
    colorscale = central_body.col)

    layout = Layout(
    ;title="Orbit Plot",
    width=1900,
    height=920,
    paper_bgcolor="#222529",
    scene = attr(xaxis = attr(autorange = false,range=[-limit,limit]),
    yaxis = attr(autorange = false,range=[-limit,limit]),
    zaxis = attr(autorange = false,range=[-limit,limit]),
    aspectratio=attr(x=1,y=1,z=1),
    aspectmode="manual"))

    p = Plot([t1,t2],layout)
    plot(p) ;

  end

  function plot_integrated(path::Integrate.Solution,central_body::Body)

    x = path.x
    return plot_integrated(x,central_body)

  end

  function plot_orbit(state::State,N::Int64=200000)

    μ = state.central_body.μ
    r_xyz, v_xyz = state.cart_vec[1:3],state.cart_vec[4:6]
    r = norm(r_xyz)
    ξ = dot(v_xyz,v_xyz)/2 - μ/r #km^2 s^-2
    a = -μ/(2ξ) #km
    T = 2π*√(a^3/μ) #s
    return plot_integrated(Integrate.propogate(state,0,T,N),state.central_body)

  end

  function plot_lamberts(
    p1::Planet,p2::Planet,date1::Real,date2::Real,N::Int64=100000)

    t1 = 86400date1; t2 = 86400date2
    p1_path = Array{Real,2}(undef, 0, 6)
    p2_path = Array{Real,2}(undef, 0, 6)
    for date in range(date1,date2,length=100)
        p1_path = [p1_path; ephemeris(p1,date).cart_vec' ]
        p2_path = [p2_path; ephemeris(p2,date).cart_vec' ]
    end
    craft_vec = [ ephemeris(p1,date1).cart_vec[1:3] ;
    solve_lamberts(p1.name,p2.name,date1,date2)[1] ]
    craft_state = State(craft_vec,Star("Sun"))
    craft_path = propogate(craft_state,t1,t2,N).x

    N = 32
    θ = collect(range(0,length=N,stop=2π))
    ϕ = collect(range(0,length=N,stop=π))
    x = cos.(θ) * sin.(ϕ)'
    y = sin.(θ) * sin.(ϕ)'
    z = repeat(cos.(ϕ)',outer=[N, 1])
    x_p,y_p,z_p = 15*rs["Sun"] .* (x,y,z)

    limit = maximum([x_p y_p z_p])
    limit = max(limit,maximum(abs.(craft_path)))
    limit = max(limit,maximum(abs.(p1_path)))
    limit = 1.1*max(limit,maximum(abs.(p2_path)))

    t_sun = surface(;x=x_p,y=y_p,z=z_p,showscale=false,colorscale = :Electric)
    t_craft = scatter3d(
    ;x=craft_path[:,1],
    y=craft_path[:,2],
    z=craft_path[:,3],
    mode="lines",
    name="Spacecraft",
    line_color="#db2100",
    line_width=4)

    t_p1 = scatter3d(
    ;x=p1_path[:,1],
    y=p1_path[:,2],
    z=p1_path[:,3],
    mode="lines",
    name=p1.name,
    line_color="#2777f7",
    line_width=3)

    t_p2 = scatter3d(
    ;x=p2_path[:,1],
    y=p2_path[:,2],
    z=p2_path[:,3],
    mode="lines",
    name=p2.name,
    line_color="#34f727",
    line_width=3)

    layout = Layout(
    ;title="Orbit Plot",
    width=1900,
    height=920,
    paper_bgcolor="#222529",
    scene = attr(xaxis = attr(autorange = false,range=[-limit,limit]),
    yaxis = attr(autorange = false,range=[-limit,limit]),
    zaxis = attr(autorange = false,range=[-limit,limit]),
    aspectratio=attr(x=1,y=1,z=1),
    aspectmode="manual"))

    p = Plot([t_sun,t_craft,t_p1,t_p2],layout)
    plot(p) ;

  end

  function plot_lamberts(
    p1_name::String,p2_name::String,date1::Real,date2::Real,N::Int64=100000)

    p1 = imd.Planet(p1_name)
    p2 = imd.Planet(p2_name)
    return plot_lamberts(p1,p2,date1,date2,N)

  end

  function emp()
  end

  function porkchop(
    p1::Planet,
    p2::Planet,
    leave::Array{Real,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      costfxn::Function=emp,
      fxnargs=())

    if "v_inf_in" in outs
      if "v_inf_out" in outs
        outputs = zeros(length(outs)+4,M,N)
      else
        outputs = zeros(length(outs)+2,M,N)
      end
    elseif "v_inf_out" in outs
      outputs = zeros(length(outs)+2,M,N)
    else
      outputs = zeros(length(outs),M,N)
    end
    i = 0
    if M == 1
      xs = [leave[1]]
    else
      xs = collect(range(leave[1],leave[2],length=M))
    end

    if N == 1
      ys = [arrive[1]]
    else
      ys = collect(range(arrive[1],arrive[2],length=N))
    end

    for t1 = xs
      i += 1 ; j = 0
      for t2 = ys
        j += 1
        p1_leave_state = ephemeris(p1,t1)
        p2_arrive_state = ephemeris(p2,t2)
        try
          outputs[:,i,j] = Array{Real,1}(lamberts(p1,p2,t1,t2,outs=outs))
        catch
          outputs[:,i,j] = NaN*ones(size(outputs[:,i,j]))
        end
      end
    end
    return costfxn != emp ? costfxn(outputs,fxnargs...) : outputs

  end

  function porkchop(
    p1::Planet,
    p2::Planet,
    leave::Array{Int64,1},
    arrive::Array{Int64,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      costfxn::Function=emp,
      fxnargs=())

    leave = Real.(leave) ; arrive = Real.(arrive)

    if costfxn != emp
      return porkchop(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      costfxn=costfxn,
      fxnargs=fxnargs)
    else
      return porkchop(p1,p2,leave,arrive;outs=outs,M=M,N=N)
    end

  end

  function porkchop(
    p1_name::String,
    p2_name::String,
    leave::Array{Real,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      costfxn::Function=emp,
      fxnargs=())

    p1 = Planet(p1_name)
    p2 = Planet(p2_name)

    if costfxn != emp
      return porkchop(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      costfxn=costfxn,
      fxnargs=fxnargs)
    else
      return porkchop(p1,p2,leave,arrive;outs=outs,M=M,N=N)
    end
  end

  function porkchop(
    p1_name::String,
    p2_name::String,
    leave::Array{Number,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      costfxn::Function=emp,
      fxnargs=())

    leave = Real.(leave) ; arrive = Real.(arrive)

    if costfxn != emp
      return porkchop(
      p1_name,
      p2_name,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      costfxn=costfxn,
      fxnargs=fxnargs)
    else
      return porkchop(p1_name,p2_name,leave,arrive;outs=outs,M=M,N=N)
    end

  end

  function porkchop_plot(
    p1::Planet,
    p2::Planet,
    leave::Array{Real,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      plt_title::String="Porkchop Plot",
      costfxn::Function=emp,
      levels::Dict{String,Array{Array{Real,1},1}}=
      Dict("v_inf" => [[13., 0.5, 15.],[15., 1., 21.]],
      "C3" => [[120., 5., 160.],[160., 10., 300.]]))

    xs = imd.JD_to_str.(collect(range(leave[1],leave[2],length=N)))
    ys = imd.JD_to_str.(collect(range(arrive[1],arrive[2],length=N)))

    if costfxn != emp
      outputs = reshape(
      porkchop(p1,p2,leave,arrive,outs=outs,N=N,costfxn=costfxn),(1,N,N))
      outs = ["cost"]
    else
      outputs = porkchop(p1,p2,leave,arrive,outs=outs,N=N)
    end

    data = Array{PlotlyBase.AbstractTrace,1}(undef,0)
    colors = :red,:blue,:green,:black,:purple,:orange
    for i = 1:length(outs)
      for level in levels[outs[i]]
        push!(data,
        contour(x=xs,
        y=ys,
        z=outputs[i,:,:],
        name=outs[i],
        contours=attr(
        start=level[1],
        size=level[2],
        coloring="none",
        showlabels=true),
        contours_end=level[3],
        line=attr(color=colors[i],width=3),
        showscale=false))
      end
    end

    layout = Layout(
    title=attr(text=plt_title,font_size=25),
    margin=attr(l=130,b=100,t=100),
    xaxis=attr(title="Date Departure",showline=true,tickformat="%Y-%m-%d"),
    yaxis=attr(title="Date Arrival",showline=true,tickformat="%Y-%m-%d"))

    display(plot(data, layout))
    return outputs

  end

  function porkchop_plot(
    p1::Planet,
    p2::Planet,
    leave::Array{Int64,1},
    arrive::Array{Int64,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      plt_title::String="Porkchop Plot",
      costfxn::Function=emp,
      levels::Dict{String,Array{Array{Real,1},1}}=
      Dict("v_inf" => [[13., 0.5, 15.],[15., 1., 21.]],
      "C3" => [[120., 5., 160.],[160., 10., 300.]]))

      leave = Real.(leave) ; arrive = Real.(arrive)

    if costfxn != emp
      return porkchop_plot(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title,
      costfxn=costfxn)
    else
      return porkchop_plot(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title)
    end

  end

  function porkchop_plot(
    p1_name::String,
    p2_name::String,
    leave::Array{Real,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      plt_title::String="Porkchop Plot",
      costfxn::Function=emp,
      levels::Dict{String,Array{Array{Real,1},1}}=
      Dict("v_inf" => [[13., 0.5, 15.],[15., 1., 21.]],
      "C3" => [[120., 5., 160.],[160., 10., 300.]]))

    p1 = Planet(p1_name)
    p2 = Planet(p2_name)

    if costfxn != emp
      return porkchop_plot(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title,
      costfxn=costfxn)
    else
      return porkchop_plot(
      p1,
      p2,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title)
    end

  end

  function porkchop_plot(
    p1_name::String,
    p2_name::String,
    leave::Array{Number,1},
    arrive::Array{Real,1};
      outs::Array{String,1}=["v_inf","C3"],
      M::Real=100,
      N::Real=100,
      plt_title::String="Porkchop Plot",
      costfxn::Function=emp,
      levels::Dict{String,Array{Array{Real,1},1}}=
      Dict("v_inf" => [[13., 0.5, 15.],[15., 1., 21.]],
      "C3" => [[120., 5., 160.],[160., 10., 300.]]))

    leave = Real.(leave) ; arrive = Real.(arrive)

    if costfxn != emp
      return porkchop_plot(
      p1_name,
      p2_name,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title,
      costfxn=costfxn)
    else
      return porkchop_plot(
      p1_name,
      p2_name,
      leave,
      arrive;
      outs=outs,
      M=M,
      N=N,
      levels=levels,
      plt_title=plt_title)
    end

  end

end
