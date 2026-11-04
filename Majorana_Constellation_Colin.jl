
# ######################################
# Majorana Constellation - Pure states #
########################################

using QuantumOptics 
using LinearAlgebra 
using Polynomials 
using GLMakie 
using Makie 
using CairoMakie 
using Colors 
using WignerSymbols 
using Cubature 

GLMakie.activate!()

function RotationArbitrary(n, θ)


    # n = axis of rotation 
    # θ = angle of rotation 

    M = Matrix{Float64}(undef, 3,3)
    M[1,1] = cos(θ) + n[1]^2*(1-cos(θ))
    M[2,2] = cos(θ) + n[2]^2*(1-cos(θ))
    M[3,3] = cos(θ) + n[3]^2*(1-cos(θ))

    M[1,2] = n[1]*n[2]*(1-cos(θ)) - n[3]*sin(θ) 
    M[1,3] = n[1]*n[3]*(1-cos(θ)) + n[2]*sin(θ) 

    M[2,1] = n[1]*n[2]*(1-cos(θ)) + n[3]*sin(θ) 
    M[2,3] = n[3]*n[2]*(1-cos(θ)) - n[1]*sin(θ) 

    M[3,1] = n[1]*n[3]*(1-cos(θ)) - n[2]*sin(θ) 
    M[3,2] = n[3]*n[2]*(1-cos(θ)) + n[1]*sin(θ) 

    return M
end

function Get_Roots(ψ::Ket, b::SpinBasis, jj::Float64)

    Maj_Poly = Vector{ComplexF64}(undef, floor(Int,2*jj+1))

    for i in 1:floor(Int,2*jj+1)

        m = (i-1) - jj

        state_m = []
        for n in 1:i-1 
            push!(state_m, 0.0 + im*0.0)
        end
        push!(state_m,1.0 + im*0.0)
        for n in i+1:floor(Int, 2*jj+1)
            push!(state_m,0.0 + im*0.0)
        end

        m_j = Ket(b,state_m)
        Maj_Poly[i] = (-1)^(jj-m)*sqrt(factorial(floor(Int,2*jj))/(factorial(floor(Int,jj - m))*factorial(floor(Int,2*jj - (jj-m)))))*ψ.data[i]#dagger(m_j)*ψ.data[i]

    end

    Roots = ComplexF64.(roots(Polynomial(Maj_Poly)))

    if length(Roots) < floor(Int, 2*jj)
        lr = floor(Int, 2*jj)-length(Roots)
        for n in 1:lr
            push!(Roots, Inf + im*Inf)
        end
    end
    #=
    for n in 1:length(Roots)
        if abs(Roots[n]) < 10^(-7)
            Roots[n] = 0.0 + im*0.0
        end
    end
    =#
    #return Polynomial(Maj_Poly)
    return Roots

end


function Stereographic_Projection(Roots::Vector)

    θ = Vector{Float64}(undef, length(Roots))
    ϕ = Vector{Float64}(undef, length(Roots))

    for i in 1:length(Roots)

        x = real(Roots[i])
        y = real(-im*Roots[i])

        if abs(x) < 10^(-12) && abs(y) < 10^(-12) 
            θ[i] = 0.0
            ϕ[i] = 0.0
        elseif abs(x) < 10^(-12) && abs(y) > 10^(-12)
            ϕ[i] = asin(y/abs(y))
            θ[i] = 2*atan(abs(y))
        elseif abs(x) > 10^(-12) && abs(y) < 10^(-12)
            ϕ[i] = acos(x/abs(x))
            θ[i] = 2*atan(abs(x))
        elseif abs(Roots[i]) > 10^10
            #println("inf here")
            θ[i] = pi 
            ϕ[i] = 0.0
        else

            if y > 0 && x > 0
                ϕ[i] = acos(abs(x)/sqrt(x^2 + y^2))
                θ[i] = 2*atan(sqrt(x^2 + y^2))
            elseif y > 0 && x < 0
                ϕ[i] = pi - acos(abs(x)/sqrt(x^2 + y^2))
                θ[i] = 2*atan(sqrt(x^2 + y^2))
            elseif y < 0 && x > 0
                ϕ[i] = - acos(abs(x)/sqrt(x^2 + y^2))
                θ[i] = 2*atan(sqrt(x^2 + y^2))
            else
                ϕ[i] = pi + acos(abs(x)/sqrt(x^2 + y^2))
                θ[i] = 2*atan(sqrt(x^2 + y^2))
            end

        end
    end

    return θ,ϕ
end

function Find_Degeneracy(roots::Vector)

    new_roots = []
    degen = []

    for i in 1:length(roots)

        if abs(roots[i]) < 10^(6)
            isthere = false 
            ind = 0
            for n in 1:length(new_roots)
                if abs(real(new_roots[n]) - real(roots[i])) < 10^(-2) && abs(real(-im*new_roots[n]) - real(-im*roots[i])) < 10^(-2)
                    isthere = true 
                    ind = n
                end
            end

            if isthere == false 
                push!(new_roots, roots[i])
                push!(degen, 1)
            elseif isthere == true 
                degen[ind] += 1
            end
        else 
            isthere = false 
            ind = 0
            for n in 1:length(new_roots)
                if abs(new_roots[n]) > 10^(6)
                    isthere = true 
                    ind = n
                end
            end

            if isthere == false 
                push!(new_roots, roots[i])
                push!(degen, 1)
            elseif isthere == true 
                degen[ind] += 1
            end
        end

    end

    return new_roots, degen

end

# Majorana constellation for pure state 

function Majorana_Constellationψ(ψ::Ket,b::SpinBasis, jj::Float64;color::Colorant=RGB(255/255,204/255,153/255), ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing)

    roots = Get_Roots(ψ,b,jj)

    println(roots)
    roots, roots_degen = Find_Degeneracy(roots)
    println(roots)
    θ,ϕ = Stereographic_Projection(roots)


    r = 1.0 

    θs = LinRange(0,pi,100)
    ϕs = LinRange(ϕmin,ϕmax, 100)

    xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
    ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
    zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

    m = Matrix{Colorant}(undef, 100,100)
    m .= color

    fig = Figure()
    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz)

    surface!(xs,ys,zs,transparency = true,color=m, alpha=0.7)
    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35, 1.35,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35, 1.35,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35, 1.35,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45,text=L"z",fontsize=25.0f0)
    text!(ax,0.0,1.45,0.0,text=L"y",fontsize=25.0f0)
    text!(ax,1.45,0.0,0.0,text=L"x",fontsize=25.0f0, align=(:right,:bottom))

    for n in 1:length(roots)
        r1 = 1.02
        scatter!([r1*sin(θ[n])*cos(ϕ[n])] , -[r1*sin(θ[n])*sin(ϕ[n])] ,-[r1*cos(θ[n])], markersize = 14.0, color=:black)
        lines!([0,r1*sin(θ[n])*cos(ϕ[n])] , -[0,r1*sin(θ[n])*sin(ϕ[n])] ,-[0,r1*cos(θ[n])], linestyle=:dash, color=:black, linewidth=1.0)
        if roots_degen[n] > 1 
            rp = 1.25
            text!([rp*sin(θ[n])*cos(ϕ[n])] , -[rp*sin(θ[n])*sin(ϕ[n])] ,-[rp*cos(θ[n])], text = "$(roots_degen[n])")
        end
        
    end

    hidespines!(ax)
    hidedecorations!(ax)
    display(fig)

    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
        println("saved !")
        display(fig)
    end

    return fig 
end

# Maj constellation for pure state, but the sphere is a heatmap 

function Majorana_Constellationψ_coloredsphere(ψ::Ket,b::SpinBasis, jj::Float64, colordata::Matrix, φs,θs; center = [-1,-1,-1], alpha=0.7,colormap=:blackbody, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing,depth_proj = 1.8)

    roots = Get_Roots(ψ,b,jj)

    println(roots)
    roots, roots_degen = Find_Degeneracy(roots)
    println(roots)
    θ,ϕ = Stereographic_Projection(roots)


    r = 1.0 


    xs = [r*sin(θs[n])*cos(φs[m]) for m in 1:length(φs), n in 1:length(θs)]
    ys = [r*sin(θs[n])*sin(φs[m]) for m in 1:length(φs), n in 1:length(θs)]
    zs = [r*cos(θs[n]) for m in 1:length(φs), n in 1:length(θs)]

    fig = Figure()
    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz)

    surface!(xs,ys,zs,transparency = true,color=colordata,colormap=colormap, alpha=alpha)


    surface!(xs[1:floor(Int, length(φs)/2),:],[-depth_proj for n in 1:floor(Int, length(φs)/2), m in 1:length(θs)],zs[1:floor(Int, length(φs)/2),:], color = colordata[1:floor(Int, length(φs)/2),:], colormap =:blackbody)
    surface!([-depth_proj for n in floor(Int, length(φs)/4):floor(Int, length(φs)/4)+floor(Int, length(φs)/2), m in 1:length(θs)],ys[floor(Int, length(φs)/4):floor(Int, length(φs)/4)+floor(Int, length(φs)/2),:],zs[floor(Int, length(φs)/4):floor(Int, length(φs)/4)+floor(Int, length(φs)/2),:], color = colordata[floor(Int, length(φs)/4):floor(Int, length(φs)/4)+floor(Int, length(φs)/2),:], colormap =:blackbody)
    surface!(xs[:,1:floor(Int, length(θs)/2)],ys[:,1:floor(Int, length(θs)/2)],[-depth_proj for n in 1:length(φs), m in 1:floor(Int, length(θs)/2)], color = colordata[:,1:floor(Int, length(θs)/2)], colormap =:blackbody)

    lines!(Point3.([[1.2,1.2,-depth_proj],[1.2,-1.2,-depth_proj],[-1.2,-1.2,-depth_proj],[-1.2,1.2,-depth_proj],[1.2,1.2,-depth_proj]]),color=:grey)
    lines!(Point3.([[1.2,-depth_proj,1.2],[1.2,-depth_proj,-1.2],[-1.2,-depth_proj,-1.2],[-1.2,-depth_proj,1.2],[1.2,-depth_proj,1.2]]),color=:grey)
    lines!(Point3.([[-depth_proj,1.2,1.2],[-depth_proj,1.2,-1.2],[-depth_proj,-1.2,-1.2],[-depth_proj,-1.2,1.2],[-depth_proj,1.2,1.2]]),color=:grey)

    #=
    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35, 1.35,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35, 1.35,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35, 1.35,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45,text=L"z",fontsize=25.0f0)
    text!(ax,0.0,1.45,0.0,text=L"y",fontsize=25.0f0)
    text!(ax,1.45,0.0,0.0,text=L"x",fontsize=25.0f0, align=(:right,:bottom))
    =#

    xa = 0.7*[1,0,0]
    ya = 0.7*[0,1,0]
    za = 0.7*[0,0,1]

    arrows3d!(Point3.([center,center,center]), Point3.([xa,ya,za]),color=:red, lengthscale=0.8)
    text!(Point3.([center+xa, center+ya, center+za]), text=[L"x", L"y", L"z"], align = [(:center,:center),(:center,:center),(:center,:center)], fontsize=25.0f0)

    for n in 1:length(roots)
        r1 = 1.02
        scatter!([r1*sin(θ[n])*cos(ϕ[n])] , [r1*sin(θ[n])*sin(ϕ[n])] ,[r1*cos(θ[n])], markersize = 14.0, color=:black)

        #lines!([0,r1*sin(θ[n])*cos(ϕ[n])] , -[0,r1*sin(θ[n])*sin(ϕ[n])] ,-[0,r1*cos(θ[n])], linestyle=:dash, color=:black, linewidth=1.0)
        if roots_degen[n] > 1 
            rp = 1.25
            text!([rp*sin(θ[n])*cos(ϕ[n])] , [rp*sin(θ[n])*sin(ϕ[n])] ,[rp*cos(θ[n])], text = "$(roots_degen[n])")
        end

        if θ[n] <= pi/2
            scatter!([sin(θ[n])*cos(ϕ[n])] , [sin(θ[n])*sin(ϕ[n])] ,[-depth_proj+0.01], markersize = 10.0, color=:black)
        end 

        if ϕ[n] >= pi/2 && ϕ[n] <= 3*pi/2
            scatter!([-depth_proj+0.01] , [sin(θ[n])*sin(ϕ[n])] ,[cos(θ[n])], markersize = 10.0, color=:black)
        end 

        if ϕ[n] <= pi 
            scatter!([sin(θ[n])*cos(ϕ[n])] , [-depth_proj+0.01] ,[cos(θ[n])], markersize = 10.0, color=:black)
        end 
    end

    hidespines!(ax)
    hidedecorations!(ax)
    display(fig)

    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
        println("saved !")
        display(fig)
    end

    return fig, ax
end


#=
jj = 2.0
b = SpinBasis(jj)
ψ = zeros(floor(Int, 2*jj+1))
ψ[1] = 1.0 
ψ[end] = 1.0
ψ = ψ/norm(ψ)
U = exp(-im*2*pi/3*sqrt(1/3)*(sigmax(b)/2 + sigmay(b)/2 + sigmaz(b)/2 ))
ψ = U*Ket(b,ψ)
color = RGB(255/255,204/255,153/255);
data_dir = "/Users/Colin/Documents/PHD/Codes/Figures/Majorana Constellation/"
filename = "/GHZ state j=2 rotated"
#Majorana_Constellationψ(ψ,b,jj,color; ϕaz = -1.0*pi/8, path_dir = data_dir, filename = filename)
Majorana_Constellationψ(ψ,b,jj,color=color, ϕaz = -1.0*pi/8)
=#
# #######################################
# Majorana Constellation - Mixed states #
#########################################

function Tlmdk(N::Int, L::Int, M::Int)

    Tlm = Array{Complex{Float64},2}(undef, N+1, N+1)

    for k in 0:N 
        for kp in 0:N 
            if N/2-k in -N/2:N/2 && N/2-kp in -N/2:N/2
                Tlm[k+1,kp+1] = sqrt((2*L+1)/(N+1))*clebschgordan(N/2,N/2-kp,L,M,N/2,N/2-k)
            else 
                Tlm[k+1,kp+1] = 0
            end
        end
    end

    return Tlm
end

function Tlm_Decomposition_h(H::Operator,b::SpinBasis, jj::Float64)
    
    h_l = Vector{Vector}(undef, floor(Int,2*jj+1))

    for L in 0:floor(Int,2*jj)

        h_l[L+1] = Vector{ComplexF64}(undef, 2*L+1)

        for M in -L:1:L

            Tlm = DenseOperator(b,Tlmdk(floor(Int,2*jj), L, M))
            h_l[L+1][M+(L+1)] = tr(dagger(Tlm)*H)
        end
    end

    return h_l
end

function Tlm_Decomposition_nh(H::Operator, b::SpinBasis, jj::Float64)

    H_S = (H + dagger(H))/2 
    H_AS = (H - dagger(H))/(2*im) 

    h_Sl = Vector{Vector}(undef, floor(Int,2*jj+1))

    for L in 0:floor(Int,2*jj)
        h_Sl[L+1] = Vector{ComplexF64}(undef, 2*L+1)
        for M in -L:1:L
            Tlm = DenseOperator(b,Tlmdk(floor(Int,2*jj), L, M))
            h_Sl[L+1][M+(L+1)] = tr(dagger(Tlm)*H_S)
        end
    end

    h_ASl = Vector{Vector}(undef, floor(Int,2*jj+1))

    for L in 0:floor(Int,2*jj)
        h_ASl[L+1] = Vector{ComplexF64}(undef, 2*L+1)
        for M in -L:1:L
            Tlm = DenseOperator(b,Tlmdk(floor(Int,2*jj), L, M))
            h_ASl[L+1][M+(L+1)] = tr(dagger(Tlm)*H_AS)
        end
    end

    return h_Sl, h_ASl

end

# transforms the roots to a quantum state 

function ConstToState(Roots::Vector)

    jj = length(Roots)/2 
    #println("test")
    #println(Roots)
    Roots_NoInf = Vector{ComplexF64}(undef, 0)
    for i in 1:length(Roots)
        if abs(Roots[i]) > 10^10

        else 
            push!(Roots_NoInf, Roots[i])
        end
    end
    #println(Roots_NoInf)
    #println(Roots)

    Poly = fromroots(Roots_NoInf)
    #println(Poly)
    #println([Poly[i-1] for i in 1:trunc(Int,2*jj+1)])
    ψ_vec = [Poly[i-1] for i in 1:trunc(Int,2*jj+1)]

    for i in 1:length(ψ_vec)
        mi = -jj + (i-1)
        Factor = (-1)^(jj - mi)*sqrt(binomial(trunc(Int, 2*jj), trunc(Int, jj - mi)))
        #println(Factor)
        ψ_vec[i] = ψ_vec[i]/Factor
    end
    ψ_vec = ψ_vec/norm(ψ_vec)


    return Ket(SpinBasis(jj),ψ_vec)
end

# colors the points 

function Colorings(L::Float64, roots::Vector, hl::Vector)

    θ,ϕ = Stereographic_Projection(roots)

    ns = Vector{Vector}(undef, floor(Int,2*L)) # all stars

    for n in 1:floor(Int,2*L)
        ns[n] = [sin(θ[n])*cos(ϕ[n]) , sin(θ[n])*sin(ϕ[n]) , cos(θ[n])]
    end

    Rs = copy(roots) # all roots

    #println(length(Rs))

    # two sets of constellation where σs[n] is antipodal to σsA[n]

    σs = Vector{Vector}(undef, floor(Int,L)) 
    Rσ = Vector{ComplexF64}(undef, floor(Int,L)) 

    σsA = Vector{Vector}(undef, floor(Int,L))
    RσA = Vector{ComplexF64}(undef, floor(Int,L)) 

    #println(ns)
    #println(Rs)

    for n in 1:floor(Int,L)

        rand_index = rand(1:length(ns))
        σs[n] = ns[rand_index]
        Rσ[n] = Rs[rand_index]

        deleteat!(ns, rand_index)
        deleteat!(Rs, rand_index)

        index_A = 0

        for m in 1:length(ns)

            if norm(σs[n] + ns[m]) < 10^(-4)
                index_A = m
            end

        end

        σsA[n] = ns[index_A]
        RσA[n] = Rs[index_A]
        deleteat!(ns, index_A)
        deleteat!(Rs, index_A)

    end

    # now we compute the state of each constellation 

    ψ_σ = ConstToState(Rσ)

    ψ_σ_L = [dagger(ψ_σ)*dagger(DenseOperator(SpinBasis(L/2),Tlmdk(floor(Int,L),floor(Int,L),floor(Int,-L + (n-1)))))*ψ_σ for n in 1:trunc(Int,2*L+1)]
    stop = false 
    n = 1 
    sign1 = 0

    while stop == false 

        Tlm = DenseOperator(SpinBasis(L/2),Tlmdk(floor(Int,L),floor(Int,L),floor(Int,-L + (n-1))))
        x = dagger(ψ_σ)*dagger(Tlm)*ψ_σ
        #println(x)
        if abs(x) > 10^(-6) 
            stop = true 
            sign1 = x
        else 
            n+=1
        end
    end

    sign2 = hl[n]
    #println("sign1 = $sign1")
    #println("hl = $hl")
    #println("ψl = $ψ_σ")
    #println("ψl exp = $ψ_σ_L")  
    if real(sign1/sign2) > 0
        return Rσ, RσA
    else 

        tamp = copy(Rσ[1])
        Rσ[1] = RσA[1]
        RσA[1] = tamp 
        #println("changed")
        return Rσ, RσA
    end
end

# Majorana constellation for Hermitian operator 

function Majorana_ConstellationH_h(H::Operator, b::SpinBasis, jj::Float64; ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing, symmetry=false,pos_star = true)
    h_l = Tlm_Decomposition_h(H,b, jj)
    radius_l = zeros(length(h_l))
    index = []

    for L in 1:floor(Int,2*jj)
        if norm(h_l[L+1]) > 10^(-5)
            radius_l[L+1] = norm(h_l[L+1]) 
            h_l[L+1] = h_l[L+1]/radius_l[L+1]
            push!(index, L+1)
        end
    end

    # index -> all L such that the constellation has a non-zero radius 
   
    colors_L = [RGB(102/255,102/255,255/255),RGB(255/255,178/255,102/255),RGB(255/255,102/255,102/255), RGB(102/255,255/255,102/255), RGB(255/255,102/255,255/255),RGB(192/255,192/255,192/255),RGB(102/255,255/255,255/255)]
    Rad_Max = maximum(radius_l)
    θs = LinRange(0,pi,100)
    ϕs = LinRange(ϕmin,ϕmax, 100)

    fig = Figure()
    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz)

    r = length(index)

    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35*r, 1.35*r,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45*r,text=L"z")
    text!(ax,0.0,1.45*r,0.0,text=L"y")
    text!(ax,1.45*r,0.0,0.0,text=L"x")

    for n in 1:length(index)

        S = index[n]-1
        r = n*1.0

        xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
        ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
        zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

        m = Matrix{Colorant}(undef, 100,100)
        if S <= length(colors_L)
            m .= colors_L[S]
        else 
            m.= RGB(102/255,102/255,255/255)
        end

        surface!(xs,ys,zs,transparency = true,color=m, alpha=0.5)

        js = 1.0*S
        base_l = SpinBasis(js)
        ψ_l = Ket(base_l, h_l[index[n]])

        roots = Get_Roots(ψ_l,base_l,js)

        Reds, Blacks = Colorings(js,roots, h_l[index[n]])

        Reds, Reds_degen = Find_Degeneracy(Reds)
        Blacks, Blacks_degen = Find_Degeneracy(Blacks)

        θreds, ϕreds = Stereographic_Projection(Reds)
        θblacks, ϕblacks = Stereographic_Projection(Blacks)

        for n in 1:length(Reds)

            r1 = r + 0.01
            scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=:red)
            if pos_star == true 
                lines!([0,r1*sin(θreds[n])*cos(ϕreds[n])] , [0,r1*sin(θreds[n])*sin(ϕreds[n])] ,[0,r1*cos(θreds[n])], linestyle=:dash, color=:red, linewidth=1.0)
            end
            if Reds_degen[n] > 1 
                rp = r + 0.1
                #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
            end 
        end

        for n in 1:length(Blacks)
            r1 = r + 0.01
            scatter!([r1*sin(θblacks[n])*cos(ϕblacks[n])] , [r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[r1*cos(θblacks[n])], markersize = 10.0, color=:black)
            if pos_star == true 
                lines!([0,r1*sin(θblacks[n])*cos(ϕblacks[n])] , [0,r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[0,r1*cos(θblacks[n])], linestyle=:dash, color=:black, linewidth=1.0)
            end
            if Blacks_degen[n] > 1 
                rp = r
                text!([rp*sin(θblacks[n])*cos(ϕblacks[n])] , [rp*sin(θblacks[n])*sin(ϕblacks[n])] ,[rp*cos(θblacks[n]) + 0.1], text = "$(Blacks_degen[n])")
            end 
        end

        if symmetry == true 

            θall = vcat(θblacks,θreds)
            ϕall = vcat(ϕblacks,ϕreds)
            r1 = r 
            for i in 1:length(θall)-1
                for j in i+1:length(ϕall)

                    lines!([r1*sin(θall[i])*cos(ϕall[i]),r1*sin(θall[j])*cos(ϕall[j])] , [r1*sin(θall[i])*sin(ϕall[i]),r1*sin(θall[j])*sin(ϕall[j])] ,[r1*cos(θall[i]),r1*cos(θall[j])], linestyle=:solid, color=:black, linewidth=1.0)

                end
            end

        end

    end
    hidedecorations!(ax)
    hidespines!(ax)

    display(fig)
    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
    end
end

# Majorana constellation for non-Hermitian operator 

function Majorana_ConstellationH_nh(H::Operator, b::SpinBasis, jj::Float64; ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing)
    h_l, h_l2 = Tlm_Decomposition_nh(H,b, jj)

    radius_l = zeros(length(h_l))
    radius_l2 = zeros(length(h_l2))
    index = []
    index2 = []

    for L in 1:floor(Int,2*jj)
        if norm(h_l[L+1]) > 10^(-5)
            radius_l[L+1] = norm(h_l[L+1]) 
            h_l[L+1] = h_l[L+1]/radius_l[L+1]
            push!(index, L+1)
        end
        if norm(h_l2[L+1]) > 10^(-5)
            radius_l2[L+1] = norm(h_l2[L+1]) 
            h_l2[L+1] = h_l2[L+1]/radius_l2[L+1]
            push!(index2, L+1)
        end
    end

    # index -> all L such that the constellation has a non-zero radius 
   
    colors_L = [RGB(102/255,102/255,255/255),RGB(255/255,178/255,102/255),RGB(255/255,102/255,102/255), RGB(102/255,255/255,102/255), RGB(255/255,102/255,255/255),RGB(192/255,192/255,192/255),RGB(102/255,255/255,255/255)]
    Rad_Max = maximum(radius_l)
    θs = LinRange(0,pi,100)
    ϕs = LinRange(ϕmin,ϕmax, 100)

    fig = Figure()


    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz, title=L"\frac{H + H^{\dagger}}{2}",titlegap =-100,protrusions=(0,0,0,0))

    r = length(index)

    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35*r, 1.35*r,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45*r,text=L"z")
    text!(ax,0.0,1.45*r,0.0,text=L"y")
    text!(ax,1.45*r,0.0,0.0,text=L"x")

    for n in 1:length(index)

        S = index[n]-1
        r = n*1.0

        xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
        ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
        zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

        m = Matrix{Colorant}(undef, 100,100)
        m .= colors_L[S]

        surface!(xs,ys,zs,transparency = true,color=m, alpha=0.7)

        js = 1.0*S
        base_l = SpinBasis(js)
        ψ_l = Ket(base_l, h_l[index[n]])

        roots = Get_Roots(ψ_l,base_l,js)

        Reds, Blacks = Colorings(js,roots, h_l[index[n]])

        Reds, Reds_degen = Find_Degeneracy(Reds)
        Blacks, Blacks_degen = Find_Degeneracy(Blacks)

        θreds, ϕreds = Stereographic_Projection(Reds)
        θblacks, ϕblacks = Stereographic_Projection(Blacks)

        for n in 1:length(Reds)

            r1 = r + 0.01
            scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=:red)
            
            if Reds_degen[n] > 1 
                rp = r + 0.1
                #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
            end 
        end

        for n in 1:length(Blacks)
            r1 = r + 0.01
            scatter!([r1*sin(θblacks[n])*cos(ϕblacks[n])] , [r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[r1*cos(θblacks[n])], markersize = 10.0, color=:black)
            if Blacks_degen[n] > 1 
                rp = r
                text!([rp*sin(θblacks[n])*cos(ϕblacks[n])] , [rp*sin(θblacks[n])*sin(ϕblacks[n])] ,[rp*cos(θblacks[n]) + 0.1], text = "$(Blacks_degen[n])")
            end 
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    ax = Axis3(fig[1,2],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz, title=L"\frac{H - H^{\dagger}}{2i}",titlegap = -100, protrusions=(0,0,0,0))

    r = length(index2)

    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35*r, 1.35*r,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45*r,text=L"z")
    text!(ax,0.0,1.45*r,0.0,text=L"y")
    text!(ax,1.45*r,0.0,0.0,text=L"x")

    for n in 1:length(index2)

        S = index2[n]-1
        r = n*1.0

        xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
        ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
        zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

        m = Matrix{Colorant}(undef, 100,100)
        m .= colors_L[S]

        surface!(xs,ys,zs,transparency = true,color=m, alpha=0.7)

        js = 1.0*S
        base_l = SpinBasis(js)
        ψ_l = Ket(base_l, h_l2[index2[n]])

        roots = Get_Roots(ψ_l,base_l,js)

        Reds, Blacks = Colorings(js,roots, h_l2[index2[n]])

        Reds, Reds_degen = Find_Degeneracy(Reds)
        Blacks, Blacks_degen = Find_Degeneracy(Blacks)

        θreds, ϕreds = Stereographic_Projection(Reds)
        θblacks, ϕblacks = Stereographic_Projection(Blacks)

        for n in 1:length(Reds)

            r1 = r + 0.01
            scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=:red)
            
            if Reds_degen[n] > 1 
                rp = r + 0.1
                #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
            end 
        end

        for n in 1:length(Blacks)
            r1 = r + 0.01
            scatter!([r1*sin(θblacks[n])*cos(ϕblacks[n])] , [r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[r1*cos(θblacks[n])], markersize = 10.0, color=:black)
            if Blacks_degen[n] > 1 
                rp = r
                text!([rp*sin(θblacks[n])*cos(ϕblacks[n])] , [rp*sin(θblacks[n])*sin(ϕblacks[n])] ,[rp*cos(θblacks[n]) + 0.1], text = "$(Blacks_degen[n])")
            end 
        end
    end
    hidedecorations!(ax)
    hidespines!(ax)


    display(fig)
    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
    end
end

# Majorana constellation for operator 

function Majorana_Constellation_H(H::Operator, b::SpinBasis, jj::Float64; ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing, symmetry=false, pos_star = true)
    
    HAS = (H-dagger(H))/(2*im)
    if norm(HAS) < 10^(-6)
        Majorana_ConstellationH_h(H,b, jj, ϕmin = ϕmin, ϕmax = ϕmax, ϕaz = ϕaz,path_dir=path_dir, filename = filename, symmetry = symmetry, pos_star = pos_star)
    else 
        Majorana_ConstellationH_nh(H,b, jj; ϕmin = ϕmin, ϕmax = ϕmax, ϕaz = ϕaz,path_dir=path_dir, filename = filename)
    end
end

# Majorana constellation of hermitian operator, but returns the figure 

function Majorana_ConstellationH_h_AddToFigure(H::Operator, b::SpinBasis, jj::Float64; ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing, symmetry=false,pos_star = true, xyzsize = 16.0f0)
    h_l = Tlm_Decomposition_h(H,b, jj)
    radius_l = zeros(length(h_l))
    index = []

    for L in 1:floor(Int,2*jj)
        if norm(h_l[L+1]) > 10^(-5)
            radius_l[L+1] = norm(h_l[L+1]) 
            h_l[L+1] = h_l[L+1]/radius_l[L+1]
            push!(index, L+1)
        end
    end

    # index -> all L such that the constellation has a non-zero radius 
   
    colors_L = [RGB(102/255,102/255,255/255),RGB(255/255,178/255,102/255),RGB(255/255,102/255,102/255), RGB(102/255,255/255,102/255), RGB(255/255,102/255,255/255),RGB(192/255,192/255,192/255),RGB(102/255,255/255,255/255)]
    Rad_Max = maximum(radius_l)
    θs = LinRange(0,pi,100)
    ϕs = LinRange(ϕmin,ϕmax, 100)

    fig = Figure()
    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz)

    r = length(index)

    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35*r, 1.35*r,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45*r,text=L"z", fontsize = xyzsize)
    text!(ax,0.0,1.45*r,0.0,text=L"y", fontsize = xyzsize)
    text!(ax,1.45*r,0.0,0.0,text=L"x", fontsize = xyzsize)

    for n in 1:length(index)

        S = index[n]-1
        r = n*1.0

        xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
        ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
        zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

        m = Matrix{Colorant}(undef, 100,100)
        m .= colors_L[S]

        surface!(xs,ys,zs,transparency = true,color=m, alpha=0.5)

        js = 1.0*S
        base_l = SpinBasis(js)
        ψ_l = Ket(base_l, h_l[index[n]])

        roots = Get_Roots(ψ_l,base_l,js)

        Reds, Blacks = Colorings(js,roots, h_l[index[n]])

        Reds, Reds_degen = Find_Degeneracy(Reds)
        Blacks, Blacks_degen = Find_Degeneracy(Blacks)

        θreds, ϕreds = Stereographic_Projection(Reds)
        θblacks, ϕblacks = Stereographic_Projection(Blacks)

        for n in 1:length(Reds)

            r1 = r + 0.01
            scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=:red)
            if pos_star == true 
                lines!([0,r1*sin(θreds[n])*cos(ϕreds[n])] , [0,r1*sin(θreds[n])*sin(ϕreds[n])] ,[0,r1*cos(θreds[n])], linestyle=:dash, color=:red, linewidth=1.0)
            end
            if Reds_degen[n] > 1 
                rp = r + 0.1
                #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
            end 
        end

        for n in 1:length(Blacks)
            r1 = r + 0.01
            scatter!([r1*sin(θblacks[n])*cos(ϕblacks[n])] , [r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[r1*cos(θblacks[n])], markersize = 10.0, color=:black)
            if pos_star == true 
                lines!([0,r1*sin(θblacks[n])*cos(ϕblacks[n])] , [0,r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[0,r1*cos(θblacks[n])], linestyle=:dash, color=:black, linewidth=1.0)
            end
            if Blacks_degen[n] > 1 
                rp = r
                text!([rp*sin(θblacks[n])*cos(ϕblacks[n])] , [rp*sin(θblacks[n])*sin(ϕblacks[n])] ,[rp*cos(θblacks[n]) + 0.1], text = "$(Blacks_degen[n])")
            end 
        end

        if symmetry == true 

            θall = vcat(θblacks,θreds)
            ϕall = vcat(ϕblacks,ϕreds)
            r1 = r 
            for i in 1:length(θall)-1
                for j in i+1:length(ϕall)

                    lines!([r1*sin(θall[i])*cos(ϕall[i]),r1*sin(θall[j])*cos(ϕall[j])] , [r1*sin(θall[i])*sin(ϕall[i]),r1*sin(θall[j])*sin(ϕall[j])] ,[r1*cos(θall[i]),r1*cos(θall[j])], linestyle=:solid, color=:black, linewidth=1.0)

                end
            end

        end

    end
    hidedecorations!(ax)
    hidespines!(ax)

    display(fig)
    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
    end
    return fig 
end

# Majorana constellation for different Hamiltonian, on the same sphere. e.g., to visualize a small unitary transformation 

function Majorana_ConstellationH_h_MultipleHamiltonian(Hs::Vector, αs::Vector, b::SpinBasis, jj::Float64; ϕmin::Float64=pi/2, ϕmax::Float64=2*pi, ϕaz::Float64=0.0, path_dir=nothing, filename=nothing, symmetry=false,pos_star = true, Num_Degen = true, all_colors = false, all_symmetry = false)
    
    
    
    h_l = [Tlm_Decomposition_h(Hs[n],b, jj) for n in 1:length(Hs)]
    radius_l = [zeros(length(h_l[n])) for n in 1:length(Hs)]
    index = [[] for n in 1:length(Hs)]
    for n in 1:length(Hs)
        for L in 1:floor(Int,2*jj)
            if norm(h_l[n][L+1]) > 10^(-5)
                radius_l[n][L+1] = norm(h_l[n][L+1]) 
                h_l[n][L+1] = h_l[n][L+1]/radius_l[n][L+1]
                push!(index[n], L+1)
            end
        end
    end
    
    for n in 1:length(Hs)-1 
        for m in 2:length(Hs)
            if length(index[n]) != length(index[m])
                println("Not the same number of constellations !")
                return "Error"
            else 
                for l in 1:length(index[n])
                    if index[n][l] != index[m][l]
                        println("Not the same irreps decompostion !")
                        return "Error"
                    end
                end
            end
        end
    end

    # index -> all L such that the constellation has a non-zero radius 
   
    colors_L = [RGB(102/255,102/255,255/255),RGB(255/255,178/255,102/255),RGB(255/255,102/255,102/255), RGB(102/255,255/255,102/255), RGB(255/255,102/255,255/255),RGB(192/255,192/255,192/255),RGB(102/255,255/255,255/255)]
    Rad_Max = maximum(radius_l[1])
    θs = LinRange(0,pi,100)
    ϕs = LinRange(ϕmin,ϕmax, 100)

    fig = Figure()
    ax = Axis3(fig[1,1],aspect=(1.0,1.0,1.0),azimuth = 0.25*pi + ϕaz)

    r = length(index[1])

    lines!([0.0 for n in 1:20],[0.0 for n in 1:20], LinRange(-1.35*r, 1.35*r,20), color=:black, linewidth=1.0)
    lines!([0.0 for n in 1:20],LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    lines!(LinRange(-1.35*r, 1.35*r,20),[0.0 for n in 1:20],[0.0 for n in 1:20],  color=:black, linewidth=1.0)
    text!(ax,0.0,0.0,1.45*r,text=L"z")
    text!(ax,0.0,1.45*r,0.0,text=L"y")
    text!(ax,1.45*r,0.0,0.0,text=L"x")

    for n in 1:length(index[1])

        S = index[1][n]-1
        r = n*1.0

        xs = [r*sin(θs[n])*cos(ϕs[m]) for n in 1:100, m in 1:100]
        ys = [r*sin(θs[n])*sin(ϕs[m]) for n in 1:100, m in 1:100]
        zs = [r*cos(θs[n]) for n in 1:100, m in 1:100]

        m = Matrix{Colorant}(undef, 100,100)
        m .= colors_L[S]

        surface!(xs,ys,zs,transparency = true,color=m, alpha=0.5)

        for nn in 1:length(Hs)

            js = 1.0*S
            base_l = SpinBasis(js)
            ψ_l = Ket(base_l, h_l[nn][index[1][n]])

            roots = Get_Roots(ψ_l,base_l,js)

            Reds, Blacks = Colorings(js,roots, h_l[nn][index[1][n]])

            Reds, Reds_degen = Find_Degeneracy(Reds)
            Blacks, Blacks_degen = Find_Degeneracy(Blacks)

            θreds, ϕreds = Stereographic_Projection(Reds)
            θblacks, ϕblacks = Stereographic_Projection(Blacks)

            for n in 1:length(Reds)

                if all_colors == true || nn == length(Hs)
                    r1 = r + 0.01
                    scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=(:red,αs[nn]))
                    if pos_star == true 
                        lines!([0,r1*sin(θreds[n])*cos(ϕreds[n])] , [0,r1*sin(θreds[n])*sin(ϕreds[n])] ,[0,r1*cos(θreds[n])], linestyle=:dash, color=(:red,αs[nn]), linewidth=1.0)
                    end
                    if Reds_degen[n] > 1 
                        rp = r + 0.1
                        #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
                    end 
                else
                    r1 = r + 0.01
                    scatter!([r1*sin(θreds[n])*cos(ϕreds[n])] , [r1*sin(θreds[n])*sin(ϕreds[n])] ,[r1*cos(θreds[n])], markersize = 10.0, color=(:black,αs[nn]))
                    if pos_star == true 
                        lines!([0,r1*sin(θreds[n])*cos(ϕreds[n])] , [0,r1*sin(θreds[n])*sin(ϕreds[n])] ,[0,r1*cos(θreds[n])], linestyle=:dash, color=(:black,αs[nn]), linewidth=1.0)
                    end
                    if Reds_degen[n] > 1 
                        rp = r + 0.1
                        #text!([rp*sin(θreds[n])*cos(ϕreds[n])] , [rp*sin(θreds[n])*sin(ϕreds[n])] ,[rp*cos(θreds[n])], text = "$(Reds_degen[n])")
                    end 

                end
            end

            for n in 1:length(Blacks)
                r1 = r + 0.01
                scatter!([r1*sin(θblacks[n])*cos(ϕblacks[n])] , [r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[r1*cos(θblacks[n])], markersize = 10.0, color=(:black,αs[nn]))
                if pos_star == true 
                    lines!([0,r1*sin(θblacks[n])*cos(ϕblacks[n])] , [0,r1*sin(θblacks[n])*sin(ϕblacks[n])] ,[0,r1*cos(θblacks[n])], linestyle=:dash, color=(:black,αs[nn]), linewidth=1.0)
                end
                if Blacks_degen[n] > 1 
                    rp = r
                    if Num_Degen == true 
                        text!([rp*sin(θblacks[n])*cos(ϕblacks[n])] , [rp*sin(θblacks[n])*sin(ϕblacks[n])] ,[rp*cos(θblacks[n]) + 0.1], text = "$(Blacks_degen[n])")
                    end
                end 
            end

            if symmetry == true 

                if all_symmetry == true && nn < length(Hs)

                    θall = vcat(θblacks,θreds)
                    ϕall = vcat(ϕblacks,ϕreds)
                    r1 = r 
                    for i in 1:length(θall)-1
                        for j in i+1:length(ϕall)

                            lines!([r1*sin(θall[i])*cos(ϕall[i]),r1*sin(θall[j])*cos(ϕall[j])] , [r1*sin(θall[i])*sin(ϕall[i]),r1*sin(θall[j])*sin(ϕall[j])] ,[r1*cos(θall[i]),r1*cos(θall[j])], linestyle=:solid, color=:black, linewidth=1.0)

                        end
                    end
                elseif nn == length(Hs)

                    θall = vcat(θblacks,θreds)
                    ϕall = vcat(ϕblacks,ϕreds)
                    r1 = r 
                    for i in 1:length(θall)-1
                        for j in i+1:length(ϕall)

                            lines!([r1*sin(θall[i])*cos(ϕall[i]),r1*sin(θall[j])*cos(ϕall[j])] , [r1*sin(θall[i])*sin(ϕall[i]),r1*sin(θall[j])*sin(ϕall[j])] ,[r1*cos(θall[i]),r1*cos(θall[j])], linestyle=:solid, color=:black, linewidth=1.0)

                        end
                    end

                end

            end

        end

    end
    hidedecorations!(ax)
    hidespines!(ax)

    display(fig)
    if path_dir != nothing 
        save(path_dir*filename*".png",fig,px_per_unit=8.0 )
    end
    return fig 
end

# Decomposition of an operator in a sum of constellation 

function Constellation_Decomposition(H::Operator, b::SpinBasis, jj::Float64, L::Int)

    N = floor(Int, 2*jj)

    
    Constellation_basis = [DenseOperator(b, Tlmdk(N,L, 0))]
    Constellation_decomp = [tr(dagger(DenseOperator(b, Tlmdk(N,L, 0)))*H)]

    for M in 1:L 
        
        hlm = tr(dagger(DenseOperator(b, Tlmdk(N,L, M)))*H)

        push!(Constellation_basis, DenseOperator(b, Tlmdk(N,L, M)) + (-1)^M*DenseOperator(b, Tlmdk(N,L, -M)))
        push!(Constellation_decomp, real(hlm))
        push!(Constellation_basis, im*(DenseOperator(b, Tlmdk(N,L, M)) - (-1)^M*DenseOperator(b, Tlmdk(N,L, -M))))
        push!(Constellation_decomp, imag(hlm))
    end

    figs = []
    norms = []
    tot_norm = sum(abs(Constellation_decomp[M]) for M in 1:length(Constellation_decomp))

    for M in 1:length(Constellation_basis)

        if norm(Constellation_decomp[M]) > 10^(-5)

            push!(figs, Majorana_ConstellationH_h_AddToFigure(Constellation_basis[M], b, jj))
            push!(norms, abs(Constellation_decomp[M])/tot_norm)
        end

    end

    return figs, norms
end