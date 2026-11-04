
##########################################
# AC states with inaccessible symmetries #
##########################################


using CairoMakie, GLMakie, Makie, Colors 
using LinearAlgebra
using Combinatorics
using LaTeXStrings
include("/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/Majorana_Constellation_Colin.jl")

# Weyl Character Formula in SU(2)
function WCF_su2(L::Any, θ::Any)

    if abs(sin(θ/2)) < 10^(-5)
        return (2*L+1)*cos((2*L+1)*θ/2)/cos(θ/2)
    else
        return sin((2*L+1)*θ/2)/sin(θ/2)
    end
end
# L = irrep 
# θ = angle of rotation of the SU(2) element 

# Multiplicity of irrep in SU(2)
function multiplicity_irrep_su2(C_G::Vector, χ_G::Vector, L::Any)

    # C_G = vector with the class elements in the same order as χ_G
    # C_G[n] = Vector with all elements of the class written in axis-angle notation

    N = length(C_G) # number of classes 

    Ntot = sum(length(C_G[n]) for n in 1:N) # number of elements 

    mult = sum(length(C_G[n])*WCF_su2(L,C_G[n][1][2])*χ_G[n] for n in 1:N)/Ntot

    return real.(round(mult,digits=5))
end
# C_G = all group elements sorted by conjugacy classes. 
# χ_G = character tables. Column = conjugacy classes (same order as C_G), Rows = irreps (one-dimensional irreps first)
# L = irrep 

# Projector 
function Projector_irrep_su2(C_G::Vector, χ_G::Vector, Π::Function)

    N = length(C_G) # number of classes 
    Ntot = sum(length(C_G[n]) for n in 1:N) # number of elements 

    P = sum(conj(χ_G[n])*sum(Π(C_G[n][m]) for m in 1:length(C_G[n])) for n in 1:N)/Ntot

    return P 
end
# Π = a given representation of SU(2). e.g., Π(U) = e^{-im θ n⋅\vec{J}} for a spin-j representation 

# Find basis of a subspace from the projector 
function get_basis(P::Matrix)

    F = eigen(P)
    println(round.(F.values,digits=4))
    indices = findall(x -> abs(x -1) < 10^(-5), F.values)

    return [F.vectors[:,i] for i in indices]
end
# P = projector 

# Representation SU(2)
function Πj(j, U)

    # U = (axis , angle)
    b = SpinBasis(j/1.0)
    J = [sigmax(b).data , sigmay(b).data , sigmaz(b).data]/2

    Jmat = collect.(J)

    return exp(-im*U[2]*sum(U[1][n]*Jmat[n] for n in 1:3))

end
# 2j+1 dimensional representation of SU(2). 

# Rotation 
function R(n, θ) # Rotation matrices

    R = zeros(3,3)

    R[1,1] = cos(θ) + n[1]^2*(1-cos(θ))
    R[2,2] = cos(θ) + n[2]^2*(1-cos(θ))
    R[3,3] = cos(θ) + n[3]^2*(1-cos(θ))

    R[1,2] = n[1]*n[2]*(1-cos(θ)) - n[3]*sin(θ)
    R[2,1] = n[1]*n[2]*(1-cos(θ)) + n[3]*sin(θ)

    R[1,3] = n[1]*n[3]*(1-cos(θ)) + n[2]*sin(θ)
    R[3,1] = n[1]*n[3]*(1-cos(θ)) - n[2]*sin(θ)

    R[2,3] = n[2]*n[3]*(1-cos(θ)) - n[1]*sin(θ)
    R[3,2] = n[2]*n[3]*(1-cos(θ)) + n[1]*sin(θ)

    return R 
end

# get all group elements 
function get_group(gens::Vector, Id::Any;order=100000)

    G = [Id]

    stop = false 

    while stop == false 

        new_elems = []

        for n in 1:length(G)
            for m in 1:length(gens)

                new = gens[m]*G[n]

                if length(findall(x -> norm(new-x)<10^(-5),vcat(new_elems,G))) == 0 
                    push!(new_elems, new)
                end

            end
        end
        if length(new_elems) == 0 
            stop = true 
        else 
            append!(G,new_elems)
        end

        if length(G) > order 
            return "not the right generators"
        end
    end

    return G 

end
# gens = generators of the group 
# Id = identity element 

# get conjugacy classes 
function get_conjugacyclass(x::Any, Group::Vector)

    class = []

    for n in 1:length(Group)

        y = Group[n] 
        conj = false 
        for m in 1:length(Group)
            gygi = Group[m]*y*inv(Group[m])
            if norm(gygi - x) < 10^(-5)
                conj = true 
            end
        end

        if conj == true 
            push!(class, y)
        end

    end

    return class 
end 
# x = an element of a group G 
# Group = all group elements 
# returns the set of element of G which are conjugate of x in G 

# classify the group in conjugacy classes 
function classify_group(Group::Vector, C_orders, C_size, E::Any;ordered=true)

    G = copy(Group)

    # first : classify 

    class = []
    class_size = []
    class_order = []

    while length(G) > 0 

        elem = G[1]
        order = 1 
        stop = false 
        while stop == false 
            elem_n = elem^order 
            if norm(elem_n - E) < 10^(-5)
                stop = true 
            else 
                order += 1
            end
        end

        conjugates = get_conjugacyclass(elem, Group)

        for n in 1:length(conjugates)
            index = findall(x -> norm(x-conjugates[n])<10^(-5), G)[1]
            deleteat!(G,index)
        end

        push!(class, conjugates)
        append!(class_order, order)
        append!(class_size, length(conjugates))
    end

    println(class_order)
    println(class_size)
    println(length(class))

    # then : get the same order as in the character tables 

    if ordered == true 
        C_elems = []

        for n in 1:length(C_orders) 

            index = findall(x -> class_order[x] - C_orders[n] == 0 && class_size[x] - C_size[n] == 0 , 1:length(class))[1]
            
            push!(C_elems, class[index])
            
            class_size[index] = 0 
            class_order[index] = 0

        end

        return C_elems
    else 
        return class
    end
    
end
# C_orders = order of the different conjuacy classes of the group
# C_size = number of element in each class 
# E = identity element 
# sorts the group by conjugacy class and sort the classes such that the order of each class and their size match the ordering of C_size and C_order 
# /!\ only works for SU(2) where the character of a unitary in the SU(2) irrep only depends on its order !!!!! /!\ 
# /!\ for SU(d), two unitary with the same order may have different character so we must be careful when we sort the table 

# from SO(3) matrix to axis-angle notation 
function SO3_to_AA(class::Vector)

    class_AA = []

    for n in 1:length(class)

        new = []

        for m in 1:length(class[n])

            M = class[n][m] 

            ang = acos(round((tr(M)-1)/2,digits=14)) 

            if abs(ang)<10^(-5) 

                push!(new, ([0.0,0,0],0.0))

            elseif abs(abs(ang)-pi) < 10^(-5)

                D = (M + Diagonal([1 , 1 , 1]))/2 
                eigens = eigen(D)
                index = findall(x -> abs(x-1) < 10^(-6), eigens.values)[1]
                vec = eigens.vectors[:,index] 
                vec = vec/norm(vec)
                push!(new, (vec, pi))

            else 

                ax = [M[3,2] - M[2,3] , M[1,3] - M[3,1] , M[2,1] - M[1,2] ]/(2*sin(ang))
                push!(new, (round.(ax,digits=12), ang))
            end
        end
        push!(class_AA, new)
    end
    return class_AA
end

# from SU(2) matrix to axis-angle notation
function CayleyKlein_to_AA(class::Vector)

    σx = [0 1;1 0]
    σy = [0 -im; im 0]
    σz = [1 0;0 -1]

    class_AA = []
    for n in 1:length(class)
        new = []
        for m in 1:length(class[n])
            M = class[n][m] 


            ang = 2*acos(round(tr(M)/2,digits=14))

            if abs(ang)<10^(-5) 
                push!(new, ([0.0,0,0],0.0))
            elseif abs(abs(ang)-2*pi) < 10^(-5)
                push!(new, ([0.0,0,1],2*pi))
            else 
                ax = [tr(M*σx) , tr(M*σy) , tr(M*σz)]/(-2*im*sin(ang/2))
                push!(new, (ax, ang))
            end
        end
        push!(class_AA, new)
    end
    return class_AA
end

# SO(3) point groups and their character tables 

D2_class = [[0.0 , pi , pi , pi] , [1 1 1 1]]
D2_all = get_group([R([0,0,1],pi),R([1 , 0 , 0],pi)], R([0,0,0],0.0))
D2_classes = classify_group(D2_all, [1 ,2 , 2 , 2], [1, 1 , 1 , 1],R([0,0,0],0.0))
D2_classes = SO3_to_AA(D2_classes)

χD2 = Matrix{Any}(undef, 4,4)
χD2[1,:] = [1 1 1 1]
χD2[2,:] = [1 1 -1 -1]
χD2[3,:] = [1 -1 1 -1]
χD2[4,:] = [1 -1 -1 1]

D3_class = [[0.0 , pi , 2*pi/3] , [1 3 2]]
D3_all = get_group([R([0,0,1],2*pi/3),R([1 , 0 , 0],pi)], R([0,0,0],0.0))
D3_classes = classify_group(D3_all, [1 ,2, 3], [1, 3 ,2],R([0,0,0],0.0))
D3_classes = SO3_to_AA(D3_classes)

χD3 = Matrix{Any}(undef, 3,3)
χD3[1,:] = [1 1 1]
χD3[2,:] = [1 -1 1]
χD3[3,:] = [2 0 -1]

T_class = [[0.0 , pi , 2*pi/3 , 4*pi/3] , [1 3 4 4]]
T_all = get_group([R([0,0,1],2*pi/3),R([-sqrt(1/6),sqrt(1/2),sqrt(1/3)],pi)], R([0,0,0],0.0))
T_classes = classify_group(T_all, [1,2,3,3], [1,3,4,4],R([0,0,0],0.0))
T_classes = SO3_to_AA(T_classes)

χT = Matrix{Any}(undef,4,4)
ω = exp(im*2*pi/3)
χT[1,:] = [1 1 1 1]
χT[2,:] = [1 1 ω ω^2]
χT[3,:] = [1 1 ω^2 ω]
χT[4,:] = [3 -1 0 0]

O_class = [[0.0 , pi , pi , 2*pi/3 , 2*pi/4] , [1 3 6 8 6]]
O_all = get_group([R([0,0,1],pi/2),R([1 0 0],pi/2)], R([0,0,0],0.0))
O_classes = classify_group(O_all, [1 , 2 , 2 , 3 , 4], [1 , 3 ,6 ,8, 6],R([0,0,0],0.0))
O_classes = SO3_to_AA(O_classes)

χO = Matrix{Any}(undef,5,5)
χO[1,:] = [1 1 1 1 1]
χO[2,:] = [1 1 -1 1 -1]
χO[3,:] = [2 2 0 -1 0]
χO[4,:] = [3 -1 -1 0 1]
χO[5,:] = [3 -1 1 0 -1]

I_class = [[0.0 , pi , 2*pi/3 , 2*pi/5 , 4*pi/5] , [1 15 20 12 12]]
Φ = (1+sqrt(5))/2 
I_all = get_group([R([0,0,1],pi/1),R([1-Φ , 0 , Φ]/sqrt(3),2*pi/3)], R([0,0,0],0.0))
I_classes = classify_group(I_all, [1 ,2 ,3 ,5 ,5], [1, 15 ,20 ,12, 12],R([0,0,0],0.0))
I_classes = SO3_to_AA(I_classes)

χI = Matrix{Any}(undef,5,5)
Φp = (1+sqrt(5))/2
Φm = (1-sqrt(5))/2
χI[1,:] = [1 1 1 1 1]
χI[2,:] = [3 -1 0 Φp Φm]
χI[3,:] = [3 -1 0 Φm Φp]
χI[4,:] = [4 0 1 -1 -1]
χI[5,:] = [5 1 -1 0 0]

# SU(2) binary point groups and their character tables 

Sx = [0 1; 1 0]/2
Sy = [0 -im; im 0]/2
Sz = [1 0; 0 -1]/2
Svec = [Sx , Sy,  Sz]
E = exp(-im*Sz*0.0)

bD2_class = [[0.0 , 2*pi , pi , pi , pi] , [1 1 2 2 2]]
bD2_gens = [exp(-im*pi*Sx),exp(-im*pi*Sy)]
bD2_all = get_group(bD2_gens,E)
bD2_classes = classify_group(bD2_all, [1 ,2, 4 ,4 ,4], [1, 1 ,2, 2, 2],E)
bD2_classes = CayleyKlein_to_AA(bD2_classes)

χ2D2 = Matrix{Any}(undef, 5,5)
χ2D2[1,:] = [1 1 1 1 1]
χ2D2[2,:] = [1 1 1 -1 -1]
χ2D2[3,:] = [1 1 -1 1 -1]
χ2D2[4,:] = [1 1 -1 -1 1]
χ2D2[5,:] = [2 -2 0 0 0]


bD3_gens = [exp(-im*2*pi/3*Sz), exp(-im*pi*Sx)]
bD3_all = get_group(bD3_gens,E)
bD3_classes = classify_group(bD3_all, [1,2,3,4,4,6],[1,1,2,3,3,2],E)
bD3_classes = CayleyKlein_to_AA(bD3_classes)

χ2D3 = Matrix{Any}(undef, 6,6)
χ2D3[1,:] = [1 1 1 1 1 1]
χ2D3[2,:] = [1 1 1 -1 -1 1]
χ2D3[3,:] = [1 -1 1 im -im -1]
χ2D3[4,:] = [1 -1 1 -im im -1]
χ2D3[5,:] = [2 2 -1 0 0 -1]
χ2D3[6,:] = [2 -2 -1 0 0 1]


bT_class = [[0.0 , 2*pi , pi , 2*pi/3 , -2*pi/3 , 4*pi/3 , -4*pi/3] , [1 1 6 4 4 4 4]]
n = [sqrt(2)/3,sqrt(2/3),1/3]
bT_gens = [exp(-im*2*pi/3*Sz) , exp(-im*2*pi/3*sum(n[m]*Svec[m] for m in 1:3))]
bT_all = get_group(bT_gens, E)
bT_classes = classify_group(bT_all,[1, 2, 4, 6, 6, 3, 3],[1, 1, 6, 4, 4, 4, 4],E)
bT_classes = CayleyKlein_to_AA(bT_classes)

χ2T = Matrix{Any}(undef,7 , 7)
ω = exp(im*2*pi/3)
χ2T[1,:] = [1 1 1 1 1 1 1]
χ2T[2,:] = [1 1 1 ω ω^2 ω^2 ω]
χ2T[3,:] = [1 1 1 ω^2 ω ω ω^2]
χ2T[4,:] = [2 -2 0 1 1 -1 -1]
χ2T[5,:] = [2 -2 0 ω^2 ω -ω^2 -ω]
χ2T[6,:] = [2 -2 0 ω^2 ω -ω -ω^2]
χ2T[7,:] = [3 3 -1 0 0 0 0]

bO_class = [[0.0 , 2*pi , pi , pi , 2*pi/3 , 4*pi/3 , pi/2 , -pi/2] , [1 1 12 6 8 8 6 6]]
bO_gens = [exp(-im*2*pi/4*Sz) , exp(-im*2*pi/4*Sx)]
bO_all = get_group(bO_gens, E)
bO_classes = classify_group(bO_all,[1,2,4,4,6,3,8,8],[1,1,12,6,8,8,6,6],E)
bO_classes = CayleyKlein_to_AA(bO_classes)

χ2O = Matrix{Any}(undef,8,8)
χ2O[1,:] = [1 1 1 1 1 1 1 1]
χ2O[2,:] = [1 1 -1 1 1 1 -1 -1]
χ2O[3,:] = [2 2 0 2 -1 -1 0 0]
χ2O[4,:] = [2 -2 0 0 1 -1 sqrt(2) -sqrt(2)]
χ2O[5,:] = [2 -2 0 0 1 -1 -sqrt(2) sqrt(2)]
χ2O[6,:] = [3 3 -1 -1 0 0 1 1]
χ2O[7,:] = [3 3 1 -1 0 0 -1 -1]
χ2O[8,:] = [4 -4 0 0 -1 1 0 0]


n=[1-Φ , 0 , Φ]/sqrt(3)
bI_gens = [exp(-im*pi*Sz) , exp(-im*2*pi/3*sum(n[i]*Svec[i] for i in 1:3))]
bI_all = get_group(bI_gens, E)
bI_classes = classify_group(bI_all,[1 , 2 , 3 , 4 , 5 , 5 , 6 , 10 , 10],[1 , 1 , 20 , 30 , 12 , 12 , 20 , 12 , 12],E)
bI_classes = CayleyKlein_to_AA(bI_classes)
χ2I = Matrix{Any}(undef,9,9)
Φp = (1+sqrt(5))/2
Φm = (1-sqrt(5))/2
χ2I[1,:] = [1 1 1 1 1 1 1 1 1]
χ2I[2,:] = [2 -2 -1 0 -Φm -Φp 1 Φp Φm]
χ2I[3,:] = [2 -2 -1 0 -Φp -Φm 1 Φm Φp]
χ2I[4,:] = [3 3 0 -1 Φm Φp 0 Φp Φm]
χ2I[5,:] = [3 3 0 -1 Φp Φm 0 Φm Φp]
χ2I[6,:] = [4 4 1 0 -1 -1 1 -1 -1]
χ2I[7,:] = [4 -4 1 0 -1 -1 -1 1 1]
χ2I[8,:] = [5 5 -1 1 0 0 -1 0 0]
χ2I[9,:] = [6 -6 0 0 1 1 0 -1 -1]


#################################################################
# Multiplicity of one-dimension irreps in H^(j) for different j #
#################################################################

data_dir = "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/G-invariant states SU2/"

######
# D2 #
######

# integer spins 

js = 1:1:10

D2_ρ1 = [multiplicity_irrep_su2(D2_classes, χD2[1,:],j) for j in js]
D2_ρ2 = [multiplicity_irrep_su2(D2_classes, χD2[2,:],j) for j in js]
D2_ρ3 = [multiplicity_irrep_su2(D2_classes, χD2[3,:],j) for j in js]
D2_ρ4 = [multiplicity_irrep_su2(D2_classes, χD2[4,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0)
lines!(js, D2_ρ1, label=L"\rho_1")
lines!(js, D2_ρ2, label=L"\rho_2")
lines!(js, D2_ρ3, label=L"\rho_3",linestyle=:dash)
lines!(js, D2_ρ4, label=L"\rho_4",linestyle=:dot)
axislegend(ax,position=:lt,labelsize=20.0f0)

#save(data_dir*"D2 results low irrep.pdf", fig)
display(fig)

# half-integer spins 

js = 0.5:0.5:20

D2_ρ1 = [multiplicity_irrep_su2(bD2_classes, χ2D2[1,:],j) for j in js]
D2_ρ2 = [multiplicity_irrep_su2(bD2_classes, χ2D2[2,:],j) for j in js]
D2_ρ3 = [multiplicity_irrep_su2(bD2_classes, χ2D2[3,:],j) for j in js]
D2_ρ4 = [multiplicity_irrep_su2(bD2_classes, χ2D2[4,:],j) for j in js]

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"L", ylabel=L"a_i")
lines!(js, D2_ρ1, label=L"\rho_1")
lines!(js, D2_ρ2, label=L"\rho_2")
lines!(js, D2_ρ3, label=L"\rho_3")
lines!(js, D2_ρ4, label=L"\rho_4")
display(fig)



#####
# D3 #
##### 

# integer spins 

js = 1:1:10

D3_ρ1 = [multiplicity_irrep_su2(D3_classes, χD3[1,:],j) for j in js]
D3_ρ2 = [multiplicity_irrep_su2(D3_classes, χD3[2,:],j) for j in js]

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"L", ylabel=L"a_i")
lines!(js, D3_ρ1, label=L"\rho_1")
lines!(js, D3_ρ2, label=L"\rho_2")
axislegend(ax,position=:lt)
display(fig)

#####
# T #
##### 

# integer spins 

js = 1:1:10

T_ρ1 = [multiplicity_irrep_su2(T_classes, χT[1,:],j) for j in js]
T_ρ2 = [multiplicity_irrep_su2(T_classes, χT[2,:],j) for j in js]
T_ρ3 = [multiplicity_irrep_su2(T_classes, χT[3,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0,yticks=([0,1,2]))
lines!(js, T_ρ1, label=L"\rho_1",linewidth=2.0)
lines!(js, T_ρ2, label=L"\rho_2",linewidth=2.0)
lines!(js, T_ρ3, label=L"\rho_3",linestyle=:dash)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)

#save(data_dir*"Tetrahedron one-dim.pdf", fig)

# integer spins 

js = 1:1:16

T_ρ1 = [multiplicity_irrep_su2(T_classes, χT[4,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0) #,yticks=([0,1,2]))
lines!(js, T_ρ1, label=L"\rho_1",linewidth=2.0)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)


# half-integer spins 

js = 0.5:1:9.5

T_ρ1 = [multiplicity_irrep_su2(bT_classes, χ2T[1,:],j) for j in js]
T_ρ2 = [multiplicity_irrep_su2(bT_classes, χ2T[2,:],j) for j in js]
T_ρ3 = [multiplicity_irrep_su2(bT_classes, χ2T[3,:],j) for j in js]

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"L", ylabel=L"a_i")
lines!(js, T_ρ1, label=L"\rho_1")
lines!(js, T_ρ2, label=L"\rho_2")
lines!(js, T_ρ3, label=L"\rho_3")
axislegend(ax,position=:lt)
display(fig)


#####
# O #
##### 

# integer spins 

js = 1:1:13

T_ρ1 = [multiplicity_irrep_su2(O_classes, χO[1,:],j) for j in js]
T_ρ2 = [multiplicity_irrep_su2(O_classes, χO[2,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0,xticks=([2,4,6,8,10,12,14]),yticks=([0,1,2]))
lines!(js, T_ρ1, label=L"\rho_1",linewidth=2.0)
lines!(js, T_ρ2, label=L"\rho_2",linestyle=:dash)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)

#save(data_dir*"Octahedron one-dim.pdf", fig)

# integer spins 2D irrep 

js = 1:1:16

T_ρ1 = [multiplicity_irrep_su2(O_classes, χO[3,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0,xticks=([2,4,6,8,10,12,14,16]))#,yticks=([0,1,2]))
lines!(js, T_ρ1, label=L"\rho_3",linewidth=2.0)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)

# integer spins 3D irrep 

js = 1:1:16

T_ρ1 = [multiplicity_irrep_su2(O_classes, χO[4,:],j) for j in js]
T_ρ2 = [multiplicity_irrep_su2(O_classes, χO[5,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0,xticks=([2,4,6,8,10,12,14,16]))#,yticks=([0,1,2]))
lines!(js, T_ρ1, label=L"\rho_4",linewidth=2.0)
lines!(js, T_ρ2, label=L"\rho_5",linewidth=2.0)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)


#####
# I #
##### 

# integer spins 

js = 1:1:13

I_ρ1 = [multiplicity_irrep_su2(I_classes, χI[1,:],j) for j in js]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0)
lines!(js, I_ρ1, label=L"\rho_1",linewidth=2.0)
axislegend(ax,position=:lt, labelsize=20.0f0)
display(fig)

#save(data_dir*"Icosahedron one-dim.pdf", fig)

# integer spins all irreps 

js = 1:1:13

I_s = [[multiplicity_irrep_su2(I_classes, χI[s,:],j) for j in js] for s in 1:5]

fig = Figure(;size=(400,300))
ax = Axis(fig[1,1], xlabel=L"j", ylabel=L"a_i",
xlabelsize = 25.0f0, ylabelsize=25.0f0,ylabelrotation=0)

for s in [1,3]
    lines!(js, I_s[s], label="$s",linewidth=2.0)
end

axislegend(ax,position=:lt, labelsize=12.0f0)
display(fig)

#save(data_dir*"Icosahedron one-dim.pdf", fig)

############################
# Some invariant operators #
############################

######
# D2 #
######

# G-invariant subspace (k=2,t=1,j=2) : j = 2 & ρ1 irrep 

D2_all = get_group([R([sqrt(2),0,1]/sqrt(3),pi),R([-1/sqrt(6) , 1/sqrt(2) , 1/sqrt(3)],pi)], R([0,0,0],0.0))
D2_classes = classify_group(D2_all, [1 ,2 , 2 , 2], [1, 1 , 1 , 1],R([0,0,0],0.0))
D2_classes = SO3_to_AA(D2_classes)

P = Projector_irrep_su2(D2_classes, χD2[1,:],x -> Πj(2,x)) # projector 
basis = get_basis(P) # basis of the subspace 
j = 2.0
s = SpinBasis(j)
P_op = DenseOperator(s,P) # projector defined as an operator in QuantumOptics.jl 
Majorana_ConstellationH_h(P_op, s ,j) # Majorana constellation of the projector 
GLMakie.activate!()
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])

fig1 = Majorana_Constellationψ(ψ1,s,j, ϕmin=pi/2,ϕaz=-pi/160) 

V1 = 1.3*[sqrt(2),0,1]/sqrt(3)
V2 = 1.3*[-sqrt(1/6),sqrt(1/2),sqrt(1/3)]
V3 = 1.3*[-sqrt(1/6),-sqrt(1/2),sqrt(1/3)]
lines!(Point3.([-V1,V1]),  color=RGB(204/255,0/255,0/255),linewidth=2.5)
lines!(Point3.([-V2,V2]),  color=RGB(204/255,0/255,0/255),linewidth=2.5)
lines!(Point3.([-V3,V3]),  color=RGB(204/255,0/255,0/255),linewidth=2.5)

Majorana_Constellationψ(ψ2,s,j, path_dir=data_dir, filename="D2 spin-2 2")

ψ_s(α,β) = (α*ψ1 + exp(-im*β)*(1-α)*ψ2)/sqrt(2*α^2 - 2*α + 1) # all states of (2,1,2)

# changing α to get from ψ1 to ψ2
αs = LinRange(0,1,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(α,0.0), s, j;path_dir = fig_dir, filename = "D2 - spin $j - ρ1 - β = 0 - α=$α" )
end
# for α=0.5, changing β to see what a phase does 
αs = LinRange(0,pi/2,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(0.5,α), s, j;path_dir = fig_dir, filename = "D2 - spin $j - ρ1 - β = $α - α=0.5" )
end

#####
# D3 #
##### 

# D3-invariant subspace (k=2,t=1,j=4) : j = 4 & ρ1 irrep 

P = Projector_irrep_su2(D3_classes, χD3[1,:],x -> Πj(4,x))
basis = get_basis(P)

j = 4.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])
Majorana_Constellationψ(ψ1,s,j)
Majorana_Constellationψ(ψ2,s,j)

ψ_s(α,β) = exp(-im*pi/3*sigmaz(s)/2)*(α*ψ1 + exp(-im*β)*(1-α)*ψ2)/sqrt(2*α^2 - 2*α + 1)

αs = LinRange(0,1,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(α,0.0), s, j;path_dir = fig_dir, filename = "D3 - spin $j - ρ1 - α=$α - β = 0.0" )
end

αs = LinRange(0,pi/2,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(0.5,α), s, j;path_dir = fig_dir, filename = "D3 - spin $j - ρ1 - α=0.5 - β = $α" )
end

# D3-invariant subspace (k=2,t=1,j=3) : j=3 & ρ2 irrep 

P = Projector_irrep_su2(D3_classes, χD3[2,:],x -> Πj(3,x))
basis = get_basis(P)

j = 3.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])
Majorana_Constellationψ(ψ1,s,j)
Majorana_Constellationψ(ψ2,s,j)

ψ_s(α,β) = (α*ψ1 + exp(-im*β)*(1-α)*ψ2)/sqrt(2*α^2 - 2*α + 1)

αs = LinRange(0,1,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(α,0.0), s, j;path_dir = fig_dir, filename = "D3 - spin $j - ρ2 - α=$α - β=0.0" )
end


αs = LinRange(0,pi/2,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(0.5,α), s, j;path_dir = fig_dir, filename = "D3 - spin $j - ρ2 - α=0.5 - β=$α" )
end

# half-integer spins 

js = 0.5:1:9.5

D3_ρ1 = [multiplicity_irrep_su2(bD3_classes, χ2D3[1,:],j) for j in js]
D3_ρ2 = [multiplicity_irrep_su2(bD3_classes, χ2D3[2,:],j) for j in js]
D3_ρ3 = [multiplicity_irrep_su2(bD3_classes, χ2D3[3,:],j) for j in js]
D3_ρ4 = [multiplicity_irrep_su2(bD3_classes, χ2D3[4,:],j) for j in js]

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"L", ylabel=L"a_i")
lines!(js, D3_ρ1, label=L"\rho_1")
lines!(js, D3_ρ2, label=L"\rho_2")
lines!(js, D3_ρ3, label=L"\rho_3")
lines!(js, D3_ρ4, label=L"\rho_4")
axislegend(ax,position=:lt)
display(fig)

# 2D3-invariant subspace (k=2, t=1, j=4.5) : spin-j=5.5 & ρ3

P = Projector_irrep_su2(bD3_classes, χ2D3[3,:],x -> Πj(4.5,x))
basis = get_basis(P)

j = 4.5
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])
Majorana_Constellationψ(ψ1,s,j)
Majorana_Constellationψ(ψ2,s,j)

#####
# T #
##### 

# T-invariant subspace (k=2,t=2,j=6), j = 6 & ρ1 irrep  

P = Projector_irrep_su2(T_classes, χT[1,:],x -> Πj(6,x))
basis = get_basis(P)

j = 6.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ1,s,j,ϕmin = 0.0,ϕmax = 2*pi)

V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3
V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end
display(fig1)
save(data_dir*"spin-6 Tetrahedron state 1.png",fig1;px_per_unit=15.0)

fig1 = Majorana_Constellationψ(ψ2,s,j,ϕmin = 0.0,ϕmax = 2*pi)

V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3
V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end

#save(data_dir*"spin-6 Tetrahedron state 2.png",fig1;px_per_unit=15.0)

ψ_s(α,β) = (α*ψ1 + exp(-im*β)*(1-α)*ψ2)/sqrt(2*α^2 - 2*α + 1)

αs = LinRange(0,1,10)
fig_dir = data_dir
for α in αs 
    Majorana_Constellationψ(ψ_s(α,0.0), s, j;path_dir = fig_dir, filename = "T - spin $j - ρ1 - α=$α - β=0.0" )
end

αs = LinRange(0,pi/2,10)
fig_dir =  data_dir
for α in αs 
    Majorana_Constellationψ(ψ_s(0.5,α), s, j;path_dir = fig_dir, filename = "T - spin $j - ρ1 - α=0.5 - β=$α" )
end

# projector, j = 8 & ρ2

P = Projector_irrep_su2(T_classes, χT[2,:],x -> Πj(8,x))
basis = get_basis(P)

j = 8.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])

ψ_s(α,β) = (α*ψ1 + exp(-im*β)*(1-α)*ψ2)/sqrt(2*α^2 - 2*α + 1)

αs = LinRange(0,1,10)
fig_dir = ""
for α in αs 
    Majorana_Constellationψ(ψ_s(α,00), s, j;path_dir = fig_dir, filename = "T ρ2 - spin $j - ρ2 - α=$α - β=0.0" )
end

# ρ2 , j=2 

P = Projector_irrep_su2(T_classes, χT[2,:],x -> Πj(2,x))
basis = get_basis(P)

j = 2.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ1,s,j,ϕmin =0.0,ϕmax = 2*pi)
arrows3d!(Point3.([[0,0,0]]), Point3.([[sqrt(2),0,1]/sqrt(3)]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)
arrows3d!(Point3.([[0,0,0]]), Point3.([[-sqrt(1/6),sqrt(1/2),1/sqrt(3)]]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)
arrows3d!(Point3.([[0,0,0]]), Point3.([[-sqrt(1/6),-sqrt(1/2),1/sqrt(3)]]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)

#save(data_dir*"D2 sub logical 0.png", fig1;px_per_unit=15)

# ρ3 , j=2 

P = Projector_irrep_su2(T_classes, χT[3,:],x -> Πj(2,x))
basis = get_basis(P)

j = 2.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ1,s,j,ϕmin = 0.0,ϕmax = 2*pi)
arrows3d!(Point3.([[0,0,0]]), Point3.([[sqrt(2),0,1]/sqrt(3)]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)
arrows3d!(Point3.([[0,0,0]]), Point3.([[-sqrt(1/6),sqrt(1/2),1/sqrt(3)]]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)
arrows3d!(Point3.([[0,0,0]]), Point3.([[-sqrt(1/6),-sqrt(1/2),1/sqrt(3)]]),color=:red, 
shaftradius=0.02,
tipradius=0.07, tiplength=0.2)

save(data_dir*"D2 sub logical 1.png", fig1;px_per_unit=15)


#####
# O #
##### 

# O-invariant state, j = 6 & ρ1 irrep

P = Projector_irrep_su2(O_classes, χO[1,:],x -> Πj(6,x))
norm(P)
basis = get_basis(P)

j = 6.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)

Jvec = [sigmax(s),sigmay(s),sigmaz(s)]/2
nnn = [1,1,0]/sqrt(2)
R_OtoT = exp(-im*pi/4*Jvec[3])*exp(-im*acos(1/sqrt(3))*sum(nnn[m]*Jvec[m] for m in 1:3))
ψ1 = R_OtoT*Ket(s, basis[1])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ1,s,j,ϕmin = 0.0,ϕmax = 2*pi)

V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3

V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end
display(fig1)

save(data_dir*"Tetrasubspace logical 0.png", fig1;px_per_unit=15)


# O-invariant state, j = 6 & ρ2 irrep

P = Projector_irrep_su2(O_classes, χO[2,:],x -> Πj(6,x))
norm(P)
basis = get_basis(P)

j = 6.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
ψ2 = R_OtoT*Ket(s, basis[1])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ2,s,j,ϕmin = 0.0,ϕmax = 2*pi)


V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3

V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end
save(data_dir*"Tetrasubspace logical 1.png", fig1;px_per_unit=15)


# O-invariant subspace, j = 8 & ρ2 irrep

P = Projector_irrep_su2(O_classes, χO[1,:],x -> Πj(12,x))
norm(P)
basis = get_basis(P)

j = 12.0
s = SpinBasis(j)
P_op = DenseOperator(s,P)
Majorana_ConstellationH_h(P_op, s ,j)
ψ1 = Ket(s, basis[1])
ψ2 = Ket(s, basis[2])
GLMakie.activate!()

fig1 = Majorana_Constellationψ(ψ1,s,j,ϕmin = 0.0,ϕmax = 2*pi)

V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3
V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end
display(fig1)
save(data_dir*"spin-6 Tetrahedron state 1.png",fig1;px_per_unit=15.0)

fig1 = Majorana_Constellationψ(ψ2,s,j,ϕmin = 0.0,ϕmax = 2*pi)

V1 = [0,0,1]
V2 = RotationArbitrary([sqrt(2) , 0 , 1]/sqrt(3) , pi)*V1
V3 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V2
V4 = RotationArbitrary([0 , 0 , 1] , 2*pi/3)*V3
V = [V1,V2,V3,V4]
edges = [1,2],[1,3],[1,4],[2,3],[3,4],[4,2]
faces = [1 2 3;
1 3 4;
1 4 2;
4 2 3]
Vtx = [V1[1] V1[2] V1[3];
V2[1] V2[2] V2[3];
V3[1] V3[2] V3[3];
V4[1] V4[2] V4[3]]
mesh!(Vtx, faces, color=RGB(255/255,153/255,153/255), alpha=0.5, transparency=true)
for n in 1:length(edges)
    lines!(Point3.([V[edges[n][1]],V[edges[n][2]]]), color=RGB(204/255,0/255,0/255))
end
#save(data_dir*"spin-6 Tetrahedron state 2.png",fig1;px_per_unit=15.0)



        ##########
        # figure # 
        ##########

CairoMakie.activate!()
fig_dir = "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/G-invariant states SU2/"

Gs = [bD2_classes , bD3_classes , bT_classes , bO_classes , bI_classes]

name_g = [L"\text{D}_2", L"\text{D}_3" , L"\text{T}" , L"\text{O}",L"\text{I}"]

Ls = 0.0:0.5:7.5

fig = Figure()

ax = Axis(fig[1,1])

lines!(Point2.([a*[0,1]+[0.5,0] for a in 0:15]), color=:black)
arrows2d!(Point2.([15*[0,1]+[0.5,0]]), Point2.([[0,0.7]]),shaftwidth=1.5, tipwidth=10)

for l in 0:7 
    lines!(Point2.([[0.5,2*l],[0.5-0.05,2*l]]), color=:black)
    text!(Point2([0.5-0.1,2*l]), text="$(l)", align=(:right,:bottom))
end
text!(Point2([0.5-0.4,8]), text=L"L",fontsize=20.0f0)

for p in 1:length(Gs)

    lines!(Point2.([a*[0,1] + [p,0] for a in 0:15]), color=:black)
    scatter!(Point2.([a*[0,1] + [p,0]  for a in 0:15]),color=:white, strokewidth=1.0, strokecolor=:black, markersize=20.0)

    text!(Point2([p,16]), text=name_g[p], align=(:center,:center), fontsize=20.0f0)

    mults = Int.([multiplicity_irrep_su2(Gs[p],[1.0 for n in 1:length(Gs[p])],L) for L in Ls])

    for m in 1:length(mults)
       
        if mults[m] > 0
            scatter!(Point2([p,m-1]), color=(:red,0.5),markersize=20.0)
            text!(Point2([p,m-1]), text="$(mults[m])", align=(:center,:center))
        else 
            scatter!(Point2([p,m-1]), color=(:green,0.5),markersize=20.0)
        end

    end
    
end

hidedecorations!(ax)
hidespines!(ax)
display(fig)

#save(fig_dir*"results low irreps.pdf", fig)

        ################
        # Beyond SU(2) #
        ################

# get the roots 

function get_roots_Cartanmatrix(N::Int)

    # get all simple roots 

    Δ = []
    for n in 1:N-1
        root = [0.0 for n in 1:N]
        root[n] = 1.0 
        root[n+1] = -1.0
        push!(Δ, root)
    end

    # Cartan matrix 

    Cartan_matrix = Matrix{Any}(undef, N-1,N-1)

    for i in 1:N-1 
        for j in 1:N-1 
            Cartan_matrix[i,j] = 2*dot(Δ[j], Δ[i])/dot(Δ[i],Δ[i])
        end
    end

    # construct all positive roots 

    stop = false 

    # all root string of height = 1
    Δ_string = []
    for n in 1:N-1 
        new_string = zeros(N-1)
        new_string[n] = 1.0 
        push!(Δ_string, new_string)
    end

    rpos_string = copy(Δ_string)

    while stop == false 

        new_strings = []

        for λ in 1:length(rpos_string) # for each positive roots already found 

            for j in 1:N-1 # for all positive simple roots 

                r = -sum(rpos_string[λ][i]*Cartan_matrix[i,j] for i in 1:N-1)

                # if r<0, we do nothing 

                # if r>0, we can construct new root, in the direction Δ[j]

                if r > 0 

                    for p in 1:r 

                        new_string = rpos_string[λ] + p*Δ_string[j]

                        if length(findall(x -> norm(x-new_string)<10^(-5), new_strings)) == 0 && length(findall(x -> norm(x-new_string)<10^(-5), rpos_string)) == 0
                            push!(new_strings,new_string)
                        end
                    end

                end
            end

        end

        if length(new_strings) == 0 
            stop = true 
        else 
            append!(rpos_string, new_strings)
        end

    end

    # now compute those roots 
    Rp = []

    for n in 1:length(rpos_string)
        push!(Rp, sum(rpos_string[n][i]*Δ[i] for i in 1:N-1))
    end

    # put negative roots 

    Rn = -Rp 

    return Δ, Rp, Rn

end 

#=
output : Tuple(3)
    -Δ the base 
    -R+ the set of positive roots 
    -R- the set of negative roots 
=#

# get the Weyl group 

function get_Weyl(N::Int, Rpos, Rneg, Δ) # first method : permutation and sign of permutation 

    Rs = copy(Rpos)
    append!(Rs,Rneg)

    println(length(Rs))

    # action of Weyl 

    function refl(a,b)
        return b - 2*dot(b,a)*a/dot(a,a)
    end

    # generators of Weyl in permutation representation

    W = [[n for n in 1:length(Rs)]]

    for a in Δ 

        w = []

        for b in Rs 

            c = refl(a,b)
            index = findall(x -> x == c , Rs)
            println(index)

            push!(w, index[1])
        end

        push!(W, w)

    end

    # product of permutation 

    function product_perm(a,b)
        c = [a[b[n]] for n in 1:length(b)]
        return c
    end

    # generate ! 

    Sd = factorial(N)

    while length(W) < Sd 

        new_perm = []

        for i in W
            for j in W
                new = product_perm(i,j)
                if length(findall(x -> x == new, new_perm)) == 0
                    push!(new_perm, new)
                end
            end
        end

        for new in new_perm 
            if length(findall(x -> x == new, W)) == 0
                push!(W,new)
            end
        end

    end

    return W, levicivita.(W)
end

#= 
output : tuple(2)
    - the group of permutation on [Rpos , Rneg]
    - the determinant of each permutation
=#

function get_Weyl_test(N::Int, Rpos, Rneg, Δ) # second method 

    Rs = copy(Rpos)
    append!(Rs,Rneg)

    println(length(Rs))

    # action of Weyl 

    function refl(a,b)
        return b - 2*dot(b,a)*a/dot(a,a)
    end

    # generators of Weyl in permutation representation

    W = [[n for n in 1:length(Rs)]]
    sign_W = [1]
    W_base = []

    for a in Δ 

        w = []

        for b in Rs 

            c = refl(a,b)
            index = findall(x -> x == c , Rs)
            println(index)

            push!(w, index[1])
        end

        push!(W, w)
        push!(W_base, w)
        push!(sign_W, -1)

    end

    # product of permutation 

    function product_perm(a,b)
        c = [a[b[n]] for n in 1:length(b)]
        return c
    end

    # generate ! 

    Sd = factorial(N)

    parity = 2

    New_layer = []

    for i in W_base
        for j in W_base
            new = product_perm(i,j)
            if length(findall(x -> x == new, New_layer)) == 0 && length(findall(x -> x == new, W)) == 0
                push!(New_layer, new)
            end
        end
    end

    for new in New_layer 
        push!(W,new)
        push!(sign_W, (-1)^parity)
    end

    Old_Layer = copy(New_layer)

    while length(W) < Sd 

        parity += 1
        
        New_layer = []

        for i in Old_Layer
            for j in W_base
                new = product_perm(i,j)
                if length(findall(x -> x == new, New_layer)) == 0 && length(findall(x -> x == new, W)) == 0
                    push!(New_layer, new)
                end
            end
        end

        for new in New_layer 
            push!(W,new)
            push!(sign_W, (-1)^parity)
        end

        Old_Layer = copy(New_layer)

    end

    return W, sign_W
end
#= 
output : tuple(2)
    - the group of permutation on [Rpos , Rneg]
    - the determinant of each permutation
=#

function Weyl_action(w::Vector, x::Vector)

    # x is an element written as l.c. of all roots 
    # w permutes the basis 

    if length(w) != length(x)
        println("big problem")
    end

    return [x[w[n]] for n in 1:length(x)]
end
#= 
output :
    - root basis decomposition of a highest weight after application of the Weyl element 
=#

# get the fundamental weights 

function get_fundamentalweights(N::Int, Δ::Vector) # get the fundamental weights 


    A = Matrix{Float64}(undef,N-1,N-1)

    for k in 1:N-1 
        for i in 1:N-1
            A[k,i] = dot(Δ[i],Δ[k])/dot(Δ[k],Δ[k])
        end
    end

    Ainv = inv(A) 

    weight_f = Vector{Vector}(undef, N-1)
    weight_f_rootbasis = Vector{Vector}(undef, N-1)

    for j in 1:N-1 

        ej = [0.0 for n in 1:N-1]
        ej[j] = 1.0
        
        decomp = 0.5*Ainv*ej

        weight_f[j] = sum(decomp[n]*Δ[n] for n in 1:N-1)
        weight_f_rootbasis[j] = decomp
    end

    return weight_f,weight_f_rootbasis
end
#= 
output : Tuple(2)
    - fundamental weight as a vector (-> diagonal matrix)
    - decomposition of the fundamental weight as a l.c. of the base 
=#

# from Dynkin label to l.c. of fundamental weights 

function Dynkin_to_highestweight(D::Vector, fw::Vector, fw_rb::Vector)

    μ = sum(D[n]*fw[n] for n in 1:length(D))

    μ_rb = sum(D[n]*fw_rb[n] for n in 1:length(D))

    return μ,μ_rb
end

#= 
output : Tuple(2)
    - highest weight as a vector (-> diagonal matrix)
    - decomposition of the highest weight as l.c. of the base 
=#

# extend Δ basis to R basis 

function extend_rootbasis(μ_rb::Vector, Δ::Vector, Rp::Vector, Rm::Vector)

    μ_ext = copy(μ_rb)
    extension = length(Rp) + length(Rm) - length(Δ)
    append!(μ_ext , [0.0 for e in 1:extension])

    return μ_ext
end

#= 
output : 
    - highest weight in complete root basis decomposition
=#

# from root basis to Vector in the Cartan 

function Base_to_h(μ::Vector, Δ::Vector, Rpos::Vector, Rneg::Vector)
    
    Rs = copy(Rpos)
    append!(Rs, Rneg)

    if length(μ) == length(Δ) # e.g. weights 

        return sum(μ[n]*Δ[n] for n in 1:length(Δ))

    elseif length(μ) == length(Rpos) # e.g. δ

        return sum(μ[n]*Rpos[n] for n in 1:length(Rpos))

    elseif length(μ) == length(Rs)

        return sum(μ[n]*Rs[n] for n in 1:length(Rs))

    end
end

#= 
output : 
    - weight as a vector (-> diagonal matrix)
=#

# Weyl dimension formula 

function Weyl_dim(D::Vector,Rpos::Vector, fw::Vector, fw_rb::Vector)

    δ = sum(r for r in Rpos)/2
    μ = Dynkin_to_highestweight(D, fw,fw_rb)[1]

    return prod(dot(a, μ+δ)/dot(a,δ) for a in Rpos)
end

#= 
output : 
    - dimension of the irrep
=#

# diagonalizes U and gets the element of the Cartan 

function U_to_H(U::Any)

    eigU = Complex.(eigvals(U))
    H = angle.(eigU) 

    # H traceless ? ok
    # if not, remove the identity otherwise its not in the Cartan 

    #H = H - sum(H[n] for n in 1:length(H))*[1.0 for n in 1:length(H)]/length(H)

    return im*H 
end

# Weyl character formula 

function Weyl_char(D::Vector, W::Vector, detW::Vector, Δ::Vector, Rpos::Vector,Rneg::Vector, fw::Vector, fw_rb::Vector,U::Any)

    h = U_to_H(U) 
    println("start")

    function inner_product(w,x,H)

        new_x = Weyl_action(w,x)

        left_side = Base_to_h(new_x,Δ,Rpos,Rneg)

        return dot(left_side,H)

    end
    
    if norm(h) < 10^(-5)

        println("identity")

        return Weyl_dim(D,Rpos,fw,fw_rb)

    elseif abs(abs(tr(U)) - length(h)) < 10^(-5)

        d = length(h)
        k = abs(imag(h[end])/(2*pi/d))

        if imag(h[end]) > 0 
            h[end] -= im*2*pi*k
        else 
            h[end] += im*2*pi*k
        end
        
        δ_rb = [0.5 for n in 1:length(Rpos)]
        append!(δ_rb, [0.0 for n in 1:length(Rneg)])

        μ = Dynkin_to_highestweight(D,fw,fw_rb)[2]

        μ = extend_rootbasis(μ,Δ,Rpos,Rneg)

        println("other element of the center")

        return Weyl_dim(D,Rpos,fw,fw_rb)*exp(inner_product(W[1],μ,h))#*U[1,1]^(sum(D[n] for n in 1:length(D)))

    else 

        if abs(sum(imag(h[n]) for n in 1:length(h)) -2*pi) < 10^(-5) 
            h[end] -= 2*pi*im 
        elseif abs(sum(imag(h[n]) for n in 1:length(h)) + 2*pi) < 10^(-5)
            h[end] += 2*pi*im
        elseif abs(sum(imag(h[n]) for n in 1:length(h)) - 4*pi) < 10^(-5)
            h[end] -= 4*pi*im
        elseif abs(sum(imag(h[n]) for n in 1:length(h)) + 4*pi) < 10^(-5)
            h[end] += 4*pi*im
        end
        
        δ_rb = [0.5 for n in 1:length(Rpos)]
        append!(δ_rb, [0.0 for n in 1:length(Rneg)])

        μ = Dynkin_to_highestweight(D,fw,fw_rb)[2]

        μ = extend_rootbasis(μ,Δ,Rpos,Rneg)

        num = sum(detW[w]*exp(inner_product(W[w],μ+δ_rb,h)) for w in 1:length(W))
        denom  = sum(detW[w]*exp(inner_product(W[w],δ_rb,h)) for w in 1:length(W))

        if abs(denom) > 0.000001

            return num/denom 

        else 

            println("L'Hospital !")

            function LHospital(o::Int)

                # all possible combinations of "o" positive roots 

                all_possibilities = Any 
                if o == 1
                    all_possibilities = [[n] for n in 1:length(Rpos)]
                elseif o == 2 
                    all_possibilities = [[n,m] for n in 1:length(Rpos), m in 1:length(Rpos)]
                elseif o == 3
                    all_possibilities = [[n,m,p] for n in 1:length(Rpos), m in 1:length(Rpos), p in 1:length(Rpos)]
                elseif o == 4 
                    all_possibilities = [[n,m,p,t] for n in 1:length(Rpos), m in 1:length(Rpos), p in 1:length(Rpos), t in 1:length(Rpos)]
                end

                # get through the possibilities 

                function denomi(vec) 
                    return sum(detW[w]*prod(inner_product(W[w],δ_rb,Rpos[vec[v]]) for v in 1:length(vec))*exp(inner_product(W[w],δ_rb,h)) for w in 1:length(W))
                end
                
                function numera(vec) 
                    return sum(detW[w]*prod(inner_product(W[w],μ+δ_rb,Rpos[vec[v]]) for v in 1:length(vec))*exp(inner_product(W[w],μ+δ_rb,h)) for w in 1:length(W))
                end

                denom = denomi.(all_possibilities)


                index = findall(x -> abs(x) > 10^(-5), denom)


                if length(index) == 0 

                    return nothing 

                else 

                    vec = all_possibilities[index[1]]

                    return numera(vec)/denomi(vec)

                end

            end

            o = 1 

            stop = false 

            while stop == false 

                Hosp = LHospital(o)

                if Hosp == nothing && o >= 4

                    println("big trouble...")

                    return "damn it"

                elseif Hosp == nothing

                    o+=1

                else 

                    println("order $o OK")

                    return Hosp
                    
                end

            end

        end
    end


end

#= 
output : 
    - character of the element h of the cartan in the given SU(d) irrep 
=#

# Multiplicity of irrep in SU(d)
function multiplicity_irrep(C_G::Vector, χ_G::Vector, χ_μ::Function, Dynkin::Vector)

    # C_G = vector with the class elements in the same order as χ_G
    # C_G[n] = Vector with all elements of the class written as a (d x d) unitary matrix

    N = length(C_G) # number of classes 
    Ntot = sum(length(C_G[n]) for n in 1:N) # number of elements 

    mult = sum(length(C_G[n])*χ_μ(C_G[n][1], Dynkin)*conj(χ_G[n]) for n in 1:N)/Ntot

    return (round(mult,digits=5))
end
# C_G = all group elements sorted by conjugacy classes. 
# χ_G = character tables. Column = conjugacy classes (same order as C_G), Rows = irreps (one-dimensional irreps first)
# L = irrep 

# test 

# step 1 : specify the Lie group 
N = 3

# step 2 : get the roots and weights 
Δ, Rp, Rn = get_roots_Cartanmatrix(N)
fw,fw_rb = get_fundamentalweights(N,Δ)

# step 3 : get the Weyl group
Weyl,detW = get_Weyl(N,Rp,Rn, Δ)

# test of the two ides for the signs of the determinant
W1,detW = get_Weyl(N,Rp,Rn, Δ)
W2,signW= get_Weyl_test(N,Rp,Rn, Δ)
println(W1)
println(detW)
println(W2)
println(signW)

Ds = [[1,0],[2,0],[1,1],[2,2],[0,3],[3,3],[4,1]]

dimensions = [Weyl_dim(d,Rp,fw,fw_rb) for d in Ds]

            #########
            # SU(3) #
            #########

Id = [1+0.0*im 0 0;0 1 0;0 0 1]
ω = exp(2*pi*im/3)
ϵ = exp(4*pi*im/9)
β = exp(2*pi*im/7)
μ1 = 0.5*(-1+sqrt(5))
μ2 = -0.5*(1+sqrt(5))

Ak(k) = collect(Diagonal([1,exp(2*pi*im/k),exp(-2*pi*im/k)]))

B = zeros(ComplexF64,3,3)
B[1,3]=-1
B[2,2]=-1
B[3,1]=-1

D = collect(Diagonal([ϵ,ϵ,ϵ*ω]))

E = zeros(ComplexF64,3,3)
E[1,2]=1
E[2,3]=1
E[3,1]=1 

F = zeros(ComplexF64,3,3)
F[1,1]=-1
F[2,3]=-ω
F[3,2]=-ω^2 

V = -im/sqrt(3)*[1 1 1; 1 ω ω^2 ;1 ω^2 ω]

W = 0.5*[-1 μ2 μ1 ; μ2 μ1 -1 ; μ1 -1 μ2]

X = -im/sqrt(3)*[1 1 ω^2 ; 1 ω ω ; ω 1 ω]

Y = collect(Diagonal([β , β^2 , β^4]))

Z = im/sqrt(7)*[β^4 - β^3 β^2 - β^5 β - β^6 ; 
β^2 - β^5 β^1 - β^6 β^4 - β^3 ;
β - β^6 β^4 - β^3 β^2 - β^5]

A = Ak(2)

C = Ak(3)

d = 3 
Δ, Rp, Rn = get_roots_Cartanmatrix(d)
fw,fw_rb = get_fundamentalweights(d,Δ)
Weyl,detW = get_Weyl(d,Rp,Rn, Δ)
χ_μ(unitary, dynkin) = Weyl_char(dynkin,Weyl,detW,Δ,Rp,Rn,fw,fw_rb,unitary)



#######################################
# G-invariant subspace multiplicity ? #
#######################################

#########################
# Test ArXiV:2605.26863 #
#########################

Sx = [0 1 0 ; 
1 0 1;
0 1 0]/sqrt(2)

Sy = [0 -im 0 ; 
im 0 -im;
0 im 0]/sqrt(2)

Sz = [1 0 0 ;
0 0 0;
0 0 -1]

Sx_t = Sx*Sz + Sz*Sx
Sy_t = Sy*Sz + Sz*Sy

gen1 = exp(-im*pi/2*Sx)
gen2 = exp(-im*pi/2*Sy)
gen3 = exp(-im*pi/2*Sx_t)
gen4 = exp(-im*pi/2*Sx_t)

Htest = [1/3 0 0 ;
0 -2/3 0;
0 0 1/3]*2
θ = 2*pi
gens5 = exp(-im*θ*Htest)

gens = [gen1, gen2,gen3,gen4, gens5]
group_all = get_group(gens,Id)
group_all = classify_group(group_all, [], [],Id,ordered=false)

# --> ZENITH generators are generators of Δ(96)

##################
# Lukin control  #
##################

function add_pulse_droid(ϕ1_ax, ϕ2_ax,θ)

    if ϕ1_ax == "x"
        ϕ1 = 0 
    elseif ϕ1_ax == "-x"
        ϕ1 = pi 
    elseif ϕ1_ax == "y"
        ϕ1 = pi/2 
    elseif ϕ1_ax == "-y"
        ϕ1 = 3*pi/2
    end

    if ϕ2_ax == "x"
        ϕ2 = 0 
    elseif ϕ2_ax == "-x"
        ϕ2 = pi 
    elseif ϕ2_ax == "y"
        ϕ2 = pi/2 
    elseif ϕ2_ax == "-y"
        ϕ2 = 3*pi/2
    end

    return exp(-im*θ*1/sqrt(2)*[0 exp(-im*ϕ1) 0;
    exp(im*ϕ1) 0 exp(-im*ϕ2);
    0 exp(im*ϕ2) 0])

end

gens = [
add_pulse_droid("-x","-x",pi), #1
add_pulse_droid("x","x",pi),
add_pulse_droid("-x","x",pi),
add_pulse_droid("x","-x",pi),

add_pulse_droid("-y","-y",pi), #5
add_pulse_droid("y","y",pi),
add_pulse_droid("-y","y",pi),
add_pulse_droid("y","-y",pi),

add_pulse_droid("-x","-x",pi/2), # 9
add_pulse_droid("x","x",pi/2),
add_pulse_droid("-x","x",pi/2),
add_pulse_droid("x","-x",pi/2),

add_pulse_droid("-y","-y",pi/2), # 13
add_pulse_droid("y","y",pi/2),
add_pulse_droid("-y","y",pi/2),
add_pulse_droid("y","-y",pi/2),

add_pulse_droid("x","y",pi/2), # 17
add_pulse_droid("-x","y",pi/2),
add_pulse_droid("x","-y",pi/2),
add_pulse_droid("-x","-y",pi/2),

] # 12 frames in Appendix B1

gens = [

add_pulse_droid("x","y",pi/4),
add_pulse_droid("-x","-y",pi/4),
add_pulse_droid("y","x",pi/4),
add_pulse_droid("-y","-x",pi/4),

] # 12 frames in Appendix B7

U = add_pulse_droid("x","y",pi/2)

Sz = [1 0 0; 0 0 0 ; 0 0 -1] 
adjoint(U)*Sz*U 

Lukin_all = get_group(gens, Id)
Δ6 = classify_group(Δ6_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ6_ρ1 = [multiplicity_irrep(Δ6, [1 for n in 1:length(Δ6)],χ_μ,D) for D in Ds]


###################
# Frobenius group #
###################

gens = [Y , E]
Frob_all = get_group(gens,Id)
Frob_classes = classify_group(Frob_all, [1 7 7 3 3], [1 3 3 7 7],Id)

χ_fro = Matrix{ComplexF64}(undef, 5,5)
χ_fro[1,:] = [1 1 1 1 1]
χ_fro[2,:] = [1 ω^2 ω 1 1]
χ_fro[3,:] = [1 ω ω^2 1 1]
χ_fro[4,:] = [3 0 0 -(1+sqrt(7))/2 (-1+sqrt(7))/2]
χ_fro[5,:] = [3 0 0 (-1+sqrt(7))/2 -(1+sqrt(7))/2]

Ds = [[0,0],
[1,0],[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Frob_ρ1 = [multiplicity_irrep(Frob_classes, χ_fro[1,:],χ_μ,D) for D in Ds]

#######
# Σ60 #
#######

gens = [E,A,W]
gens = [E,A*W]
gens = [E*A,W]

Σ60_all = get_group(gens, Id)
Σ60 = classify_group(Σ60_all,[1 2 3 5 5],[1 15 20 12 12],Id)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ60_ρ1 = [multiplicity_irrep(Σ60, χI[1,:],χ_μ,D) for D in Ds]

multiplicity_irrep(Σ60, χI[1,:],χ_μ,[1,1])

χ_μ(Σ60[end-3][1] , [1,1])

Σ60[end-3][1]
U_to_H(Σ60[end-3][1])

#####################
# Δ(6n^2) n=1-> 6 #
#####################

gens = [Ak(1),E,B]

Δ6_all = get_group(gens, Id)
Δ6 = classify_group(Δ6_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ6_ρ1 = [multiplicity_irrep(Δ6, [1 for n in 1:length(Δ6)],χ_μ,D) for D in Ds]


#####################
# Δ(6n^2) n=2-> Δ24 #
#####################

gens = [Ak(2),E,B]

Δ24_all = get_group(gens, Id)
Δ24 = classify_group(Δ24_all,[1 2 2 3 4],[1 3 6 8 6],Id)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ24_ρ1 = [multiplicity_irrep(Δ24, χO[1,:],χ_μ,D) for D in Ds]


#####################
# Δ(6n^2) n=3-> Δ54 #
#####################

gens = [Ak(3),E,B]

Δ54_all = get_group(gens, Id)
Δ54 = classify_group(Δ54_all,[1 1 1 1 2 3],[1 1 3 6 6 9],Id,ordered=false)

Ds = [[1,0],[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ54_ρ1 = [multiplicity_irrep(Δ54, [1.0 for n in 1:10] ,χ_μ,D) for D in Ds]

#####################
# Δ(6n^2) n=4-> Δ96 #
#####################

gens = [Ak(4),E,B]

Δ96_all = get_group(gens, Id)
Δ96 = classify_group(Δ96_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],[2,0],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ96_ρ1 = [multiplicity_irrep(Δ96, [1.0 for n in 1:10],χ_μ,D) for D in Ds]


#####################
# Δ(6n^2) n=5-> Δ150 #
#####################

gens = [Ak(5),E,B]

Δ150_all = get_group(gens, Id)
Δ150 = classify_group(Δ150_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ150_ρ1 = [multiplicity_irrep(Δ150, [1.0 for n in 1:length(Δ150)],χ_μ,D) for D in Ds]


###########
# Σ(36x3) #
###########

gens = [C,E,V]

gens = [C,V]
Σ36x3_all = get_group(gens, Id)
Σ36x3 = classify_group(Σ36x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ36x3_ρ1 = [multiplicity_irrep(Σ36x3, [1.0 for n in 1:14],χ_μ,D) for D in Ds]


#######################
# Δ(3xn^2) n=2 -> Δ12 #
#######################

gens = [Ak(2),E]


Δ12_all = get_group(gens, Id)

Δ12 = classify_group(Δ12_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ12_ρ1 = [multiplicity_irrep(Δ12, [1.0 for n in 1:10],χ_μ,D) for D in Ds]


#######################
# Δ(3xn^2) n=3 -> Δ27 #
#######################

gens = [Ak(3),E]

Δ27_all = get_group(gens, Id)

Δ27 = classify_group(Δ27_all,[],[],Id,ordered=false)

display(Δ27[1][1])
display(Δ27[end-1][1])
display(Δ27[end][1])
Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ27_ρ1 = [multiplicity_irrep(Δ27, [1.0 for n in 1:11],χ_μ,D) for D in Ds]

#######################
# Δ(3xn^2) n=4 -> Δ48 #
#######################

gens = [Ak(4),E]

Δ48_all = get_group(gens, Id)
Δ48 = classify_group(Δ48_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ48_ρ1 = [multiplicity_irrep(Δ48, [1.0 for n in 1:8],χ_μ,D) for D in Ds]

######################
# Δ(3xn^2) n=5 -> 75 #
######################

gens = [Ak(5),E]

Δ75_all = get_group(gens, Id)
Δ75 = classify_group(Δ75_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ75_ρ1 = [multiplicity_irrep(Δ75, [1.0 for n in 1:11],χ_μ,D) for D in Ds]

#######################
# Δ(3xn^2) n=6 -> 108 #
#######################

gens = [Ak(6),E]

Δ108_all = get_group(gens, Id)
Δ108 = classify_group(Δ108_all,[],[],Id,ordered=false)
display(Δ108[1][1])
display(Δ108[end-1][1])
display(Δ108[end][1])

Ds = [[1,0],
[1,1],
[2,0],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ108_ρ1 = [multiplicity_irrep(Δ108, [1.0 for n in 1:20],χ_μ,D) for D in Ds]

#######################
# Δ(3xn^2) n=7 -> 147 #
#######################

gens = [Ak(7),E]

Δ147_all = get_group(gens, Id)
Δ147 = classify_group(Δ147_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ147_ρ1 = [multiplicity_irrep(Δ147, [1.0 for n in 1:20],χ_μ,D) for D in Ds]

#######################
# Δ(3xn^2) n=8 -> 192 #
#######################

gens = [Ak(8),E]

Δ192_all = get_group(gens, Id)
Δ192 = classify_group(Δ192_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Δ192_ρ1 = [multiplicity_irrep(Δ192, [1.0 for n in 1:24],χ_μ,D) for D in Ds]


########
# Σ168 # 168 elements ! 
########

gens = [Y,E,Z]
gens = [Y,E]

Σ168_all = get_group(gens, Id)
Σ168 = classify_group(Σ168_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ168_ρ1 = [multiplicity_irrep(Σ168, [1.0 for n in 1:6],χ_μ,D) for D in Ds]


########
# Σ168 # subgroup ? 
########

gens = [Z,Y]

Σ168_all = get_group(gens, Id)
Σ168 = classify_group(Σ168_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ168_ρ1 = [multiplicity_irrep(Σ168, [1.0 for n in 1:6],χ_μ,D) for D in Ds]


###########
# Σ(72x3) # 216 element !
###########

gens = [C,E,V,X]

Σ72x3_all = get_group(gens, Id)
Σ72x3 = classify_group(Σ72x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ72x3_ρ1 = [multiplicity_irrep(Σ72x3, [1.0 for n in 1:16],χ_μ,D) for D in Ds]


###########
# Σ(72x3) # 216 element ! subgroup ? 
###########

z1 = (1/sqrt(5-2*sqrt(3)))*( (sqrt(3)-1/2) + (1-sqrt(3)/2)*im  )
angle(z1)

acos((sqrt(3)-1/2)/sqrt(5-2*sqrt(3)))
gens = [V,X]

(2+sqrt(3))/(sqrt(3)-2)
(1+2*sqrt(3))/(1-2*sqrt(3))

V^4
X^4
(V*X)^4
(V*X^2)^4

a_v = -(sqrt(234)/36 + sqrt(3)/(18*sqrt(78))) + 72*im/(18*sqrt(78))
b_v = -0.5/sqrt(18)*(sqrt(4-sqrt(3))*(0.5 + sqrt(3) + im*(1 + sqrt(3)/2)  ) + sqrt(4+sqrt(3))/sqrt(13)*( 4*sqrt(3) - 7/2 + im*( 2 + sqrt(3)/2 )  )  )

H_V = -pi/2*sqrt(3)/2*[2/9 a_v 0;
conj(a_v) 5/18 b_v;
0 conj(b_v) -0.5]

H_X = pi/2*[sqrt(3)/9 sqrt(78)/9 0;
sqrt(78)/9 -sqrt(3)/9 0;
0 0 0]

U_V = exp(im*H_V)
U_X = exp(im*H_X)

Σ72x3_all = get_group([U_V,U_X], Id)

Σ72x3_all = get_group(gens, Id)

Σ72x3 = classify_group(Σ72x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ72x3_ρ1 = [multiplicity_irrep(Σ72x3, [1.0 for n in 1:17],χ_μ,D) for D in Ds]


############
# Σ(216x3) # 
############

gens = [C,E,V,D]

Σ216x3_all = get_group(gens, Id)
Σ216x3 = classify_group(Σ216x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ216x3_ρ1 = [multiplicity_irrep(Σ216x3, [1.0 for n in 1:24],χ_μ,D) for D in Ds]

############
# Σ(216x3) # 648 element ! subgroup ?
############

gens = [V,D]

Σ216x3_all = get_group(gens, Id)
Σ216x3 = classify_group(Σ216x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ216x3_ρ1 = [multiplicity_irrep(Σ216x3, [1.0 for n in 1:27],χ_μ,D) for D in Ds]


############
# Σ(360x3) # 
############

gens = [A,E,W,F]

Σ360x3_all = get_group(gens, Id)
Σ360x3 = classify_group(Σ360x3_all,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ360x3_ρ1 = [multiplicity_irrep(Σ360x3, [1.0 for n in 1:24],χ_μ,D) for D in Ds]


############
# Σ(360x3) # subgroup ? 
############

gens = [A,E,W,F]

Σ360x3_all_sub = get_group(gens, Id)
Σ360x3_sub  = classify_group(Σ360x3_all_sub ,[],[],Id,ordered=false)

Ds = [[1,0],
[1,1],
[2,2],[0,3],
[3,3],[4,1],
[4,4],[5,2],[6,0]]

Σ360x3_sub_ρ1 = [multiplicity_irrep(Σ360x3_sub , [1.0 for n in 1:length(Σ360x3_sub)],χ_μ,D) for D in Ds]

########################################
# Figure inaccessible symmetries SU(3) #
########################################

Gs = [Frob_classes, Δ6 , Δ24 , Δ54 , Δ96 , Δ150 , Δ12 , Δ27 , Δ48 , Δ75 , Δ108 , Δ147 , Δ192 , Σ60 , Σ168 , Σ36x3 , Σ72x3 , Σ216x3 , Σ360x3]
names = [L"T_7" , L"\Delta(6)",  L"\Delta(24)", L"\Delta(54)", L"\Delta(96)", L"\Delta(150)", L"\Delta(12)", L"\Delta(27)", L"\Delta(48)", L"\Delta(75)", L"\Delta(108)", L"\Delta(147)", L"\Delta(192)",L"\Sigma(60)", L"\Sigma(168)", L"\Sigma(36\times 3)", L"\Sigma(72\times 3)",L"\Sigma(216\times 3)",L"\Sigma(360\times 3)"]
names_file = ["T_7" , "Δ(6)",  "Δ(24)", "Δ(54)", "Δ(96)", "Δ(150)", "Δ(12)", "Δ(27)", "Δ(48)", "Δ(75)", "Δ(108)", "Δ(147)", "Δ(192)","Σ(60)", "Σ(168)", "Σ(36x3)", "Σ(72x3)","Σ(216x3)","Σ(360x3)"]

fig_dir = "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/G-invariant  states SU3/"

for g in 1:length(Gs)
        
    # choose group # 

    G = Gs[g]
    name = names[g]
    namefile = names_file[g]

    fig = Figure()
    ax = Axis(fig[1,1],title=name,titlesize=25.0f0)

    a = [1,0]
    b = [cos(pi/3),sin(pi/3)]
    # axes 
    lines!([n for n in 0:8], [0.0 for n in 0:8], color=:black)
    lines!([n for n in 0:8], [tan(pi/3)*n for n in 0:8], color=:black)

    O = -(a+b)/2

    arrows!(Point2(O), Point2(a), shaftwidth=1.5, tipwidth=10,tiplength=10)
    arrows!(Point2(O), Point2(b), shaftwidth=1.5, tipwidth=10,tiplength=10)

    text!(Point2(O+a/2), text=L"\mu_1", fontsize=22.0f0, align=(:center, :top))
    text!(Point2(O+b/2), text=L"\mu_2", fontsize=22.0f0, align=(:right, :bottom))
    # grid 

    for n in 0:8
        lines!(Point2.([n*a+m*b for m in 0:8]),color=:grey, linewidth=0.2)
    end
    for m in 0:8
        lines!(Point2.([n*a+m*b for n in 0:8]),color=:grey, linewidth=0.2)
    end
    for n in 1:8 
        lines!(Point2.([n*a,n*b]),color=:grey, linewidth=0.2)
    end

    # where to look 

    Dynk_op = [0,0],[1,1],[2,2],[3,0],[0,3],[3,3],[4,1],[1,4],[4,4],[5,2],[6,0],[2,5],[0,6]
    colors = [:orange,:blue, :orange,:orange, :orange,:purple,:purple,:purple, :brown,:brown, :brown,:brown, :brown]

    for d in 1:length(Dynk_op)
        irrep = Dynk_op[d][1]*a+Dynk_op[d][2]*b 
        scatter!(Point2(irrep), markersize=60.0, color=(colors[d],0.3))
    end

    # points 

    for n in 0:8
        for m in 0:8
            irrep = n*a+m*b 
            scatter!(Point2(irrep),color=:white, strokewidth=1.0, strokecolor=:black, markersize=20.0)
        end
    end

    Dynkins = Matrix{Vector}(undef, 8,8)
    for i in 1:8
        for j in 1:8
            Dynkins[i,j] = [i-1,j-1]
        end
    end

    # calculate multiplicities 

    Mult = Matrix{Any}(undef, 8,8)
    for i in 1:8
        for j in 1:8
            Mult[i,j] = multiplicity_irrep(G, [1.0 for n in 1:length(G)],χ_μ,Dynkins[i,j])
        end
    end
    Mult = Int.(round.(Mult))
    mmax = maximum(Mult)
    mmin = 0

    CR = range(colorant"green", stop=colorant"red", length=mmax+1)

    for i in 1:8
        for j in 1:8
            irrep = Dynkins[i,j][1]*a+Dynkins[i,j][2]*b 
            if Mult[i,j] > 0 
                scatter!(Point2(irrep), markersize=20.0,color=(:red, 0.5))
                text!(Point2(irrep), text="$(Mult[i,j])",align=(:center,:center))
            else 
                scatter!(Point2(irrep), markersize=20.0,color=(:green, 0.5))
            end 
            
        end
    end
    colors_l = [:blue, :orange, :purple, :brown]

    legend_marker = [MarkerElement(marker=:dot,markersize=27.0,color=(colors_l[n],0.4)) for n in 1:4]
    legend_name = [L"1-\text{local}",L"2-\text{local}",L"3-\text{local}",L"4-\text{local}"]

    axislegend(ax,legend_marker, legend_name, position=:lt, framevisible=false)

    xlims!(nothing, 6.5)
    ylims!(nothing, 6.08)
    hidedecorations!(ax)
    hidespines!(ax)
    #display(fig)

    save(fig_dir*namefile*".pdf",fig)
end


#####################################
# Cayley graph and adjacency matrix #
#####################################

function get_GmodZ(G::Vector, Z::Vector)

    O = length(G) 
    GmodZ = Vector{Matrix}(undef, 0)

    for i in 1:O 

        element = G[i] 
        isthere = false 

        for n in 1:length(Z)
            index = findall(x -> norm(Z[n]*element - x) < 10^(-6), GmodZ)
            if length(index) > 0 
                isthere = true 
            end
        end

        if isthere == false 
            push!(GmodZ, element)
        end

    end

    # check 
    OmodZ = length(GmodZ) 
    if abs(OmodZ - O/3) > 0 
        println(abs(OmodZ - O/3))
        return "problem"
    else 
        return GmodZ 
    end

end

function get_adjacencymatrix(G::Vector, Γ::Vector)
    
    O = length(G)

    AdjMat = zeros(O,O)

    for i in 1:O
        for j in 1:O 

            p = G[j]*adjoint(G[i])

            index = findall(x -> norm(p-x) < 10^(-6) , Γ)
            if length(index) > 0 
                AdjMat[i,j] = index[1]
            end

        end
    end

    return AdjMat 
end

function get_adjacencymatrix_modZ(GmodZ::Vector, Γ::Vector,Z::Vector)
    
    O = length(GmodZ)

    AdjMat = zeros(O,O)

    for i in 1:O
        for j in 1:O 

            p = GmodZ[j]*adjoint(GmodZ[i])

            for n in 1:length(Z)

                index = findall(x -> norm(Z[n]*p-x) < 10^(-6) , Γ)

                if length(index) > 0 
                    AdjMat[i,j] = index[1]
                end

            end

        end
    end

    # sanity check 

    for i in 1:O 

        p = sum(AdjMat[i,n] for n in 1:O)
        test = abs(p - sum(n for n in 1:length(Γ)))
        if test > 0 
            println("problem")
        end

        p = sum(AdjMat[n,i] for n in 1:O)
        test = abs(p - sum(n for n in 1:length(Γ)))
        if test > 0 
            println("problem")
        end
    end

    return AdjMat 
end

using JLD 

data_dir = "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/Adjacency Matrix SU3/"

ξ(d,n) = exp(2*pi*im*n/d)
Z3 = [ξ(3,n) for n in 0:2] # center

##############
# Σ(72x3)/Z3 #
##############

gens = [V,X]

Σ72x3_all = get_group(gens, Id)
#Σ72x3 = classify_group(Σ72x3_all,[],[],Id,ordered=false)

Σ72x3modZ = get_GmodZ(Σ72x3_all,Z3)

AdjMat_Σ72x3modZ = get_adjacencymatrix_modZ(Σ72x3modZ, gens, Z3)

jldopen(data_dir*"Adjacency matrix : Σ72x3 - V,X - modZ.jld", "w") do file 
    write(file, "data", AdjMat_Σ72x3modZ)
end

##############
# Σ(36x3)/Z3 #
##############

gens = [V,E] #[E,V]
(V*E)^4

Σ36x3_all = get_group(gens, Id)
#Σ36x3 = classify_group(Σ36x3_all,[],[],Id,ordered=false)

Σ36x3modZ = get_GmodZ(Σ36x3_all,Z3)

AdjMat_Σ36x3modZ = get_adjacencymatrix_modZ(Σ36x3modZ, gens, Z3)

jldopen(data_dir*"Adjacency matrix : Σ36x3 - V,C - modZ.jld", "w") do file 
    write(file, "data", AdjMat_Σ36x3modZ)
end

##########
# Σ(168) #
##########

gens = [Y,Z]

Σ168_all = get_group(gens, Id)
#Σ168 = classify_group(Σ168_all,[],[],Id,ordered=false)

#Σ168modZ = get_GmodZ(Σ168_all,Z3)

AdjMat_Σ168 = get_adjacencymatrix(Σ168_all, gens)
#AdjMat_Σ168 = get_adjacencymatrix_modZ(Σ168_all, gens,Z3)

jldopen(data_dir*"Adjacency matrix : Σ168 - Y,Z.jld", "w") do file 
    write(file, "data", AdjMat_Σ168)
end


####################
# Σ(168) subgroups #
####################

d = [1 0 0;
0 -1 0;
0 0 -1]

a = [0 1 0;
0 0 1;
1 0 0]

b = -[0 0 1;
0 1 0;
1 0 0]

Zp = im/sqrt(7)*[β^2-β^5 β-β^6 β^4-β^3;
β-β^6 β^4-β^3 β^2-β^5;
β^4-β^3 β^2-β^5 β-β^6]

gens = [S,T] # generator of tetrahedron (S)^2, T^3
gens = [S*T,U] # generator of Octahedron (S x Up)^4, U^2

gens=[Y,Zp] # gens 

A3 = im/sqrt(7)*[β^2-β^5 β-β^6 β^4-β^3;
β-β^6 β^4-β^3 β^2-β^5;
β^4-β^3 β^2-β^5 β-β^6]

B3 = im/sqrt(7)*[β^3-β^6 β^3-β β-1;
β^2-1 β^6-β^5 β^6-β^2;
β^5-β^4 β^4-1 β^5-β^3]

w1 = B3*(A3*B3)^2*B3*A3*B3*A3

N = (28*(-342125 - 349668*β + 283769*β^2 + 9406*β^3 - 501928*β^4 + 287955*β^6))^(1/3)
x = (−64 - 15*β - 56*β^3 - 46*β^4 + 5*β^5 - 27*β^6)/N
y = (73*β + 156*β^2 + 46*β^3 + 12*β^4 + 137*β^5 + 115*β^6)/N
z = (15 + 34*β − 35*β^2 − 23*β^3 + 41*β^4 − 46*β^6)/N

w2 = [x y z;
z x y;
y z x]

adjoint(wtransfo)*a^2*b*wtransfo

wtransfo = w2*w1


gens = [A3,B3]
Σ168_all = get_group(gens, Id)

for n in 1:length(Σ168_all)

    if norm(adjoint(wtransfo)*a*wtransfo - Σ168_all[n]) < 10^(-5)
        println("ok")
    end

end

gens = [adjoint(wtransfo)*d*wtransfo,adjoint(wtransfo)*a*wtransfo] # tetra 
gens = [adjoint(wtransfo)*d*a*wtransfo,adjoint(wtransfo)*a^2*b*wtransfo] # octa

round.(adjoint(wtransfo)*d*a*wtransfo;digits=10)
round.(adjoint(wtransfo)*a^2*b*wtransfo;digits=10)
(a^2*b)^2
det(d*a)
round.(Σ168_all[4];digits=10)

Σ168_all = get_group(gens, Id)
for i in 1:length(Σ168_all)

    Σ168_all[i]  =round.(Σ168_all[i];digits=10)
    println("$i")
    display(Σ168_all[i])

end
#Σ168 = classify_group(Σ168_all,[],[],Id,ordered=false)

AdjMat_Σ168 = get_adjacencymatrix(Σ168_all, gens)

jldopen(data_dir*"Adjacency matrix : Σ168 - Y,Z.jld", "w") do file 
    write(file, "data", AdjMat_Σ168)
end


############
# Σ(27)/Z3 #
############

gens = [C,E]

Σ27_all = get_group(gens, Id)
#27 = classify_group(Σ27_all,[],[],Id,ordered=false)

Σ27modZ = get_GmodZ(Σ27_all,Z3)

AdjMat_Σ27modZ = get_adjacencymatrix_modZ(Σ27modZ, gens, Z3)

jldopen(data_dir*"Adjacency matrix : Δ27 - C,E.jld", "w") do file 
    write(file, "data", AdjMat_Σ27modZ)
end

############
# Σ(54)/Z3 #
############

gens = [C,E,B]

Σ54_all = get_group(gens, Id)

U = eigvecs(Ak(3))
adjoint(B)*adjoint(Ak(3))*B
adjoint(U)*Ak(3)*U
adjoint(U)*adjoint(B)*adjoint(Ak(3))*B*U


Σ54modZ = get_GmodZ(Σ54_all,Z3)

AdjMat_Σ54modZ = get_adjacencymatrix_modZ(Σ54modZ, gens, Z3)

jldopen(data_dir*"Adjacency matrix : Δ54 - A(3),E,B.jld", "w") do file 
    write(file, "data", AdjMat_Σ54modZ)
end

#####
# T #
#####

Sp = (1/3)*[-1 2 2;
2 -1 2;
2 2 -1]

Tp = [ω^2 0 0;
0 1 0;
0 0 ω]

gens = [Sp, Tp]

T_all = get_group(gens, Id)
#27 = classify_group(Σ27_all,[],[],Id,ordered=false)

#Σ27modZ = get_GmodZ(Σ27_all,Z3)

AdjMat_T = get_adjacencymatrix(T_all, gens)

jldopen(data_dir*"Adjacency matrix : T - Sp,Tp.jld", "w") do file 
    write(file, "data", AdjMat_T)
end


#####
# O #
#####

S = (1/3)*[-1 2 2;
2 -1 2;
2 2 -1]

T = [1 0 0;
0 ω^2 0;
0 0 ω]

Up = -[1 0 0 ;
0 0 ω^2;
0 ω 0]

U = -[1 0 0;
0 0 1;
0 1 0]

gens = [S*T, U]

O_all = get_group(gens, Id)
#27 = classify_group(Σ27_all,[],[],Id,ordered=false)

#Σ27modZ = get_GmodZ(Σ27_all,Z3)

AdjMat_O = get_adjacencymatrix(O_all, gens)

jldopen(data_dir*"Adjacency matrix : O - ST,U.jld", "w") do file 
    write(file, "data", AdjMat_O)
end


###################
# CodeSpace SU(3) #
###################

χ_Σ36 = zeros(ComplexF64,14,14)
Σ36_orders = [1 3 3 3 3 2 6 6 4 12 12 4 12 12]
Σ36_nbrs = [1 1 1 12 12 9 9 9 9 9 9 9 9 9]
χ_Σ36[1,:] = [1 1 1 1 1 1 1 1 1 1 1 1 1 1]
χ_Σ36[2,:] = [1 1 1 1 1 -1 -1 -1 im im im -im -im -im]
χ_Σ36[3,:] = [1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1]
χ_Σ36[4,:] = [1 1 1 1 1 -1 -1 -1 -im -im -im im im im]

χ_Σ72 = zeros(ComplexF64,16,16)
Σ72_orders = [1 3 3 3 2 6 6 4 12 12 4 12 12 12 12 4]
Σ72_nbrs = [1 1 1 24 9 9 9 18 18 18 18 18 18 18 18 18]
χ_Σ72[1,:] = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
χ_Σ72[2,:] = [1 1 1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 -1]
χ_Σ72[3,:] = [1 1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1]
χ_Σ72[4,:] = [1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 1 1 1]

χ_Σ216 = zeros(ComplexF64,24,24)
Σ216_orders = [1 3 3 3 2 6 6 4 12 12 3 3 9 9 9 9 9 9 18 18 18 18 18 18]
Σ216_nbrs = [1 1 1 24 9 9 9 54 54 54 72 72 12 12 12 12 12 12 36 36 36 36 36 36]
χ_Σ216[1,:] = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
χ_Σ216[2,:] = [1 1 1 1 1 1 1 1 1 1 ω ω^2 ω ω ω ω^2 ω^2 ω^2 ω ω ω ω^2 ω^2 ω^2]
χ_Σ216[3,:] = [1 1 1 1 1 1 1 1 1 1 ω^2 ω ω^2 ω^2 ω^2 ω ω ω ω^2 ω^2 ω^2 ω ω ω]

gens_24 = [Ak(2),E,B]
Δ24_all = get_group(gens_24, Id)
Δ24 = classify_group(Δ24_all,[1 2 2 3 4],[1 3 6 8 6],Id)

gens_168 = [Y,E,Z]
Σ168_all = get_group(gens_168, Id)
Σ168 = classify_group(Σ168_all,[],[],Id,ordered=false)

C1 = [Id]
C2 = ω*C1
C3 = ω^2*C1

C4 = append!([C,C^2,E,E^2],ω*[C,C^2,E,E^2],ω^2*[C,C^2,E,E^2])

C5 = append!([C*E, C^2*E,C*E^2, C^2*E^2],ω*[C*E, C^2*E,C*E^2, C^2*E^2],ω^2*[C*E, C^2*E,C*E^2, C^2*E^2])

C6 = [V^2,E*V^2,E^2*V^2,C*V^2,C^2*V^2,ω*C*E^2*V^2,ω*C^2*E*V^2,ω^2*C*E*V^2,ω^2*C^2*E^2*V^2]
C7 = ω*C6 
C8 = ω^2*C6

C9 = [V,ω*C*V,ω*C^2*V,ω*E*V,ω*E^2*V,C^2*E*V,C*E^2*V,ω*C*E*V,ω*C^2*E^2*V]
C10 = ω*C9 
C11 = ω^2*C9

C12 = [V^3,ω^2*C*V^3,ω^2*C^2*V^3,ω^2*E*V^3,ω^2*E^2*V^3,ω^2*C^2*E*V^3,ω^2*C*E^2*V^3,C*E*V^3,C^2*E^2*V^3]
C13 = ω*C12 
C14 = ω^2*C12

Σ36x3 = [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14]

D1 = C1 
D2 = ω*D1 
D3 = ω^2*D1

D4 = append!(C4,C5)

D5 = C6 
D6 = ω*D5 
D7 = ω^2*D5 

D8 = append!(C9,C12)
D9 = ω*D8 
D10 = ω^2*D8 

D11 = [X , ω^2*V^2*X , ω^2*E^2*V^2*E , ω^2*C^2*E*V^2*X , ω*C^2*X , ω*E*V^2*X , C^2*V^2*X , ω*C*X , ω^2*C*V^2*X , ω^2*C*E^2*V^2*X , ω*E*X , ω*C*E*V^2*X , ω*C^2*E*X , C^2*E^2*V^2*X , C*E^2*X , ω^2*C^2*E^2*X, C*E*X , E^2*X]
D12 = ω*D11 
D13 = ω^2*D11 

D14 = [V*X , C*V^3*X , ω*C*E*V^3*X , ω^2*C^2*E*V^3*X , ω^2*E*V*X , ω^2*C^2*V*X , C*E^2*V^3*X , ω*C^2*E^2*V^3*X , ω^2*C*E*V*X , ω*C*E^2*V*X , ω*E^2*V^3*X , ω^2*C*V*X , C^2*E*V*X , E^2*V*X , ω*C^2*V^3*X , ω*V^3*X , C^2*E^2*V*X , ω^2*E*V^3*X]
D15 = ω*D14 
D16 = ω^2*D14

Σ72x3 = [D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16]

E1 = D1 
E2 = D2 
E3 = D3 

E4 = D4

E5 = D5 
E6 = D6 
E7 = D7 

E8 = append!(D8,D11,D16)
E9 = ω*E8 
E10 = ω^2*E8 

Es = [E1[1] ,E2[1],E3[1],E4[1],E5[1],E6[1],E7[1],E8[1],E9[1] ,E10[1] ,
 E*D , E*D^2 , 
 D , ω*D , ω^2*D , D^2 , ω*D^2 , ω^2*D^2 ,
  V^2*D , ω*V^2*D,ω^2*V^2*D , V*D^2,ω*V*D^2, ω^2*V*D^2]

function sort_classes(C_sorted, C_unsorted)

    Conjugacy_classes = Vector{Vector}(undef, 0)

    for n in 1:length(C_sorted)

        element = C_sorted[n]

        m = 1 
        stop = false 

        while stop == false && m <= length(C_unsorted)

            index = findall(x -> norm(x - element)<10^(-4), C_unsorted[m])

            if length(index) == 0 
                m+=1
            else 
                stop = true 
            end
        end

        if m > length(C_unsorted) 
            return "trouble"
        else 
            push!(Conjugacy_classes,C_unsorted[m])
        end
    end

    return Conjugacy_classes
end

Σ216x3 = sort_classes(Es ,Σ216x3 )

Ns = [n for n in 1:16]

fig = Figure(;size=(600,400))

ax = Axis(fig[1,1], title=L"\Delta(24)",
 xticklabelsvisible=false , titlesize=25.0f0,
 ylabel=L"a_i",ylabelsize=25.0f0,ylabelrotation=0,
 yticklabelsize=15.0f0, titlegap = 13.0)
lines!(Ns , real.(round.([multiplicity_irrep(Δ24, [1.0 for n in 1:5], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1",linewidth=2.0)
axislegend(ax, position=:lt,labelsize=17.0f0)

ax = Axis(fig[2,1], title=L"\Sigma(168)" , titlesize=25.0f0, xlabel=L"N",xlabelsize=20.0f0,
 ylabel=L"a_i",ylabelsize=25.0f0,ylabelrotation=0,
 yticklabelsize=17.0f0, titlegap = 13.0, yticks=([0,1,2]))
lines!(Ns , real.(round.([multiplicity_irrep(Σ168, [1.0 for n in 1:6], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1",linewidth=2.0)
axislegend(ax, position=:lt,labelsize=17.0f0)

ax = Axis(fig[1,2],title=L"\Sigma(36\times 3)" ,
 xticklabelsvisible=false , titlesize=25.0f0,
 yticklabelsize=15.0f0, titlegap = 13.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3, χ_Σ36[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3",linewidth=2.0,linestyle=:dash)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_4",linewidth=2.0,linestyle=:dot)
axislegend(ax, position=:lt,labelsize=17.0f0, backgroundcolor=(:white, 0.4))

ax = Axis(fig[2,2],title=L"\Sigma(72\times 3)"  , titlesize=25.0f0, xlabel=L"N",xlabelsize=20.0f0,
 yticklabelsize=17.0f0, titlegap = 13.0)

lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3, χ_Σ72[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3",linewidth=2.0,linestyle=:dash)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_4",linewidth=2.0,linestyle=:dot)
axislegend(ax, position=:lt,labelsize=17.0f0, backgroundcolor=(:white, 0.4))

display(fig)
save(data_dir*"codespace qutrit.pdf",fig)


fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(168)",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.(round.([multiplicity_irrep(Σ168, [1.0 for n in 1:6], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1")

display(fig)


fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N)", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Delta(24))",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.(round.([multiplicity_irrep(Δ24, [1.0 for n in 1:5], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1")

display(fig)


fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(72\times 3)",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3, χ_Σ72[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1")
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2")
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3")
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_4")

axislegend(ax, position=:lt)
display(fig)

fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(36\times 3)",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.([multiplicity_irrep(Σ36x3, χ_Σ36[1,:], χ_μ,[n,0]) for n in Ns]), label=L"\rho_1")
lines!(Ns , real.([multiplicity_irrep(Σ36x3,  χ_Σ36[2,:], χ_μ,[n,0]) for n in Ns]), label=L"\rho_2")
lines!(Ns , real.([multiplicity_irrep(Σ36x3,  χ_Σ36[3,:], χ_μ,[n,0]) for n in Ns]), label=L"\rho_3")
lines!(Ns , real.([multiplicity_irrep(Σ36x3,  χ_Σ36[4,:], χ_μ,[n,0]) for n in Ns]), label=L"\rho_4")

axislegend(ax, position=:lt)
display(fig)


fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(36\times 3)",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.(([multiplicity_irrep(Σ36x3, χ_Σ36[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1")
lines!(Ns , real.(([multiplicity_irrep(Σ36x3,  χ_Σ36[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2")
lines!(Ns , real.(([multiplicity_irrep(Σ36x3,  χ_Σ36[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3")
lines!(Ns , real.(([multiplicity_irrep(Σ36x3,  χ_Σ36[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_4")

axislegend(ax, position=:lt)
display(fig)


fig = Figure()

ax = Axis(fig[1,1], xlabel=L"N", ylabel = L"a_i",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(216\times 3)",
ylabelrotation=0,titlesize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(Ns , real.(([multiplicity_irrep(Σ216x3, χ_Σ216[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1")
lines!(Ns , real.(([multiplicity_irrep(Σ216x3,  χ_Σ216[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2")
lines!(Ns , real.(([multiplicity_irrep(Σ216x3,  χ_Σ216[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3")

axislegend(ax, position=:lt)
display(fig)


[multiplicity_irrep(Σ216x3, χ_Σ216[1,:], χ_μ,[n,0]) for n in Ns]

############
# Appendix #
############

Ns = [n for n in 1:14]
CairoMakie.activate!()
fig = Figure(;size=(650,400))

ax = Axis(fig[1,1],title=L"\Sigma(36\times 3)" , ylabel = L"a_i",ylabelrotation=0,ylabelsize=25.0f0 , titlesize=25.0f0,
 yticklabelsize=15.0f0, titlegap = 13.0, xticks=([4,8,12,16]))
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3, χ_Σ36[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_1",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_2",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_3",color=:red,linewidth=2.0,linestyle=:dash)
lines!(Ns , real.(round.([multiplicity_irrep(Σ36x3,  χ_Σ36[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\rho_4",linewidth=2.0,linestyle=:dot)
axislegend(ax, position=:lt,labelsize=17.0f0, backgroundcolor=(:white, 0.4))

ax = Axis(fig[1,2],title=L"\Sigma(72\times 3)"  , titlesize=25.0f0, xlabel=L"N",xlabelsize=20.0f0,
 yticklabelsize=15.0f0, titlegap = 13.0, xticks=([4,8,12,16]))

lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3, χ_Σ72[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\sigma_1",linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\sigma_2",color=:red,linewidth=2.0)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\sigma_3",linewidth=2.0,linestyle=:dash)
lines!(Ns , real.(round.([multiplicity_irrep(Σ72x3,  χ_Σ72[4,:], χ_μ,[n,0]) for n in Ns])), label=L"\sigma_4",linewidth=2.0,linestyle=:dot)
axislegend(ax, position=:lt,labelsize=17.0f0, backgroundcolor=(:white, 0.4))

ax = Axis(fig[1,3], xlabel=L"N",
xlabelsize=25.0f0, ylabelsize=25.0f0, title=L"\Sigma(216\times 3)",
ylabelrotation=0,titlesize=25.0f0,
yticklabelsize=15.0f0, xticks=([4,8,12,16]),yticks=([0,1,2]))

lines!(Ns , real.(([multiplicity_irrep(Σ216x3, χ_Σ216[1,:], χ_μ,[n,0]) for n in Ns])), label=L"\tau_1",linewidth=2.0)
lines!(Ns , real.(([multiplicity_irrep(Σ216x3,  χ_Σ216[2,:], χ_μ,[n,0]) for n in Ns])), label=L"\tau_2",linewidth=2.0,linestyle=:dash)
lines!(Ns , real.(([multiplicity_irrep(Σ216x3,  χ_Σ216[3,:], χ_μ,[n,0]) for n in Ns])), label=L"\tau_3",linewidth=2.0,linestyle=:dot)
axislegend(ax, position=:lt,labelsize=17.0f0, backgroundcolor=(:white, 0.4))

display(fig)
save(data_dir*"codespace qutrit 2.pdf",fig)


#########
# SU(4) #
#########

ω = exp(2*pi/3*im)
β = exp(2*pi/7*im)
p = β + β^2 + β^4
q = β^3+β^5+β^6
s = β^2+β^5
t = β^3+β^4
u = β+β^6

F1 = [1 0 0 0;
0 1 0 0;
0 0 ω 0;
0 0 0 ω^2]

F2 = sqrt(1/3)*[1 0 0 sqrt(2);
0 -1 sqrt(2) 0;
0 sqrt(2) 1 0;
sqrt(2) 0 0 -1]

F3 = [sqrt(3)/2 1/2 0 0;
1/2 -sqrt(3)/2 0 0 ;
0 0 0 1;
0 0 1 0]

F2p = (1/3)*[3 0 0 0;
0 -1 2 2;
0 2 -1 2;
0 2 2 -1]

F3p = (1/4)*[-1 sqrt(15) 0 0;
sqrt(15) 1 0 0;
0 0 0 4;
0 0 4 0]

F4 = [0 1.0 0 0;
1 0 0 0;
0 0 0 -1;
0 0 -1 0]

S = [1 0 0 0;
0 β 0 0;
0 0 β^4 0;
0 0 0 β^2]

T = [1 0 0 0;
0 0 1 0;
0 0 0 1;
0 1 0 0]

W = 1/(im*sqrt(7))*[p^2 1 1 1;
1 -q -p -p;
1 -p -q -p;
1 -p -p -q]

R_ = sqrt(1/7)*[1 1 1 1;
2 s t u;
2 t u s;
2 u s t]

C = [1 0 0 0;
0 1 0 0;
0 0 ω 0;
0 0 0 ω^2]

D = [ω 0 0 0;
0 ω 0 0
0 0 ω 0;
0 0 0 1]

V = 1/(im*sqrt(3))*[im*sqrt(3) 0 0 0;
0 1 1 1;
0 1 ω ω^2;
0 1 ω^2 ω]

F = [0 0 -1 0;
0 1 0 0;
-1 0 0 0;
0 0 0 -1]

Fp = (1+im)/sqrt(2)*[1 0 0 0;
0 1 0 0;
0 0 0 1;
0 0 1 0]

Fpp = [0 1 0 0;
-1 0 0 0;
0 0 0 1;
0 0 -1 0]

A = (1+im)/sqrt(2)*[1 0 0 0;
0 im 0 0;
0 0 im 0;
0 0 0 1]

B = (1+im)/sqrt(2)*[1 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 -1]

Sp = (1+im)/sqrt(2)*[im 0 0 0;
0 im 0 0;
0 0 1 0;
0 0 0 1]

Tp = (1+im)/2*[-im 0 0 im;
0 1 1 0;
1 0 0 1;
0 -im im 0]

Rp = sqrt(1/2)*[1 im 0 0 ;
im 1 0 0;
0 0 im 1;
0 0 -1 -im]

A1 = [1 0 0 0;
0 1 0 0;
0 0 -1 0;
0 0 0 -1]

A2 = [0 0 1 0;
0 0 0 1;
1 0 0 0;
0 1 0 0]

A3 = [0 1 0 0;
1 0 0 0;
0 0 0 1;
0 0 1 0]

A4 = [0 0 1 0;
0 0 0 1;
1 0 0 0;
0 1 0 0]

S_su2 = [-1+im -1+im;
1+im -1-im]/2 

U_su2 = [1+im 0;
0 1-im]/sqrt(2)

V_su2 = [im/2 (1-sqrt(5))/4-im*(1+sqrt(5)/4);
-(1-sqrt(5))/4-im*(1+sqrt(5)/4) -im/2]

x1 = [1 1;
im -im]/sqrt(2)

x2 = [im im;
-1 1]/sqrt(2)

x3 = [-1 -1;
-1 1]/sqrt(2)

x4 = [im 1;
1 im]/sqrt(2)

x5 = [1 -1;
-im -im]/sqrt(2)

x6 = [im -im;
1 1]/sqrt(2)

function tensor_product_2x2(x,y)

    xy = Matrix{Any}(undef, 4,4)

    xy[1,1] = x[1,1]*y[1,1]
    xy[1,2] = x[1,1]*y[1,2]
    xy[2,1] = x[1,1]*y[2,1]
    xy[2,2] = x[1,1]*y[2,2]

    xy[1,3] = x[1,2]*y[1,1]
    xy[1,4] = x[1,2]*y[1,2]
    xy[2,3] = x[1,2]*y[2,1]
    xy[2,4] = x[1,2]*y[2,2]

    xy[3,1] = x[2,1]*y[1,1]
    xy[3,2] = x[2,1]*y[1,2]
    xy[4,1] = x[2,1]*y[2,1]
    xy[4,2] = x[2,1]*y[2,2]

    xy[3,3] = x[2,2]*y[1,1]
    xy[3,4] = x[2,2]*y[1,2]
    xy[4,3] = x[2,2]*y[2,1]
    xy[4,4] = x[2,2]*y[2,2]

    return xy 

end


G1 = [ComplexF64.(F1) , ComplexF64.(F2) , ComplexF64.(F3)]

G2 = [ComplexF64.(F1) , ComplexF64.(F2p) , ComplexF64.(F3p)]

G4 = [S,T,W]

G5 = [S,T,R_]

G7 = [F1 , F2 , F3,Fpp]

G8 = [F1 , F2p , F3p ,Fp]

G11 = [tensor_product_2x2(x1,x2), tensor_product_2x2(x1,transpose(x2)) , tensor_product_2x2(x3,x4) , tensor_product_2x2(x5,x6)]

G10_left = [S_su2, U_su2^2]
G10_right = [S_su2, U_su2^2]

G14_left = [S_su2, U_su2]
G14_right = [S_su2, U_su2]

K = [A1,A2,A3,A4]

G22 = [A1,A2,A3,A4,Tp]

G23 = [A1,A2,A3,A4,Tp,Rp^2]

G24= [A1,A2,A3,A4,T,R_]

d = 4

Δ,Rp,Rn = get_roots_Cartanmatrix(d)

fw,fw_rb = get_fundamentalweights(d,Δ)

Weyl,detW = get_Weyl(d,Rp,Rn, Δ)

χ_μ(unitary, dynkin) = Weyl_char(dynkin,Weyl,detW,Δ,Rp,Rn,fw,fw_rb,unitary)

Id = [1.0+0.0*im 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 1]

# dimensions : 

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

for d in 1:length(Ds)
    println(round(Weyl_dim(Ds[d], Rp,fw,fw_rb),digits=2))
end


######
# G1 #
######

G1_all = get_group(G1, Id)
G1_class = classify_group(G1_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G1_ρ1 = [multiplicity_irrep(G1_class, [1.0 for n in 1:20],χ_μ,D) for D in Ds]

######
# G2 #
######

G2_all = get_group(G2, Id)
G2_class = classify_group(G2_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G2_ρ1 = [multiplicity_irrep(G2_class, [1.0 for n in 1:20],χ_μ,D) for D in Ds]

######
# G4 #
######

G4_all = get_group(G4, Id)

G4_class = classify_group(G4_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G4_ρ1 = [multiplicity_irrep(G4_class, [1.0 for n in 1:20],χ_μ,D) for D in Ds]

######
# G5 #
######

G5_all = get_group(G5, Id)
G5_class = classify_group(G5_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G5_ρ1 = [multiplicity_irrep(G5_class, [1.0 for n in 1:length(G5_class)],χ_μ,D) for D in Ds]

######
# G7 #
######

G7_all = get_group(G7, Id)
G7_class = classify_group(G7_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G7_ρ1 = [multiplicity_irrep(G7_class, [1.0 for n in 1:24],χ_μ,D) for D in Ds]

######
# G8 #
######

G8_all = get_group(G8, Id)
G8_class = classify_group(G8_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G8_ρ1 = [multiplicity_irrep(G8_class, [1.0 for n in 1:28],χ_μ,D) for D in Ds]

######
# G10 #
######

G10_all_left = get_group(G10_left, [1 0+0.0*im; 0 1],order=288*3)
G10_all_right = get_group(G10_right, [1 0+0.0*im; 0 1],order=288*3)

G10_all = Vector{typeof(Id)}(undef,0)
for n in 1:length(G10_all_left)
    for m in 1:length(G10_all_right)
        push!(G10_all, tensor_product_2x2(G10_all_left[n],G10_all_right[m]))
    end
end

function check_duplicate(G)

    new= []
    for g in G 
        index = findall(x-> norm(x-g)<10^(-5), new)
        if length(index) == 0
            push!(new, g)
        end
    end
    return new
end

G10_all = check_duplicate(G10_all)

G10_class = classify_group(G10_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G10_ρ1 = [multiplicity_irrep(G10_class, [1.0 for n in 1:29],χ_μ,D) for D in Ds]


######
# G14 #
######

G14_all_left = get_group(G14_left, [1 0+0.0*im; 0 1],order=288*3)
G14_all_right = get_group(G14_right, [1 0+0.0*im; 0 1],order=288*3)

G14_all = Vector{typeof(Id)}(undef,0)
for n in 1:length(G14_all_left)
    for m in 1:length(G14_all_right)
        push!(G14_all, tensor_product_2x2(G14_all_left[n],G14_all_right[m]))
    end
end

function check_duplicate(G)

    new= []
    for g in G 
        index = findall(x-> norm(x-g)<10^(-5), new)
        if length(index) == 0
            push!(new, g)
        end
    end
    return new
end

G14_all = check_duplicate(G14_all)

G14_class = classify_group(G14_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G14_ρ1 = [multiplicity_irrep(G14_class, [1.0 for n in 1:length(G14_class)],χ_μ,D) for D in Ds]

######
# G11 #
######

G11_all = get_group(G11, Id,order=288*3)

G11_class = classify_group(G11_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G11_ρ1 = [multiplicity_irrep(G11_class, [1.0 for n in 1:28],χ_μ,D) for D in Ds]



#####
# K #
#####

K_all = get_group(K, Id)
K_class = classify_group(K_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

K_ρ1 = [multiplicity_irrep(K_class, [1.0 for n in 1:28],χ_μ,D) for D in Ds]

######
# G22 #
######

G22_all = get_group(G22, Id)
G22_class = classify_group(G22_all,[],[],Id,ordered=false)

Ds = [[1,0,1],[2,0,2],[2,1,0],[0,1,2],[0,2,0]]

G22_ρ1 = [multiplicity_irrep(G22_class, [1.0 for n in 1:28],χ_μ,D) for D in Ds]


####################
# Young's Tableaux #
####################


function Dynkin2Young(D::Vector, d::Int; D_i=0, color=false)


    function DrawSquare(x,y,a)

        Pts = Point2.([[x,y] , [x,y+a],[x+a,y+a],[x+a,y],[x,y] ])

        lines!(Pts, color=:black,linewidth=3.0)
    end

    function DrawColoredSquare(x,y,a,col)

        Pts = Point2.([[x,y] , [x,y+a],[x+a,y+a],[x+a,y],[x,y] ])

        vertices = [x y;
        x y+a;
        x+a y+a;
        x+a y]
        faces=[1 2 3;
        3 4 1]

        mesh!(vertices, faces,color=(col,0.8))

    end

    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())

    if D_i > 0 

        for n in 1:D_i 
            for m in 1:d 

                DrawSquare(n,-m,-1)

            end
        end
    end

    all_colors = [
        RGB(153/255,204/255,255/255),
        RGB(255/255,153/255,153/255),
        RGB(153/255,255/255,153/255),
        RGB(255/255,204/255,153/255),
        RGB(204/255,153/255,255/255)
    ]

    for n in 1:d-1 

        if n == 1
            start = D_i 
        else
            start = D_i + sum(D[m] for m in 1:n-1)
        end

        for m in 1:D[n]

            for k in 1:d-n 
                if color == false 
                    
                elseif color == true 
                    DrawColoredSquare(start+m,-k,-1,all_colors[k])
                else
                    DrawColoredSquare(start+m,-k,-1,color[k])
                end
            end

        end
    end


    for n in 1:d-1 

        if n == 1
            start = D_i 
        else
            start = D_i + sum(D[m] for m in 1:n-1)
        end

        for m in 1:D[n]

            for k in 1:d-n 

                DrawSquare(start+m,-k,-1)

            end

        end
    end


    hidespines!(ax)
    hidedecorations!(ax)

    return fig 
end

function Dynkin2Young_TP(D1::Vector,D2::Vector, d::Int;color=false)

    function DrawSquare(x,y,a)

        Pts = Point2.([[x,y] , [x,y+a],[x+a,y+a],[x+a,y],[x,y] ])

        lines!(Pts, color=:black,linewidth=3.0)
    end

    function DrawColoredSquare(x,y,a,col)

        Pts = Point2.([[x,y] , [x,y+a],[x+a,y+a],[x+a,y],[x,y] ])

        vertices = [x y;
        x y+a;
        x+a y+a;
        x+a y]
        faces=[1 2 3;
        3 4 1]

        mesh!(vertices, faces,color=(col,0.8))

    end

    function Draw_Oplus(x,y,a)
        center = [x,y] + [a,a]/2
        circ = Point2.([center + 0.3*a*[cos(n),sin(n)] for n in LinRange(0,2*pi,100)])

        lines!(circ, color=:black, linewidth=2.0)

        lines!(Point2.([center+0.3*[a,0] , center-0.3*[a,0]]),color=:black, linewidth=2.0)
        lines!(Point2.([center+0.3*[0,a] , center-0.3*[0,a]]),color=:black, linewidth=2.0)

    end

    function Draw_Otimes(x,y,a)
        center = [x,y] + [a,a]/2
        circ = Point2.([center + 0.3*a*[cos(n),sin(n)] for n in LinRange(0,2*pi,100)])

        lines!(circ, color=:black, linewidth=2.0)

        lines!(Point2.([center+0.3*[a,a]/sqrt(2) , center-0.3*[a,a]/sqrt(2)]),color=:black, linewidth=2.0)
        lines!(Point2.([center+0.3*[a,-a]/sqrt(2) , center-0.3*[a,-a]/sqrt(2)]),color=:black, linewidth=2.0)

    end

    all_colors = [
        RGB(153/255,204/255,255/255),
        RGB(255/255,153/255,153/255),
        RGB(153/255,255/255,153/255),
        RGB(255/255,204/255,153/255),
        RGB(204/255,153/255,255/255)
    ]

    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())


    for n in 1:d-1 

        if n == 1
            start = 0 
        else
            start = sum(D1[m] for m in 1:n-1)
        end

        for m in 1:D1[n]

            for k in 1:d-n 

                DrawSquare(start+m,-k,-1)

            end

        end
    end

    Draw_Otimes(sum(D1[m] for m in 1:d-1),-2,1)
    
    
    for n in 1:d-1 

        if n == 1
            start = sum(D1[m] for m in 1:d-1) + 1
        else
            start = sum(D1[m] for m in 1:d-1) + 1 + sum(D2[m] for m in 1:n-1)
        end

        for m in 1:D2[n]

            for k in 1:d-n 

                if color == false
                elseif color == true
                    DrawColoredSquare(start+m,-k,-1,all_colors[k])
                else
                    DrawColoredSquare(start+m,-k,-1,color[k])
                end

            end

        end
    end


    for n in 1:d-1 

        if n == 1
            start = sum(D1[m] for m in 1:d-1) + 1
        else
            start = sum(D1[m] for m in 1:d-1) + 1 + sum(D2[m] for m in 1:n-1)
        end

        for m in 1:D2[n]

            for k in 1:d-n 

                DrawSquare(start+m,-k,-1)

            end

        end
    end

    hidespines!(ax)
    hidedecorations!(ax)

    return fig 
end

all_colors = [
        RGB(153/255,204/255,255/255),
        RGB(255/255,153/255,153/255),
        RGB(255/255,255/255,255/255),
        RGB(255/255,204/255,153/255),
        RGB(204/255,153/255,255/255)
    ]


fig=Dynkin2Young([1,0,0],4,D_i=0,color=all_colors[2:4])
save(fig_dir*"(0,...0,1).svg",fig)

Dynkin2Young_TP([1,1],[1,1],3)

CairoMakie.activate!()
fig=  Dynkin2Young_TP([1,0,0,1],[1,0,0,1],5, color=all_colors)
save(fig_dir*"(1,...,1) TP (0,...,1).svg",fig)



#################################
# G-anticoherence, some figures #
#################################

        # spectrum of U_t(su(3))

CairoMakie.activate!()

fig = Figure()
ax = Axis(fig[1,1],titlesize=25.0f0)

a = [1,0]
b = [cos(pi/3),sin(pi/3)]

# axes 
lines!([n for n in 0:10], [0.0 for n in 0:10], color=:black)
lines!([n for n in 0:10], [tan(pi/3)*n for n in 0:10], color=:black)

O = -(a+b)/2

arrows!(Point2(O), Point2(a), shaftwidth=1.5, tipwidth=10,tiplength=10)
arrows!(Point2(O), Point2(b), shaftwidth=1.5, tipwidth=10,tiplength=10)

text!(Point2(O+a/2), text=L"p", fontsize=22.0f0, align=(:center, :top))
text!(Point2(O+b/2), text=L"q", fontsize=22.0f0, align=(:right, :bottom))

# grid 

for n in 0:10
    lines!(Point2.([n*a+m*b for m in 0:10]),color=:grey, linewidth=0.2)
end
for m in 0:10
    lines!(Point2.([n*a+m*b for n in 0:10]),color=:grey, linewidth=0.2)
end
for n in 1:10
    lines!(Point2.([n*a,n*b]),color=:grey, linewidth=0.2)
end


# convex hull 

for t in 1:7 

    r1 = [1,-2]
    r2 = [-2,1]

    Dyn = [t,t]
    Dyn_1 = Dyn + t/2*r1
    Dyn_2 = Dyn + t/2*r2

    irrep = Dyn[1]*a + Dyn

    lines!(Point2.([Dyn_1[1]*a + Dyn_1[2]*b,Dyn[1]*a + Dyn[2]*b,Dyn_2[1]*a + Dyn_2[2]*b]), linestyle=:dash, color=:black)

    text!(Point2(Dyn_2[1]*a + Dyn_2[2]*b  + 0.2*r2), text="t=$t",align=(:right,:bottom),fontsize=15.0f0)
end

# points 

for n in 0:10
    for m in 0:10
        irrep = n*a+m*b 
        scatter!(Point2(irrep),color=:white, strokewidth=1.0, strokecolor=:black, markersize=20.0)
    end
end

# where to look 

Dynk_op = [0,0],[1,1],[2,2],[3,0],[0,3],[3,3],[4,1],[1,4],[4,4],[5,2],[6,0],[2,5],[0,6] ,[7,1],[1,7],[6,3],[3,6],[5,5] , [9,0],[0,9], [8,2],[2,8],[7,4],[4,7],[6,6] , [10,1],[1,10],[9,3],[3,9],[8,5],[5,8],[7,7]
#colors = [:orange,:blue, :orange,:orange, :orange,:purple,:purple,:purple, :brown,:brown, :brown,:brown, :brown, :green, :green, :green, :green, :green]
colors = [:red for n in 1:length(Dynk_op)]
for d in 1:length(Dynk_op)
    irrep = Dynk_op[d][1]*a+Dynk_op[d][2]*b 
    scatter!(Point2(irrep), markersize=20.0, color=(colors[d],0.3))
end

xlims!(nothing, 8.5)
ylims!(nothing, 8.2)
hidedecorations!(ax)
hidespines!(ax)

display(fig)

#save("/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/convex hull.pdf",fig)