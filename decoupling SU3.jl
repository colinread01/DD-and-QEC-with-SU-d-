

###################
# Qutrit sequence #
###################


using QuantumOptics 
using LinearAlgebra 
using JLD
using Colors
using GLMakie 

function GellMann()

    basis = GenericBasis(3)

    λ1 = zeros(ComplexF64,3,3)
    λ1[1,2] = 1.0 
    λ1[2,1] = 1.0

    λ2 = zeros(ComplexF64,3,3)
    λ2[1,2] = -im
    λ2[2,1] = im

    λ3 = zeros(ComplexF64,3,3)
    λ3[1,1] = 1
    λ3[2,2] = -1

    λ4 = zeros(ComplexF64,3,3)
    λ4[1,3] = 1.0 
    λ4[3,1] = 1.0

    λ5 = zeros(ComplexF64,3,3)
    λ5[1,3] = -im
    λ5[3,1] = im

    λ6 = zeros(ComplexF64,3,3)
    λ6[2,3] = 1
    λ6[3,2] = 1

    λ7 = zeros(ComplexF64,3,3)
    λ7[2,3] = -im
    λ7[3,2] = im

    λ8 = zeros(ComplexF64,3,3)
    λ8[1,1] = 1/sqrt(3)
    λ8[2,2] = 1/sqrt(3)
    λ8[3,3] = -2/sqrt(3)

    return DenseOperator(basis,λ1), DenseOperator(basis,λ2) , DenseOperator(basis,λ3), DenseOperator(basis,λ4) , DenseOperator(basis,λ5) , DenseOperator(basis,λ6) , DenseOperator(basis,λ7) , DenseOperator(basis,λ8)

end

function TwoBody(i,j,N,Oi,Oj,Id)

    if i == 1 

        H  = Oi 
        for n in 2:j-1 
            H = tensor(H,Id)
        end
        H = tensor(H,Oj)
        for n in j+1:N 
            H = tensor(H,Id)
        end

    else 

        H = Id 
        for n in 2:i-1 
            H  = tensor(H,Id)
        end
        H  = tensor(H,Oi)
        for n in i+1:j-1 
            H = tensor(H,Id)
        end
        H = tensor(H,Oj)
        for n in j+1:N 
            H = tensor(H,Id)
        end

    end

    return H 

end

function OneBody(i,N,Oi,Id)

    
    if i == 1 

        H  = Oi 
        for k in 2:N 
            H = tensor(H,Id)
        end

    else 

        H = Id 
        for k in 2:i-1 
            H  = tensor(H,Id)
        end
        H  = tensor(H,Oi)
        for k in i+1:N
            H = tensor(H,Id)
        end

    end

    return H 

end

function ThreeBody(i,j,k,N,Oi,Oj,Ok,Id)

    if i == 1 

        H  = Oi 
        for n in 2:j-1 
            H = tensor(H,Id)
        end
        H = tensor(H,Oj)
        for n in j+1:k-1
            H = tensor(H,Id)
        end
        H = tensor(H,Ok)
        for n in k+1:N
            H = tensor(H,Id)
        end
    else 

        H = Id 
        for n in 2:i-1 
            H  = tensor(H,Id)
        end
        H  = tensor(H,Oi)
        for n in i+1:j-1 
            H = tensor(H,Id)
        end
        H = tensor(H,Oj)
        for n in j+1:k-1
            H = tensor(H,Id)
        end
        H = tensor(H,Ok)
        for n in k+1:N 
            H = tensor(H,Id)
        end

    end

    return H 

end

function print_sequence(seq::Vector, colors::Vector, Ncol::Int)


    function single_pulse(color, pos_x , pos_y , height,τ)

        points = [[pos_x-τ , pos_y] , [pos_x-τ , pos_y+height],[pos_x , pos_y+height], [pos_x , pos_y]]

        lines!(Point2.(points), color=color, linewidth=2.0)
        band!(Point2.([[pos_x-τ,pos_y ],[pos_x,pos_y]]),Point2.(points[2:3]), color=(color,0.5))

    end


    height = 0.55
    τ = 0.2

    Nrow = floor(Int,length(seq)/Ncol)

    leftover = ceil(length(seq)/Ncol) - floor(Int,length(seq)/Ncol)

    fig = Figure(;size=(800,100*ceil(length(seq)/Ncol)))

    ax = Axis(fig[1,1])
    hidespines!(ax)
    hidedecorations!(ax)

    for m in 1:Nrow

        pos_y = -(m-1)

        for n in 1:Ncol

            single_pulse(colors[seq[Ncol*(m-1)+n]] , n , pos_y , height , τ)

        end

        lines!(Point2.([[0-τ,pos_y],[Ncol+0.5,pos_y]]), color=:black, linewidth=1.5)
        
        if m < Nrow || leftover > 0

            lines!(Point2.([[Ncol+0.5,pos_y],[Ncol+1.5,pos_y]]), color=:black, linewidth=1.5, linestyle=:dash)

        else

            arrows2d!(Point2([Ncol+0.5,pos_y]) , Point2([1.0,0]) , color=:black , shaftwidth = 0.00 , shaftlength=0.2 , tipwidth = 0.1 , tiplength=0.1,normalize=false,lengthscale=1.0)

        end

    end

    if leftover > 0

        pos_y = -Nrow

        for n in 1:length(seq) - Nrow*Ncol

            
            single_pulse(colors[seq[Ncol*Nrow+n]] , n , pos_y , height , τ)
        end

        lines!(Point2.([[Ncol+0.5,pos_y],[Ncol+1,pos_y]]), color=:black, linewidth=1.5, linestyle=:dash)
        arrows2d!([Ncol+0.5,pos_y] , [0.5,0] , color=:black , shaftwidth = 0.1 , shaftlength=1.0,tipwidth = 0.1 , tiplength=0.5, normalize=false,lengthscale=1.0)

    end 

    ylims!(-(ceil(length(seq)/Ncol)-1)-0.2,height+0.2)

    return fig 
end

function Hierholzer_cycle_random(start_vertex::Int, AdjVec::Vector)

    curr_vertex = start_vertex 

    curr_path = [curr_vertex]

    stop = false 

    while stop == false
        if length(AdjVec[curr_vertex]) == 0
            stop = true 
        else
            old_vertex = curr_vertex
            rand_inx = rand(1:length(AdjVec[curr_vertex]))
            curr_vertex = AdjVec[curr_vertex][rand_inx]
            push!(curr_path, curr_vertex)
            splice!(AdjVec[old_vertex],rand_inx)
        end
    end

    return curr_path[2:end], AdjVec
end

function Hierholzer_random(AdjMat::Matrix)

    N = floor(Int,sqrt(length(AdjMat)))
    AdjVec = Vector{Vector}(undef, N )

    for i in 1:N 
        AdjVec[i] = []
        for j in 1:N 
            if AdjMat[i,j] > 10^(-6)
                push!(AdjVec[i], j)
            end
        end
    end

    # AdjVec tells you to which vertex each vertex is connected (but does not care about the type of vertex)

    start = 1 
    circuit = [1]
    loop, AdjVec = Hierholzer_cycle_random(start, AdjVec)
    circuit = vcat(circuit, loop)
    bigstop = false 
    #println(circuit)

    while bigstop == false 

        edges_left = 0 
        for i in 1:N 
            edges_left += length(AdjVec[i])
        end

        if edges_left == 0 
            bigstop = true 
        end

        if bigstop == false 

            stop = false 

            n = rand(1:length(circuit))

            while stop == false 

                start = circuit[n]

                if length(AdjVec[start]) != 0 

                    stop = true 
                    loop, AdjVec = Hierholzer_cycle_random(start, AdjVec)

                    circuit1 = circuit[1:n]
                    circuit2 = circuit[(n+1):end]
                    circuit = vcat(circuit1, loop)
                    circuit = vcat(circuit, circuit2)
                else 
                    if n == length(circuit) 
                        n = 1
                    else 
                        n+=1
                    end
                end
            end

        end
    end

    return circuit 
    
end

function path_to_sequence(AdjMat::Matrix, path::Vector)

    A = copy(AdjMat)

    seq = []

    for i in 1:length(path)-1
        push!(seq, Int(A[path[i],path[i+1]]))
        A[path[i],path[i+1]] = 0
    end

    if length(findall(x -> x == 0,seq))==0 
        return seq 
    else     
        println("problem")    
        return "X"
    end

end

function Hamiltonian_disorder(N,λ,Id)

    n = [rand()-0.5 for n in 1:8]
    n = n/norm(n)

    Op = sum(n[k]*λ[k] for k in 1:8)
    H = (rand()-0.5)*OneBody(1,N,Op,Id)

    for m in 2:N 

        n = [rand()-0.5 for n in 1:8]
        n = n/norm(n)
        Op = sum(n[k]*λ[k] for k in 1:8)

        H += (rand()-0.5)*OneBody(m,N,Op,Id)
    end

    return H 
end

function Hamiltonian_interaction(N,λ,Id)

    H = 0.0*OneBody(1,N,Id,Id)

    for k in 1:N-1 
        for p in k+1:N 

            J = rand(8,8)
            J[:,:] .-= 0.5
            Inv = collect(Diagonal([1 for n in 1:8]))
            J = J - tr(J)*Inv/8
            println(tr(J))
            J = J/norm(J)

            H += (rand()-0.5)*sum(sum(J[i,j]*TwoBody(k,p,N, λ[i], λ[j],Id) for i in 1:8) for j in 1:8)
        end
    end

    return H 


end

function propagate_ideal(N,H, Seq, Gens,τ)

    U0 = exp(dense(-im*H*τ))
    U = exp(dense(-im*0.0*H))

    for n in 1:length(Seq)
        pulse = Gens[Seq[n]]
        for m in 2:N 
            pulse = tensor(pulse,Gens[Seq[n]])
        end
        U = pulse*U0*U
    end

    return log10(sqrt(abs(1 - (1/3^N)*abs(tr(U)))))

end

function propagate_finite(N,H, Seq, Gens,τ0,τp,Id)

    U0 = exp(dense(-im*H*τ0))
    U = exp(dense(-im*0.0*H))

    for n in 1:length(Seq)
        pulse = exp(-im*τp*(sum(OneBody(i,N,Gens[Seq[n]],Id) for i in 1:N)/τp + H))
        U = pulse*U0*U
    end

    return log10(sqrt(abs(1 - (1/3^N)*abs(tr(U)))))

end

function propagate_nodd(N,H, M, τ)

    U = exp(dense(-im*H*M*τ))

    return log10(sqrt(1 - (1/3^N)*abs(tr(U))))
end

function propagate_errors(Gens_err, Seq, Gens, Ufinal) 

    U = exp(dense(-im*0.0*Gens_err[Seq[1][1]]))

    for n in 1:length(Seq)

        pulse = exp(dense(-im*(Gens[Seq[n][1]] + norm(Gens[Seq[n][1]],Inf)*Gens_err[Seq[n][1]])))

        for m in 2:length(Seq[n])
            
            newpulse = exp(dense(-im*(Gens[Seq[n][m]] + norm(Gens[Seq[n][m]],Inf)*Gens_err[Seq[n][m]])))

            pulse = newpulse*pulse

        end

        U = pulse*U

    end
    
    return log10(sqrt(abs(1 - (1/3)*abs(tr(U*dagger(Ufinal))))))

end

function Hamiltonian_disorder_NV(N,λ,Id)

    H = 0.0*OneBody(1,N,Id,Id)

    for n in 1:N 
        a = rand()-0.5 
        H+= a*OneBody(n,N,(λ[3]+sqrt(3)*λ[8]),Id)
    end

    return H 

end

function Hamiltonian_disorder_NV_quad(N,λ,Id)

    H = 0.0*OneBody(1,N,Id,Id)

    for n in 1:N 
        a = rand()-0.5 
        H+= a*OneBody(n,N,(λ[3]+sqrt(3)*λ[8])^2,Id)
    end

    return H 

end

function Hamiltonian_interaction_NV(N,λ,Id)

    H = 0.0*OneBody(1,N,Id,Id)

    for k in 1:N-1 
        for p in k+1:N 

            J = rand()-0.5 

            H += J*TwoBody(k,p,N,λ[1],λ[1],Id)+ J*TwoBody(k,p,N,λ[2],λ[2],Id)+ J*TwoBody(k,p,N,λ[6],λ[6],Id)+ J*TwoBody(k,p,N,λ[7],λ[7],Id)- J*TwoBody(k,p,N,λ[3]+sqrt(3)*λ[8],λ[3]+sqrt(3)*λ[8],Id)
        end
    end

    return H 


end

# nested sequences 

function propagate_nested(N,H, Seq, Gens,τ,Ufinal)

    U0 = exp(dense(-im*H*τ))
    U = exp(dense(-im*0.0*H))

    for n in 1:length(Seq)

        pulse = Gens[Seq[n][1]]
        for m in 2:N 
            pulse = tensor(pulse,Gens[Seq[n][1]])
        end

        for m in 2:length(Seq[n])
            
            newpulse = Gens[Seq[n][m]]
            for j in 2:N 
                newpulse = tensor(newpulse,Gens[Seq[n][m]])
            end
            pulse = newpulse*pulse
        end

        U = pulse*U0*U
    end

    Utrg = copy(Ufinal)
    for m in 2:N 
        Utrg = tensor(Utrg,Ufinal)
    end
    return log10(sqrt(1 - (1/3^N)*abs(tr(U*dagger(Utrg)))))

end

# zero field sequences 

function Hamiltonian_interaction_NV_ZF(N,λ,Id)

    H = 0.0*OneBody(1,N,Id,Id)

    for k in 1:N-1 
        for p in k+1:N 

            J = rand()-0.5
            θ = rand()*pi
            ϕ = rand()*2*pi

            H += J*(1-3*cos(θ)^2)/4*(TwoBody(k,p,N,λ[1],λ[1],Id)+ TwoBody(k,p,N,λ[2],λ[2],Id)+ TwoBody(k,p,N,λ[6],λ[6],Id)+ TwoBody(k,p,N,λ[7],λ[7],Id)- TwoBody(k,p,N,λ[3]+sqrt(3)*λ[8],λ[3]+sqrt(3)*λ[8],Id))

            H += 3*J*sin(θ)^2/2*(exp(im*ϕ)*(TwoBody(k,p,N,λ[1]+im*λ[2],λ[6]+im*λ[7],Id) + TwoBody(k,p,N,λ[6]+im*λ[7],λ[1]+im*λ[2],Id)) + exp(-im*ϕ)*(TwoBody(k,p,N,λ[1]-im*λ[2],λ[6]-im*λ[7],Id) + TwoBody(k,p,N,λ[6]-im*λ[7],λ[1]-im*λ[2],Id)))

        end
    end

    return H 

end

function propagate_cyclic_ideal(N,H, Seq, Gens,τ,Ufinal, Nrepet) # Gens = generating pulses 

    U0 = exp(dense(-im*H*τ))
    U = exp(dense(-im*0.0*H))

    for n in 1:length(Seq)

        pulse = Gens[Seq[n][1]]
        for m in 2:N 
            pulse = tensor(pulse,Gens[Seq[n][1]])
        end

        for m in 2:length(Seq[n])
            
            newpulse = Gens[Seq[n][m]]
            for j in 2:N 
                newpulse = tensor(newpulse,Gens[Seq[n][m]])
            end
            pulse = newpulse*pulse
        end

        U = pulse*U0*U
    end
    
    for n in 1:length(Seq)

        pulse = Gens[Seq[end+1-n][1]]
        for m in 2:N 
            pulse = tensor(pulse,Gens[Seq[end+1-n][1]])
        end

        for m in 2:length(Seq[end+1-n])
            
            newpulse = Gens[Seq[end+1-n][m]]

            for j in 2:N 
                newpulse = tensor(newpulse,Gens[Seq[end+1-n][m]])
            end

            pulse = newpulse*pulse

        end

        U = U0*dagger(pulse)*U
    end

    U = U^Nrepet
    
    Utrg = copy(Ufinal)
    for m in 2:N 
        Utrg = tensor(Utrg,Ufinal)
    end

    return log10(sqrt(abs(1 - (1/3^N)*abs(tr(U*dagger(Utrg))))))

end

function propagate_cyclic_finite(N,H, Seq, Gens,τ0,τp, Id,Nrepet) # gens = generating Hamiltonians 

    U0 = exp(dense(-im*H*τ0))
    U = exp(dense(-im*0.0*H))

    for n in 1:length(Seq)

        pulse = exp(dense(-im*τp*(sum(OneBody(i,N,Gens[Seq[n][1]]/τp,Id) for i in 1:N) + H)))

        for m in 2:length(Seq[n])
            
            newpulse = exp(dense(-im*τp*(sum(OneBody(i,N,Gens[Seq[n][m]]/τp,Id) for i in 1:N) + H)))

            pulse = newpulse*pulse
        end

        U = pulse*U0*U

    end
    
    for n in 1:length(Seq)

        pulse = exp(dense(-im*τp*(-sum(OneBody(i,N,Gens[Seq[end+1-n][end]]/τp,Id) for i in 1:N) + H)))

        for m in 2:length(Seq[end+1-n])
            
            newpulse = exp(dense(-im*τp*(-sum(OneBody(i,N,Gens[Seq[end+1-n][end+1-m]]/τp,Id) for i in 1:N) + H)))
 
            pulse = newpulse*pulse
        end

        U = U0*pulse*U
    end

    U = U^Nrepet
    
    Utrg = copy(Id)
    for m in 2:N 
        Utrg = tensor(Utrg,Id)
    end
    
    return log10(sqrt(abs(1 - (1/3^N)*abs(tr(U*dagger(Utrg))))))

end

function propagate_pulseerrors(N,Gens_err, Seq, Gens, Ufinal) 

    U = exp(dense(-im*0.0*Gens_err[Seq[1][1]]))

    for n in 1:length(Seq)

        pulse = exp(dense(-im*(sum(OneBody(i,N,Gens[Seq[n][1]],Id) for i in 1:N) + norm(Gens[Seq[n][1]],Inf)*Gens_err[Seq[n][1]])))

        for m in 2:length(Seq[n])
            
            newpulse = exp(dense(-im*(sum(OneBody(i,N,Gens[Seq[n][m]],Id) for i in 1:N) + norm(Gens[Seq[n][m]],Inf)*Gens_err[Seq[n][m]])))

            pulse = newpulse*pulse

        end

        U = pulse*U

    end
    
    Utrg = copy(Ufinal)
    for m in 2:N 
        Utrg = tensor(Utrg,Ufinal)
    end
    
    return log10(sqrt(abs(1 - (1/3^N)*abs(tr(U*dagger(Utrg))))))

end
GLMakie.activate!()

##############
# generators #
##############

# All those matrices in Eq.(41)

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

#####################
# qutrit decoupling #
#####################

data_dir = "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/Adjacency Matrix SU3/"
fig_dir =  "/home/colin/Documents/JULIA/Codes persos/Irreps decomposition/decoupling SU3/" 

# the adjacency matrix are calculated in the file "G-invariant states.jl"
# the Eulerian path used in the paper are provided 

Δ27_mat = load(data_dir*"Adjacency matrix : Δ27 - C,E.jld", "data")
#=Δ27_sequence = Hierholzer_random(Δ27_mat)
Δ27_sequence = path_to_sequence(Δ27_mat,Δ27_sequence)
jldopen(data_dir*"path Δ27.jld", "w") do file 
    write(file, "data", Δ27_sequence)
end=#
Δ27_sequence = load(data_dir*"path Δ27.jld", "data")

Δ54_mat = load(data_dir*"Adjacency matrix : Δ54 - A(3),E,B.jld", "data")
#=Δ54_sequence = Hierholzer_random(Δ54_mat)
Δ54_sequence = path_to_sequence(Δ54_mat,Δ54_sequence)
jldopen(data_dir*"path Δ54.jld", "w") do file 
    write(file, "data", Δ54_sequence)
end=#
Δ54_sequence = load(data_dir*"path Δ54.jld", "data")

Σ36x3_mat = load(data_dir*"Adjacency matrix : Σ36x3 - V,C - modZ.jld", "data")
#Σ36_sequence = Hierholzer_random(Σ36x3_mat)
#=Σ36_sequence = path_to_sequence(Σ36x3_mat,Σ36_sequence)
jldopen(data_dir*"path Σ36.jld", "w") do file 
    write(file, "data", Σ36_sequence)
end=#
Σ36_sequence = load(data_dir*"path Σ36.jld", "data")

Σ72x3_mat = load(data_dir*"Adjacency matrix : Σ72x3 - V,X - modZ.jld", "data")
#=Σ72_sequence = Hierholzer_random(Σ72x3_mat)
Σ72_sequence = path_to_sequence(Σ72x3_mat,Σ72_sequence)
jldopen(data_dir*"path Σ72.jld", "w") do file 
    write(file, "data", Σ72_sequence)
end=#
Σ72_sequence = load(data_dir*"path Σ72.jld", "data")

#=
Σ168_mat = load(data_dir*"Adjacency matrix : Σ168 - Y,Z.jld", "data")
Σ168_sequence = Hierholzer_random(Σ168_mat)
Σ168_sequence = path_to_sequence(Σ168_mat,Σ168_sequence)
jldopen(data_dir*"path Σ168.jld", "w") do file 
    write(file, "data", Σ168_sequence)
end=#
Σ168_sequence = load(data_dir*"path Σ168.jld", "data")


color_C = RGB(255/255,102/255,102/255)
color_E = RGB(102/255,178/255,255/255)
color_V = RGB(178/255,105/255,255/255)
color_V = RGB(178/255,105/255,255/255)
color_X = RGB(255/255,178/255,102/255)
color_O1 = RGB(76/255,153/255,0/255)
color_O2 = RGB(153/255,0/255,0/255)

CairoMakie.activate!()

# Δ27 

fig = print_sequence(Δ27_sequence,[color_C,color_E], 18)
display(fig)
#save(fig_dir*"sequence Δ27.pdf", fig)

# Σ36

fig = print_sequence(Σ36_sequence,[color_V,color_C], 36)
display(fig)
#save(fig_dir*"sequence Σ36.pdf", fig)

# Σ72 

fig = print_sequence(Σ72_sequence,[color_V,color_X], 36)
display(fig)
#save(fig_dir*"sequence Σ72.pdf", fig)

# Klein-4 

K4 = [1 , 2 , 1 , 2]
fig = print_sequence(K4,[color_X,color_V], 4)
display(fig)
#save(fig_dir*"sequence K4.pdf", fig)

# Q8

K4 = [1 ,1 , 1, 2 , 1 , 1 , 1, 2]
fig = print_sequence(K4,[color_X,color_V], 8)
display(fig)
#save(fig_dir*"sequence Q8.pdf", fig)

# <<E>>

EE = [1 , 1 , 1]
fig = print_sequence(EE,[color_E], 3)
display(fig)
#save(fig_dir*"sequence E.pdf", fig)

# ZENITH

ZENITH = [5,2,5,3,5,1,5,3,5,3,5,6]
colors_Z = [RGB(86/255,180/255,233/255) , RGB(0,114/255,178/255) , RGB(230/255,159/255,0) , RGB(213/255,94/255,0),RGB(116/255,218/255,98/255),:black]
fig = print_sequence(ZENITH,colors_Z, 12)
display(fig)

# Δ27h

seq = [1,1,2,1,1,2,1,1,2]
fig = print_sequence(seq,[color_E,color_C], 9)
display(fig)
#save(fig_dir*"sequence Δ27h.pdf", fig)

# Δ27h

seq = [1,1,2,1,1,2,1,1,2,2,2,1,2,2,1,2,2,1]
fig = print_sequence(seq,[color_E,color_C], 9)
display(fig)
#save(fig_dir*"sequence Δ27e.pdf", fig)


# O 

O_seq = load(data_dir*"path O.jld", "data")
fig = print_sequence(O_seq,[color_O1,color_O2], 8)
display(fig)
#save(fig_dir*"sequence O.pdf", fig)

##################################################
# Generators of the 3 groups used in Section 4.3 #
##################################################

b = GenericBasis(3)

Δ27_gens = DenseOperator(b,C),DenseOperator(b,E)

Δ54_gens = DenseOperator(b,Ak(3)),DenseOperator(b,E),DenseOperator(b,B)

Σ36x3_gens = DenseOperator(b,V),DenseOperator(b,C)

Σ72x3_gens = DenseOperator(b,V),DenseOperator(b,X)

Σ168_gens = DenseOperator(b,Y),DenseOperator(b,Z)

Δ27_gens_h = DenseOperator(b,-2*pi/3*[0 0 0 ; 0 1 0 ; 0 0 -1]) , DenseOperator(b, -2*pi/3*sqrt(1/3)*[0 -im im; im 0 -im ; -im im 0])

Σ36x3_gens_h = DenseOperator(b,-pi/2*sqrt(3)/3*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]) , DenseOperator(b,-2*pi/3*[0 0 0 ; 0 1 0 ; 0 0 -1])

U = eigvecs(X)
Σ72x3_gens_h = DenseOperator(b,-pi/2*sqrt(3)/3*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]) , DenseOperator(b, U*[pi/2 0 0 ; 0 -pi/2 0; 0 0 0]*adjoint(U)) 

##################################
# Ideal pulses --> unused figure #
##################################

λ = GellMann()

N = 3 
τ = 1.0
nsample = 1000

dis = LinRange(-3,0,50)
inter = LinRange(-3,0,50)

data_Δ27 = zeros(length(dis),length(inter))
data_Σ36x3 = zeros(length(dis),length(inter))
data_Σ72x3 = zeros(length(dis),length(inter))
data_nodd = zeros(length(dis),length(inter))

Δ27_sequence_extended = vcat(Δ27_sequence,Δ27_sequence)
Δ27_sequence_extended = vcat(Δ27_sequence_extended,Δ27_sequence_extended)
Δ27_sequence_extended = vcat(Δ27_sequence_extended,Δ27_sequence_extended)
Σ36_sequence_extended = vcat(Σ36_sequence,Σ36_sequence)


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder(N,λ,one(b))
    H_dip = Hamiltonian_interaction(N,λ,one(b))

    for i in 1:length(dis)

        Threads.@threads for j in 1:length(inter)

            H = exp10(dis[i])*H_dis + exp10(inter[j])*H_dip 
            
            data_Δ27[i,j]+=propagate_ideal(N,H,Δ27_sequence_extended,Δ27_gens, τ)/nsample

            data_Σ36x3[i,j]+=propagate_ideal(N,H,Σ36_sequence_extended,Σ36x3_gens, τ)/nsample
            data_Σ72x3[i,j]+=propagate_ideal(N,H,Σ72_sequence,Σ72x3_gens, τ)/nsample

            data_nodd[i,j]+=propagate_nodd(N,H,length(Σ72_sequence),τ)/nsample

        end

        println("n=$n i=$i")
    end
end
=#


#=
jldopen(fig_dir*"/disorder and interaction : Δ27.jld", "w") do file 
    write(file, "data", data_Δ27)
end

jldopen(fig_dir*"/disorder and interaction : Σ36x3.jld", "w") do file 
    write(file, "data", data_Σ36x3)
end
jldopen(fig_dir*"/disorder and interaction : Σ72x3.jld", "w") do file 
    write(file, "data", data_Σ72x3)
end

jldopen(fig_dir*"/disorder and interaction : nodd.jld", "w") do file 
    write(file, "data", data_nodd)
end
=#

data_Δ27 = load(fig_dir*"/disorder and interaction : Δ27.jld", "data")
data_Σ36x3 = load(fig_dir*"/disorder and interaction : Σ36x3.jld", "data")
data_Σ72x3 = load(fig_dir*"/disorder and interaction : Σ72x3.jld", "data")
data_nodd = load(fig_dir*"/disorder and interaction : nodd.jld", "data")

colors = [RGB(204/255,153/255,255/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(153/255 , 204/255 ,255/255) ,
RGB(192/255,192/255,192/255)]

datas = [data_Δ27 , data_Σ36x3 , data_Σ72x3 , data_nodd]
lab = [L"\Delta(27)/\mathbb{Z}_3"  , L"\Sigma(36\times 3)/\mathbb{Z}_3" , L"\Sigma(72\times 3)/\mathbb{Z}_3", L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau\Delta", ylabel = L"\tau \Gamma" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1,0], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"1"]),yticks = ([-3,-2,-1,0], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"1"]),
zticks = ([-4 ,-2 , 0], [L"10^{-4}", L"10^{-2}", L"1"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50)

for n in 1:4
    cmat = Matrix{Colorant}(undef, 50,50)
    cmat[:,:].= colors[n]
    surface!(ax, dis, inter, datas[n], color=cmat, alpha=0.6)
end

midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in 1:4

        lines!(ax, [dis[midlines[n]]], inter, datas[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, dis, [inter[midlines[n]]], datas[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end

lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in 1:4]
Legend(fig[2,1], lab_mark, [lab[k] for k in 1:4],framevisible=false, 
nbanks=2, tellwidth=false, tellheight=true,
labelsize=27.0f0)
display(fig)

#save(fig_dir*"/disorder and interaction : Σ72x3 - Σ36x3 - Δ27.png", fig ; px_per_unit = 15.0)


###########################
# Finite pulses --> Fig 5 #
###########################

λ = GellMann()

Id = one(b)

N = 3 

nsample = 500

τ0 = LinRange(-3,0,50)
τp = LinRange(-3,0,50)

data_Δ27_dis = zeros(length(τ0),length(τp))
data_Σ36x3_dis = zeros(length(τ0),length(τp))
data_Σ72x3_dis = zeros(length(τ0),length(τp))
data_nodd_dis = zeros(length(τ0),length(τp))

data_Δ27_dip = zeros(length(τ0),length(τp))
data_Σ36x3_dip = zeros(length(τ0),length(τp))
data_Σ72x3_dip = zeros(length(τ0),length(τp))
data_nodd_dip = zeros(length(τ0),length(τp))

Δ27_sequence_extended = vcat(Δ27_sequence,Δ27_sequence)
Δ27_sequence_extended = vcat(Δ27_sequence_extended,Δ27_sequence_extended)
Δ27_sequence_extended = vcat(Δ27_sequence_extended,Δ27_sequence_extended)
Σ36_sequence_extended = vcat(Σ36_sequence,Σ36_sequence)

#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder(N,λ,one(b))
    H_dis = H_dis/norm(H_dis)
    H_dip = Hamiltonian_interaction(N,λ,one(b))
    H_dip = H_dip/norm(H_dip)

    for i in 1:length(τ0)

        Threads.@threads for j in 1:length(τp)

            t0 = exp10(τ0[i])
            tp = exp10(τp[j])

            data_Δ27_dis[i,j]+=propagate_finite(N,H_dis,Δ27_sequence_extended,Δ27_gens_h,t0,tp,Id )/nsample
            data_Σ36x3_dis[i,j]+=propagate_finite(N,H_dis,Σ36_sequence_extended,Σ36x3_gens_h,t0,tp,Id)/nsample
            data_Σ72x3_dis[i,j]+=propagate_finite(N,H_dis,Σ72_sequence,Σ72x3_gens_h,t0,tp,Id)/nsample
            
            data_nodd_dis[i,j]+=propagate_nodd(N,H_dis,144,t0+tp)/nsample

            data_Δ27_dip[i,j]+=propagate_finite(N,H_dip,Δ27_sequence_extended,Δ27_gens_h,t0,tp,Id )/nsample
            data_Σ36x3_dip[i,j]+=propagate_finite(N,H_dip,Σ36_sequence_extended,Σ36x3_gens_h,t0,tp,Id)/nsample
            data_Σ72x3_dip[i,j]+=propagate_finite(N,H_dip,Σ72_sequence,Σ72x3_gens_h,t0,tp,Id)/nsample
            
            data_nodd_dip[i,j]+=propagate_nodd(N,H_dip,144,t0+tp)/nsample

        end

        println("n=$n i=$i")
    end
end
=#


#=
jldopen(fig_dir*"/disorder and interaction finite dis: Δ27.jld", "w") do file 
    write(file, "data", data_Δ27_dis)
end
jldopen(fig_dir*"/disorder and interaction finite dis: Σ36x3.jld", "w") do file 
    write(file, "data", data_Σ36x3_dis)
end
jldopen(fig_dir*"/disorder and interaction finite dis: Σ72x3.jld", "w") do file 
    write(file, "data", data_Σ72x3_dis)
end
jldopen(fig_dir*"/disorder and interaction finite dis: nodd.jld", "w") do file 
    write(file, "data", data_nodd_dis)
end

jldopen(fig_dir*"/disorder and interaction finite dip: Δ27.jld", "w") do file 
    write(file, "data", data_Δ27_dip)
end
jldopen(fig_dir*"/disorder and interaction finite dip: Σ36x3.jld", "w") do file 
    write(file, "data", data_Σ36x3_dip)
end
jldopen(fig_dir*"/disorder and interaction finite dip: Σ72x3.jld", "w") do file 
    write(file, "data", data_Σ72x3_dip)
end
jldopen(fig_dir*"/disorder and interaction finite dip: nodd.jld", "w") do file 
    write(file, "data", data_nodd_dip)
end
=#

data_Δ27_dis = load(fig_dir*"/disorder and interaction finite dis: Δ27.jld", "data")
data_Σ36x3_dis = load(fig_dir*"/disorder and interaction finite dis: Σ36x3.jld", "data")
data_Σ72x3_dis = load(fig_dir*"/disorder and interaction finite dis: Σ72x3.jld", "data")
data_nodd_dis = load(fig_dir*"/disorder and interaction finite dis: nodd.jld", "data")

data_Δ27_dip = load(fig_dir*"/disorder and interaction finite dip: Δ27.jld", "data")
data_Σ36x3_dip = load(fig_dir*"/disorder and interaction finite dip: Σ36x3.jld", "data")
data_Σ72x3_dip = load(fig_dir*"/disorder and interaction finite dip: Σ72x3.jld", "data")
data_nodd_dip = load(fig_dir*"/disorder and interaction finite dip: nodd.jld", "data")

colors = [RGB(204/255,153/255,255/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(153/255 , 204/255 ,255/255) , 
RGB(192/255,192/255,192/255)]

datadis = [data_Δ27_dis , data_Σ36x3_dis , data_Σ72x3_dis , data_nodd_dis]
datadip = [data_Δ27_dip , data_Σ36x3_dip , data_Σ72x3_dip , data_nodd_dip]

lab = [L"\Delta(27)/\mathbb{Z}_3" , L"\Sigma(36\times 3)/\mathbb{Z}_3" , L"\Sigma(72\times 3)/\mathbb{Z}_3" , L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau_0 \Delta", ylabel = L"\tau_p \Delta" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-6,-4,-2 , 0], [L"10^{-6}",L"10^{-4}",L"10^{-2}", L"1"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50, title=L"\Gamma = 0",titlesize=30.0f0)


for n in 1:4
    cmat = Matrix{Colorant}(undef, 50,50)
    cmat[:,:].= colors[n]
    surface!(ax, τ0, τp, datadis[n], color=cmat, alpha=0.8)
end

midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in 1:4

        lines!(ax, [τ0[midlines[n]]], τp, datadis[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, τ0, [τp[midlines[n]]], datadis[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end

lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in 1:4]
Legend(fig[2,1], lab_mark, [lab[k] for k in 1:4],framevisible=false, 
nbanks=2, tellwidth=false, tellheight=true,
labelsize=27.0f0)
display(fig)

#save(fig_dir*"/disorder and interaction : Σ72x3 - Σ36x3 - Δ27 finite dis.png", fig ; px_per_unit = 15.0)


fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau_0 \Gamma", ylabel = L"\tau_p \Gamma" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-4,-2 , 0], [L"10^{-4}",L"10^{-2}", L"1"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50, title=L"\Delta = 0",titlesize=30.0f0)


for n in 1:4
    cmat = Matrix{Colorant}(undef, 50,50)
    cmat[:,:].= colors[n]
    surface!(ax, τ0, τp, datadip[n], color=cmat, alpha=0.7)
end

midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in 1:4

        lines!(ax, [τ0[midlines[n]]], τp, datadip[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, τ0, [τp[midlines[n]]], datadip[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end

lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in 1:4]
Legend(fig[2,1], lab_mark, [lab[k] for k in 1:4],framevisible=false, 
nbanks=2, tellwidth=false, tellheight=true,
labelsize=27.0f0)
display(fig)

#save(fig_dir*"/disorder and interaction : Σ72x3 - Σ36x3 - Δ27 finite dip.png", fig ; px_per_unit = 15.0)

##########################
# pulse errors --> Fig 5 #
##########################

Id = one(b)
λ = GellMann()

N = 1

nsample = 1000

logϵ = LinRange(-3,-1,500)
ϵ = exp10.(logϵ)

Ufinal = Id

data_Δ27 = zeros(length(logϵ))
data_Σ36x3 = zeros(length(logϵ))
data_Σ72x3 = zeros(length(logϵ))

data_Δ272 = zeros(length(logϵ))
data_Σ36x32 = zeros(length(logϵ))
data_Σ72x32 = zeros(length(logϵ))

colors = [RGB(204/255,153/255,255/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(153/255 , 204/255 ,255/255) ]

#=
for n in 1:nsample

    error1 = sum((rand()-0.5)*λ[j] for j in 1:8) # a very random error
    error1= error1/norm(error1,Inf)

    error2 = sum((rand()-0.5)*λ[j] for j in 1:8) # a very random error
    error2= error2/norm(error2,Inf)

    Threads.@threads for e in 1:length(logϵ)

        data_Δ27[e] += propagate_errors([ϵ[e]*error1/norm(Δ27_gens_h[1],Inf),ϵ[e]*error2/norm(Δ27_gens_h[2],Inf)],Δ27_sequence_extended,Δ27_gens_h,Ufinal)/nsample
        data_Σ36x3[e] += propagate_errors([ϵ[e]*error1/norm(Σ36x3_gens_h[1],Inf),ϵ[e]*error2/norm(Σ36x3_gens_h[2],Inf)],Σ36_sequence_extended,Σ36x3_gens_h,Ufinal)/nsample
        data_Σ72x3[e] += propagate_errors([ϵ[e]*error1/norm(Σ72x3_gens_h[1],Inf),ϵ[e]*error2/norm(Σ72x3_gens_h[2],Inf)],Σ72_sequence,Σ72x3_gens_h,Ufinal)/nsample

        data_Δ272[e] += propagate_errors([ϵ[e]*error1,ϵ[e]*error2],Δ27_sequence_extended,Δ27_gens_h,Ufinal)/nsample
        data_Σ36x32[e] += propagate_errors([ϵ[e]*error1,ϵ[e]*error2],Σ36_sequence_extended,Σ36x3_gens_h,Ufinal)/nsample
        data_Σ72x32[e] += propagate_errors([ϵ[e]*error1,ϵ[e]*error2],Σ72_sequence,Σ72x3_gens_h,Ufinal)/nsample

    end

end
=# 

GLMakie.activate!()
fig = Figure(;size=(370,370))

ax = Axis(fig[1,1], ylabel=L"1-\overline{D}", xlabel=L"\epsilon",
xlabelsize=30.0f0, ylabelsize=30.0f0,
xticks=([0,0.05,0.1],[L"0.0",L"0.05",L"0.1"]),
yticks=([0.8,0.85,0.9,0.95,1.0],[L"0.8",L"0.85",L"0.9",L"0.95",L"1"]),
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(data_Δ27), color=colors[1], linewidth=4.0)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(data_Σ36x3), color=colors[2], linewidth=4.0)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(data_Σ72x3), color=colors[3], linewidth=4.0)

elems1 = [LineElement(color=:black, linewidth=3.0,linestyle=:solid),LineElement(color=:black, linewidth=3.0,linestyle=(:dash,:dense))]
elems2 = [PolyElement(color=colors[i], strokewidth=2.0, strokecolor=:black) for i in 1:3]
leg1 = [L"\text{not scaled}" , L"\text{scaled}"]
leg2 = [L"\Delta(27)/\mathbb{Z}_3" , L"\Sigma(36\times 3)/\mathbb{Z}_3" , L"\Sigma(72\times 3)/\mathbb{Z}_3"]

axislegend(ax,[ elems2],[leg2],[""],tellwidth=false,tellheight=false , nbanks=1,position=:lb,labelsize=20.0f0,framevisible=false)

ylims!(0.8,nothing)
ax_inset = Axis(fig[1, 1],
width=Relative(0.35),
height=Relative(0.3),
halign=0.90,
valign=0.97,xticks=([-3,-1],[L"10^{-3}",L"10^{-1}"]),
yticks = ([-5,-3,-1],[L"10^{-5}",L"10^{-3}",L"10^{-1}"]),
xaxisposition=:bottom,yaxisposition=:left,xticklabelsize=18.0f0,yticklabelsize=18.0f0,
xlabel=L"\epsilon", ylabel=L"\bar{D}",ylabelrotation=0,
ylabelsize=20.0f0,xlabelsize=20.0f0,xlabelpadding=-20)

lines!(ax_inset,logϵ,data_Δ27, color=colors[1], linewidth=3.0)
lines!(ax_inset,logϵ,data_Σ36x3, color=colors[2], linewidth=3.0)
lines!(ax_inset,logϵ,data_Σ72x3, color=colors[3], linewidth=3.0)

display(fig)

#save(fig_dir*"/disorder and interaction : Σ72x3 - Σ36x3 - Δ27 errors.png", fig ; px_per_unit = 20.0)


##############################
# NV - centers : Section 4.4 #
##############################

W0 = collect(Diagonal([1 , 1 , 1]))

W1 = zeros(3,3)
W1[1,1] = -1
W1[2,3] = 1
W1[3,2] = 1
det(W1)

W2 = zeros(3,3)
W2[2,2] = -1
W2[1,3] = 1
W2[3,1] = 1

W3 = zeros(3,3)
W3[3,3] = -1
W3[1,2] = 1
W3[2,1] = 1

W4 = zeros(3,3)
W4[1,2] = 1
W4[2,3] = 1
W4[3,1] = 1

W5 = zeros(3,3)
W5[1,3] = 1
W5[3,2] = 1
W5[2,1] = 1

##########################################
# Σ(72 x 3) = Σ(36 x 3) {1,X} --> Fig 10 #
##########################################

GLMakie.activate!()

U = eigvecs(X)
round.(adjoint(U)*X*U,digits=13) # ok 

UVU = adjoint(U)*V*U
UCU = adjoint(U)*C*U

Σ36x3_gens = DenseOperator(b,V),DenseOperator(b,C)

Σ36x3_gens_U_W0 = DenseOperator(b,UVU),DenseOperator(b,UCU)
Σ36x3_gens_U_W1 = DenseOperator(b,adjoint(W1)*UVU*W1),DenseOperator(b,adjoint(W1)*UCU*W1)
Σ36x3_gens_U_W2 = DenseOperator(b,adjoint(W2)*UVU*W2),DenseOperator(b,adjoint(W2)*UCU*W2)
Σ36x3_gens_U_W3 = DenseOperator(b,adjoint(W3)*UVU*W3),DenseOperator(b,adjoint(W3)*UCU*W3)
Σ36x3_gens_U_W4 = DenseOperator(b,adjoint(W4)*UVU*W4),DenseOperator(b,adjoint(W4)*UCU*W4)
Σ36x3_gens_U_W5 = DenseOperator(b,adjoint(W5)*UVU*W5),DenseOperator(b,adjoint(W5)*UCU*W5)

Σ72x3_gens = DenseOperator(b,V),DenseOperator(b,X)

λ = GellMann()

N = 3 
τ = 1.0
nsample = 500

dis = LinRange(-3,0,50)
inter = LinRange(-3,0,50)

data_Δ27 = zeros(length(dis),length(inter))

data_Σ36x3 = zeros(length(dis),length(inter))

data_Σ36x3_U1 = zeros(length(dis),length(inter))
data_Σ36x3_U2 = zeros(length(dis),length(inter))
data_Σ36x3_U3 = zeros(length(dis),length(inter))
data_Σ36x3_U4 = zeros(length(dis),length(inter))
data_Σ36x3_U5 = zeros(length(dis),length(inter))
data_Σ36x3_U6 = zeros(length(dis),length(inter))

data_Σ72x3 = zeros(length(dis),length(inter))

data_nodd = zeros(length(dis),length(inter))
#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis,Inf)
    H_dip = Hamiltonian_interaction_NV(N,λ,one(b))
    H_dip = H_dip/norm(H_dip,Inf)

    for i in 1:length(dis)

        Threads.@threads for j in 1:length(inter)

            H = exp10(dis[i])*H_dis + exp10(inter[j])*H_dip 
            
            #data_Δ27[i,j]+=propagate_ideal(N,H,Δ27_seq,Δ27_gens, τ)/nsample

            #data_Σ36x3[i,j]+=propagate_ideal(N,H,Σ36x3_seq,Σ36x3_gens, τ)/nsample

            data_Σ36x3_U1[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W0, τ)/nsample
            data_Σ36x3_U2[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W1, τ)/nsample
            data_Σ36x3_U3[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W2, τ)/nsample
            data_Σ36x3_U4[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W3, τ)/nsample
            data_Σ36x3_U5[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W4, τ)/nsample
            data_Σ36x3_U6[i,j]+=propagate_ideal(N,H,Σ36_sequence,Σ36x3_gens_U_W5, τ)/nsample

            #data_Σ72x3[i,j]+=propagate_ideal(N,H,Σ72x3_seq,Σ72x3_gens, τ)/nsample

            
            data_nodd[i,j]+=propagate_nodd(N,H,length(Σ36_sequence),τ)/nsample
        end

        println("n=$n i=$i")
    end
end

jldopen(fig_dir*"/disorder and interaction NV : Δ27.jld", "w") do file 
    write(file, "data", data_Δ27)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3.jld", "w") do file 
    write(file, "data", data_Σ36x3)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U1.jld", "w") do file 
    write(file, "data", data_Σ36x3_U1)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U2.jld", "w") do file 
    write(file, "data", data_Σ36x3_U2)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U3.jld", "w") do file 
    write(file, "data", data_Σ36x3_U3)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U4.jld", "w") do file 
    write(file, "data", data_Σ36x3_U4)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U5.jld", "w") do file 
    write(file, "data", data_Σ36x3_U5)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ36x3 U6.jld", "w") do file 
    write(file, "data", data_Σ36x3_U6)
end

jldopen(fig_dir*"/disorder and interaction NV : Σ72x3.jld", "w") do file 
    write(file, "data", data_Σ72x3)
end

jldopen(fig_dir*"/disorder and interaction NV : nodd 36x2.jld", "w") do file 
    write(file, "data", data_nodd)
end
=#

colors = [RGB(178/255,102/255,255/255) , 
RGB(255/255,204/255,153/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(153/255 , 153/255 ,255/255) , 
RGB(153/255 , 255/255 ,153/255) ,
RGB(255/255,255/255,153/255) ,
RGB(255/255,102/255,178/255) ,
RGB(192/255,192/255,192/255)]


data_Σ36x3_U1 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U1.jld", "data")

data_Σ36x3_U2 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U2.jld", "data")

data_Σ36x3_U3 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U3.jld", "data")

data_Σ36x3_U4 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U4.jld", "data")

data_Σ36x3_U5 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U5.jld", "data")

data_Σ36x3_U6 = load(fig_dir*"/disorder and interaction NV : Σ36x3 U6.jld", "data")

data_nodd = load(fig_dir*"/disorder and interaction NV : nodd 36x2.jld", "data")

datas = [data_Σ36x3_U1 , data_Σ36x3_U2 , data_Σ36x3_U3 , data_Σ36x3_U4 , data_Σ36x3_U5 , data_Σ36x3_U6 , data_Σ72x3, data_nodd]
lab = [L"U" , L"U(23),\, U(132)" ,L"U,\,U(13)" ,L"U(12),\,U(123)",L"U(123)",L"U(132)", L"\Sigma(72\times 3)/\mathbb{Z}_3", L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau_0\Delta", ylabel = L"\tau_0 J" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-5 ,-3 , -1], [L"10^{-5}", L"10^{-3}", L"10^{-1}"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50)

Ns = [2 , 3 , 4 , 8]


for n in Ns
    cmat = Matrix{Colorant}(undef, 50,50)
    cmat[:,:].= colors[n]
    surface!(ax, dis, inter, datas[n], color=cmat, alpha=0.6)
end
midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in Ns

        lines!(ax, [dis[midlines[n]]], inter, datas[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, dis, [inter[midlines[n]]], datas[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end
lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in Ns]
Legend(fig[2,1], lab_mark, [lab[k] for k in Ns],framevisible=false, 
nbanks=4, tellwidth=false, tellheight=true,
labelsize=20.0f0)
display(fig)

#save(fig_dir*"/disorder and interaction NV : UΣ36x3U.png", fig;px_per_unit=15)


###############################
# nested sequences --> Fig 12 #
###############################

gens = DenseOperator(b, E) , DenseOperator(b, X) , DenseOperator(b, V)

Seq1 = [[1] , [1] , [1,2], [1] , [1] , [1,3],[1] , [1] , [1,2],[1] , [1] , [1,3]]

Seq2 = [[1] , [1] , [1,2],[1] , [1] , [1,2],[1] , [1] , [1,2],[1] , [1] , [1,3],[1] , [1] , [1,2],[1] , [1] , [1,2],[1] , [1] , [1,2],[1] , [1] , [1,3]]

Ufinal_seq1 = DenseOperator(b, V*X*V*X)
Ufinal_seq2 = DenseOperator(b, V*X^3*V*X^3)

λ = GellMann()

N = 3 
τ = 1.0
nsample = 100

dis = LinRange(-3,0,50)
inter = LinRange(-3,0,50)

data_nodd = zeros(length(dis),length(inter))
data_seq1 = zeros(length(dis),length(inter))
data_seq2 = zeros(length(dis),length(inter))


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis)
    H_dip = Hamiltonian_interaction_NV(N,λ,one(b))
    H_dip = H_dip/norm(H_dip)

    for i in 1:length(dis)

        for j in 1:length(inter)
            
            H = exp10(dis[i])*H_dis + exp10(inter[j])*H_dip 

            data_seq1[i,j]+=propagate_nested(N,H,Seq1,gens, τ,Ufinal_seq1)/nsample

            data_seq2[i,j]+=propagate_nested(N,H,Seq2,gens, τ,Ufinal_seq2)/nsample
            
            data_nodd[i,j]+=propagate_nodd(N,H,length(Seq1),τ)/nsample
        end

        println("n=$n i=$i")
    end
end

jldopen(fig_dir*"/disorder and interaction NV : nested Klein 4-group.jld", "w") do file 
    write(file, "data", data_seq1)
end


jldopen(fig_dir*"/disorder and interaction NV : nested Quaternion.jld", "w") do file 
    write(file, "data", data_seq2)
end

jldopen(fig_dir*"/disorder and interaction NV : nodd 4τ.jld", "w") do file 
    write(file, "data", data_nodd)
end
=#


data_nodd = load(fig_dir*"/disorder and interaction NV : nodd 4τ.jld", "data")
data_seq1 = load(fig_dir*"/disorder and interaction NV : nested Klein 4-group.jld", "data")
data_seq2 = load(fig_dir*"/disorder and interaction NV : nested Quaternion.jld", "data")

GLMakie.activate!()

colors = [RGB(204/255,153/255,255/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(192/255,192/255,192/255)]

datas = [data_seq1 , data_seq2 , data_nodd]
lab = [L"K_4[E\,]" , L"Q_8[E\,]" , L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau_0\Delta", ylabel = L"\tau_0 J" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-6,-4 ,-2 , 0], [L"10^{-6}",L"10^{-4}", L"10^{-2}", L"1"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50)

Ns = [1 ,2 , 3]
for n in Ns
    cmat = Matrix{Colorant}(undef, 50,50)
    cmat[:,:].= colors[n]
    surface!(ax, dis, inter, datas[n], color=cmat, alpha=0.6)
end

midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in Ns

        lines!(ax, [dis[midlines[n]]], inter, datas[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, dis, [inter[midlines[n]]], datas[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end

lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in Ns]
Legend(fig[2,1], lab_mark, [lab[k] for k in Ns],framevisible=false, 
nbanks=3, tellwidth=false, tellheight=true,
labelsize=27.0f0)

display(fig)


#save(fig_dir*"/disorder and interaction NV : nested sequences.png", fig;px_per_unit=15.0)


#########################################
# nested sequences ZF --> section 4.4.4 #
#########################################

# Fig 14.A 

U = [(1+sqrt(3))/sqrt(6+2*sqrt(3)) (1-sqrt(3))/sqrt(6-2*sqrt(3)) 0 ; 
1/sqrt(6+2*sqrt(3)) 1/sqrt(6-2*sqrt(3)) -exp(4*pi*im/3)/sqrt(2) ;
exp(2*pi*im/3)/sqrt(6+2*sqrt(3)) exp(2*pi*im/3)/sqrt(6-2*sqrt(3)) 1/sqrt(2)]

gens = DenseOperator(b, E) , DenseOperator(b, X) ,
 DenseOperator(b, V) , DenseOperator(b, C) 

gens_ham = DenseOperator(b, -2*pi/3*sqrt(1/3)*[0 -im im; im 0 -im ; -im im 0]) , DenseOperator(b, U*[pi/2 0 0 ; 0 -pi/2 0; 0 0 0]*adjoint(U)) , 
DenseOperator(b,-pi/2*sqrt(3)/3*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]) , DenseOperator(b,-2*pi/3*[0 0 0 ; 0 1 0 ; 0 0 -1])
 
Sx = [0 1 0;
1 0 1;
0 1 0]/sqrt(2)

Sy = [0 -im 0;
im 0 -im;
0 im 0]/sqrt(2)

Sz = [1 0 0;
0 0 0;
0 0 -1]

gens_zenith = DenseOperator(b, exp(-im*pi/2*Sx)) , DenseOperator(b, exp(-im*pi/2*Sy)) ,
 DenseOperator(b, exp(-im*pi/2*(Sx*Sz+Sz*Sx))) , DenseOperator(b, exp(-im*pi/2*(Sy*Sz+Sz*Sy))) , DenseOperator(b, exp(-im*pi/2*Sx)*exp(-im*pi/2*Sy)*exp(im*pi/2*Sx)), DenseOperator(b,exp(-im*0.0*Sz))

gens_zenith_ham = DenseOperator(b, pi/2*Sx) , DenseOperator(b, pi/2*Sy) ,
 DenseOperator(b,pi/2*(Sx*Sz+Sz*Sx)) , DenseOperator(b, pi/2*(Sy*Sz+Sz*Sy)) , DenseOperator(b, pi/2*Sz), DenseOperator(b,0.0*Sz)

K4_E = [[1] , [1] , [1,2], [1] , [1] , [1,3],[1] , [1] , [1,2],[1] , [1] , [1,3]]

K4_Δ27_hamiltonian = [[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4,2], 
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4,3] , 
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4,2] ,
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4,3]]

K4_Δ27_eulerian = [[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4], [4] , [4] , [1] , [4] ,[4] ,[1] ,[4] , [4] ,[1,2], 
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4], [4] , [4] , [1] , [4] ,[4] ,[1] ,[4] , [4] ,[1,3] , 
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4], [4] , [4] , [1] , [4] ,[4] ,[1] ,[4] , [4] ,[1,2] ,
[1] , [1] , [4] , [1] ,[1] ,[4] ,[1] , [1] ,[4], [4] , [4] , [1] , [4] ,[4] ,[1] ,[4] , [4] ,[1,3]]

ZENITH = [[5],[2],[5],[3],[5],[1],[5],[3],[5],[3],[5],[6]]

Id = one(b)
λ = GellMann()

N = 3
τ = 1.0
nsample = 100

dis = LinRange(-3,-0.0,50)
inter = LinRange(-3,-0.0,50)

data_nodd = zeros(length(dis),length(inter))
data_seq1 = zeros(length(dis),length(inter))
data_seq2 = zeros(length(dis),length(inter))
data_seq3 = zeros(length(dis),length(inter))
data_seq4 = zeros(length(dis),length(inter))


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dip = Hamiltonian_interaction_NV_ZF(N,λ,one(b))

    for i in 1:length(dis)

        Threads.@threads for j in 1:length(inter)
            
            H = exp10(dis[i])*H_dis/norm(H_dis,Inf) + exp10(inter[j])*H_dip/norm(H_dip,Inf)

            data_seq1[i,j]+=propagate_cyclic_ideal(N,H,K4_E,gens, τ,Id,6)/nsample

            data_seq2[i,j]+=propagate_cyclic_ideal(N,H,K4_Δ27_hamiltonian,gens, τ,Id,2)/nsample

            data_seq3[i,j]+=propagate_cyclic_ideal(N,H,K4_Δ27_eulerian,gens, τ,Id,1)/nsample

            data_seq4[i,j]+=propagate_cyclic_ideal(N,H,ZENITH,gens_zenith,τ,Id,6)/nsample
            
            data_nodd[i,j]+=propagate_nodd(N,H,length(K4_Δ27_eulerian)*2,τ)/nsample

        end

        println("n=$n i=$i")
    end
end

jldopen(fig_dir*"/disorder and interaction NV ZF : K4_E.jld", "w") do file 
    write(file, "data", data_seq1)
end

jldopen(fig_dir*"/disorder and interaction NV ZF : K4_Δham.jld", "w") do file 
    write(file, "data", data_seq2)
end

jldopen(fig_dir*"/disorder and interaction NV ZF : K4_Δeul.jld", "w") do file 
    write(file, "data", data_seq3)
end

jldopen(fig_dir*"/disorder and interaction NV ZF : Zenith.jld", "w") do file 
    write(file, "data", data_seq4)
end

jldopen(fig_dir*"/disorder and interaction NV ZF : nodd.jld", "w") do file 
    write(file, "data", data_nodd)
end
=#


data_nodd = load(fig_dir*"/disorder and interaction NV ZF : nodd.jld", "data")
data_seq1 = load(fig_dir*"/disorder and interaction NV ZF : K4_E.jld", "data")
data_seq2 = load(fig_dir*"/disorder and interaction NV ZF : K4_Δham.jld", "data")
data_seq3 = load(fig_dir*"/disorder and interaction NV ZF : K4_Δeul.jld", "data")
data_seq4 = load(fig_dir*"/disorder and interaction NV ZF : Zenith.jld", "data")

GLMakie.activate!()

colors = [RGB(204/255,153/255,255/255) , 
RGB(255/255 , 102/255 , 102/255) , RGB(153/255 , 204/255 , 255/255) ,RGB(255/255 , 204/255 , 153/255) ,
RGB(192/255,192/255,192/255)]

datas = [data_seq1 , data_seq2 ,data_seq3 , data_seq4 , data_nodd]
lab = [L"K_4[\langle E\,\rangle]" , L"K_4[\Delta(27)h\,]",L"K_4[\Delta(27)e\,]",L"\mathrm{ZENITH}" , L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau\Delta", ylabel = L"\tau J" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-6,-4,-2 ,0], [L"10^{-6}",L"10^{-4}", L"10^{-2}", L"1"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50)
alphas = [0.6,0.7,0.6,0.8,0.6]
for n in 1:5
    cmat = Matrix{Any}(undef, 50,50)
    for i in 1:50 
        for j in 1:50 
            if datas[n][i,j] > -6
                cmat[i,j] = colors[n]
            else 
                cmat[i,j] = :transparent
            end
        end
    end
    surface!(ax, dis, inter, datas[n], color=cmat, alpha=0.6,transparency=false)
end

midlines = [1 , 10 , 20 , 30 , 40 , 50]

for n in 1:length(midlines)
    for m in 1:5

        c1 = Vector{Any}(undef,50)
        for i in 1:50 
            if datas[m][midlines[n],i] > -6 
                c1[i] = colors[m] 
            else 
                c1[i] = :transparent
            end
        end
        c2 = Vector{Any}(undef,50)
        for i in 1:50 
            if datas[m][i,midlines[n]] > -6 
                c2[i] = colors[m] 
            else 
                c2[i] = :transparent
            end
        end
        lines!(ax, [dis[midlines[n]]], inter, datas[m][midlines[n],:], color=c1, linewidth=2.0)
        lines!(ax, dis, [inter[midlines[n]]], datas[m][:,midlines[n]], color=c2, linewidth=2.0)

    end
end

lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in 1:5]
Legend(fig[2,1], lab_mark, [lab[k] for k in 1:5],framevisible=false, 
nbanks=3, tellwidth=false, tellheight=true,
labelsize=27.0f0)

display(fig)

#save(fig_dir*"/disorder and interaction NV ZF : nested sequences.png", fig;px_per_unit=15.0)

# finite pulses -> Fig 15 D

Id = one(b)
λ = GellMann()

N = 3

nsample = 30

τ0_log = LinRange(-3,-1,50)
τp_log = LinRange(-3,-1,50)

τ0 = exp10.(τ0_log)
τp = exp10.(τp_log)

data_seq2 = zeros(length(τ0),length(disinter),length(τp))
data_seq3 = zeros(length(τ0),length(disinter),length(τp))
data_seq4 = zeros(length(τ0),length(disinter),length(τp))

data_seq2_bis = zeros(length(τ0),length(disinter),length(τp))
data_seq3_bis = zeros(length(τ0),length(disinter),length(τp))
data_seq4_bis = zeros(length(τ0),length(disinter),length(τp))


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis,Inf)
    H_dip = Hamiltonian_interaction_NV_ZF(N,λ,one(b))
    H_dip = H_dip/norm(H_dip,Inf)

    for p in 1:length(τp)

        Threads.@threads for i in 1:length(τ0)

             for j in 1:length(disinter)
                
                H = H_dis 

                data_seq2[i,j,p]+=propagate_cyclic_finite(N,H,K4_Δ27_hamiltonian,gens_ham, τ0[i],τp[p],Id,2)/nsample

                data_seq3[i,j,p]+=propagate_cyclic_finite(N,H,K4_Δ27_eulerian,gens_ham, τ0[i],τp[p],Id,1)/nsample

                data_seq4[i,j,p]+=propagate_cyclic_finite(N,H,ZENITH,gens_zenith_ham,τ0[i],τp[p],Id,6)/nsample
                
                H = H_dip 


                data_seq2_bis[i,j,p]+=propagate_cyclic_finite(N,H,K4_Δ27_hamiltonian,gens_ham, τ0[i],τp[p],Id,2)/nsample

                data_seq3_bis[i,j,p]+=propagate_cyclic_finite(N,H,K4_Δ27_eulerian,gens_ham, τ0[i],τp[p],Id,1)/nsample

                data_seq4_bis[i,j,p]+=propagate_cyclic_finite(N,H,ZENITH,gens_zenith_ham,τ0[i],τp[p],Id,6)/nsample

            end
            
        end
        println("n=$n p=$p")
    end
end

jldopen(fig_dir*"/disorder and interaction NV ZF finite : K4_Δham.jld", "w") do file 
    write(file, "data", data_seq2)
end

jldopen(fig_dir*"/disorder and interaction NV ZF finite : K4_Δeul.jld", "w") do file 
    write(file, "data", data_seq3)
end

jldopen(fig_dir*"/disorder and interaction NV ZF finite : Zenith.jld", "w") do file 
    write(file, "data", data_seq4)
end

jldopen(fig_dir*"/disorder and interaction NV ZF finite bis : K4_Δham.jld", "w") do file 
    write(file, "data", data_seq2_bis)
end

jldopen(fig_dir*"/disorder and interaction NV ZF finite bis : K4_Δeul.jld", "w") do file 
    write(file, "data", data_seq3_bis)
end


jldopen(fig_dir*"/disorder and interaction NV ZF finite bis : Zenith.jld", "w") do file 
    write(file, "data", data_seq4_bis)
end
=#

data_seq2 = load(fig_dir*"/disorder and interaction NV ZF finite : K4_Δham.jld", "data")
data_seq3 = load(fig_dir*"/disorder and interaction NV ZF finite : K4_Δeul.jld", "data")
data_seq4 = load(fig_dir*"/disorder and interaction NV ZF finite : Zenith.jld", "data")

data_seq2_bis = load(fig_dir*"/disorder and interaction NV ZF finite bis : K4_Δham.jld", "data")
data_seq3_bis = load(fig_dir*"/disorder and interaction NV ZF finite bis : K4_Δeul.jld", "data")
data_seq4_bis = load(fig_dir*"/disorder and interaction NV ZF finite bis : Zenith.jld", "data")

datas = [data_seq2 , data_seq3 ,data_seq4]
datas_bis = [data_seq2_bis , data_seq3_bis ,data_seq4_bis]

colors = [ RGB(255/255 , 102/255 , 102/255) , RGB(153/255 , 204/255 , 255/255) , RGB(255/255 , 204/255 , 153/255)]
GLMakie.activate!()
fig = Figure(;size=(500,470))

ax = Axis3(fig[1,1], xlabel = L"\tau_0 \,\Delta", ylabel = L"\tau_p \, \Delta" ,title = L"J = 0", zlabel=L"\overline{D}",
xticks = ([-3,-2,-1], [L"10^{-3}",L"10^{-2}",L"10^{-1}"]),
yticks = ([-3,-2,-1], [L"10^{-3}",L"10^{-2}",L"10^{-1}"]), 
zticks = ([-2,-1], [L"10^{-2}",L"10^{-1}"]), 
xlabelsize=25.0f0,ylabelsize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0,
titlesize = 30.0f0,titlegap = -50,
zlabelsize=30.0f0,zticklabelsize=20.0f0,protrusions=(80,30,0,0),
zlabelrotation=0, zlabeloffset=80
)
alphas = [0.7,0.6,0.8]
for n in 1:3

    mat = Matrix{Any}(undef, length(noise), length(τp))
    for i in 1:length(noise), j in 1:length(τp)
        if datas[n][i,1,j] <-8 
            mat[i,j] = :transparent
        else 
            mat[i,j] = colors[n]
        end
    end

    surface!(τ0_log, τp_log , [datas[n][i,1,j] for i in 1:length(τ0_log), j in 1:length(τp)], 
    color=mat, alpha = alphas[n] , transparency = false)

    for m in [1,10,20,30,40,50]

        lines!(τ0_log, [τp_log[m] for p in 1:length(τ0_log)] , [datas[n][i,1,m] for i in 1:length(τ0_log)] , color = colors[n] , linestyle=(:solid))
        lines!([τ0_log[m] for p in 1:length(τp_log)], τp_log , [datas[n][m,1,i] for i in 1:length(τp_log)] , color = colors[n] , linestyle=(:solid))

    end
end
display(fig)

#save(fig_dir*"/disorder and interaction NV ZF : nested sequences finite pulses Δ.png", fig;px_per_unit=15.0)

fig = Figure(;size=(500,470))

ax = Axis3(fig[1,1], xlabel = L"\tau_0 \,J", ylabel = L"\tau_p \,J" ,title = L"\Delta = 0", zlabel=L"\overline{D}",
xticks = ([-3,-2,-1], [L"10^{-3}",L"10^{-2}",L"10^{-1}"]),
yticks = ([-3,-2,-1], [L"10^{-3}",L"10^{-2}",L"10^{-1}"]), 
zticks = ([-2,-1], [L"10^{-2}",L"10^{-1}"]), 
xlabelsize=25.0f0,ylabelsize=25.0f0,
xticklabelsize=20.0f0,yticklabelsize=20.0f0,
titlesize = 30.0f0,titlegap = -50,
zlabelsize=30.0f0,zticklabelsize=20.0f0,protrusions=(80,30,0,0),
zlabelrotation=0, zlabeloffset=80
)

alphas = [0.7,0.6,0.8]

for n in 1:3

    mat = Matrix{Any}(undef, length(noise), length(τp))
    for i in 1:length(noise), j in 1:length(τp)
        if datas_bis[n][i,1,j] < -8 
            mat[i,j] = :transparent
        else 
            mat[i,j] = colors[n]
        end
    end
    surface!(τ0_log, τp_log , [datas_bis[n][i,1,j] for i in 1:length(τ0_log), j in 1:length(τp)], 
    color=mat, alpha = alphas[n]  , transparency = false)

    for m in [1,10,20,30,40,50]

        lines!(τ0_log, [τp_log[m] for p in 1:length(τ0_log)] , [datas_bis[n][i,1,m] for i in 1:length(τ0_log)] , color = colors[n] , linestyle=(:solid))
        lines!([τ0_log[m] for p in 1:length(τp_log)], τp_log , [datas_bis[n][m,1,i] for i in 1:length(τp_log)] , color = colors[n] , linestyle=(:solid))

    end
end

display(fig)

#save(fig_dir*"/disorder and interaction NV ZF : nested sequences finite pulses J.png", fig;px_per_unit=15.0)


# pulse errors --> Fig 6 B

Id = one(b)
λ = GellMann()

N = 1

nsample = 1000

logϵ = LinRange(-4,-1,500)
ϵ = exp10.(logϵ)

Ufinal = DenseOperator(b,(V*X)^2)

data_seq1_e1 = zeros(length(logϵ))
data_seq1_e2 = zeros(length(logϵ))
data_seq1_e3 = zeros(length(logϵ))
data_seq1_e4 = zeros(length(logϵ))

data_seq2_e1 = zeros(length(logϵ))
data_seq2_e2 = zeros(length(logϵ))
data_seq2_e3 = zeros(length(logϵ))
data_seq2_e4 = zeros(length(logϵ))


#=
for n in 1:nsample

    error_tot = sum(OneBody(i,N,sum((rand()-0.5)*λ[j] for j in 1:8),Id) for i in 1:N) # a very random error
    error_tot = error_tot/norm(error_tot,Inf)

    Threads.@threads for e in 1:length(logϵ)

        data_seq1_e1[e] += propagate_pulseerrors(N,[ϵ[e]*error_tot,0*error_tot,0*error_tot,0*error_tot],K4_Δ27_hamiltonian,gens_ham,Ufinal)/nsample
        data_seq1_e2[e] += propagate_pulseerrors(N,[0*error_tot,ϵ[e]*error_tot,0*error_tot,0*error_tot],K4_Δ27_hamiltonian,gens_ham,Ufinal)/nsample
        data_seq1_e3[e] += propagate_pulseerrors(N,[0*error_tot,0*error_tot,ϵ[e]*error_tot,0*error_tot],K4_Δ27_hamiltonian,gens_ham,Ufinal)/nsample
        data_seq1_e4[e] += propagate_pulseerrors(N,[0*error_tot,0*error_tot,0*error_tot,ϵ[e]*error_tot],K4_Δ27_hamiltonian,gens_ham,Ufinal)/nsample

        data_seq2_e1[e] += propagate_pulseerrors(N,[ϵ[e]*error_tot,0*error_tot,0*error_tot,0*error_tot],K4_Δ27_eulerian,gens_ham,Ufinal)/nsample
        data_seq2_e2[e] += propagate_pulseerrors(N,[0*error_tot,ϵ[e]*error_tot,0*error_tot,0*error_tot],K4_Δ27_eulerian,gens_ham,Ufinal)/nsample
        data_seq2_e3[e] += propagate_pulseerrors(N,[0*error_tot,0*error_tot,ϵ[e]*error_tot,0*error_tot],K4_Δ27_eulerian,gens_ham,Ufinal)/nsample
        data_seq2_e4[e] += propagate_pulseerrors(N,[0*error_tot,0*error_tot,0*error_tot,ϵ[e]*error_tot],K4_Δ27_eulerian,gens_ham,Ufinal)/nsample

    end
end

datas_seq1 = [data_seq1_e1,data_seq1_e2,data_seq1_e3,data_seq1_e4]
datas_seq2 = [data_seq2_e1,data_seq2_e2,data_seq2_e3,data_seq2_e4]

jldopen(fig_dir*"/pulse errors : K4_Δham.jld", "w") do file 
    write(file, "data", datas_seq1)
end
jldopen(fig_dir*"/pulse errors : K4_Δeul.jld", "w") do file 
    write(file, "data", datas_seq2)
end
=#


datas_seq1 = load(fig_dir*"/pulse errors : K4_Δham.jld", "data")
datas_seq2 = load(fig_dir*"/pulse errors : K4_Δeul.jld", "data")

fig = Figure(;size=(370,370))
ax = Axis(fig[1,1], ylabel=L"1-\overline{D}", xlabel=L"\epsilon",
xlabelsize=30.0f0, ylabelsize=30.0f0,
xticks=([0,0.05,0.1],[L"0",L"0.05",L"0.1"]),
yticks=([0.7,0.8,0.9,1],[L"0.7",L"0.8",L"0.9",L"1"]),
xticklabelsize=20.0f0,yticklabelsize=20.0f0)

lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq1[1]), color=color_E, linewidth=3.0)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq1[2]), color=color_X, linewidth=3.0)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq1[3]), color=color_V, linewidth=3.0)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq1[4]), color=color_C, linewidth=3.0)


lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq2[1]), color=color_E, linewidth=3.0, linestyle=:dash)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq2[2]), color=color_X, linewidth=3.0, linestyle=:dash)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq2[3]), color=color_V, linewidth=3.0, linestyle=:dash)
lines!(ϵ,[1 for n in 1:length(ϵ)] - exp10.(datas_seq2[4]), color=color_C, linewidth=3.0, linestyle=:dash)

elems1 = [LineElement(color=:black, linewidth=3.0,linestyle=:solid),LineElement(color=:black, linewidth=3.0,linestyle=(:dash,:dense))]
elems2 = [PolyElement(color=color_E, strokewidth=2.0, strokecolor=:black),
PolyElement(color=color_X, strokewidth=2.0, strokecolor=:black),
PolyElement(color=color_V, strokewidth=2.0, strokecolor=:black),
PolyElement(color=color_C, strokewidth=2.0, strokecolor=:black)]
leg1 = ["Hamiltonian" , "Eulerian"]
leg2 = ["E", "X", "V", "C"]
axislegend(ax,[elems1, elems2],[leg1,leg2],["",""],tellwidth=false,tellheight=false , nbanks=4,position=:lb,labelsize=16.0f0)

display(fig)

#save(fig_dir*"/disorder and interaction NV ZF : nested sequences errors.png", fig;px_per_unit=15.0)



# Cayley graph Δ27 --> Fig 6 A

function R2(θ)
    M = zeros(2,2)
    M[1,1] = cos(θ)
    M[2,2] = cos(θ)
    M[1,2] = -sin(θ)
    M[2,1] = sin(θ)

    return M
end

t = [[0,1.0]]
for n in 1:2 
    push!(t, R2(2*pi/3)*t[end])
end

for m in 1:2
    for n in 1:3
        push!(t, (m+1)*t[n])
    end
end

for n in 1:3 
    t[3+n] = R2(14/180*pi)*t[3+n]
    t[6+n] = R2(28/180*pi)*t[6+n]
end

l1 = [[1,2],[2,3],[3,1],[4,5],[5,6],[6,4],[7,8],[8,9],[9,7]]
l2 = [[1,4],[4,7],[7,1],[2,5],[5,8],[8,2],[3,6],[6,9],[9,3]]

using CairoMakie
CairoMakie.activate!()
fig = Figure()

ax = Axis(fig[1,1], aspect=DataAspect())
for n in 1:length(l1)
    p1 = t[l1[n][1]]
    p2 = t[l1[n][2]] 

    u = p2-p1

    lines!(Point2.([p1,p2]), color=color_E,linewidth=3.0)
    arrows2d!(Point2.([p1]) , Point2.([u/2]), color=color_E,tailwidth=0)
end
for n in 1:length(l2)
    p1 = t[l2[n][1]]
    p2 = t[l2[n][2]] 

    u = p2-p1

    lines!(Point2.([p1,p2]), color=color_C,linewidth=3.0)
    arrows2d!(Point2.([p1]) , Point2.([u/2]), color=color_C,tailwidth=0)
end

scatter!(Point2.(t), color=:black, markersize=15.0,depth_shift=-1)

hidespines!()
hidedecorations!()
display(fig)
#save(fig_dir*"/Cayley graph Δ27.pdf",fig)



###############################################
# Orientations of Σ(72 x 3) --> Section 4.4.5 #
###############################################

# orientation 1 

H_X_1 = pi/2*(1/9)*[sqrt(3) sqrt(78) 0;
sqrt(78) -sqrt(3) 0;
0 0 0]

a_v = -(sqrt(234)/36 + sqrt(3)/(18*sqrt(78))) + 72*im/(18*sqrt(78))
b_ϕ(ϕ) = -sqrt(3)/2*(0.5+sqrt(3)+im*(1+sqrt(3)/2))/sqrt(5+2*sqrt(3))*exp(im*ϕ)
b_ϕp(ϕ) = -sqrt(3)/2*exp(im*ϕ)

H_V_1(ϕ) = -pi/2*sqrt(3)/2*[2/9 a_v 0;
conj(a_v) 5/18 b_ϕ(ϕ);
0 conj(b_ϕ(ϕ)) -0.5]

H_V_1p(ϕ) = -pi/2*sqrt(3)/2*[2/9 a_v 0;
conj(a_v) 5/18 b_ϕp(ϕ);
0 conj(b_ϕp(ϕ)) -0.5]

test_1 = -pi/2*[1/3*1/sqrt(3) 1/sqrt(3)*(-5/3*sqrt(2/13)+sqrt(6/13)*im) 0;
1/sqrt(3)*(-5/3*sqrt(2/13)-sqrt(6/13)*im) -1/3*1/sqrt(3)+sqrt(3)/4 -3/4;
0 -3/4 -sqrt(3)/4]

H_V_1p(0.0) - test_1

U_X_1 = exp(im*H_X_1)
U_V_1(ϕ) = exp(im*H_V_1(ϕ))

# orientation 2

H_X_2 = pi/2*[-1 0 0;
0 5/(14+sqrt(3)) -sqrt(5)*sqrt(9+sqrt(3))/(14+sqrt(3));
0 -sqrt(5)*sqrt(9+sqrt(3))/(14+sqrt(3)) (9+sqrt(3))/(14+sqrt(3))]

θspecial = atan(sqrt((9+sqrt(3))/5))

c_ϕ(ϕ) = sqrt(14+sqrt(3))/(2*sqrt(6))*(-1/sqrt(5)+2*im/sqrt(5))*exp(im*ϕ)
d_v = -(3-sqrt(3))/8*sin(2*θspecial) -sqrt(3+sqrt(3))/(2*sqrt(6)*sqrt(5)*sqrt(5+2*sqrt(3)))*(8*im+cos(2*θspecial))

H_V_2(ϕ) = -pi/2*sqrt(3)/2*[(1+sqrt(3))/4 c_ϕ(ϕ) 0;
conj(c_ϕ(ϕ)) -(17+7*sqrt(3))/(4*(14+sqrt(3))) d_v;
0 conj(d_v) -2*sqrt(3)/(14+sqrt(3)) ]

U_X_2 = exp(im*H_X_2)
U_V_2(ϕ) = exp(im*H_V_2(ϕ))


# representation of the pulses 

using CairoMakie
CairoMakie.activate!()

fig = Figure(;size=(300,300))
ax = Axis(fig[1,1])

c_x = [cos(i) for i in LinRange(0,2*pi, 100)]
c_y = [sin(i) for i in LinRange(0,2*pi, 100)]

lines!(c_x,c_y, color=:black)
text!(Point2([1.2,0]), text=L"x", fontsize=19.0f0)
text!(Point2([0,1.2]), text=L"y", fontsize=19.0f0)

band!(Point2f.([[0,0] for i in LinRange(0,acos(-5/sqrt(52)), 100)]) , Point2f.([[cos(i),sin(i)] for i in LinRange(0,acos(-5/sqrt(52)), 100)]) , color=(:red,0.4))
lines!([0,0],[-1.2,1.2], color=:black)
lines!([-1.2,1.2],[0,0], color=:black)
arrows2d!(Point2([0,0]) , Point2([-5/sqrt(52),3*sqrt(3)/sqrt(52)]), color=:red)

lines!(Point2.([[-5/sqrt(52),3*sqrt(3)/sqrt(52)] , [0,3*sqrt(3)/sqrt(52)]]), color=:black, linestyle=:dot)
text!(Point2([0,3*sqrt(3)/sqrt(52)]), text=L"\frac{3\sqrt{3}}{\sqrt{52}}", align=(:left, :center), offset = (5,-2))

lines!(Point2.([[-5/sqrt(52),3*sqrt(3)/sqrt(52)] , [-5/sqrt(52),0]]), color=:black, linestyle=:dot)
text!(Point2([-5/sqrt(52),0]), text=L"-\frac{5}{\sqrt{52}}", align=(:center, :top), offset = (0,-5))

hidedecorations!(ax)
hidespines!(ax)
xlims!(-1.2,1.5)
ylims!(-1.2,1.5)
display(fig)

#save(fig_dir*"/phase 1-0 - orientation A.pdf", fig)

fig = Figure(;size=(300,300))
ax = Axis(fig[1,1])

c_x = [cos(i) for i in LinRange(0,2*pi, 100)]
c_y = [sin(i) for i in LinRange(0,2*pi, 100)]

lines!(c_x,c_y, color=:black)
text!(Point2([1.2,0]), text=L"x", fontsize=19.0f0)
text!(Point2([0,1.2]), text=L"y", fontsize=19.0f0)

band!(Point2f.([[0,0] for i in LinRange(0,acos(-5/sqrt(52)), 100)]) , Point2f.([[cos(i),sin(i)] for i in LinRange(0,acos(-5/sqrt(52)), 100)]) , color=(:red,0.4))
lines!([0,0],[-1.2,1.2], color=:black)
lines!([-1.2,1.2],[0,0], color=:black)
arrows2d!(Point2([0,0]) , Point2([-5/sqrt(52),3*sqrt(3)/sqrt(52)]), color=:red)

lines!(Point2.([[-5/sqrt(52),3*sqrt(3)/sqrt(52)] , [0,3*sqrt(3)/sqrt(52)]]), color=:black, linestyle=:dot)
text!(Point2([0,3*sqrt(3)/sqrt(52)]), text=L"\frac{3\sqrt{3}}{\sqrt{52}}", align=(:left, :center), offset = (5,-2))

lines!(Point2.([[-5/sqrt(52),3*sqrt(3)/sqrt(52)] , [-5/sqrt(52),0]]), color=:black, linestyle=:dot)
text!(Point2([-5/sqrt(52),0]), text=L"-\frac{5}{\sqrt{52}}", align=(:center, :top), offset = (0,-5))

hidedecorations!(ax)
hidespines!(ax)
xlims!(-1.2,1.5)
ylims!(-1.2,1.5)
display(fig)

#save(fig_dir*"/phase 1-0 - orientation A.pdf", fig)


fig = Figure()

fig = Figure(;size=(300,300))
ax = Axis(fig[1,1])

c_x = [cos(i) for i in LinRange(0,2*pi, 100)]
c_y = [sin(i) for i in LinRange(0,2*pi, 100)]

lines!(c_x,c_y, color=:black)
text!(Point2([1.2,0]), text=L"x", fontsize=19.0f0)
text!(Point2([0,1.2]), text=L"y", fontsize=19.0f0)

lines!([0,0],[-1.2,1.2], color=:black)
lines!([-1.2,1.2],[0,0], color=:black)
arrows2d!(Point2([0,0]) , Point2([1,0]), color=:blue)

hidedecorations!(ax)
hidespines!(ax)
xlims!(-1.2,1.5)
ylims!(-1.2,1.5)
display(fig)

#save(fig_dir*"/phase -1-0 - orientation A.pdf", fig)

GLMakie.activate!()

ϕa_V = acos(5/(2*sqrt(13)))
ϕb_V = 0

fig = Figure(;size=(300,300))
ax = Axis(fig[1,1])

c_x = [cos(i) for i in LinRange(0,2*pi, 100)]
c_y = [sin(i) for i in LinRange(0,2*pi, 100)]

lines!(c_x,c_y, color=:black)
text!(Point2([1.2,0]), text=L"x", fontsize=19.0f0)
text!(Point2([0,1.2]), text=L"y", fontsize=19.0f0)


band!(Point2f.([[0,0] for i in LinRange(0,2*pi, 100)]) , Point2f.([[cos(i),sin(i)] for i in LinRange(0,2*pi, 100)]) , color = 1:100,colormap = :hsv, alpha=0.9)
lines!([0,0],[-1.2,1.2], color=:black)
lines!([-1.2,1.2],[0,0], color=:black)

arrows2d!(Point2([0,0]) , Point2([cos(ϕa_V),sin(ϕa_V)]), color=:black)

arrows2d!(Point2([0,0]) , Point2([1,0]), color=:black)

hidedecorations!(ax)
hidespines!(ax)
xlims!(-1.2,1.5)
ylims!(-1.2,1.5)
display(fig)

#save(fig_dir*"/phase X- orientation 1.png", fig,px_per_unit=15)

#####################################
# Σ(168) = Δ(24) x <<Y>> --> fig 11 #
#####################################

using DoubleFloats
# Weyl element 

W1 = zeros(3,3)
W1[1,1] = -1.0+0*im
W1[2,3] = 1
W1[3,2] = 1
det(W1)

W2 = zeros(3,3)
W2[2,2] = -1.0+0*im
W2[1,3] = 1
W2[3,1] = 1

W3 = zeros(3,3)
W3[3,3] = -1.0+0*im
W3[1,2] = 1
W3[2,1] = 1

W4 = zeros(3,3)
W4[1,2] = 1.0+0*im
W4[2,3] = 1
W4[3,1] = 1

W5 = zeros(3,3)
W5[1,3] = 1.0+0*im
W5[3,2] = 1
W5[2,1] = 1

# generator Σ(168)

β = exp(ComplexDF64(2.0*im)/7.0*pi)

A3 = im*[β^2-β^5 β-β^6 β^4-β^3;
β-β^6 β^4-β^3 β^2-β^5;
β^4-β^3 β^2-β^5 β-β^6]/sqrt(7)

B3 = im*[β^3-β^6 β^3-β β-1.0;
β^2-1.0 β^6-β^5 β^6-β^2;
β^5-β^4 β^4-1.0 β^5-β^3]/sqrt(7)

det(A3)
det(B3)

# unitary transfo 

w1 = B3*(A3*B3)*(A3*B3)*B3*A3*B3*A3

N = (28*(-342125 - 349668*β + 283769*β^2 + 9406*β^3 - 501928*β^4 + 287955*β^6))^(1/3)

x = (−64.0 - 15*β - 56*β^3 - 46*β^4 + 5*β^5 - 27*β^6)/N

y = (73.0*β + 156*β^2 + 46*β^3 + 12*β^4 + 137*β^5 + 115*β^6)/N

z = (15.0 + 34*β − 35*β^2 − 23*β^3 + 41*β^4 − 46*β^6)/N

w2 = [x y z;
z x y;
y z x]

w = w2*w1

det(w2)
log10(sqrt(abs(1 - abs(tr(w2*adjoint(w2)))/3)))

# generator Δ(24)

M1 = [0.0+0.0*im 1.0 0;
0 0 -1.0;
-1.0 0 0]

M2 = [-1.0+0.0*im 0 0;
0 0 -1.0;
0 -1.0 0]

O_mat = load(data_dir*"Adjacency matrix : O - ST,U.jld", "data")
#O_path = Hierholzer_random(O_mat)
#O_seq = path_to_sequence(O_mat,O_path)
#jldopen(data_dir*"path O.jld", "w") do file 
#    write(file, "data", O_seq)
#end

O_seq = load(data_dir*"path O.jld", "data")

O_gens = adjoint(w)*M1*w,  adjoint(w)*M2*w

wM1 = adjoint(w)*M1*w
wM2 = adjoint(w)*M2*w

O_gens_U_W0 = DenseOperator(b,wM1),DenseOperator(b,wM2)
O_gens_U_W1 = DenseOperator(b,adjoint(W1)*wM1*W1),DenseOperator(b,adjoint(W1)*wM2*W1)
O_gens_U_W2 = DenseOperator(b,adjoint(W2)*wM1*W2),DenseOperator(b,adjoint(W2)*wM2*W2)
O_gens_U_W3 = DenseOperator(b,adjoint(W3)*wM1*W3),DenseOperator(b,adjoint(W3)*wM2*W3)
O_gens_U_W4 = DenseOperator(b,adjoint(W4)*wM1*W4),DenseOperator(b,adjoint(W4)*wM2*W4)
O_gens_U_W5 = DenseOperator(b,adjoint(W5)*wM1*W5),DenseOperator(b,adjoint(W5)*wM2*W5)

# data 

λ = GellMann()

N = 3 
τ = 1.0
nsample = 100

dis = LinRange(-3,0,25)
inter = LinRange(-3,0,25)

data_O_U1 = zeros(length(dis),length(inter))
data_O_U2 = zeros(length(dis),length(inter))
data_O_U3 = zeros(length(dis),length(inter))
data_O_U4 = zeros(length(dis),length(inter))
data_O_U5 = zeros(length(dis),length(inter))
data_O_U6 = zeros(length(dis),length(inter))

data_nodd = zeros(length(dis),length(inter))


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis,Inf)
    H_dip = Hamiltonian_interaction_NV(N,λ,one(b))
    H_dip = H_dip/norm(H_dip,Inf)

    for i in 1:length(dis)

        Threads.@threads for j in 1:length(inter)

            H = exp10(dis[i])*H_dis + exp10(inter[j])*H_dip

            data_O_U1[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W0, τ)/nsample
            data_O_U2[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W1, τ)/nsample
            data_O_U3[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W2, τ)/nsample
            data_O_U4[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W3, τ)/nsample
            data_O_U5[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W4, τ)/nsample
            data_O_U6[i,j]+=propagate_ideal(N,H,O_seq,O_gens_U_W5, τ)/nsample

            data_nodd[i,j]+=propagate_nodd(N,H,length(O_seq),τ)/nsample

        end

        println("n=$n i=$i")
    end
end

jldopen(fig_dir*"/disorder and interaction NV : O U1.jld", "w") do file 
    write(file, "data", data_O_U1)
end

jldopen(fig_dir*"/disorder and interaction NV : O U2.jld", "w") do file 
    write(file, "data", data_O_U2)
end

jldopen(fig_dir*"/disorder and interaction NV : O U3.jld", "w") do file 
    write(file, "data", data_O_U3)
end

jldopen(fig_dir*"/disorder and interaction NV : O U4.jld", "w") do file 
    write(file, "data", data_O_U4)
end

jldopen(fig_dir*"/disorder and interaction NV : O U5.jld", "w") do file 
    write(file, "data", data_O_U5)
end

jldopen(fig_dir*"/disorder and interaction NV : O U6.jld", "w") do file 
    write(file, "data", data_O_U6)
end

jldopen(fig_dir*"/disorder and interaction NV : nodd.jld", "w") do file 
    write(file, "data", data_nodd)
end
=#


data_O_U1 = load(fig_dir*"/disorder and interaction NV : O U1.jld", "data")
data_O_U2 = load(fig_dir*"/disorder and interaction NV : O U2.jld", "data")
data_O_U3 = load(fig_dir*"/disorder and interaction NV : O U3.jld", "data")
data_O_U4 = load(fig_dir*"/disorder and interaction NV : O U4.jld", "data")
data_O_U5 = load(fig_dir*"/disorder and interaction NV : O U5.jld", "data")
data_O_U6 = load(fig_dir*"/disorder and interaction NV : O U6.jld", "data")
data_nodd = load(fig_dir*"/disorder and interaction NV : nodd.jld", "data")

colors = [RGB(178/255,102/255,255/255) , 
RGB(255/255,204/255,153/255) , 
RGB(255/255 , 102/255 , 102/255) , 
RGB(153/255 , 153/255 ,255/255) , 
RGB(153/255 , 255/255 ,153/255) ,
RGB(255/255,255/255,153/255) ,
RGB(192/255,192/255,192/255)]
GLMakie.activate!()

datas = [data_O_U1 , data_O_U2 , data_O_U3 , data_O_U4 , data_O_U5 , data_O_U6 , data_nodd]
lab = [L"U" , L"U(23),\, U(132)" ,L"U,\,U(13)" ,L"U(12),\,U(123)",L"U(123)",L"U(132)", L"\mathrm{NoDD}"]

fig = Figure()
ax = Axis3(fig[1,1], xlabel = L"\tau_0\Delta", ylabel = L"\tau_0 J" , zlabel=L"\overline{D}",
xlabelsize=25.0f0, ylabelsize=25.0f0,zlabelsize=25.0f0,
zlabelrotation=0,zlabeloffset=85,
xticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),yticks = ([-3,-2,-1], [L"10^{-3}", L"10^{-2}", L"10^{-1}"]),
zticks = ([-5 ,-3, -1], [L"10^{-5}", L"10^{-3}", L"10^{-1}"]),
xticklabelsize=22.0f0,yticklabelsize=22.0f0,
zticklabelsize=22.0f0,zlabelalign=(:left, :bottom),
ylabeloffset=50)

Ns = [2 , 3 , 4 , 7]

for n in Ns
    cmat = Matrix{Colorant}(undef, 25,25)
    cmat[:,:].= colors[n]
    surface!(ax, dis, inter, datas[n], color=cmat, alpha=0.6)
end
midlines = [1 , 5 , 10 , 15 , 20 , 25]

for n in 1:length(midlines)
    for m in Ns

        lines!(ax, [dis[midlines[n]]], inter, datas[m][midlines[n],:], color=colors[m], linewidth=2.0)
        lines!(ax, dis, [inter[midlines[n]]], datas[m][:,midlines[n]], color=colors[m], linewidth=2.0)

    end
end
lab_mark = [[PolyElement(color=(colors[n],0.7),strokewidth=2.0,strokecolor=colors[n]), LineElement(color=colors[n],linewidth=2.0)] for n in Ns]
Legend(fig[2,1], lab_mark, [lab[k] for k in Ns],framevisible=false, 
nbanks=2, tellwidth=false, tellheight=true,
labelsize=25.0f0,colgap=40,rowgap=10)
display(fig)

#save(fig_dir*"/disorder and interaction NV : UOU.png", fig;px_per_unit = 15.0)


###################################
# Comparison with DROID --> Fig13 #
###################################

U = eigvecs(X)
round.(adjoint(U)*X*U,digits=13) # ok 

UVU = adjoint(U)*V*U
UCU = adjoint(U)*C*U

# best Σ36  

Σ36x3_gensP = DenseOperator(b,adjoint(W5)*UVU*W5),DenseOperator(b,adjoint(W5)*UCU*W5)
Σ36x3_gensH =  DenseOperator(b,-pi/2*sqrt(3)/3* adjoint(W5)*adjoint(U)*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]*U*W5) , DenseOperator(b,-2*pi/3*adjoint(W5)*adjoint(U)*[0 0 0 ; 0 1 0 ; 0 0 -1]*U*W5)
Σ36seq = load(data_dir*"path Σ36.jld", "data")

Σ36sequence = Vector{Vector}(undef, length(Σ36seq))
for n in 1:length(Σ36seq)
    Σ36sequence[n] = [Σ36seq[n]]
end
τp36 = [1 , 4/3]

# one Σ72, not the best one 

Σ72x3_gensP = DenseOperator(b,V),DenseOperator(b,X)
Σ72x3_gensH = DenseOperator(b,-pi/2*sqrt(3)/3*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]),DenseOperator(b, U*[pi/2 0 0 ; 0 -pi/2 0; 0 0 0]*adjoint(U)) 
Σ72seq = load(data_dir*"path Σ72.jld", "data")

Σ72sequence = Vector{Vector}(undef, length(Σ72seq))
for n in 1:length(Σ72seq)
    Σ72sequence[n] = [Σ72seq[n]]
end
τp72 = [1 , 1]

# nested sequence 

K4_gensP = DenseOperator(b,E),DenseOperator(b,X),DenseOperator(b,V)
K4_gensH = DenseOperator(b, -2*pi/3*sqrt(1/3)*[0 -im im; im 0 -im ; -im im 0]) , DenseOperator(b, U*[pi/2 0 0 ; 0 -pi/2 0; 0 0 0]*adjoint(U)) , 
DenseOperator(b,-pi/2*sqrt(3)/3*[-1 -1 -1 ; -1 1/2 1/2 ; -1 1/2 1/2]) 

Seq_nested = [[1] , [1] , [1,2], [1] , [1] , [1,3],[1] , [1] , [1,2],[1] , [1] , [1,3]]

Ufinal_seq = DenseOperator(b, V*X*V*X)

# best O

M1 = [0.0+0.0*im 1.0 0;
0 0 -1.0;
-1.0 0 0]

eigen(M1).values
H_M1 = (eigen(M1).vectors)*[2*pi/3 0 0 ; 0 -2*pi/3 0 ; 0 0 0]*adjoint(eigen(M1).vectors)

norm(exp(-im*H_M1) - M1)

M2 = [-1.0+0.0*im 0 0;
0 0 -1.0;
0 -1.0 0]

eigen(M2).values
H_M2 = (eigen(M2).vectors)*[pi 0 0 ; 0 -pi 0 ; 0 0 0]*adjoint(eigen(M2).vectors)

exp(-im*H_M2) - M2

w1 = B3*(A3*B3)^2*B3*A3*B3*A3

N = (28*(-342125 - 349668*β + 283769*β^2 + 9406*β^3 - 501928*β^4 + 287955*β^6))^(1/3)
x = (−64 - 15*β - 56*β^3 - 46*β^4 + 5*β^5 - 27*β^6)/N
y = (73*β + 156*β^2 + 46*β^3 + 12*β^4 + 137*β^5 + 115*β^6)/N
z = (15 + 34*β − 35*β^2 − 23*β^3 + 41*β^4 − 46*β^6)/N

w2 = [x y z;
z x y;
y z x]

w = ComplexF64.(w2*w1 )

wM1 = adjoint(w)*M1*w
wM2 = adjoint(w)*M2*w

O_gensP = DenseOperator(b,adjoint(W1)*wM1*W1),DenseOperator(b,adjoint(W1)*wM2*W1)
O_gensH = DenseOperator(b,adjoint(W1)*adjoint(w)*H_M1*w*W1),DenseOperator(b,adjoint(W1)*adjoint(w)*H_M2*w*W1)

Oseq = load(data_dir*"path O.jld", "data")
Osequence = Vector{Vector}(undef, length(Oseq))
for n in 1:length(Oseq)
    Osequence[n] = [Oseq[n]]
end
τpO = [4/3 , 2]

# DROID-C3PO 

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

    return DenseOperator(b,exp(-im*θ*1/sqrt(2)*[0 exp(-im*ϕ1) 0;
    exp(im*ϕ1) 0 exp(-im*ϕ2);
    0 exp(im*ϕ2) 0]))

end

function add_pulse_droid_ham(ϕ1_ax, ϕ2_ax,θ)

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

    return DenseOperator(b,θ*1/sqrt(2)*[0 exp(-im*ϕ1) 0;
    exp(im*ϕ1) 0 exp(-im*ϕ2);
    0 exp(im*ϕ2) 0])

end

DC3_gensP = [
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

]

DC3_gensH = [add_pulse_droid_ham("-x","-x",pi), #1
add_pulse_droid_ham("x","x",pi),
add_pulse_droid_ham("-x","x",pi),
add_pulse_droid_ham("x","-x",pi),

add_pulse_droid_ham("-y","-y",pi), #5
add_pulse_droid_ham("y","y",pi),
add_pulse_droid_ham("-y","y",pi),
add_pulse_droid_ham("y","-y",pi),

add_pulse_droid_ham("-x","-x",pi/2), # 9
add_pulse_droid_ham("x","x",pi/2),
add_pulse_droid_ham("-x","x",pi/2),
add_pulse_droid_ham("x","-x",pi/2),

add_pulse_droid_ham("-y","-y",pi/2), # 13
add_pulse_droid_ham("y","y",pi/2),
add_pulse_droid_ham("-y","y",pi/2),
add_pulse_droid_ham("y","-y",pi/2),

]


DROIDC3PO = [[1],[9,14],[1],[9,14],[2],[10,12],[5],[13,16],[6],[14,15],[6],
[14,9],[8],[16,12],[8],[16,12],[8],[16,13],[3],[11,10],[4],[12,9],[3],
[11,15],[1],[9,14],[2],[10,13],[2],[10,12],[5],[13,16],[5],[13,16],[6],
[14,9],[8],[16,12],[7],[15,11],[8],[16,13],[3],[11,10],[3],[11,10],[3],
[3],[4],[10,12],[4],[10,12],[4],[14,16],[8],[12,15],[7],[11,16],[8],
[9,13],[5],[16,14],[6],[16,14],[6],[11,10],[2],[14,10],[2],[13,9],[1],
[15,12],[4],[9,11],[3],[10,12],[4],[14,16],[8],[11,16],[8],[11,16],[8],
[9,13],[5],[15,13],[5],[16,14],[6],[11,10],[2],[13,9],[1],[13,9],[1],[7]] 


τp_DC3 = [2,2,2,2,
2,2,2,2,
1,1,1,1,
1,1,1,1]

function propagate_finite_err(N,H,Seq,Gens,Gens_err,ϵ,τ0,τps,Id,Mrepet,Ufinal)

    U = Id
    for n in 2:N 
        U = tensor(U,Id)
    end

    U0 = exp(dense(-im*τ0*H))

    for n in 1:length(Seq)

        U = U0*U
        for m in 1:length(Seq[n])

            pulse = exp(dense(-im*τps[Seq[n][m]]*(sum(OneBody(i,N,Gens[Seq[n][m]]+(rand()-0.5)*ϵ*Gens_err[Seq[n][m]],Id) for i in 1:N)/τps[Seq[n][m]] + H)))

            U = pulse*U

        end
    end

    U = U^Mrepet

    Uf = Ufinal 
    for i in 2:N 
        Uf = tensor(Uf,Ufinal)
    end

    return log10(sqrt(abs(1 - abs(tr(U*dagger(Uf)))/(3^N))))

end

function propagate_finite_err_subsamp(N,H,Seq,Gens,Gens_err,ϵ,τ0,τps,Id,Mrepet,Ufinal)

    U = Id
    for n in 2:N 
        U = tensor(U,Id)
    end

    U0 = exp(dense(-im*τ0*H))

    for n in 1:length(Seq)

        U = U0*U
        for m in 1:length(Seq[n])

            pulse = exp(dense(-im*τps[Seq[n][m]]*(sum(OneBody(i,N,Gens[Seq[n][m]]+(rand()-0.5)*ϵ*Gens_err[Seq[n][m]],Id) for i in 1:N)/τps[Seq[n][m]] + H)))

            U = pulse*U

        end
    end

    Uf = Ufinal 
    for i in 2:N 
        Uf = tensor(Uf,Ufinal)
    end

    return [1-sqrt(abs(1 - abs(tr(U^n*dagger(Uf)^n))/(3^N))) for n in 1:Mrepet]

end

function propagate_finite_err_subsampDROID(N,H,Seq,Gens,Gens_err,ϵ,τ0,τps,Id,Mrepet,Ufinal)

    U = Id
    for n in 2:N 
        U = tensor(U,Id)
    end

    U0 = exp(dense(-im*τ0*H))
    U02 = exp(dense(-im*τ0/2*H))

    for n in 1:length(Seq)

        if n == 1
            U = U02*U
        else 
            U = U0*U
        end

        for m in 1:length(Seq[n])

            pulse = exp(dense(-im*τps[Seq[n][m]]*(sum(OneBody(i,N,Gens[Seq[n][m]]+(rand()-0.5)*ϵ*Gens_err[Seq[n][m]],Id) for i in 1:N)/τps[Seq[n][m]] + H)))

            U = pulse*U

        end
    end

    U = U02*U 

    Uf = Ufinal 
    for i in 2:N 
        Uf = tensor(Uf,Ufinal)
    end

    return [1-sqrt(abs(1 - abs(tr(U^n*dagger(Uf)^n))/(3^N))) for n in 1:Mrepet]

end

function propagate_finite_err_subsampCYCLIC(N,H,Seq,Gens,Gens_err,ϵ,τ0,τps,Id,Mrepet,Ufinal)

    U = Id
    for n in 2:N 
        U = tensor(U,Id)
    end

    U0 = exp(dense(-im*τ0*H))

    for n in 1:length(Seq)

        U = U0*U

        for m in 1:length(Seq[n])

            pulse = exp(dense(-im*τps[Seq[n][m]]*(sum(OneBody(i,N,Gens[Seq[n][m]]+(rand()-0.5)*ϵ*Gens_err[Seq[n][m]],Id) for i in 1:N)/τps[Seq[n][m]] + H)))

            U = pulse*U

        end
    end

    for n in 1:length(Seq)

        for m in 1:length(Seq[end+1-n])

            pulse = exp(dense(-im*τps[Seq[end+1-n][end+1-m]]*(sum(OneBody(i,N,-Gens[Seq[end+1-n][end+1-m]]+(rand()-0.5)*ϵ*Gens_err[Seq[end+1-n][end+1-m]],Id) for i in 1:N)/τps[Seq[end+1-n][end+1-m]] + H)))

            U = pulse*U

        end

        U = U0*U
    end
    Uf = Ufinal 
    for i in 2:N 
        Uf = tensor(Uf,Ufinal)
    end

    return [1-sqrt(abs(1 - abs(tr(U^n*dagger(Uf)^n))/(3^N))) for n in 1:Mrepet]

end

λ = GellMann()
Id = one(b)
N = 3 
τ = 1.0
nsample = 20

τ0 = LinRange(-3,0,20)
τp = LinRange(-3,0,20)

data_DC3 = [zeros(length(τ0),length(τp)) for i in 1:3]
data_nodd = [zeros(length(τ0),length(τp)) for i in 1:3]

DC3_gens_err = [one(b) for i in 1:length(DC3_gensH)]

Ufinal = DenseOperator(b, [-1 0 0;0 1 0; 0 0 -1])

#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis)
    H_disq = Hamiltonian_disorder_NV_quad(N,λ,one(b))
    H_disq = H_disq/norm(H_disq)
    H_dip = Hamiltonian_interaction_NV(N,λ,one(b))
    H_dip = H_dip/norm(H_dip)   


    for i in 1:length(τ0)

        Threads.@threads for j in 1:length(τp)

            data_DC3[1][i,j] += propagate_finite_err(N,H_dis,DROIDC3PO,DC3_gensH,DC3_gens_err,0.0,exp10(τ0[i]),exp10(τp[j])*τp_DC3,Id,1,Ufinal)/nsample
            data_DC3[2][i,j] += propagate_finite_err(N,H_disq,DROIDC3PO,DC3_gensH,DC3_gens_err,0.0,exp10(τ0[i]),exp10(τp[j])*τp_DC3,Id,1,Ufinal)/nsample
            data_DC3[3][i,j] += propagate_finite_err(N,H_dip,DROIDC3PO,DC3_gensH,DC3_gens_err,0.0,exp10(τ0[i]),exp10(τp[j])*τp_DC3,Id,1,Ufinal)/nsample

            data_nodd[1][i,j]+=propagate_nodd(N,H_dis,96,(exp10(τ0[i]) + 2*exp10(τp[j])))/nsample
            data_nodd[2][i,j]+=propagate_nodd(N,H_disq,96,(exp10(τ0[i]) + 2*exp10(τp[j])))/nsample
            data_nodd[3][i,j]+=propagate_nodd(N,H_dip,96,(exp10(τ0[i]) + 2*exp10(τp[j])))/nsample

        end

        println("n=$n i=$i")

    end
end

fig = Figure()
ax = Axis3(fig[1,1])
surface!(τ0,τp,data_DC3[1])
surface!(τ0,τp,data_nodd[1])
display(fig)

fig = Figure()
ax = Axis3(fig[1,1])
surface!(τ0,τp,data_DC3[2])
surface!(τ0,τp,data_nodd[2])
display(fig)

fig = Figure()
ax = Axis3(fig[1,1])
surface!(τ0,τp,data_DC3[3])
surface!(τ0,τp,data_nodd[3])
display(fig)
=#

λ = GellMann()
Id = one(b)
N = 4
nsample = 750

MHz = 1.0 
kHz = 10^(-3)*MHz 
μs = 1/MHz 
ns = 10^(-3)*μs 

τpi2 = 8*ns
τ0 = 25*ns

exps = [(2*pi*4*MHz,2*pi*MHz,2*pi*40*kHz,0.01) , (2*pi*4*MHz,2*pi*MHz,2*pi*40*kHz*10,0.01), (2*pi*4*MHz/10,2*pi*MHz/10,2*pi*40*kHz,0.01), (2*pi*4*MHz/10,2*pi*MHz/10,2*pi*40*kHz/10,0.01) ,

(2*pi*4*MHz,2*pi*MHz,2*pi*40*kHz,0.05) , (2*pi*4*MHz,2*pi*MHz,2*pi*40*kHz*10,0.05), (2*pi*4*MHz/10,2*pi*MHz/10,2*pi*40*kHz,0.05), (2*pi*4*MHz/10,2*pi*MHz/10,2*pi*40*kHz/10,0.05) ]

nexps = 8

n_DC3 = 7 
n_Σ72 = 6
n_Σ36 = 12 
n_Σ72_cyclic = 3
n_Σ36_cyclic = 6
n_O = 14
n_O_cyclic = 7

t_DC3 = 96*τ0 + 96*2*τpi2 
t_72 = 144*τ0 + 144*τpi2 
t_36 = 72*τ0 + 36*τpi2 + 36*(4/3)*τpi2/2
t_72_cyclic =  2*t_72
t_36_cyclic = 2*t_36
t_O = 48*τ0 + 24*2*τpi2 + 24*4/3*τpi2
t_O_cyclic =2*t_O

data_DC3 = [zeros(n_DC3) for p in 1:nexps]
data_Σ72 = [zeros(n_Σ72) for p in 1:nexps]
data_Σ36 = [zeros(n_Σ36) for p in 1:nexps]
data_Σ72_cyclic = [zeros(n_Σ72_cyclic) for p in 1:nexps]
data_Σ36_cyclic = [zeros(n_Σ36_cyclic) for p in 1:nexps]
data_O = [zeros(n_O) for p in 1:nexps]
data_O_cyclic = [zeros(n_O_cyclic) for p in 1:nexps]


#=
for n in 1:nsample 

    H_dis = Hamiltonian_disorder_NV(N,λ,one(b))
    H_dis = H_dis/norm(H_dis)
    H_disq = Hamiltonian_disorder_NV_quad(N,λ,one(b))
    H_disq = H_disq/norm(H_disq)
    H_dip = Hamiltonian_interaction_NV(N,λ,one(b))
    H_dip = H_dip/norm(H_dip)  
            
    errors_DC3 = [(rand()-0.5)*DC3_gensH[i] for i in 1:16]
    for i in 1:16
        errors_DC3[i] = τp_DC3[i]*errors_DC3[i]/norm(errors_DC3[i])
    end
    errors_72 = [sum(λ[j]*(rand()-0.5) for j in 1:8) for i in 1:2]
    for i in 1:2
        errors_72[i] = τp72[i]*errors_72[i]/norm(errors_72[i])
    end

    errors_36 = [sum(λ[j]*(rand()-0.5) for j in 1:8) for i in 1:2]
    for i in 1:2
        errors_36[i] = τp36[i]*errors_36[i]/norm(errors_36[i])
    end
    errors_O = [sum(λ[j]*(rand()-0.5) for j in 1:8) for i in 1:2]
    for i in 1:2
        errors_O[i] = τpO[i]*errors_O[i]/norm(errors_O[i])
    end

    for m in 1:nexps

        H = exps[m][1]*H_dis + exps[m][2]*H_disq + exps[m][3]*H_dip 

        data_DC3[m] += propagate_finite_err_subsampDROID(N,H,DROIDC3PO,DC3_gensH,errors_DC3,exps[m][4],τ0,τpi2*τp_DC3,Id, n_DC3,Ufinal)/nsample

        data_Σ72[m] += propagate_finite_err_subsamp(N,H,Σ72sequence,Σ72x3_gensH,errors_72,exps[m][4],τ0,τpi2*τp72,Id, n_Σ72,Id)/nsample
        data_Σ36[m] += propagate_finite_err_subsamp(N,H,Σ36sequence,Σ36x3_gensH,errors_36,exps[m][4],τ0,τpi2*τp36,Id, n_Σ36,Id)/nsample

        data_Σ72_cyclic[m] += propagate_finite_err_subsampCYCLIC(N,H,Σ72sequence,Σ72x3_gensH,errors_72,exps[m][4],τ0,τpi2*τp72,Id, n_Σ72_cyclic,Id)/nsample
        data_Σ36_cyclic[m] += propagate_finite_err_subsampCYCLIC(N,H,Σ36sequence,Σ36x3_gensH,errors_36,exps[m][4],τ0,τpi2*τp36,Id, n_Σ36_cyclic,Id)/nsample

        data_O[m] += propagate_finite_err_subsamp(N,H,Osequence,O_gensH,errors_O,exps[m][4],τ0,τpi2*τpO,Id, n_O,Id)/nsample
        data_O_cyclic[m] += propagate_finite_err_subsampCYCLIC(N,H,Osequence,O_gensH,errors_O,exps[m][4],τ0,τpi2*τpO,Id, n_O_cyclic,Id)/nsample

    end

    println("$n/$nsample")
end


jldopen(fig_dir*"Comparison SOTA DROID.jld", "w") do file 
    write(file, "data", data_DC3)
end
jldopen(fig_dir*"Comparison SOTA 72.jld", "w") do file 
    write(file, "data", data_Σ72)
end
jldopen(fig_dir*"Comparison SOTA 36.jld", "w") do file 
    write(file, "data", data_Σ36)
end
jldopen(fig_dir*"Comparison SOTA O.jld", "w") do file 
    write(file, "data", data_O)
end
jldopen(fig_dir*"Comparison SOTA 72c.jld", "w") do file 
    write(file, "data", data_Σ72_cyclic,)
end
jldopen(fig_dir*"Comparison SOTA 36c.jld", "w") do file 
    write(file, "data", data_Σ36_cyclic)
end
jldopen(fig_dir*"Comparison SOTA Oc.jld", "w") do file 
    write(file, "data", data_O_cyclic)
end
=#
data_DC3 = load(fig_dir*"Comparison SOTA DROID.jld", "data")
data_Σ72 = load(fig_dir*"Comparison SOTA 72.jld", "data")
data_Σ36 = load(fig_dir*"Comparison SOTA 36.jld", "data")
data_O = load(fig_dir*"Comparison SOTA O.jld", "data")
data_Σ72_cyclic = load(fig_dir*"Comparison SOTA 72c.jld", "data")
data_Σ36_cyclic = load(fig_dir*"Comparison SOTA 36c.jld", "data")
data_O_cyclic = load(fig_dir*"Comparison SOTA Oc.jld", "data")

ms = [ :rtriangle,:star4,:diamond,:dot,:ltriangle,:star5, :rect]
ls = [ :solid, :solid,:solid,:solid, (:dot,:dense),(:dot,:dense), (:dot,:dense)]
colors = [ RGB(153/255 , 204/255 ,255/255) ,RGB(255/255 , 102/255 , 102/255) , RGB(0/255,153/255,0/255),RGB(255/255,178/255,102/255) , RGB(204/255,155/255,255/255), RGB(102/255,102/255,255/255) ,RGB(160/255,160/255,160/255)  ]

name = [L"\Sigma(72\times 3)",L"\Sigma(36\times 3)",L"\Delta(24)",L"\text{DROID}",L"\Sigma(72\times 3)*",L"\Sigma(36\times 3)*",L"\Delta(24)*"]

datas = [data_Σ72,data_Σ36, data_O,data_DC3,data_Σ72_cyclic,data_Σ36_cyclic , data_O_cyclic]
ts = [[t_72*n for n in 1:n_Σ72] , [t_36*n for n in 1:n_Σ36] ,[t_O*n for n in 1:n_O],[t_DC3*n for n in 1:n_DC3] ,  [t_72_cyclic*n for n in 1:n_Σ72_cyclic] , [t_36_cyclic*n for n in 1:n_Σ36_cyclic] , [t_O_cyclic*n for n in 1:n_O_cyclic]]

fig = Figure(;size=(700,400))

for n in 1:4 

    ax = Axis(fig[1,n], xticklabelsvisible=false)

    for p in 1:7
        scatterlines!(ts[p] , datas[p][n] , marker = ms[p] , strokecolor = colors[p] , strokewidth = 2.2 , markercolor = :transparent, color=colors[p], linestyle=ls[p])
    end


    ax = Axis(fig[2,n], xticks=([5,15,25]), xlabel="time [μs]")

    for p in 1:7
        scatterlines!(ts[p] , datas[p][n+4] , marker = ms[p] , strokecolor = colors[p] , strokewidth = 2.0 , markercolor = :transparent, color=colors[p], linestyle=ls[p])
    end
    

end

Legend(fig[3,1:4], [[LineElement(color=colors[p], linestyle=ls[p]),MarkerElement(markersize=15.0, marker = ms[p] , strokecolor = colors[p] , strokewidth = 2.2 , markercolor = :transparent)] for p in 1:7],
 name,patchsize=(35,20),nbanks=4,framevisible=false,tellheight=true)

display(fig)
CairoMakie.activate!()
#save(fig_dir*"DROID comparison.pdf",fig)
