using DelimitedFiles
using CSV
using Statistics
using LinearAlgebra
using IterativeSolvers

function POD_SVD(Unm)
    Nx = size(Unm,1)
    Nt = size(Unm,2)
    Ux0 = mean(Unm,dims=1)
    Ut0 = mean(Unm,dims=2)
    Unm = (Unm .- Ut0)
    F = svd(Unm,full=true)

    ϕ = F.U #N_mode=Nx
    λ = @. F.S^2 /Nt
    a = Unm'*ϕ

    return Ut0,ϕ,a,λ
end

function gappy_reconstruct(Nx1,Nx2,Ut0_mf,ϕ_mf,k,Ux)    
    mask = zeros(Nx1+Nx2)
    mask[Nx1+1:end] .= 1

    ϕ = similar(ϕ_mf)
    for i in eachindex(mask)
        ϕ[i,:] = mask[i] .* ϕ_mf[i,:]
    end

    gappy = zeros(Nx1+Nx2)
    gappy[Nx1+1:end] = Ux - Ut0_mf[Nx1+1:end]
    
    ϕ = ϕ_mf[Nx1+1:end,1:k]
    K = ϕ'*ϕ
    f = ϕ'*gappy[Nx1+1:end]
    ξ = gmres(K,f)

    # reconstruction with i known basis function
    Ux_pred = zeros(Nx1+Nx2)
    for i=1:k, j in eachindex(Ux_pred)
        Ux_pred[j] += ξ[i]*ϕ_mf[j,i]
    end
    Ux_pred = Ux_pred + Ut0_mf

    ϕ_bar = ϕ_mf[1:Nx1,1:k]
    ϕ_til = ϕ_mf[Nx1+1:end,1:k]

    ϕ_inv = pinv(ϕ_til,rtol = sqrt(eps(real(float(one(eltype(ϕ_til)))))))
    return Ux_pred,ϕ_bar,ϕ_til,ϕ_inv,ξ
end

function evaluate_error(U_real,U_pred)
    e = similar(U_real)
    for i in eachindex(e)
        e[i] = (U_real[i] - U_pred[i])^2
    end
    error = sum(e)
    error1 = error/sum(U_real.^2)
    rmse = sqrt(error/length(U_real))
    return error1,e,rmse
end

function error_of_modes(ua,up,u_test,u_true)
    snapshot = vcat(ua,up)
    Nx_a = size(ua,1)
    Nx_p = size(up,1)

    Ut0_mf,ϕ_mf,a_mf,λ_mf = POD_SVD(snapshot)
    n=length(λ_mf)
    mode_error = zeros(n)
    rmse = zeros(n)
    λ_min = zeros(n)
    c1 = zeros(n)
    c2 = zeros(n)
    c3 = zeros(n)
    c = zeros(n)
    c4 = zeros(n)
    for k in eachindex(λ_mf)
        Ux_pred,ϕ_bar,ϕ_til,ϕ_inv,ξ = gappy_reconstruct(Nx_a,Nx_p,Ut0_mf,ϕ_mf,k,u_test)
        an_bar = (Ux_pred[1:Nx_a]-Ut0_mf[1:Nx_a])'*ϕ_mf[1:Nx_a,:]
        an_til = (Ux_pred[Nx_a+1:end]-Ut0_mf[Nx_a+1:end])'*ϕ_mf[Nx_a+1:end,:]
        an = an_bar+an_til
        
        mode_error[k],e_point,rmse[k] = evaluate_error(u_true,Ux_pred)

        e1 = (ϕ_bar*(I(k)-ϕ_inv*ϕ_til)*an[1:k]) 
        c1[k] = sqrt(sum(e1 .^2)/Nx_a)
        if k==n
            e2=zeros(Nx_a)
            e3=zeros(Nx_a)
            c2[k]=0.0
            c3[k]=0.0
        else
            Ur_bar = ϕ_mf[1:Nx_a,k+1:n]*a_mf[5,k+1:n]
            Ur_til = ϕ_mf[Nx_a+1:end,k+1:n]*a_mf[5,k+1:n]
            
            e2 = (-ϕ_bar*ϕ_inv*Ur_til) 
            e3 = (Ur_bar)
  
            c2[k] = sqrt(sum(e2 .^2)/Nx_a)
            c3[k] = sqrt(sum(e3 .^2)/Nx_a)
        end

        c[k] = c1[k]+c2[k]+c3[k]
        c4[k] = sqrt(sum((e1+e2+e3) .^2)/Nx_a)
        F = svd(ϕ_til)
        λ_min[k] = minimum(F.S)
    end

    return mode_error,λ_min,c1,c2,c3,c,rmse,c4
end
#--------------------------------------
# main program
#--------------------------------------
ns = 11
tf = 1.0

solution = readdlm("solution_p.txt")
up = solution[:,2:ns+1]

solution = readdlm("solution_a.txt")
ua = solution[:,2:ns+1]

Ut0_p,ϕ_p,an_p,Ds_p = POD_SVD(up)
Ut0_a,ϕ_a,an_a,Ds_a = POD_SVD(ua)

u_test = up[:,5]
u_true = ua[:,5]

k = 10
snapshot = vcat(ua,up)
Ut0_mf,ϕ_mf,a_mf,λ_mf = POD_SVD(snapshot)
Nx_a = size(ua,1)
Nx_p = size(up,1)
u_pred,ϕ_bar,ϕ_til,ϕ_inv,ξ = gappy_reconstruct(Nx_a,Nx_p,Ut0_mf,ϕ_mf,k,u_test)
mode_error,λ_min,c1,c2,c3,c,rmse,c4= error_of_modes(ua,up,u_test,u_true)
