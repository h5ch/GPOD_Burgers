using DelimitedFiles
using CSV
using CairoMakie


#-----------------------------------------------------------------------------#
# fig1
#-----------------------------------------------------------------------------#
ns = 11
tf = 1.0

solution = readdlm("solution_p.txt")
xp = solution[:,1]
up = solution[:,2:ns+1]

solution = readdlm("solution_a.txt")
xa = solution[:,1]
ua = solution[:,2:ns+1]

fig = Figure(resolution=(500,400),fonts = (; regular = "Times New Roman"),fontsize=16)
kwargs = (;xminorgridvisible = false, xgridvisible = false, yminorgridvisible = false, ygridvisible = false, xtickalign=1, ytickalign=1)
ax = Axis(fig[1, 1],xlabel="x",ylabel="u";kwargs...)
lines!(ax, xp, up[:,1], label="Numerical", color=:blue)
lines!(ax, xa, ua[:,1], label="Analytical", color=:red, linestyle=:dash)
for i = 2:ns
    lines!(ax, xp, up[:,i], color=:blue)
end
for i = 2:ns
    lines!(ax, xa, ua[:,i], color=:red, linestyle=:dash)
end
axislegend()
fig
save("Burgers_solution.pdf", fig; pt_per_unit=2)

#-----------------------------------------------------------------------------#
# fig2
#-----------------------------------------------------------------------------#
fig2 = Figure(resolution=(500,400),fonts = (; regular = "Times New Roman"),fontsize=16)
ax = Axis(fig2[1,1],xlabel="x",ylabel="u";kwargs...)
lines!(ax,x,up[:,5],label="Numerical",color =:blue,linestyle = :dot)
lines!(ax,x,ua[:,5],label="Analytical",color =:red)
lines!(ax,x,u_pred[1:Nx_a],label="Reconstruction",color =:green,linestyle = :dash)
axislegend()
fig2
save("Reconstruction.pdf", fig2; pt_per_unit=2)

#-----------------------------------------------------------------------------#
# fig3
#-----------------------------------------------------------------------------#
fig3 = Figure(resolution=(500,400),fonts = (; regular = "Times New Roman"),fontsize=16)
kwargs = (;xminorgridvisible = false, xgridvisible = false, yminorgridvisible = false, ygridvisible = false, xtickalign=1, ytickalign=1)
x_mode = 11
ax = Axis(fig3[1,1],xlabel="k";kwargs...)
scatterlines!(ax,1:x_mode,c1,  color=:black,label="c1")
scatterlines!(ax,1:x_mode,c2,  color=:green,label="c2",marker=:diamond)
scatterlines!(ax,1:x_mode,c3,  color=:purple,label="c3",linestyle=:dash)
scatterlines!(ax,1:x_mode,c4,  color=:blue,label="c1+c2+c3",linestyle=:dot,marker=:diamond)
scatter!(     ax,1:x_mode,rmse,color=:red,label="RMSE",marker=:cross)
axislegend(;position=:rt,framevisible = false)
fig3
save("Error.pdf", fig3; pt_per_unit=2)

#-----------------------------------------------------------------------------#
# fig4
#-----------------------------------------------------------------------------#
fig4 = Figure(resolution=(1500, 400),fonts = (; regular = "Times New Roman"),fontsize=16)
kwargs = (;xminorgridvisible = false, xgridvisible = false, yminorgridvisible = false, ygridvisible = false, xtickalign=1, ytickalign=1)
ax1 = Axis(fig4[1,1],xlabel="k";kwargs...)
ax2 = Axis(fig4[1,2],xlabel="k";kwargs...)
ax3 = Axis(fig4[1,3],xlabel="k";kwargs...)
lines!(ax1,x,ϕ_a[:,1],label="HF Mode 1",color=:red)
lines!(ax2,x,ϕ_a[:,2],label="HF Mode 2",color=:red)
lines!(ax3,x,ϕ_a[:,3],label="HF Mode 3",color=:red)
lines!(ax1,x,ϕ_p[:,1],label="LF Mode 1",color=:blue,linestyle = :dash)
lines!(ax2,x,ϕ_p[:,2],label="LF Mode 2",color=:blue,linestyle = :dash)
lines!(ax3,x,ϕ_p[:,3],label="LF Mode 3",color=:blue,linestyle = :dash)
axislegend(ax1,position = :lt)
axislegend(ax2,position = :lt)
axislegend(ax3,position = :lt)
fig4
save("Modes.pdf", fig4; pt_per_unit=2)