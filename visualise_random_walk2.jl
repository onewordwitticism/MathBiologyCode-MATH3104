using GLMakie
import GLMakie: scatter!, lines!
using Printf
#using Compat

const Ntraj=150
const N=1000
const T=15.0
const dt=T/N
const D=0.1
const L=6.0

#############################

fig = Figure()

axtop = Axis(fig[1:4, 1], limits = ((-L, L), (-L, L))) #, title = "position of random walkers")

sl_t    = Slider(fig[5, 1], range = range(0, N, step=1), startvalue=0)

axbottom = Axis(fig[6, 1]) #, title = "position of random walkers")

alltraj=Vector{Vector{Vector{Float64}}}(undef, Ntraj)
msd=Vector{Float64}(undef, N)
for i=1:Ntraj
   trajx=cumsum(sqrt(2*D*dt)*randn(N))
   trajy=cumsum(sqrt(2*D*dt)*randn(N))
   traj=[[trajx[j], trajy[j]] for j=1:N]
   alltraj[i]=traj
end
for k=1:N
   msd[k]=sum(alltraj[j][k][1]^2+alltraj[j][k][2]^2 for j=1:Ntraj)/Ntraj
end

graphs = lift(sl_t.value) do tind
   walkerx=zeros(Ntraj)
   walkery=zeros(Ntraj)
   if tind>0
      for i=1:Ntraj
         walkerx[i]=alltraj[i][tind][1]
         walkery[i]=alltraj[i][tind][2]
      end
   end
   (x=walkerx, y=walkery, myt=tind*dt, mymsd=sum(walkerx.^2+walkery.^2)/Ntraj)
end


walkerposx=@lift($graphs.x)
walkerposy=@lift($graphs.y)
mymsd=@lift($graphs.mymsd)
myt=@lift($graphs.myt)

col=rand(1:10, Ntraj)
scatter!(axtop, walkerposx, walkerposy, markersize = 8,  color = col, colormap = :tab10, colorrange = (1, 10))

lines!(axbottom, collect(0:N)*dt, [0.0;msd])
scatter!(axbottom, myt, mymsd,  markersize = 8, color = :red )
display(fig)

#mysol=getsol(0.02, 1.7)
#testing a random change for git puurposes