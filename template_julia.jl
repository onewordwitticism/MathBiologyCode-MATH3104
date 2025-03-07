using Plots
using LinearAlgebra

const a, b = 0.0, 10.0
const Dcoef=0.1
const N=100
const dt, dx=0.1, (b-a)/N

const xval=a:dx:(b-dx)
const M= diff(diff([[zeros(1, N-1) 1.0]; I(N) ; [1.0 zeros(1, N - 1)]], dims = 1), dims = 1) / dx^2

function main()
  ρ=zeros(N); ρ[1+N ÷ 2]=100.0/dx
  t=0.0
  while t<10.0
    t=t+dt
    ρold=ρ
    ρ=(I(N)-Dcoef*dt*M)\ρold

    display(plot(xval, ρ, ylim=(0.0, 100.0), show=true))
    sleep(0.1)
  end

end

main()