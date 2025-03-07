a=0; b=10; 	Dcoef=0.1;
N=100; deltat=0.1;
xval=linspace(a, b, N);
deltax=xval(2)-xval(1)
M=diff(diff([[zeros(1, N-1), 1]; eye(N); [1, zeros(1, N-1)]]))/deltax^2;
h=figure;
	
rho=zeros(N, 1); rho(N/2)=100/deltax;
t=0;
while t<20
  t=t+deltat;
  rhoold=rho;
  rho=(eye(N)-Dcoef*deltat*M)\rhoold;   
  figure(h);
  plot(xval, rho);
  yl=ylim;
  ylim([0, 100]);
  drawnow
end

