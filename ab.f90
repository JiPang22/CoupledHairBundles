program abc
implicit none
integer j,jmax,n,nmax
real X1,X2,Xa1,Xa2,V1,V2,X1dot,X2dot,Xa1dot,Xa2dot,V1dot,V2dot,Fmax1,Fmax2,gam,u1,u2,z1,z2
real delta,fHB1,fHB2,lambda,lambdaa,kgs,ksp,D,m,dt,Emedia,k,dk,S,A,Po1,Po2,sumt,t
parameter(S=0.65,jmax=100,nmax=20000,dk=0.005,A=exp(10.),Fmax1=48.,Fmax2=51.,gam=1.3,m=1.)
parameter(dt=0.01,kgs=0.75,ksp=0.6,delta=4.53,lambda=2.8,lambdaa=10.,D=60.9)
real noise(nmax)
open(1,file='aa')
open(2,file='ab')

call random_seed()
do n = 1, nmax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
noise(2 * n - 1) = z1
noise(2 * n) = z2
end do

!do j=0,jmax
!k=j*dk
k=0.1
t=0.
V1=1.
V2=1.
X1=-1.
X2=1.
Xa1=1.
Xa2=1.
sumt=0.

do n=1,nmax
t=n*dt

Po1=(1.+A*exp(-(X1-Xa1)/delta))**(-1.)
Po2=(1.+A*exp(-(X2-Xa2)/delta))**(-1.)

fHB1=-lambda*V1-kgs*(X1-Xa1-D*Po1)
fHB2=-lambda*V2-kgs*(X2-Xa2-D*Po2)

V1dot=-gam*X1+(k/m)*(X2-X1)+fHB1/m!+noise(n)
V2dot=-gam*X2+(k/m)*(X1-X2)+fHB2/m!+noise(n)

X1dot=V1
X2dot=V2
Xa1dot=(kgs/lambda)*(X1-Xa1-D*Po1)-(Fmax1/lambdaa)*(1.-S*Po1)!+noise(n)
Xa2dot=(kgs/lambda)*(X2-Xa2-D*Po2)-(Fmax2/lambdaa)*(1.-S*Po2)!+noise(n)

V1=V1+V1dot*dt
V2=V2+V2dot*dt
X1=X1+X1dot*dt
X2=X2+X2dot*dt
Xa1=Xa1+Xa1dot*dt
Xa2=Xa2+Xa2dot*dt

sumt=sumt+(X1-X2)
!write(1,*) t,(1./2.)*m*V1**2+(1./2.)*ksp*X1**2
!write(2,*) t,(1./2.)*m*V2**2+(1./2.)*ksp*X2**2
write(1,*) t,X1
write(2,*) t,X2
enddo

!Emedia=(1./2.)*k*(((1./(dt*nmax))*sumt)**2)
!write(1,*) k,Emedia
!enddo
end
