program molecular_dynamics_phase_transition_three_dim
implicit none
include "fftw3.f"
integer ( kind = 4 ), parameter :: N = 1000
real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0
real ( kind = 8 ) pi,a,M,epsilon,wpd,L,Lx,Ly,Lz,Vxmax,Vymax,Vzmax,svx,svy,svz,nbar,r
real ( kind = 8 ) Gamma,Temp,K,fx,fy,fz,KE,PE,TE,B,alpha,epsa
real ( kind = 8 ) t,tmax,dt,xdiff,ydiff,zdiff,scale,tau,sumvx,sumvy,sumvz,virial,f,Pressure
integer ( kind = 8 ) i,j,k,N,p
real ( kind = 8 ) x(N),y(N),z(N),vx(N),vy(N),vz(N),ax(N),ay(N),az(N),ux(N),uy(N),uz(N),avx(N),avy(N),avz(N),v0x(N),v0y(N),v0z(N)
integer,parameter :: seed = 99999999
!pi = 4.0*atan(1.0)
call srand(seed)

!====================== User Inputs =============================

!Coupling Parameter...
Gamma = 250.0
Temp = 1.0/Gamma

!Screening Parameter...
K = 1.0

!Mass of the Dust Particles...
M = 1.0

!Free Space Permittivity...
epsilon = 1.0 !8.8541878*10.0**(-12)

!Total Number of Particles...
!N = 1000

!Areal Density of Particles...
!nbar = float(N)/(2.0*Lx*2.0*Ly)
nbar = 1.0/pi !3.0/(4.0*pi)

!Interparticle Distance...
a = 1.0!sqrt(1.0/(pi*nbar))

!Dust Plasma Frequency...
wpd = sqrt(nbar/(M*a*epsilon))   ! 2-D system -- 'a' is present in the denominator...

!Normalized System Size...
!L = (float(N)/(8.0*nbar))**(1.0/3.0)
!Lx = 7.105!L/a
!Ly = 7.105!L/a
!Lz = 7.105!L/a
L = sqrt(float(N)/(4.0*nbar))
Lx = L/a
Ly = L/a

!Initial Normalized Maximum Velocity of the Particles
Vxmax = 1.0*sqrt(2.0)/a
Vymax = 1.0*sqrt(2.0)/a
!Vzmax = 1.0*sqrt(2.0)/a

!Normalized Final Time and Time Steps...
tmax = 2000.0!/sqrt(2.0)
dt = 0.01!/sqrt(2.0)

svy = 0.0
!svz = 0.0
virial = 0.0

dx = 0.01
dy = 0.01
x0 = -Lx
y0 = -Ly
Nx = int(2.0*Lx/dx)
Ny = int(2.0*Ly/dy)

!====================== Output Filenames ===========================

open(unit=1,file='Initial_Configuration.dat',status='unknown')
open(unit=2,file='Average_Velocity.dat',status='unknown')
open(unit=10,file='Energy.dat',status='unknown')

!====================== Definition of initial state =================

!Definition of the initial random positions and velocities in -Lx to Lx and -Ly to Ly rectangular box... 
do i = 1, N, 1
x(i) = (rand())*2.0*Lx - Lx 
y(i) = (rand())*2.0*Ly - Ly 
!z(i) = (rand())*2.0*Lz - Lz
vx(i) = (rand())*Vxmax - Vxmax/2.0
svx = svx + vx(i) 
vy(i) = (rand())*Vymax - Vymax/2.0
svy = svy + vy(i)                                               ! Center of Mass has Zero y Velocity...
!vz(i) = (rand())*Vzmax - Vzmax/2.0
!svz = svz + vz(i)                                               ! Center of Mass has Zero z Velocity...
enddo 

!Definitions of corrected initial velocities...
do i = 1,N,1
vx(i) = vx(i) - svx/float(N)
vy(i) = vy(i) - svy/float(N)
!vz(i) = vz(i) - svz/float(N)
enddo

!Calculating the initial accelerations of the particles...
do xg = 1,Nx,1
do yg = 1,Ny,1
rho(xg,yg) = 0.0
enddo
enddo

do i = 1,N,1
xlc(i) = (x(i)-x0)/dx + 1.0
ylc(i) = (y(i)-y0)/dy + 1.0
xg = int(xlc(i))
yg = int(ylc(i))
dxg = xlc(i) - float(xg)
dyg = ylc(i) - float(yg)
rho(xg,yg) = rho(xg,yg) + (1.0 - dxg) * (1.0 - dyg)
rho(xg+1,yg) = rho(xg+1,yg) + dxg * (1.0 - dyg)
rho(xg,yg+1) = rho(xg,yg+1) + (1.0 - dxg) * dyg
rho(xg+1,yg+1) = rho(xg+1,yg+1) + dxg * dyg
enddo

do xg = 1,Nx,1
do yg = 1,Ny,1
rho_dum(xg,yg) = rho(xg,yg)
enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

do xg = 1,Nx/2+1
  do yg = 1,Ny/2
  kx = 2.0d0*pi*float(xg-1)/Lx
  ky = 2.0d0*pi*float(yg-1)/Ly
    if (xg == 1 .and. yg == 1) then
    phi_k(xg,yg) = (0.0d0,0.0d0)
    else
    phi_k(xg,yg) = rhok(xg,yg)/( kx*kx + ky*ky )
    endif
  phi_k_dum(xg,yg) = phi_k(xg,yg) 
  enddo
  do yg = Ny/2+1,Ny
  kx = 2.0d0*pi*float(xg-1)/Lx
  ky = 2.0d0*pi*float((yg-1)-Ny)/Ly
  phi_k(xg,yg) = rhok(xg,yg)/( kx*kx + ky*ky )
  phi_k_dum(xg,yg) = phi_k(xg,yg) 
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, phi_k_dum, phi, FFTW_ESTxgMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

do xg = 1, Nx
  do yg = 1, Ny
  phi(xg,yg) = phi(xg,yg)/(float(Nx)*float(Ny))
  !write(10,*) x(xg),y(yg),rho(xg,yg),phi(xg,yg)
  fx(xg,yg) = ( phi(xg+1,yg) - phi(xg-1,yg) ) / (2.0*dx)
  fy(xg,yg) = ( phi(xg,yg+1) - phi(xg,yg-1) ) / (2.0*dy)
  end do
end do

do i = 1,N,1
xlc(i) = (x(i)-x0)/dx + 1.0
ylc(i) = (y(i)-y0)/dy + 1.0
xg = int(xlc(i))
yg = int(ylc(i))
dxg = xlc(i) - float(xg)
dyg = ylc(i) - float(yg)
ax(i) = fx(xg,yg)*(1.0d0-dxg)*(1.0d0-dyg) + fx(xg+1,yg)*dxg*(1.0d0-dyg) + fx(xg,yg+1)*(1.0d0-dxg)*dyg + fx(xg+1,yg+1)*dxg*dyg
ay(i) = fy(xg,yg)*(1.0d0-dxg)*(1.0d0-dyg) + fy(xg+1,yg)*dxg*(1.0d0-dyg) + fy(xg,yg+1)*(1.0d0-dxg)*dyg + fy(xg+1,yg+1)*dxg*dyg
enddo


do i = 1,N,1
!write(1,*) 0,x(i),y(i),z(i),vx(i),vy(i),vz(i),ax(i),ay(i),az(i)
write(1,*) 0,x(i),y(i),vx(i),vy(i),ax(i),ay(i)
enddo

!====== MD Time Evolution using velocity verlet method started =========

do t = 0.010,tmax,dt
sumvx = 0.0
sumvy = 0.0
!sumvz = 0.0
virial = 0.0
Pressure = 0.0
PE = 0.0
KE = 0.0

!Calculating velocity after time dt/2.0 and position after time dt... 
do i = 1,N,1
ux(i) = vx(i) + dt*ax(i)/2.0
x(i) = x(i) + dt*ux(i)
!ux(i) = ux(i) - abs(int(x(i)/Lx))*ux(i)*2.0                  ! Perfectly Reflecting Boundary
!x(i) = x(i) - (int(x(i)/Lx))*(abs(x(i))-Lx)*2.0              ! Perfectly Reflecting Boundary
x(i) = x(i) - (int(x(i)/Lx))*2.0*Lx                         ! Periodic Boundary Condition
uy(i) = vy(i) + dt*ay(i)/2.0
y(i) = y(i) + dt*uy(i)
!uy(i) = uy(i) - abs(int(y(i)/Ly))*uy(i)*2.0                  ! Perfectly Reflecting Boundary
!y(i) = y(i) - (int(y(i)/Ly))*(abs(y(i))-Ly)*2.0              ! Perfectly Reflecting Boundary
y(i) = y(i) - (int(y(i)/Ly))*2.0*Ly                         ! Periodic Boundary Condition
!uz(i) = vz(i) + dt*az(i)/2.0
!z(i) = z(i) + dt*uz(i)
!z(i) = z(i) - (int(z(i)/Lz))*2.0*Lz                         ! Periodic Boundary Condition
enddo ! i

!Calculating the acceleration after time dt...

do xg = 1,Nx,1
do yg = 1,Ny,1
rho(xg,yg) = 0.0
enddo
enddo

do i = 1,N,1
xlc(i) = (x(i)-x0)/dx + 1.0
ylc(i) = (y(i)-y0)/dy + 1.0
xg = int(xlc(i))
yg = int(ylc(i))
dxg = xlc(i) - float(xg)
dyg = ylc(i) - float(yg)
rho(xg,yg) = rho(xg,yg) + (1.0 - dxg) * (1.0 - dyg)
rho(xg+1,yg) = rho(xg+1,yg) + dxg * (1.0 - dyg)
rho(xg,yg+1) = rho(xg,yg+1) + (1.0 - dxg) * dyg
rho(xg+1,yg+1) = rho(xg+1,yg+1) + dxg * dyg
enddo

do xg = 1,Nx,1
do yg = 1,Ny,1
rho_dum(xg,yg) = rho(xg,yg)
enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

do xg = 1,Nx/2+1
  do yg = 1,Ny/2
  kx = 2.0d0*pi*float(xg-1)/Lx
  ky = 2.0d0*pi*float(yg-1)/Ly
    if (xg == 1 .and. yg == 1) then
    phi_k(xg,yg) = (0.0d0,0.0d0)
    else
    phi_k(xg,yg) = rhok(xg,yg)/( kx*kx + ky*ky )
    endif
  phi_k_dum(xg,yg) = phi_k(xg,yg) 
  enddo
  do yg = Ny/2+1,Ny
  kx = 2.0d0*pi*float(xg-1)/Lx
  ky = 2.0d0*pi*float((yg-1)-Ny)/Ly
  phi_k(xg,yg) = rhok(xg,yg)/( kx*kx + ky*ky )
  phi_k_dum(xg,yg) = phi_k(xg,yg) 
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, phi_k_dum, phi, FFTW_ESTxgMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

do xg = 1, Nx
  do yg = 1, Ny
  phi(xg,yg) = phi(xg,yg)/(float(Nx)*float(Ny))
  !write(10,*) x(xg),y(yg),rho(xg,yg),phi(xg,yg)
  fx(xg,yg) = ( phi(xg+1,yg) - phi(xg-1,yg) ) / (2.0*dx)
  fy(xg,yg) = ( phi(xg,yg+1) - phi(xg,yg-1) ) / (2.0*dy)
  end do
end do

do i = 1,N,1
xlc(i) = (x(i)-x0)/dx + 1.0
ylc(i) = (y(i)-y0)/dy + 1.0
xg = int(xlc(i))
yg = int(ylc(i))
dxg = xlc(i) - float(xg)
dyg = ylc(i) - float(yg)
ax(i) = fx(xg,yg)*(1.0d0-dxg)*(1.0d0-dyg) + fx(xg+1,yg)*dxg*(1.0d0-dyg) + fx(xg,yg+1)*(1.0d0-dxg)*dyg + fx(xg+1,yg+1)*dxg*dyg
ay(i) = fy(xg,yg)*(1.0d0-dxg)*(1.0d0-dyg) + fy(xg+1,yg)*dxg*(1.0d0-dyg) + fy(xg,yg+1)*(1.0d0-dxg)*dyg + fy(xg+1,yg+1)*dxg*dyg
enddo


!do i = 1,N
!ax(i) = 0.0
!ay(i) = 0.0
!!az(i) = 0.0
!do j = 1,N
!if (i .ne. j) then                                          ! .ne. is used, PE should be halved...
!xdiff = (x(i)-x(j)) - nint((x(i)-x(j))/(2.0*Lx))*2.0*Lx     ! Minimum Image Convension Introduced...
!ydiff = (y(i)-y(j)) - nint((y(i)-y(j))/(2.0*Ly))*2.0*Ly     ! Minimum Image Convension Introduced...
!!zdiff = (z(i)-z(j)) - nint((z(i)-z(j))/(2.0*Lz))*2.0*Lz     ! Minimum Image Convension Introduced...
!!r = sqrt((xdiff)**2 + (ydiff)**2 + (zdiff)**2)
!r = sqrt((xdiff)**2 + (ydiff)**2)
!f = (1.0+k*r)*exp(-k*r)/r**3
!fx = (xdiff)*(1.0+k*r)*exp(-k*r)/r**3 
!ax(i) = ax(i) + fx
!fy = (ydiff)*(1.0+k*r)*exp(-k*r)/r**3 
!ay(i) = ay(i) + fy
!!fz = (zdiff)*(1.0+k*r)*exp(-k*r)/r**3 
!!az(i) = az(i) + fz
!virial = virial + r**2 * f
!PE = PE + (exp(-k*r)/(2.0*r))                               ! Calculation of Potential Energy...
!endif
!enddo ! j
!enddo ! i

!Calculating the velocity after time dt...
do i = 1,N,1
sumvx = sumvx + vx(i)                                       ! Check for average x-velocity...
sumvy = sumvy + vy(i)                                       ! Check for average y-velocity...
!sumvz = sumvz + vz(i)                                       ! Check for average y-velocity...
vx(i) = ux(i) + ax(i)*dt/2.0
vy(i) = uy(i) + ay(i)*dt/2.0
!vz(i) = uz(i) + az(i)*dt/2.0
!KE = KE + (vx(i)**2 + vy(i)**2 + vz(i)**2)/2.0              ! Calculation of Kinetic Energy...
KE = KE + (vx(i)**2 + vy(i)**2)/2.0 
enddo ! i

!Writing the instantaneous position, velocity, acceleration in files...
do i = 1,N,1
p = int(t/dt)
if (p  .ge. tmax/(4.0*dt) .and. mod(float(p),1.0/dt) == 0.0) then
!write(p,*) t,x(i),y(i),z(i),vx(i),vy(i),vz(i),ax(i),ay(i),az(i)
write(p,*) t,x(i),y(i),vx(i),vy(i),ax(i),ay(i)
endif
enddo

!Total Energy and Pressure...
!Pressure = (2.0*KE + virial/2.0)/(3.0*2.0*Lx*2.0*Ly*2.0*Lz) 
!Pressure = (2.0*KE + virial/2.0)/(2.0*Lx*2.0*Ly) ! CALCULATE THE PRESSURE OF INDIVIDUAL BOX FOR PERFECT 2D SYSTEM...
TE = KE + PE
write(10,*) t,TE/float(N),2.0*KE/float(N),PE/float(N)!,Pressure


!====================  Thermostat =======================

tau = 10.0*dt
if (t .le. tmax/2.0) then
do i = 1,N,1
scale = sqrt(1.0 + (dt/tau)*((Temp/(KE/(float(N)/2.0) ))-1.0))
!scale = sqrt(1.0 + (dt/tau)*((Temp/(2.0*KE/(3.0*float(N)) ))-1.0))
vx(i) = scale*vx(i)
vy(i) = scale*vy(i)
!vz(i) = scale*vz(i)
enddo
else
do i = 1,N,1
vx(i) = vx(i)
vy(i) = vy(i)
!vz(i) = vz(i)
enddo
endif

enddo ! t

 close(1)
 close(10)

end
