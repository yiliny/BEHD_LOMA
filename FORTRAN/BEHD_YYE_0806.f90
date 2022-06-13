!!!!! Yilin YE @ Jun. 2022
!!!!! EHD + Brownian motion
!!!!! Based on 1412.0162
!!!!! ONLY ONE VARIABLE «z»

!!!!! Here WAS a file based on the technique shown in the article 
!!!!! J. Phys. Chem. B 2014, 118, 6466−6474
!!!!! See ""Support Information"" for the Multi-dimensional case

!!!!! Here is a file based on the Euler–Maruyama method.
!!!!! without different gamma values for frictions and noises.


Program main
    use MATHS
    implicit none    
    integer :: Nmax=100000 !! Define the maximum number of steps
    !real*8,parameter :: pi=3.14159265358979d0 !! constant
    !complex*16,parameter :: Im=(0.d0,1.0d0)    !! imaginary unit
    integer i,j,k,l !! loop variables
    real*8 :: time_begin,time_end !! define variables to record time consumed
    real*8,external :: rectgauss    
    
    include "para.h"    
    !real*8 :: p12,v14,v24,v34 !! old intermediate variables.
    !real*8 :: coefa(6),coefb(6),coefc(6) !! coefficients for equations of motion
    real*8 :: masse(3),amplitude(3,3),fext(3)
    !real*8 :: spuriousforce !! spuriousforce, & 3 components z/x/Θ
    real*8,allocatable :: position(:,:),velocity(:,:),force(:) !! define arrays for position/velocity/force
    allocate(position(3,Nmax)); allocate(velocity(3,Nmax)); allocate(force(Nmax))
    
    call cpu_time(time_begin)
    open(unit=31,file="random_number_test.txt")
    open(unit=32,file="positions.txt")
    open(unit=33,file="velocitys.txt")
    open(unit=34,file="forces.txt")
    !open(unit=35,file="inter_Euler_Maruyama.txt")
    open(unit=36,file="Mass_Matrix_Inverse.txt")
    open(unit=371,file="gamma0.txt")
    open(unit=372,file="gamma1.txt")
    open(unit=373,file="gamma1v.txt")

    61 format(2x,'Time Step',15x,'Noise z',15x,'Noise x',15x,'Noise Θ')
    62 format(2x,'Time Step',12x,'Position z/∆',12x,'Position x',12x,'Angle Θ')
    63 format(2x,'Time Step',12x,'Veloctiy z',15x,'Velocity x',15x,'Velocity Θ')
    64 format(2x,'Time Step',12x,'Force z',12x,'Noise z',12x,'Noise x',12x,'Noise Θ')
    !65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
    66 format(2x,'Time Step',12x,'M_zz',12x,'M_xx',12x,'M_xΘ',12x,'M_Θx',12x,'M_ΘΘ')
    67 format(2x,'Time Step',10x,'γv_zz',10x,'γv_xx',10x,'γv_ΘΘ',10x,'γv_zx',10x,'γv_zΘ',10x,'γv_xΘ')
    write(31,61); write(32,62); write(33,63); write(34,64); !write(35,65)
    write(36,66); write(371,67); write(372,67); write(373,67)
    
    71 format(f10.1,3(5X,ES18.8))
    72 format(f10.1,3(5X,ES18.8)) !! format for file(32)=position_z.txt
    73 format(f10.1,3(5X,ES20.6)) !! format for file(33)=velocity_z.txt
    74 format(f10.1,4(5X,ES15.6)) !! format for file(34)=force_z.txt
    !75 format(f10.1,4(5X,ES20.10)) !! format for file(35)=intermediates.txt
    76 format(f10.1,5(5X,ES15.6))
    77 format(f8.1,6(5X,ES12.4))
    
    99 format('It spent',f8.3,' seconds for the whole program.')
    

    rayon = 1.5d-6; rhosty = 1.06e3; rhosol = 1.00e3; grav = 9.80665
    mass = pi*rayon**2*rhosty; mz = mass; mx = mass; mt = mass*rayon**2/2; masse(1)=mz; masse(2)=mx; masse(3)=mt
    temp = 298.0d0; k_B = 1.38064852d-23; beta = 1.d0/(k_B*temp)
    clight = sqrt(2.d0*grav*rayon*(1.d0-rhosol/rhosty))

    kappa = 1.0d-4; eps = 0.1d0; xi = 0.1d0; kxi = kappa*xi; kxe = kappa*xi*eps
    !coefa(1)=xi;  coefa(2)=21.d0*kxi/4.d0; coefa(3)=-kxi/4.d0; coefa(4)=kxi/2.d0; coefa(5)=-15.d0*kxi/8.d0; coefa(6)=1.d0
    !coefb(1)=2.d0/3.d0*eps*xi;  coefb(2)=19.d0*kxe/24.d0; coefb(3)=-kxe/6.d0; coefb(4)=kxe/12.d0; coefb(5)=-kxe/12.d0; coefb(6)=0.d0
    !coefc(1)=4.d0/3.d0*eps*xi;  coefc(2)=19.d0*kxe/12.d0; coefc(3)=-kxe/3.d0; coefc(4)=kxe/6.0d0; coefc(5)=-kxe/6.0d0; coefc(6)=0.d0
    dt = 2.0d-3 !! time gap for each step

    !! Initiation 
    position=0.0d0; velocity=0.0d0; force=0.0d0; amplitude=0.d0
    gmaeff=0.0d0; Minv=0.0d0; fext=0.d0
    
    position(1,1)=1.0d0; position(2,1)=0.d0; position(3,1)=0.d0
    !p12=0.0d0; v14=0.0d0; v24=0.0d0; v34=0.0d0

    do i=1,Nmax-1
        noise = normaldist(0.0d0,1.0d0,3)
        !do j=1,3
        !    noise(j,1) = noise(j,1)*rectGauss(noise(j,1),2.d0)
        !end do
        write(31,71) 1.d0*i,noise(1,1),noise(2,1),noise(3,1)

        do j=1,3
            do k=1,3
                do l=1,3
                    call updategamma(gmaeff(j,k,l),j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                end do
                call updateMinverse(Minv(j,k),position(1,i),j,k)
            end do
            amplitude(j,j) = sqrt(2.d0*(gmaeff(1,j,j) - gmaeff(2,j,j))/(beta*masse(j)))
        end do


        call updateforce(fext(1),position(1,i))
        velocity(:,i+1) = velocity(:,i) + &
            &fext(:) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),velocity(:,i)) + MATMUL(amplitude(:,:),noise(:,1))
        position(:,i+1) = position(:,i) + velocity(:,i) * dt
        !position(1,i+1) = position(1,1)

        !call updategamma(gamma,position(i),velocity(i))
        !call updateintm
        !call updateMinverse(Minv,position(i))
        !v14 = sqrt(intma)*velocity(i) + sqrt((1.0d0-intma)*Minv/beta)*noise(1,1)
        
        write(32,72) 1.0d0*i,position(1,i),position(2,i),position(3,i)
        write(33,73) 1.0d0*i,velocity(1,i),velocity(2,i),velocity(3,i)
        write(34,74) 1.0d0*i,force(i), amplitude(1,1)*noise(1,1), amplitude(2,2)*noise(2,1), amplitude(3,3)*noise(3,1)
        !write(35,75) 1.0d0*i,gamma,intma,intmb,Minv
        write(36,76) 1.0d0*i,Minv(1,1),Minv(2,2),Minv(1,2),Minv(2,1),Minv(2,2)
        write(371,77) 1.0d0*i,gmaeff(1,1,1),gmaeff(1,2,2),gmaeff(1,3,3),gmaeff(1,1,2),gmaeff(1,1,3),gmaeff(1,2,3)
        write(372,77) 1.0d0*i,gmaeff(2,1,1),gmaeff(2,2,2),gmaeff(2,3,3),gmaeff(2,1,2),gmaeff(2,1,3),gmaeff(2,2,3)
        write(373,77) 1.0d0*i,gmaeff(3,1,1),gmaeff(3,2,2),gmaeff(3,3,3),gmaeff(3,1,2),gmaeff(3,1,3),gmaeff(3,2,3)

    end do


    deallocate(velocity); deallocate(position); deallocate(force)
    close(31); close(32); close(33); close(34); close(35); close(36); close(371); close(372); close(373)
    call cpu_time(time_end)
    write(*,99) time_end-time_begin

end Program main


!*** Here is the module to furnish NORMAL DISTRIBUTION, namely the white noises.
!!! Reference: https://www.cxyzjd.com/article/za36mize/78948490
!!! Reference: https://en.wikipedia.org/wiki/Box–Muller_transfo
MODULE MATHS
    implicit none
    real(kind=8),parameter :: pi=4.0d0*atan(1.0d0),twopi=2.0d0*pi
    CONTAINS
    function normaldist(mean,std,n) result(r)
        implicit none
        real(kind=8),intent(in) :: mean,std
        integer,intent(in) :: n
        real(kind=8) :: r(n,2)
        real(kind=8),dimension(n,2) :: zeta
        !call random_seed()
        call random_number(zeta)
        r(:,1) = dsqrt(-2.0d0*log(zeta(:,1)))*cos(twopi*zeta(:,2))
        r(:,2) = dsqrt(-2.0d0*log(zeta(:,1)))*sin(twopi*zeta(:,2))
        r = mean + std * r
    end function normaldist
END MODULE MATHS


!*** Here is the function to calculate the gamma (matrix)
real*8 function gammavalue(no,i,j,z,vz,vx,vt)
    implicit none
    include "para.h"
    integer :: no,i,j !! no: gamma No.?; i/j: matrix index
    real*8 :: z,vz,vx,vt
    real*8 :: zroot
    zroot = sqrt(z)

    select case(no)
    case(1)
        if (i.eq.j) then
            if (i.eq.1) gammavalue = xi/(z*zroot)
            if (i.eq.2) gammavalue = (2.d0*xi*eps)/(3.d0*zroot)
            if (i.eq.3) gammavalue = (4.d0*xi*eps)/(3.d0*zroot)
        else
            gammavalue = 0.d0
        end if
    case(2)
        if ((i-1)*(j-1).eq.0) gammavalue = 0.0d0
        if ((i.eq.1).and.(j.eq.1)) gammavalue = (15.d0*kxi*xi)/(8.d0*z**4)
        if ((i.eq.2).and.(j.eq.2)) gammavalue = (kxi*xi*eps**2)/(18.d0*z**3)
        if ((i.eq.3).and.(j.eq.3)) gammavalue = (2.d0*kxi*xi*eps**2)/(9.d0*z**3)
        if (i*j.eq.6) gammavalue = -(kxi*xi*eps**2)/(9.d0*z**3)
    case(3)
        if ((i.eq.1).and.(j.eq.1)) gammavalue = (21.d0*kxi*vz)/(4.d0*zroot*z**4)
        if ((i.eq.2).and.(j.eq.2)) gammavalue = (kxi*vz*(6.d0+19.d0*eps))/(24.d0*zroot*z**3)
        if ((i.eq.3).and.(j.eq.3)) gammavalue = (kxi*vz*(3.d0+19.d0*eps))/(12.d0*zroot*z**3)
        if (i*j.eq.2) gammavalue = (kxi*((eps+3.d0)*vt-3.d0*vx))/(12.d0*zroot*z**3)
        if (i*j.eq.3) gammavalue = (kxi*((3.d0-eps)*vx-3.d0*vt))/(12.d0*zroot*z**3)
        if (i*j.eq.6) gammavalue = -(kxi*vz*(eps+1.d0))/(4.d0*zroot)
    end select
    return
end function gammavalue

!*** Here is the function to calculate the force depending on z.
real*8 function forcevalue(z)
    !This function would furnish the force on each direction
    !z: vertical position; vz: velocity vector
    implicit none
    real*8 :: z
    real*8 :: zroot,gzzz,gzxx,gztt,coulombmax
    include "para.h"
    
    zroot = sqrt(z)
    gzzz = (21.d0*kxi)/(4.d0*zroot*z**4)
    gzxx = -(kxi)/(4*zroot*z**3)
    gztt = -(kxi)/(4*zroot*z**3)
    coulombmax = 1e-15*0

    forcevalue = -grav*mz*(1-rhosol/rhosty) + (gzzz+gzxx+gztt)/(beta*mz) + coulombmax/(z**2)*exp(-z)

    return
end function forcevalue

!*** Here is the function to calculate inverse mass matrix
real*8 function Minverse(z,i,j)
    implicit none
    include "para.h"
    real*8 :: z,zroot
    integer :: i,j !! matrix index
    zroot = sqrt(z)

    select case(i)
    case(1)
        if (j.eq.1) then
            Minverse = (1.d0/mz) * (1.d0 + (15.d0*kxi)/(8.d0*zroot*z**2))
        else
            Minverse = 0.d0
        end if
    case(2)
        if (j.eq.1) Minverse = 0.d0
        if (j.eq.2) Minverse = (1.d0/mx) * (1.d0 + kxe/(12.d0*zroot*z**2))
        if (j.eq.3) Minverse = -kxe/(12.d0*zroot*z**2*mt)
    case(3)
        if (j.eq.1) Minverse = 0.d0
        if (j.eq.2) Minverse = -kxe/(6.d0*zroot*z**2*mx)
        if (j.eq.3) Minverse = (1.d0/mt) * (1.d0 + kxe/(6.d0*zroot*z**2*mt))
    end select
    return
end function Minverse


!*** Here is the function to cut off Gauss white noise.
real*8 function rectGauss(num,value)
    implicit none
    real*8 num,value
    if (abs(value).le.num) then
        rectGauss=1.d0
    else
        rectGauss=0.d0
    end if
    return
end function rectGauss


!*** Here is the procedure to update gamma.
subroutine updategamma(valeur,no,i,j,z,vz,vx,vt)
    implicit none
    real*8 :: valeur,z,vz,vx,vt
    integer :: no,i,j
    real*8,external :: gammavalue
    include "para.h"
    valeur = gammavalue(no,i,j,z,vz,vx,vt)
end subroutine


!*** Here is the procedure to update force along z direction.
subroutine updateforce(valeur,z)
    implicit none
    real*8 :: valeur,z
    real*8,external :: forcevalue
    !include "para.h"
    valeur = forcevalue(z)
end subroutine


!*** Here is the procedure to update the inverse mass matrix
subroutine updateMinverse(valeur,z,i,j)
    implicit none
    real*8 :: valeur,z
    integer :: i,j
    real*8,external :: Minverse
    valeur = Minverse(z,i,j)
end subroutine


!*** Here is the procedure to update intermediate variables.
subroutine updateintm
    implicit none
    include "para.h"
    real*8 :: gt2
    !intma = exp(-gamma*dt)
    !gt2 = gamma*dt/2.d0
    !intmb = sqrt(tanh(gt2)/gt2)
end subroutine


