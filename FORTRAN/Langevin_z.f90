!!!!! Yilin YE @ Mar. 2022
!!!!! EHD + Brownian motion
!!!!! Based on 1412.0162
!!!!! ONLY ONE VARIABLE «z»

!!!!! Here is a file based on the technique shown in the article 
!!!!! J. Phys. Chem. B 2014, 118, 6466−6474
!!!!! See ""Support Information"" for the Multi-dimensional case


Program main
    use MATHS
    implicit none    
    integer :: Nmax=40000
    !real*8,parameter :: pi=3.14159265358979d0 !!常量赋值
    complex*16,parameter :: Im=(0.d0,1.0d0)    !!虚数单位
    integer i,j,k ! 循环变量；
    real*8 :: time_begin,time_end ! 定义程序开始、结束时间
    !real*8,external :: normaldist
    
    include "para.h"    
    real*8 :: p12,v14,v24,v34,grandN(1,2),coefa(6)
    real*8,allocatable :: position(:),velocity(:),force(:)
    allocate(position(Nmax)); allocate(velocity(Nmax)); allocate(force(Nmax))
    call cpu_time(time_begin)
    open(unit=31,file="random_number_test.txt")
    open(unit=32,file="position_z.txt")
    open(unit=33,file="velocity_z.txt")
    open(unit=34,file="force_z.txt")
    open(unit=35,file="intermediates.txt")

    62 format(2x,'Time Step',15x,'Position z',15x,'p12')
    63 format(2x,'Time Step',12x,'Veloctiy z',15x,'v14',15x,'v24',15x,'v34')
    64 format(2x,'Time Step',12x,'Force z',12x,'grandN+',12x,'grandN-')
    65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
    write(32,62); write(33,63); write(34,64); write(35,65)
    
    72 format(f10.1,2(5X,ES18.8)) !! format for file(32)=position_z.txt
    73 format(f10.1,4(5X,ES15.6)) !! format for file(33)=velocity_z.txt
    74 format(f10.1,3(5X,ES15.6)) !! format for file(34)=force_z.txt
    75 format(f10.1,4(5X,ES20.10)) !! format for file(35)=intermediates.txt
    99 format('It spent',f8.3,' seconds for the whole program.')
    

    rayon = 1.5d-6; rhosty = 1.06e3; rhosol = 1.00e3; grav = 9.80665
    mass = pi*rayon**2*rhosty; mz = mass; mt = mass*rayon**2/2
    temp = 298.0d0; k_B = 1.38064852d-23; beta = 1.d0/(k_B*temp)
    clight = sqrt(2.d0*grav*rayon*(1.d0-rhosol/rhosty))

    kappa = 0.1d0; eps = 0.1d0; xi = 10.0d0; kxi = kappa*xi; kxe = kappa*xi*eps
    coefa(1) = xi;  coefa(2) = 21*kxi/4; coefa(3) = -kxi/4; coefa(4) = kxi/2; coefa(5) = -15*kxi/8; coefa(6) = 1.d0
    dt = 2.0d-3

    !! Initiation 
    position=0.0d0; velocity=0.0d0; force=0.0d0
    gamma=0.0d0; Minv=0.0d0; intma=0.0d0; intmb=0.0d0
    
    position(1)=1.0d-7
    p12=0.0d0; v14=0.0d0; v24=0.0d0; v34=0.0d0

    do i=1,Nmax-1
        grandN = normaldist(0.0d0,1.0d0,1)
        !write(31,*) "grandN = ",grandN
        !write(31,*) "N1 = ",grandN(1,1)
        !write(31,*) "N2 = ",grandN(1,2)
        !write(31,*)

        call updategamma(gamma,position(i),velocity(i))
        call updateintm
        !call updateforce(force(i),position(i),velocity(i))
        call updateMinverse(Minv,position(i))
        v14 = sqrt(intma)*velocity(i) + sqrt((1.0d0-intma)*Minv/beta)*grandN(1,1)
        
        call updategamma(gamma,position(i),v14)
        call updateintm
        call updateforce(force(i),position(i))
        v24 = v14 + dt/2*intmb*Minv*force(i)
        
        call updategamma(gamma,position(i),v24)
        call updateintm
        p12 = position(i) + dt/2*intmb*v24
        
        call updategamma(gamma,p12,v24)
        call updateintm
        position(i+1) = p12 + dt/2*intmb*v24
        
        call updategamma(gamma,position(i+1),v24)
        call updateintm
        call updateforce(force(i+1),position(i+1))
        call updateMinverse(Minv,position(i+1))
        v34 = v24 + dt/2*intmb*Minv*force(i+1)
        
        call updategamma(gamma,position(i+1),v34)
        call updateintm
        velocity(i+1) = sqrt(intma)*v34 + sqrt((1-intma)*Minv/beta)*grandN(1,2)

        write(32,72) 1.0d0*i,position(i),p12
        write(33,73) 1.0d0*i,velocity(i),v14,v24,v34
        write(34,74) 1.0d0*i,force(i),grandN(1,1),grandN(1,2)
        write(35,75) 1.0d0*i,gamma,intma,intmb,Minv
    end do


    deallocate(velocity); deallocate(position); deallocate(force)
    close(31); close(32); close(33); close(34); close(35)
    call cpu_time(time_end)
    write(*,99) time_end-time_begin

end Program main


!*** Here is the module to furnish NORMAL DISTRIBUTION, namely the white noises.
!!! Reference: https://www.cxyzjd.com/article/za36mize/78948490
!!! Reference: https://en.wikipedia.org/wiki/Box–Muller_transfo
MODULE MATHS
    implicit none
    real(kind=8),parameter :: pi=4.0d0*atan(1.d0),twopi=2.0d0*pi
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
real*8 function gammavalue(z,vz)
    implicit none
    include "para.h"
    real*8 :: z,vz
    real*8 :: zroot
    zroot = sqrt(z)
    gammavalue = xi/zroot**3 + ((15*xi)/(8*z**4)+(21*vz)/(4*zroot**9))*kxi
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
    gzzz = 21*kxi/(4*zroot**9)
    gzxx = -kxi/(4*zroot**7)
    gztt = -kxi/(4*zroot**7)
    coulombmax = 1e-15*0

    forcevalue = -grav*mz*(1-rhosol/rhosty) + (gzzz+gzxx+gztt)/(beta*mz) + coulombmax/(z**2)*exp(-z)
    return
end function forcevalue

!*** Here is the function to calculate inverse mass matrix
real*8 function Minverse(z)
    implicit none
    include "para.h"
    real*8 :: z,zroot 
    zroot = sqrt(z)

    Minverse = 1/mz * (1 + (15*kxi)/(8*zroot**5))
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
subroutine updategamma(valeur,z,vz)
    implicit none
    real*8 :: valeur,z,vz
    real*8,external :: gammavalue
    include "para.h"
    valeur = gammavalue(z,vz)
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
subroutine updateMinverse(valeur,z)
    implicit none
    real*8 :: valeur,z
    real*8,external :: Minverse
    valeur = Minverse(z)
end subroutine


!*** Here is the procedure to update intermediate variables.
subroutine updateintm
    implicit none
    include "para.h"
    real*8 :: gt2
    intma = exp(-gamma*dt)
    gt2 = gamma*dt/2.d0
    intmb = sqrt(tanh(gt2)/gt2)
end subroutine


