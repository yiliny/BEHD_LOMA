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
    integer :: Nmax = 5000*12 !! Define the maximum number of steps
    integer :: dtmax = 1000*20 !! Define the maximum ∆t to calculate MSD
    !real*8,parameter :: pi=3.14159265358979d0 !! constant
    integer i,j,k,l !! loop variables
    real*8 :: time_begin,time_end !! define variables to record time consumed
    real*8,external :: rectgauss    
    !complex*16,parameter :: Im=(0.d0,1.0d0)    !! imaginary unit
    include "para.h"
    real*8 :: masse(3),amplitude(3,3),fext(3),tratio,msdanax,msdanat
    real*8,allocatable :: position(:,:),velocity(:,:),force(:) !! define arrays for position/velocity/force
    real*8,allocatable :: msdx(:),msdt(:),sumx(:),sumt(:),Dcoefx(:),Dcoeft(:)
    
    
    !! Must decalre all allocatable variables in advance.
    allocate(position(3,Nmax)); allocate(velocity(3,Nmax)); allocate(force(Nmax))
    allocate(msdx(Nmax)); allocate(msdt(Nmax)); allocate(sumx(dtmax)); allocate(sumt(dtmax));
    allocate(Dcoefx(Nmax)); allocate(Dcoeft(Nmax))
    

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
    open(unit=381,file="MSD.txt")


    61 format(2x,'Time Step',15x,'Noise z',15x,'Noise x',15x,'Noise Θ')
    62 format(2x,'Time Step',12x,'Position z/∆',12x,'Position x',12x,'Angle Θ')
    63 format(2x,'Time Step',12x,'Veloctiy z',15x,'Velocity x',15x,'Velocity Θ')
    64 format(2x,'Time Step',12x,'Force z',12x,'Noise z',12x,'Noise x',12x,'Noise Θ')
    !65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
    66 format(2x,'Time Step',12x,'M_zz',12x,'M_xx',12x,'M_xΘ',12x,'M_Θx',12x,'M_ΘΘ')
    67 format(2x,'Time Step',10x,'γv_zz',10x,'γv_xx',10x,'γv_ΘΘ',10x,'γv_zx',10x,'γv_zΘ',10x,'γv_xΘ')
    68 format(2x,'Time Step',9x,'ana. MSD_x',9x,'<∆x**2>',9x,'D_x',9x,'ana. MSD_t',9x,'<∆Θ**2>',9x,'D_t')
    write(31,61); write(32,62); write(33,63); write(34,64); !write(35,65)
    write(36,66); write(371,67); write(372,67); write(373,67); write(381,68)
    
    
    71 format(f10.4,3(5X,ES18.8)) !! format for file(31)=random_number_test.txt
    72 format(f10.4,3(5X,ES18.8)) !! format for file(32)=position_z.txt
    73 format(f10.4,3(5X,ES20.6)) !! format for file(33)=velocity_z.txt
    74 format(f10.4,4(5X,ES15.6)) !! format for file(34)=force_z.txt
    !75 format(f10.1,4(5X,ES20.10)) !! format for file(35)=intermediates.txt
    76 format(f10.4,5(5X,ES15.6)) !! format for file(36)=Mass_Matrix_Inverse.txt
    77 format(f8.4,6(5X,ES12.4)) !! format for file(371/372/373)
    78 format(f10.4,6(5X,ES10.3)) !! format for file(38-)
    99 format('It spent',f8.2,' seconds for the whole program.')
    

    rayon = 1.5d-6 !! unit m
    rhosty = 1.06e3; rhosol = 1.00e3 !! unit kg/m**3
    grav = 9.80665d0 !! unit m/s**2
    mass = pi*rayon**2*rhosty !! unit kg/m
    mz = mass; mx = mass; mt = mass*rayon**2/2; masse(1)=mz; masse(2)=mx; masse(3)=mt
    temp = 298.0d0 !! unit K
    k_B = 1.38064852d-23 !! unit J/K = 
    beta = 1.d0/(k_B*temp) !! unit J**-1 = 
    clight = sqrt(2.d0*grav*rayon*(1.0d0-rhosol/rhosty)) !! unit m/s
    kappa = 1.0d-4; eps = 0.1d0; xi = 0.1d0; kxi = kappa*xi; kxe = kappa*xi*eps !! non-dimensional parameters w/o units.
    !coefa(1)=xi;  coefa(2)=21.d0*kxi/4.d0; coefa(3)=-kxi/4.d0; coefa(4)=kxi/2.d0; coefa(5)=-15.d0*kxi/8.d0; coefa(6)=1.d0
    !coefb(1)=2.d0/3.d0*eps*xi;  coefb(2)=19.d0*kxe/24.d0; coefb(3)=-kxe/6.d0; coefb(4)=kxe/12.d0; coefb(5)=-kxe/12.d0; coefb(6)=0.d0
    !coefc(1)=4.d0/3.d0*eps*xi;  coefc(2)=19.d0*kxe/12.d0; coefc(3)=-kxe/3.d0; coefc(4)=kxe/6.0d0; coefc(5)=-kxe/6.0d0; coefc(6)=0.d0
    

    !! Time gap for each step, t = T * r * sqrt(2*eps) / clight. Here we pose dt ~ ∆T
    tratio = 1.0d0 / (20.0d0*10)
    dt = 1.0d-3*clight/(rayon*sqrt(2.0d0*eps)) * tratio


    !! Initiation 
    position=0.0d0; velocity=0.0d0; force=0.0d0; amplitude=0.0d0
    gmaeff=0.0d0; Minv=0.0d0; fext=0.0d0
    position(2,1)=0.0d0; position(3,1)=0.0d0
    position(1,1)=1.0d0 !! Define initial height ∆(0)

    do i=1,Nmax-1
        !noise = normaldist(0.0d0,1.0d0,3)
        noise = normaldist(0.0d0,dt,3)
        !do j=1,3
        !    noise(j,1) = noise(j,1)*rectGauss(noise(j,1),2.d0)
        !end do
        write(31,71) i*tratio,noise(1,1),noise(2,1),noise(3,1)

        !! Update all elements in the effective friction matrix, and mass matrix
        do j=1,3
            do k=1,3
                do l=1,3
                    call updategamma(gmaeff(j,k,l),j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                end do
                call updateMinverse(Minv(j,k),position(1,i),j,k)
            end do
            amplitude(j,j) = sqrt( 2.0d0 * (gmaeff(1,j,j) - gmaeff(2,j,j)) / (beta * masse(j)) )
        end do
        amplitude(1,1) = amplitude(1,1) * sqrt(2.0d0 / (clight**2 * eps))
        amplitude(2,2) = amplitude(2,2) * clight
        amplitude(3,3) = amplitude(3,3) * rayon / clight


        call updateforce(force(i),position(1,i))
        fext(1) = force(i)


        velocity(:,i+1) = velocity(:,i) + MATMUL(amplitude(:,:),noise(:,1)) + &
        & dt * ( fext(:) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),velocity(:,i)) )
        !velocity(1,i) = velocity(1,i) / clight * sqrt(2.0d0 / eps)
        !velocity(2,i) = velocity(2,i) / clight
        !velocity(3,i) = velocity(3,i) / clight * rayon

        position(:,i+1) = position(:,i) + velocity(:,i) * dt
        !position(1,i+1) = position(1,1)


        Dcoefx(i) = (gmaeff(1,2,2) - gmaeff(2,2,2))/(beta*mx*gmaeff(1,2,2)**2)
        Dcoeft(i) = (gmaeff(1,3,3) - gmaeff(2,3,3))/(beta*mt*gmaeff(1,3,3)**2)


        write(32,72) i*tratio,position(1,i),position(2,i),position(3,i)
        write(33,73) i*tratio,velocity(1,i),velocity(2,i),velocity(3,i)
        write(34,74) i*tratio,force(i), amplitude(1,1)*noise(1,1), amplitude(2,2)*noise(2,1), amplitude(3,3)*noise(3,1)
        !write(35,75) i*tratio,gamma,intma,intmb,Minv
        write(36,76) i*tratio,Minv(1,1),Minv(2,2),Minv(1,2),Minv(2,1),Minv(2,2)
        write(371,77) i*tratio,gmaeff(1,1,1),gmaeff(1,2,2),gmaeff(1,3,3),gmaeff(1,1,2),gmaeff(1,1,3),gmaeff(1,2,3)
        write(372,77) i*tratio,gmaeff(2,1,1),gmaeff(2,2,2),gmaeff(2,3,3),gmaeff(2,1,2),gmaeff(2,1,3),gmaeff(2,2,3)
        write(373,77) i*tratio,gmaeff(3,1,1),gmaeff(3,2,2),gmaeff(3,3,3),gmaeff(3,1,2),gmaeff(3,1,3),gmaeff(3,2,3)
    end do


    !! MSD, Mean Square Displacement
    
    do i=1,dtmax !! ∆t?
        sumx(i) = 0.0d0; sumt(i) = 0.0d0
        do j=1,Nmax-dtmax
            msdx(j) = (position(2,j) - position(2,i+j))**2
            sumx(i) = sumx(i) + msdx(j)

            msdt(j) = (position(3,j) - position(3,i+j))**2
            sumt(i) = sumt(i) + msdt(j)
        end do

        msdanax = 1.0d0/(beta*mx*gmaeff(1,2,2)) * ((exp(-i*gmaeff(1,2,2)))/gmaeff(1,2,2) + i)
        msdanat = 0.0d0

        write(381,78) i*tratio, msdanax, sumx(i)/j, Dcoefx(i)*i*1.0d0, msdanat, sumt(i)/j, Dcoeft(i)*i*1.0d0
    end do


    deallocate(velocity); deallocate(position); deallocate(force) 
    deallocate(msdx); deallocate(msdt); deallocate(sumx); deallocate(sumt); deallocate(Dcoefx); deallocate(Dcoeft)
    close(31); close(32); close(33); close(34); close(35); close(36); close(371); close(372); close(373); close(381)
    call cpu_time(time_end)
    write(*,99) time_end-time_begin
    write(*,*) 'clight = ',clight

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
            if (i.eq.2) gammavalue = (2.0d0*xi*eps)/(3.0d0*zroot)
            if (i.eq.3) gammavalue = (4.0d0*xi*eps)/(3.0d0*zroot)
        else
            gammavalue = 0.0d0
        end if
    case(2)
        if ((i-1)*(j-1).eq.0) gammavalue = 0.0d0
        if ((i.eq.1).and.(j.eq.1)) gammavalue = (15.0d0*kxi*xi)/(8.0d0*z**4)
        if ((i.eq.2).and.(j.eq.2)) gammavalue = (kxi*xi*eps**2)/(18.0d0*z**3)
        if ((i.eq.3).and.(j.eq.3)) gammavalue = (2.0d0*kxi*xi*eps**2)/(9.0d0*z**3)
        if (i*j.eq.6) gammavalue = -(kxi*xi*eps**2)/(9.0d0*z**3)
    case(3)
        if ((i.eq.1).and.(j.eq.1)) gammavalue = (21.0d0*kxi*vz)/(4.0d0*zroot*z**4)
        if ((i.eq.2).and.(j.eq.2)) gammavalue = (kxi*vz*(6.0d0+19.0d0*eps))/(24.0d0*zroot*z**3)
        if ((i.eq.3).and.(j.eq.3)) gammavalue = (kxi*vz*(3.0d0+19.0d0*eps))/(12.0d0*zroot*z**3)
        if (i*j.eq.2) gammavalue = (kxi*((3.0d0+eps)*vt-3.0d0*vx))/(12.0d0*zroot*z**3)
        if (i*j.eq.3) gammavalue = (kxi*((3.0d0-eps)*vx-3.0d0*vt))/(12.0d0*zroot*z**3)
        if (i*j.eq.6) gammavalue = -(kxi*vz*(eps+1.0d0))/(4.0d0*zroot)
    end select
    return
end function gammavalue

!*** Here is the function to calculate the force depending on z.
real*8 function forcevalue(z)
    !This function would furnish the force on each direction
    !z: vertical position; vz: velocity vector
    implicit none
    real*8 :: z
    real*8 :: zroot,gzzz,gzxx,gztt!,coulombmax
    include "para.h"
    
    zroot = sqrt(z)
    gzzz = (21.0d0*kxi)/(4.0d0*zroot*z**4)
    gzxx = -(kxi)/(4.0d0*zroot*z**3)
    gztt = -(kxi)/(4.0d0*zroot*z**3)
    !coulombmax = 1.0e-15*0.0d0

    !! Attention to the unit! As for the force, it should be [kg·m·s^{-2}]
    !! However, since m = rho * π * r**2 = [kg/m], the force unit would be [kg/s^2]
    forcevalue = ( - grav*(1.0d0-rhosol/rhosty) + (gzzz+gzxx+gztt)/(mz*beta) ) !* 2.0d0 * rayon / clight**2
    !+ coulombmax/(z**2)*exp(-z)
    

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
            Minverse = (1.0d0/mz) * (1.0d0 + (15.0d0*kxi)/(8.0d0*zroot*z**2))
        else
            Minverse = 0.0d0
        end if
    case(2)
        if (j.eq.1) Minverse = 0.0d0
        if (j.eq.2) Minverse = (1.0d0/mx) * (1.0d0 + kxe/(12.0d0*zroot*z**2))
        if (j.eq.3) Minverse = -kxe/(12.0d0*zroot*z**2*mt)
    case(3)
        if (j.eq.1) Minverse = 0.0d0
        if (j.eq.2) Minverse = -kxe/(6.0d0*zroot*z**2*mx)
        if (j.eq.3) Minverse = (1.0d0/mt) * (1.0d0 + kxe/(6.0d0*zroot*z**2))
    end select
    return
end function Minverse


!*** Here is the function to cut off Gauss white noise.
real*8 function rectGauss(num,value)
    implicit none
    real*8 num,value
    if (abs(value).le.num) then
        rectGauss=1.0d0
    else
        rectGauss=0.0d0
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
