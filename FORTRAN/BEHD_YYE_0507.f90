!!!!! Yilin YE @ Jun. 2022
!!!!! EHD + Brownian motion
!!!!! discretisation by Euler–Maruyama method
!!!!! effective friction matrix
!!!!! modified noise correlator amplitude


Program main
    use MATHS
    implicit none    
    integer :: Nmax = 8000*16 !! Define the maximum number of steps
    integer :: dtmax = 4000*10 !! Define the maximum ∆t to calculate MSD
    integer i,j,k,l,count_kappa,count_delta !! loop variables
    real*8 :: time_begin,time_end,temps_debut,temps_fin !! define variables to record time consumed
    real*8,external :: rectgauss,Minverse,gammavalue,forcevalue

    include "para.h"
    include "loop.h"
    !real*8,parameter :: pi=3.14159265358979d0 !! constant
    !complex*16,parameter :: Im=(0.d0,1.0d0)    !! imaginary unit
    character*128 ::  rue_total, no_kappa, no_delta
    integer :: date(8)
    
    real*8 :: masse(3),amplitude(3,3),fext(3),tratio,msdanax,msdanat
    real*8,allocatable :: position(:,:),velocity(:,:),force(:) !! define arrays for position/velocity/force
    real*8,allocatable :: msdx(:),msdt(:),sumx(:),sumt(:),Dcoefx(:),Dcoeft(:)
    
    
    call cpu_time(temps_debut)
    call cpu_time(time_begin)
    
    61 format(2x,'Time Step',15x,'Noise z',15x,'Noise x',15x,'Noise Θ')
    62 format(2x,'Time Step',12x,'Position z/∆',12x,'Position x',12x,'Angle Θ')
    63 format(2x,'Time Step',12x,'Veloctiy z',15x,'Velocity x',15x,'Velocity Θ')
    64 format(2x,'Time Step',12x,'Force z',12x,'Noise z',12x,'Noise x',12x,'Noise Θ')
    !65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
    66 format(2x,'Time Step',12x,'M_zz',12x,'M_xx',12x,'M_xΘ',12x,'M_Θx',12x,'M_ΘΘ')
    67 format(2x,'Time Step',10x,'γv_zz',10x,'γv_xx',10x,'γv_ΘΘ',10x,'γv_zx',10x,'γv_zΘ',10x,'γv_xΘ')
    68 format(2x,'Time Step',9x,'ana. MSD_x',9x,'<∆x**2>',9x,'D_x',9x,'ana. MSD_t',9x,'<∆Θ**2>',9x,'D_t')
    
    71 format(f10.4,3(5X,ES18.8)) !! format for file(31)=random_number_test.txt
    72 format(f10.4,3(5X,ES18.8)) !! format for file(32)=position_z.txt
    73 format(f10.4,3(5X,ES20.6)) !! format for file(33)=velocity_z.txt
    74 format(f10.4,4(5X,ES15.6)) !! format for file(34)=force_z.txt
    !75 format(f10.1,4(5X,ES20.10)) !! format for file(35)=intermediates.txt
    76 format(f10.4,5(5X,ES15.6)) !! format for file(36)=Mass_Matrix_Inverse.txt
    77 format(f8.4,6(5X,ES12.4)) !! format for file(371/372/373)
    78 format(f10.4,6(5X,ES10.3)) !! format for file(38-)
    
    95 format('It spent',f8.2,' seconds for the whole initiation.')
    96 format('!!! Here is the case. κ = ',f10.3,' ; and ∆_0 = ',f10.3)
    97 format(f8.2,' seconds spent for Euler method.')
    98 format(f8.2,' seconds spent for MSD calculations.')
    99 format('It spent',f8.2,' seconds for the whole program.')
    

    i=1 !! Suppose i the index to select input file. 1, input.txt; 2, input_real.txt.
    call obtaindata(i)
    !goto 1106

    1107 continue
    !! Parameters read from input.txt
    !grav = 9.80665d0 !! unit m/s**2
    !temp = 298.0d0 !! unit K
    !rayon = 1.5d-6 !! unit m
    !rhosty = 1.06e3 !! unit kg/m**3
    !rhosol = 1.00e3 !! unit kg/m**3

    !! Constant Parameters
    mass = pi*rayon**2*rhosty !! unit kg/m
    mz = mass; mx = mass; mt = mass*rayon**2/2; masse(1)=mz; masse(2)=mx; masse(3)=mt
    k_B = 1.38064852d-23 !! unit J/K = 
    beta = 1.0d0/(k_B*temp) !! unit J**-1 = 
    clight = sqrt(2.d0*grav*rayon*(1.0d0-rhosol/rhosty)) !! unit m/s, the maximum speed
    if (i.eq.1) then
        eps = 0.1d0 !! dimensionless parameter
        xi = 1.0d0 !! dimensionless parameter
    end if

    !! Time gap for each step, t = T * r * sqrt(2*eps) / clight. Here we pose dt ~ ∆T
    tratio = 1.0d0 / (20.0d0*10)
    dt = 2.0d-6*clight/(rayon*sqrt(2.0d0*eps)) * tratio
    1108 continue


    count_kappa = 0
    count_delta = 0
    call date_and_time(values=date)
    open(unit=30,file="output_record.txt")
    call cpu_time(time_end)
    write(30,95) time_end-time_begin !! It spent ? seconds for the whole initiation.
    write(30,*)
    write(30,*) "Below we simulate the Brownian motion of a cyclindrical particle with parameters:"
    write(30,*) "Radius = ",rayon,"µm;   ","Temperature = ",temp,"K;   "
    write(30,*); write(30,*)


    !write(*,*) 'clight = ',clight

    1201 continue
    select case(loop_kappa)
    case(0)
        kappa = 1.0d-2
    case(1)
        kappa = min_kappa + gap_kappa * count_kappa
    end select
    kxi = kappa*xi; kxe = kappa*xi*eps !! non-dimensional parameters w/o units.
    
    1202 continue
    select case(loop_delta)
    case(0)
        ini_height = 1.0d0
        !position(1,1) = 1.0d0
    case(1)
        ini_height = min_delta + gap_delta * count_delta
        !position(1,1) = min_delta + gap_delta * count_delta
    end select

    
    call cpu_time(time_begin)
    write(30,96) kappa, ini_height

    !! Must decalre all allocatable variables in advance.
    allocate(position(3,Nmax)); allocate(velocity(3,Nmax)); allocate(force(Nmax))
    allocate(msdx(Nmax)); allocate(msdt(Nmax)); allocate(sumx(dtmax)); allocate(sumt(dtmax));
    allocate(Dcoefx(Nmax)); allocate(Dcoeft(Nmax))
    !! Initiation 
    position=0.0d0; velocity=0.0d0; force=0.0d0; amplitude=0.0d0
    gmaeff=0.0d0; Minv=0.0d0; fext=0.0d0
    !position(2,1)=0.0d0; position(3,1)=0.0d0
    !position(1,1)=1.0d0 !! Define initial height ∆(0)
    position(1,1) = ini_height

    !rue_total  = "/Results/"//char(date(1))//char(date(2))//char(date(3))//"/"//&
    !&"kappa="//char(count_kappa)//"delta="//char(count_delta)
    write(no_kappa,*) count_kappa; no_kappa = ADJUSTL(no_kappa); write(*,*) "no_kappa = ", no_kappa
    write(no_delta,*) count_delta; no_delta = ADJUSTL(no_delta); write(*,*) "no_delta = ", no_delta
    rue_total = "k_"//trim(no_kappa)//"-"//"d_"//trim(no_delta)//"-"
    write(*,*) "rue_total = ", rue_total

    !! add the below path to each file
    !TRIM(rue_total)+
    open(unit=31,file=TRIM(rue_total)//"random_number_test.txt"); write(31,96) kappa, ini_height
    open(unit=32,file=TRIM(rue_total)//"positions.txt"); write(32,96) kappa, ini_height
    open(unit=33,file=TRIM(rue_total)//"velocitys.txt"); write(33,96) kappa, ini_height
    open(unit=34,file=TRIM(rue_total)//"forces.txt"); write(34,96) kappa, ini_height
    !open(unit=35,file="inter_Euler_Maruyama.txt")
    open(unit=36,file=TRIM(rue_total)//"Mass_Matrix_Inverse.txt"); write(36,96) kappa, ini_height
    open(unit=371,file=TRIM(rue_total)//"gamma0.txt"); write(371,96) kappa, ini_height
    open(unit=372,file=TRIM(rue_total)//"gamma1.txt"); write(372,96) kappa, ini_height
    open(unit=373,file=TRIM(rue_total)//"gamma1v.txt"); write(373,96) kappa, ini_height
    open(unit=381,file=TRIM(rue_total)//"MSD.txt"); write(381,96) kappa, ini_height
    write(31,61); write(32,62); write(33,63); write(34,64); !write(35,65)
    write(36,66); write(371,67); write(372,67); write(373,67); write(381,68)
    
    

    1101 continue
    do i=1,Nmax-1
        noise = normaldist(0.0d0,dt,3)
        write(31,71) i*tratio,noise(1,1),noise(2,1),noise(3,1)

        !! Update all elements in the effective friction matrix, and mass matrix
        do j=1,3
            do k=1,3
                do l=1,3
                    !call updategamma(gmaeff(j,k,l),j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                    gmaeff(j,k,l) = gammavalue(j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                end do
                !call updateMinverse(Minv(j,k),position(1,i),j,k)
                Minv(j,k) = Minverse(position(1,i),j,k)
            end do
            !amplitude(j,j) = sqrt( 2.0d0 * (gmaeff(1,j,j) - gmaeff(2,j,j)) / (beta * masse(j)) )
        end do
        !amplitude(1,1) = amplitude(1,1) * sqrt(2.0d0 / (clight**2 * eps))
        !amplitude(1,1) = sqrt(2.0d0/(clight**2 * eps)) * sqrt(2.0d0*(gmaeff(1,1,1) - gmaeff(2,1,1))/(beta*mz))
        amplitude(1,1) = sqrt(4.0d0 * mz * (gmaeff(1,1,1) - gmaeff(2,1,1)) /(clight**2 * eps * beta)) * Minv(1,1)
        !amplitude(2,2) = amplitude(2,2) * clight
        !amplitude(2,2) = 1.0d0/clight * sqrt(2.0d0*(gmaeff(1,2,2) - gmaeff(2,2,2))/(beta*mx))
        amplitude(2,2) = 1.0d0/clight * sqrt(2.0d0 * mx * (gmaeff(1,2,2) - gmaeff(2,2,2))/beta) * Minv(2,2)
        !amplitude(3,3) = amplitude(3,3) * rayon / clight
        !amplitude(3,3) = rayon/clight * sqrt(2.0d0*(gmaeff(1,3,3) - gmaeff(2,3,3))/(beta*mt))
        amplitude(3,3) = rayon/clight * sqrt(2.0d0 * mt * (gmaeff(1,3,3) - gmaeff(2,3,3))/beta) * Minv(3,3)


        !call updateforce(force(i),position(1,i))
        force(i) = forcevalue(position(1,i))
        !fext(1) = forcevalue(position(1,i))
        fext(1) = force(i)!; fext(2)=0.0d0; fext(3)=0.0d0


        velocity(:,i+1) = velocity(:,i) + MATMUL(amplitude(:,:),noise(:,1)) + &
        & dt * ( MATMUL(Minv(:,:),fext(:)) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),velocity(:,i)) )
        !velocity(1,i) = velocity(1,i) / clight * sqrt(2.0d0 / eps)
        !velocity(2,i) = velocity(2,i) / clight
        !velocity(3,i) = velocity(3,i) / clight * rayon

        position(:,i+1) = position(:,i) + velocity(:,i) * dt
        !position(1,i+1) = position(1,1)


        !Dcoefx(i) = (gmaeff(1,2,2) - gmaeff(2,2,2))/(beta*mx*gmaeff(1,2,2)**2)
        Dcoefx(i) = (1.0d0 - gmaeff(2,2,2)/gmaeff(1,2,2))/(beta*mx*gmaeff(1,2,2)) * (1.0/clight)**2
        !Dcoeft(i) = (gmaeff(1,3,3) - gmaeff(2,3,3))/(beta*mt*gmaeff(1,3,3)**2)
        Dcoeft(i) = (1.0d0 - gmaeff(2,3,3)/gmaeff(1,3,3))/(beta*mt*gmaeff(1,3,3)) * (rayon/clight)**2


        write(32,72) i*tratio,position(1,i),position(2,i),position(3,i)
        write(33,73) i*tratio,velocity(1,i),velocity(2,i),velocity(3,i)
        write(34,74) i*tratio, Minv(1,1)*force(i), amplitude(1,1)*noise(1,1), amplitude(2,2)*noise(2,1), amplitude(3,3)*noise(3,1)
        !write(35,75) i*tratio,gamma,intma,intmb,Minv
        write(36,76) i*tratio,Minv(1,1),Minv(2,2),Minv(1,2),Minv(2,1),Minv(2,2)
        write(371,77) i*tratio,gmaeff(1,1,1),gmaeff(1,2,2),gmaeff(1,3,3),gmaeff(1,1,2),gmaeff(1,1,3),gmaeff(1,2,3)
        write(372,77) i*tratio,gmaeff(2,1,1),gmaeff(2,2,2),gmaeff(2,3,3),gmaeff(2,1,2),gmaeff(2,1,3),gmaeff(2,2,3)
        write(373,77) i*tratio,gmaeff(3,1,1),gmaeff(3,2,2),gmaeff(3,3,3),gmaeff(3,1,2),gmaeff(3,1,3),gmaeff(3,2,3)
    end do
    1102 continue
    call cpu_time(time_end)
    write(30,97) time_end-time_begin !! Print the tine consumed for Euler method.
    

    call cpu_time(time_begin)
    1103 continue
    goto 1104
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
        msdanat = 1.0d0/(beta*mt*gmaeff(1,3,3)) * ((exp(-i*gmaeff(1,3,3)))/gmaeff(1,3,3) + i)

        write(381,78) i*tratio, msdanax, sumx(i)/j, Dcoefx(i)*i*1.0d0, msdanat, sumt(i)/j, Dcoeft(i)*i*1.0d0
    end do
    1104 continue
    call cpu_time(time_end)
    write(30,98) time_end-time_begin !! Print the time consumed for MSD calculations

    
    1105 continue
    deallocate(velocity); deallocate(position); deallocate(force) 
    deallocate(msdx); deallocate(msdt); deallocate(sumx); deallocate(sumt); deallocate(Dcoefx); deallocate(Dcoeft)
    close(31); close(32); close(33); close(34); close(35); close(36); close(371); close(372); close(373); close(381)
    1106 continue 
    write(30,*) !! Print a white rank in the file.
    
    if ((loop_delta.eq.1).and.(ini_height.lt.max_delta)) then
        count_delta = count_delta + 1
        goto 1202
    end if
    if ((loop_kappa.eq.1).and.(kappa.lt.max_kappa)) then
        count_kappa = count_kappa + 1
        goto 1201
    end if
    
    

    call cpu_time(temps_fin)
    write(30,*); write(30,*); write(30,*)
    write(30,99) temps_fin-temps_debut
    write(*,99) temps_fin-temps_debut
    close(30)

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


!*** Here is the subroutine to obtain parameters
subroutine obtaindata(i)
    implicit none
    include "para.h"
    include "loop.h"
    integer :: i,status1
    character*128 :: msg,ruein

    select case(i)
    case(1)
        ruein="input.txt"
        open(unit=201,file=TRIM(ruein),form='formatted',status='old',action='read',iostat=status1,iomsg=msg)
        read(201,*) grav !! gravity acceleration
        read(201,*) temp !! ambient temperature
        read(201,*) rayon !! particle radius
        read(201,*) rhosty !! particle density
        read(201,*) rhosol !! solvent density
        read(201,*) l_D !! Debye length
        read(201,*) l_B !! Boltzmann length
        read(201,*) loop_kappa; read(201,*) min_kappa; read(201,*) max_kappa; read(201,*) gap_kappa
        read(201,*) loop_delta; read(201,*) min_delta; read(201,*) max_delta; read(201,*) gap_delta
        close(201)
    case(2)
        ruein="input_real.txt"
        !! Space left for possible parameters in reality such as viscosity.
    case default
        ruein="not_exist"
    end select
    
    !write(*,*) grav,temp,rayon,rhosty,rhosol
    !write(*,*) l_D, l_B
    !write(*,*) loop_kappa,min_kappa,max_kappa,gap_kappa
    !write(*,*) loop_delta,min_delta,max_delta,gap_delta
end subroutine



!*** Here is the function to calculate the gamma (effective friction matrix)
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

!*** Here is the function to calculate the force depending on z/∆.
real*8 function forcevalue(z)
    !This function would furnish the force on each direction
    !z: vertical position; vz: velocity vector
    implicit none
    real*8 :: z
    real*8 :: zroot,gzzz,gzxx,gztt,coulombmax
    include "para.h"
    
    zroot = sqrt(z)
    gzzz = (21.0d0*kxi)/(4.0d0*zroot*z**4)
    gzxx = -(kxi)/(4.0d0*zroot*z**3)
    gztt = -(kxi)/(4.0d0*zroot*z**3)
    coulombmax = 4.0d0 !! namely the "B" shown by Maxime, PHYSICAL REVIEW RESEARCH 3, L032011 (2021)

    
    !! Attention to the unit! As for the force, it should be [kg·m·s^{-2}]
    !! However, since m = rho * π * r**2 = [kg/m], the force unit would be [kg/s^2]
    !! Finally, we take all parameters dimensionless.
    forcevalue =  - mz*grav*(1.0d0-rhosol/rhosty) * 2.0d0 * rayon / clight**2 + &
    & 1.0d0 * (gzzz+gzxx+gztt)/(beta) * 2.0d0 / (clight**2 * eps) + &
    & coulombmax/(beta*l_D) * exp(-z*rayon*eps/l_D) * 2.0d0 * rayon / clight**2
    
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


!*** Here is the procedure to update force along z/∆ direction.
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
