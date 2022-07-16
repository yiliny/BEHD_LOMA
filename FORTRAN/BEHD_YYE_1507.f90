!!!!! Yilin YE @ Jun. 2022
!!!!! EHD + Brownian motion
!!!!! discretisation by Euler–Maruyama method
!!!!! effective friction matrix
!!!!! modified noise correlator amplitude


Program main
    use MATHS
    implicit none    
    integer :: Nmax = 10000*1 !! Define the maximum number of steps. Note ∆t=0.1ms,Nmax=1e5,T->100
    integer :: dtmax = 8000*10 !! Define the maximum ∆t to calculate MSD
    integer i,j,k,l,count_kappa,count_delta,ratio_sim,ratio_msd,count_sim,count_msd,crash,ratio_pdf,num_particle !! loop variables
    real*8 :: time_begin,time_end,temps_debut,temps_fin,time_sim !! define variables to record time consumed
    real*8,external :: rectgauss,Minverse,gammavalue,forcevalue

    include "para.h"
    include "loop.h"
    !real*8,parameter :: pi=3.14159265358979d0 !! constant
    !complex*16,parameter :: Im=(0.d0,1.0d0)    !! imaginary unit
    character*128 ::  rue_total, no_kappa, no_delta !! Defin path name for outputs
    integer :: date(8)
    
    real*8 :: masse(3),amplitude(3,3),fext(3) !! Define matrices for Euler method, used in MATMUL()
    real*8 :: tratio !! Define the ratio for modified dimensionless time gap ∆T
    real*8 :: itm_p(3,2),itm_v(3,2),itm_f !! Define intermediates to facilitate the codes
    real*8 :: msdx_ana,msdt_ana,msdc_ana !! Define three analytical MSD, <∆r_x **2> & <∆r_Θ **2> & <∆r_x * ∆r_Θ>
    real*8,allocatable :: position(:,:),velocity(:,:),force(:) !! Define arrays for position/velocity/force
    integer,allocatable :: pdfz(:) !! Define array for PDF (Probability Distribution Function), count the number in the limited height
    real*8,allocatable :: msdx_num(:),msdt_num(:),sumx(:),sumt(:),Dcoefx(:),Dcoeft(:)
    
    
    call cpu_time(temps_debut)
    call cpu_time(time_begin)
    
    61 format(2x,'Time',15x,'Noise z',15x,'Noise x',15x,'Noise Θ')
    62 format(2x,'Time',12x,'Position z/∆',12x,'Position x',12x,'Angle Θ')
    63 format(2x,'Time',12x,'Veloctiy z',15x,'Velocity x',15x,'Velocity Θ')
    64 format(2x,'Time',12x,'Force z',12x,'Noise z',12x,'Noise x',12x,'Noise Θ')
    !65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
    66 format(2x,'Time',12x,'M_zz',12x,'M_xx',12x,'M_xΘ',12x,'M_Θx',12x,'M_ΘΘ')
    67 format(2x,'Time',10x,'γv_zz',10x,'γv_xx',10x,'γv_ΘΘ',10x,'γv_zx',10x,'γv_zΘ',10x,'γv_xΘ')
    681 format(2x,'Time',9x,'ana. MSD_x',9x,'<∆x**2>',9x,'D_x',9x,'ana. MSD_t',9x,'<∆Θ**2>',9x,'D_t')
    682 format(2x,'Index',12x,'Height min',12x,'Height max',12x,'Ratio',12x,'Gibbs-Boltzmann')
    
    71 format(f16.9,3(5X,ES18.8)) !! format for file(31)=random_number_test.txt
    72 format(f16.9,3(5X,ES18.8)) !! format for file(32)=position_z.txt
    73 format(f16.9,3(5X,ES20.6)) !! format for file(33)=velocity_z.txt
    74 format(f16.9,4(5X,ES15.6)) !! format for file(34)=force_z.txt
    !75 format(f10.1,4(5X,ES20.10)) !! format for file(35)=intermediates.txt
    76 format(f16.9,5(5X,ES15.6)) !! format for file(36)=Mass_Matrix_Inverse.txt
    77 format(f8.4,6(5X,ES12.4)) !! format for file(371/372/373)
    781 format(f16.9,6(5X,ES10.3)) !! format for file(381)=MSD.txt
    782 format(I8,2(10X,f12.6),12X,f12.7,12X,ES15.6) !! format for file(382)=PDFz.txt
    
    95 format('It spent',f8.2,' seconds for the whole initiation.')
    96 format('!!! Here is the case. κ = ',f10.3,' ; and ∆_0 = ',f10.3)
    97 format(f8.2,' seconds spent for Euler method.')
    981 format(f8.2,' seconds spent for MSD calculations.')
    982 format(f8.2,' seconds spent for PDF output.')
    99 format('It spent',f8.2,' seconds for the whole program.')
    

    i=1 !! Suppose i the index to select input file. 1, input.txt; 2, input_real.txt.
    call obtaindata(i) !! Read all necessary parameters.
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
    k_B = 1.38064852d-23 !! unit J/K
    beta = 1.0d0/(k_B*temp) !! unit J**-1 = kg**-1 * m**-2 * s**2
    clight = sqrt(2.d0*grav*rayon*(1.0d0-rhosol/rhosty)) !! unit m/s, the maximum speed
    csquare = clight**2
    if (i.eq.1) then
        eps = 0.1d0 !! dimensionless parameter, delta = ∆ * r * eps
        xi = 1.0d0 !! dimensionless parameter
    end if

    !! Time gap for each step, t = T * r * sqrt(2*eps) / clight. Here we pose dt ~ ∆T
    tratio = 1.0d0 / (20.0d0*10) !! Define the time discretisation ratio
    dt = 1.0d-4 * clight / (rayon * sqrt(2.0d0 * eps)) * tratio
    1108 continue


    ratio_sim = 40; ratio_msd = 40 !! Here are two variables to accelerate calculations for Euler method & MSD.
    ratio_pdf = 800 !! Here is the ratio to determine how many data would be taken part for PDF calculation.
    
    call date_and_time(values=date)
    open(unit=30,file="output_record.txt")
    call cpu_time(time_end) !! Record the end time for the whole initiation.
    write(30,95) time_end-time_begin !! It spent ? seconds for the whole initiation.
    write(30,*)
    write(30,*) "Below we simulate the Brownian motion of a cyclindrical particle with parameters:"
    write(30,*) "Radius = ",rayon," (m);   ","Temperature = ",temp," (K);   "
    write(30,*); write(30,*)
    !write(*,*) 'clight = ',clight

    count_kappa = 0 !! Reset 
    1201 continue
    select case(loop_kappa) !! Determine if loop or not for kappa
    case(0)
        kappa = 1.0d-2 * 0.0d0
    case(1)
        kappa = min_kappa + gap_kappa * count_kappa
    end select
    kxi = kappa*xi; kxe = kappa*xi*eps !! non-dimensional parameters w/o units.
    
    count_delta = 0 !! Reset
    1202 continue
    select case(loop_delta) !! Determine if loop of not for ∆
    case(0)
        ini_height = 1.0d0
        !position(1,1) = 1.0d0
    case(1)
        ini_height = min_delta + gap_delta * count_delta
        !position(1,1) = min_delta + gap_delta * count_delta
    end select

    
    call cpu_time(time_begin)
    write(30,96) kappa, ini_height

    1207 continue !! Restart if the simulation crash (the particle drops into the surface)
    count_sim = 0; count_msd = 0 !! Here are two variables to record the number of big loops for Euler method.
    crash = 0 !! Define the variable to determine whether the particle drops into the surface

    !! Must decalre all allocatable variables in advance.
    allocate(position(3,Nmax+1)); allocate(velocity(3,Nmax+1)); allocate(force(Nmax+1))
    allocate(pdfz(ratio_pdf))
    allocate(msdx_num(Nmax)); allocate(msdt_num(Nmax)); allocate(sumx(dtmax)); allocate(sumt(dtmax));
    allocate(Dcoefx(Nmax)); allocate(Dcoeft(Nmax))
    !! Initiation 
    position=0.0d0; velocity=0.0d0; force=0.0d0; amplitude=0.0d0
    gmaeff=0.0d0; Minv=0.0d0; fext=0.0d0
    itm_p = 0.0d0; itm_v = 0.0d0; itm_f = 0.0d0; pdfz = 0.0d0
    !position(2,1)=0.0d0; position(3,1)=0.0d0
    !position(1,1)=1.0d0 !! Define initial height ∆(0)
    position(1,1) = ini_height; itm_p(1,2) = position(1,1)

    !rue_total  = "/Results/"//char(date(1))//char(date(2))//char(date(3))//"/"//&
    !&"kappa="//char(count_kappa)//"delta="//char(count_delta)
    write(no_kappa,*) count_kappa; no_kappa = ADJUSTL(no_kappa) !! Obtain the ordinal number of kappa as the string w/o spaces.
    write(no_delta,*) count_delta; no_delta = ADJUSTL(no_delta) !! Obtain the ordinal number of delta as the string w/o spaces.
    !write(*,*) "no_kappa = ", no_kappa
    !write(*,*) "no_delta = ", no_delta
    write(*,96) kappa, ini_height
    rue_total = "k_"//trim(no_kappa)//"-"//"d_"//trim(no_delta)//"-"
    write(*,*) "rue_total = ", rue_total

    !! Add the below path to each file
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
    open(unit=382,file=TRIM(rue_total)//"PDFz.txt"); write(382,96) kappa, ini_height
    write(31,61); write(32,62); write(33,63); write(34,64); !write(35,65)
    write(36,66); write(371,67); write(372,67); write(373,67); write(381,681); write(382,682)
    
    

    1101 continue
    do i = 1 , ratio_sim * (Nmax - 1)
        noise = normaldist(0.0d0,dt,3)
        !! itm_p(:,1) refers to the position at present; itm_p(:,2) refers to the next step value. Same for itm_v
        itm_p(:,1) = itm_p(:,2); itm_v(:,1) = itm_v(:,2) !! Get the last values for recursion
        !! Note, we change variables: velocity(:,i) -> itm_v(:,1); position(:,i) -> itm_p(:,1)
        !! Note, we change variables: velocity(:,i+1) -> itm_v(:,2); position(:,i+1) -> itm_p(:,2)
        !! Then we can advance towards much longer time with limited storage space...
        

        !! Update all elements in the effective friction matrix, and mass matrix
        do j=1,3
            do k=1,3
                do l=1,3
                    !gmaeff(j,k,l) = gammavalue(j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                    gmaeff(j,k,l) = gammavalue(j,k,l,itm_p(1,1),itm_v(1,1),itm_v(2,1),itm_v(3,1))
                end do
                !Minv(j,k) = Minverse(position(1,i),j,k)
                Minv(j,k) = Minverse(itm_p(1,1),j,k)
            end do
        end do
        amplitude(1,1) = Minv(1,1) * sqrt(4.0d0/beta * mz * (gmaeff(1,1,1) - gmaeff(2,1,1)) / (rayon * csquare * eps))
        amplitude(2,2) = Minv(2,2) * sqrt(2.0d0/beta * mx * (gmaeff(1,2,2) - gmaeff(2,2,2)) / (rayon * csquare))
        amplitude(3,3) = Minv(3,3) * sqrt(2.0d0/beta * mt * (gmaeff(1,3,3) - gmaeff(2,3,3)) * rayon / csquare)


        !force(i) = forcevalue(position(1,i))
        itm_f = forcevalue(itm_p(1,1))
        !fext(1) = forcevalue(position(1,i))
        !fext(1) = force(i)!; fext(2)=0.0d0; fext(3)=0.0d0
        fext(1) = itm_f


        !velocity(:,i+1) = velocity(:,i) + MATMUL(amplitude(:,:),noise(:,1)) + &
        !& dt * ( MATMUL(Minv(:,:),fext(:)) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),velocity(:,i)) )
        itm_v(:,2) = itm_v(:,1) + MATMUL(amplitude(:,:),noise(:,1)) + &
        & dt * ( MATMUL(Minv(:,:),fext(:)) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),itm_v(:,1)) )

        !position(:,i+1) = position(:,i) + velocity(:,i) * dt
        itm_p(:,2) = itm_p(:,1) + itm_v(:,1) * dt
        !position(1,i+1) = position(1,1)
        if (itm_p(1,2).le.0.0d0) then !! Extra treatment while the particle falls down.
            crash = 1
            write(30,*) "   # The simulation crashed once... Start again!"
            write(*,*) "   # The simulation crashed once... Start again!"
            goto 1105 !! End the current Euler loop directly since it drops.
        end if

        !! Only take part of data based on « ratio_sim »
        if (MOD(i,ratio_sim).eq.0) then
            count_sim = count_sim + 1 !! count_sim = 0 in the beginning.
            time_sim = i * dt
            position(:,count_sim + 1) = itm_p(:,2)
            velocity(:,count_sim + 1) = itm_v(:,2)
            force(count_sim) = itm_f

            j = ceiling(ratio_pdf * itm_p(1,2)) !! Determine the current height belongs to which interval
            pdfz(j) = pdfz(j) + 1

            !Dcoefx(i) = (gmaeff(1,2,2) - gmaeff(2,2,2))/(beta*mx*gmaeff(1,2,2)**2)
            Dcoefx(count_sim) = (1.0d0 - gmaeff(2,2,2)/gmaeff(1,2,2))/(beta*mx*gmaeff(1,2,2)) * (1.0d0/clight)**2
            !Dcoeft(i) = (gmaeff(1,3,3) - gmaeff(2,3,3))/(beta*mt*gmaeff(1,3,3)**2)
            Dcoeft(count_sim) = (1.0d0 - gmaeff(2,3,3)/gmaeff(1,3,3))/(beta*mt*gmaeff(1,3,3)) * (rayon/clight)**2

            !! We only keep two sub-output files for shorter running time.
            !write(31,71) time_sim,noise(1,1),noise(2,1),noise(3,1)
            write(32,72) time_sim,position(1,count_sim),position(2,count_sim),position(3,count_sim)
            write(33,73) time_sim,velocity(1,count_sim),velocity(2,count_sim),velocity(3,count_sim)
            !write(34,74) time_sim,force(count_sim),amplitude(1,1)*noise(1,1),amplitude(2,2)*noise(2,1),amplitude(3,3)*noise(3,1)
            !write(35,75) i*tratio,gamma,intma,intmb,Minv
            !write(36,76) time_sim,Minv(1,1),Minv(2,2),Minv(1,2),Minv(2,1),Minv(2,2)
            !write(371,77) time_sim,gmaeff(1,1,1),gmaeff(1,2,2),gmaeff(1,3,3),gmaeff(1,1,2),gmaeff(1,1,3),gmaeff(1,2,3)
            !write(372,77) time_sim,gmaeff(2,1,1),gmaeff(2,2,2),gmaeff(2,3,3),gmaeff(2,1,2),gmaeff(2,1,3),gmaeff(2,2,3)
            !write(373,77) time_sim,gmaeff(3,1,1),gmaeff(3,2,2),gmaeff(3,3,3),gmaeff(3,1,2),gmaeff(3,1,3),gmaeff(3,2,3)
        end if
    end do
    1102 continue
    call cpu_time(time_end)
    write(30,97) time_end-time_begin !! Print the time consumed for Euler method.


    call cpu_time(time_begin)
    do k=1,ratio_pdf
        itm_f = (2.0d0 * k - 1.0d0)/(2.0d0 * ratio_pdf) * rayon * eps !! Just calculate the average height for k-th data. Do not want to define a new varaible.....
        write(382,782) k, 1.0d0*(k-1)/ratio_pdf, 1.0d0*k/ratio_pdf, 1.0d0*pdfz(k)/Nmax, &
        !& 4.0d-7/rayon*exp(-itm_f/l_D) + beta*mz*rayon*grav*itm_f
        & exp(-beta*mz*(1.0d0-rhosol/rhosty)*rayon*grav*itm_f)/100.0d0
        !&-beta*mz*rayon*grav*itm_f*1.0d0 - 4.0d-7*exp(-itm_f/l_D)!*0.0d0 + itm_f
    end do
    call cpu_time(time_end)
    write(30,982) time_end-time_begin !! Print the time consumed for PDF output.
    

    call cpu_time(time_begin)
    1103 continue
    goto 1104
    !! MSD, Mean Square Displacement
    do i=1,dtmax !! ∆t?
        sumx(i) = 0.0d0; sumt(i) = 0.0d0
        do j=1,Nmax-dtmax
            msdx_num(j) = (position(2,j) - position(2,i+j))**2
            sumx(i) = sumx(i) + msdx_num(j)

            msdt_num(j) = (position(3,j) - position(3,i+j))**2
            sumt(i) = sumt(i) + msdt_num(j)
        end do

        msdx_ana = 1.0d0/(beta * mx * gmaeff(1,2,2)) * ((exp(-time_sim * gmaeff(1,2,2)))/gmaeff(1,2,2) + time_sim)
        msdt_ana = 1.0d0/(beta * mt * gmaeff(1,3,3)) * ((exp(-time_sim * gmaeff(1,3,3)))/gmaeff(1,3,3) + time_sim)
        msdc_ana = 0.0d0


        write(381,781) time_sim, msdx_ana, sumx(i)/j, Dcoefx(i)*i*1.0d0, msdt_ana, sumt(i)/j, Dcoeft(i)*i*1.0d0
    end do
    1104 continue
    call cpu_time(time_end)
    write(30,981) time_end-time_begin !! Print the time consumed for MSD calculations

    
    1105 continue !! Delete all intermediate data and close all sub- output files.
    deallocate(velocity); deallocate(position); deallocate(force); deallocate(pdfz)
    deallocate(msdx_num); deallocate(msdt_num); deallocate(sumx); deallocate(sumt); deallocate(Dcoefx); deallocate(Dcoeft)
    close(31); close(32); close(33); close(34); close(35); close(36); close(371); close(372); close(373); close(381); close(382)
    if (crash.eq.1) goto 1207 !! Restart the Euler loop until the particle always keeps a positive height.
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
    write(30,*); write(30,*); write(30,*) !! Blank rank for the output file.
    write(30,99) temps_fin-temps_debut !! Write the total time consumed into the output file.
    write(*,99) temps_fin-temps_debut !! Write the total time consumed to the screen.
    close(30) !! Close the total output file.

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
        !read(201,*) l_B !! Boltzmann length
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
    real*8 :: zroot,gabb!,gzzz,gzxx,gztt,
    real*8 :: coulombmax
    include "para.h"
    
    zroot = sqrt(z)
    !gzzz = (21.0d0*kxi)/(4.0d0*zroot*z**4)
    !gzxx = -(kxi)/(4.0d0*zroot*z**3)
    !!gztt = -(kxi)/(4.0d0*zroot*z**3)
    !gztt = gzxx
    gabb = kxi/(2.0d0*zroot*z**3) * (21.0d0/(2.0d0*z) - 1.0d0)
    coulombmax = 0.4d-6 * 1.0d0 !! namely the "B" shown by Maxime, PHYSICAL REVIEW RESEARCH 3, L032011 (2021)

    
    !! Attention to the unit! As for the force, it should be [kg·m·s^{-2}]
    !! However, since m = rho * π * r**2 = [kg/m], the force unit would be [kg/s^2]
    !! Finally, we take all parameters dimensionless.
    forcevalue = ( - mz*grav*(1.0d0-rhosol/rhosty) * rayon + &
    & 1.0d0/(beta * rayon) * ( gabb/eps + coulombmax/l_D * exp(-z*rayon*eps/l_D) ) ) * 2.0d0 / csquare
    !forcevalue =  - mz*grav*(1.0d0-rhosol/rhosty) * 2.0d0 * rayon / csquare + &
    !& 1.0d0/beta * (gzzz+gzxx+gztt)/(rayon*eps) * (2.0d0 * rayon) / (rayon * csquare) + &
    !& coulombmax/(rayon * beta * l_D) * exp(-z*rayon*eps/l_D) * 2.0d0 / csquare
    
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
