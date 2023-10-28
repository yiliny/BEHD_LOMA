!!!!! Yilin YE @ Jul. 2023
!!!!! EHD + Brownian motion
!!!!! Fokker-Planck method 
!!!!! UPDATE effective friction matrix
!!!!! modified noise correlator amplitude


Program main
    use MATHS
    implicit none    
    !integer :: Nmax = 10000 !! Define the maximum number of steps. Note ∆t=0.1ms,Nmax=1e5,T->100
    !integer :: dtmax = 8000 !! Define the maximum ∆t to calculate MSD
    integer i,j,k,l,count_kappa,count_delta,count_sim,count_msd,crash,count_particle !! Loop variables
    integer :: ratio_sim,ratio_msd,ratio_pdfz,pdft_max,ratio_pdfx,ratio_vz,ratio_vx,outp_max !! Constant
    real*8 :: time_begin,time_end,temps_debut,temps_fin,time_sim !! define variables to record time consumed
    real*8 :: zmax,xmax,pdf_trapz,pdf_trapx,vzmax,vxmax,pdf_trapvz,pdf_trapvx
    real*8 :: reff,aint,tvl,wi,bint !! 20230228 add variable to simplify calculations
    integer :: node !! 20230228 add for Gauss-Legendre integral
    real*8,external :: rectgauss,Minverse,gammavalue,spurious_force,repulsion,Laguerre,GLpoint,GBdist,GBvites,Mass_eff

    include "para.h"
    include "loop.h"
    !real*8,parameter :: pi=3.14159265358979d0 !! constant
    !complex*16,parameter :: Im=(0.d0,1.0d0)    !! imaginary unit
    character*128 ::  rue_total, no_kappa, no_delta !! Defin path name for outputs
    integer :: date(8)
    
    real*8 :: masse(3),amplitude(3,3),fext(3),vext !! Define matrices for Euler method, used in MATMUL()
    real*8 :: tratio !! Define the ratio for modified dimensionless time gap ∆T
    real*8 :: itm_q(3,2),itm_v(3,2),itm_f(2),itm_p(3,2),zsqrt!,itm_g(6) !! Define intermediates to facilitate the codes. 20230911, Position->itm_q & Momenta->itm_p
    real*8 :: msdz_num!,msdx_ana,msdx_num,msdt_num !! Define three analytical MSD, <∆r_x **2> & <∆r_Θ **2> & <∆r_x * ∆r_Θ>
    real*8,allocatable :: position(:,:),velocity(:,:),force(:),momenta(:,:) !! Define arrays for position/velocity/force
    !! Define array for PDF (Probability Distribution Function), count the number in the limited height
    integer,allocatable :: pdfzt(:,:),pdfxt(:,:),pdfvz(:,:),pdfvx(:,:)!,pdfz(:) !! pdfz for the single particle's time average. pdft for particle ensemble average at different time. 20228022 add pdfx. 20230301 remove pdfz, rename pdft -> pdfzt & pdfx -> pdfxt
    real*8,allocatable :: dltseq(:),zseq(:),Upot(:),Pzana(:) !! 20230228 add for G-B distribution
    real*8,allocatable :: vseq(:),Pvana(:) !! 20230301 add for velocity distribution
    real*8,allocatable :: sumz(:),sumx(:),sumt(:)!,Dcoefx(:),Dcoeft(:)!,msdx_num(:),msdt_num(:),;; 20220819 night add sumz; 20221003 sumt ~ <∆r_x * ∆r_Θ>
    
    
    call cpu_time(temps_debut)
    call cpu_time(time_begin)
    
    continue !! Here we list all format declarations.
        60 format('#',I5)
        61 format(2x,'Time',8x,'Position ∆',8x,'γ_z',8x,'Veloctiy ∆',8x,'Fz total',8x,'electro î',8x,'spurious',8x,'Noise z')
        62 format('#',2x,'Time',8x,'Position ∆',12x,'Position X',12x,'Angle Θ')
        !63 format('#',2x,'Time',8x,'Veloctiy ∆',15x,'Velocity X',15x,'Velocity Θ')
        64 format(2x,'Time',8x,'mass',12x,'Force z',12x,'Noise z',12x,'Noise x',12x,'Noise Θ')
        !!!!65 format(2x,'Time Step',12x,'γ coeff',12x,'a=exp(-γ.∆t)',12x,'b=sqrt(tanh(gt2)/gt2)',12x,'M inverse')
        !66 format(2x,'Time',12x,'M_zz',12x,'M_xx',12x,'M_xΘ',12x,'M_Θx',12x,'M_ΘΘ')
        !67 format(2x,'Time',10x,'γv_zz',10x,'γv_xx',10x,'γv_ΘΘ',10x,'γv_zx',10x,'γv_zΘ',10x,'γv_xΘ')
        !681 format('Time',9x,'sumx',9x,'<∆x^2>',9x,'ana.MSDx',9x,'sumt',9x,'<∆x∆Θ>',9x,'ana.MSD_xt',9x,'sumz',9x,'<∆z^2>')
        6812 format('Time',9x,'sumz',9x,'<∆z^2>')
        !682 format('# Index.. ',2x,'Height',12x,'Ratio',12x,'G-B theory')
        683 format('# Index.. ',2x,'Height',12x,'Particle Number')
        684 format('# Index.. ',2x,'Parallel Displacement',12x,'Particle Number')
        685 format(2x,'Height ∆',12x,'Potential U (kT)',12x,'G-B theory',12x,'P normalized')
        686 format('# Indexx.. ',2x,'Velocity',12x,'Particle Number')
        687 format(2x,'Speed',12x,'raw GB theory',12x,'P normalized')
    
        71 format(f11.5,8(3X,ES12.4)) !! format for file(31) = debug.txt
        72 format(f11.5,3(",",4X,ES14.4)) !! format for file(32) = position_z.txt
        !73 format(f11.5,3(5X,ES14.4)) !! format for file(33) = velocity_z.txt
        74 format(f11.5,5(5X,ES15.6)) !! format for file(34) = force_z.txt
        !!!!75 format(f10.1,4(5X,ES20.10)) !! format for file(35) = intermediates.txt
        !76 format(f16.9,5(5X,ES15.6)) !! format for file(36) = Mass_Matrix_Inverse.txt
        !77 format(f8.4,6(5X,ES12.4)) !! format for file(371/372/373)
        !781 format(f10.3,8(2X,ES13.5)) !! format for file(381) = MSD.txt
        7812 format(f15.5,2(2X,ES13.5))
        !782 format(f12.6,",",4X,f12.7,",",4X,ES15.6) !! format for file(382) = PDFz.txt
        785 format(f12.6,",",3(10X,ES15.7,",")) !! format for file(385) = GBhdist.txt
        787 format(f12.6,",",2(10X,ES12.5,","))
    
        94 format('BEHD @ YYE. Sim. done on ',I4,2('-',I2),'   ',I2,'h',I2)
        95 format('It spent',f8.2,' second for the whole initiation.')
        961 format('### Sim. w/ κ =',f9.3,' ; and ∆_0 =',f6.2,'; T =',ES12.2)
        962 format(' ## r =',ES10.2,' ; rhosty =',ES10.2,' ; rhosol =',ES10.2)
        963 format(' ## k =',f6.2,' ; ∆_eq =',f6.2,' ; num_particle =',I8)
        964 format(' ## Nmax =',I8,' ; dtmax =',I8)
        !97 format(f8.2,' seconds spent for this Euler loop.')
        97 format('    This Euler loop spent',f8.2,' second.')
        !980 format(f8.2,' seconds spent for this MSD calculations.')
        981 format(f8.2,' second spent for this MSD output.')
        982 format(f8.2,' second spent for this PDF output.')
        99 format('It spent',f10.2,' second for the whole program.')
    continue
    

    i=1 !! Suppose i the index to select input file. 1, input.txt; 2, input_real.txt.
    call obtaindata(i) !! Read all necessary parameters.
    !goto 1106

    continue !! Constant Parameters
        rhoa = (1.0d0-rhosol/rhosty)
        clight = sqrt(2.d0*grav*rayon*rhoa) !! unit m/s, the maximum speed, 1.3e-3
        csquare = clight**2

        mass = pi*rayon**2*rhosty !! unit kg/m
        mz = mass; mx = mass; mt = mass*rayon**2/2; masse(1)=mz; masse(2)=mx; masse(3)=mt
        
        if (i.eq.1) then
            eps = 0.1d0 !! dimensionless parameter, delta = ∆ * r * eps
            xi = 1.0d0 !! dimensionless parameter
        end if
        reff = rayon * eps

        k_B = 1.38064852d-23 !! unit J/K
        !beta = 1.0d0/(k_B*temp) !! unit J**-1 = kg**-1 * m**-2 * s**2
        beta = (mass*rhoa*grav*reff)/(k_B*temp) !! unit = m**-1
        
        !! Time gap for each step, t = T * r * sqrt(2*eps) / clight. Here we pose dt ~ ∆T
        tratio = 1.0d0 / (20.0d0*10) !! Define the time discretisation ratio
        dt = 0.50d-4 * clight / (rayon * sqrt(2.0d0 * eps)) !* tratio

        !! Define loop parameters
        ratio_sim = 1; ratio_msd = 400 !! Here are two variables to accelerate calculations for Euler method & MSD.
        zmax = 3.0d0 !! Define the maximum height considered for PDF
        xmax = 1.0d0 !! 20220822 added. Define the maximum absolute displacement considered for PDF along x.
        vzmax = 1.2d0; vxmax = 0.9d0 !! 20230301 added for velocity distribution
        ratio_pdfz = ceiling(zmax*50) !! Here is the ratio to determine how many data would be taken part for PDF calculation.
        ratio_pdfx = ceiling(2*xmax*50)
        ratio_vz = ceiling(vzmax*50); ratio_vx = ceiling(vxmax*50) !! 20230301 added for velocity distribution
        pdft_max = Nmax !10000 !! Define the maximum time section number. We ignore the further development after this value.
        !! Below two parameters are moved inside 'loop.h'
        !num_particle = 200 !! Define the number of particles. We would repeat the same simulation to obtain the average property.
        !zfix = 2 !! Determine whether we fix the vertical coordinate. 0, fix at the initial height all the time; 1, not fixed with constant initial height; 2, not fixed with random height
        outp_max = 20 !! Max particle number in position.txt file output
    continue
    
    
    call date_and_time(values=date)
    open(unit=30,file="output_record.txt")
    call cpu_time(time_end) !! Record the end time for the whole initiation.
    write(30,94) date(1), date(2), date(3), date(5), date(6); write(30,*); write(30,*); write(30,*)
    write(30,95) time_end-time_begin !! It spent ? seconds for the whole initiation.
    write(30,*)
    !write(30,*) "Below we simulate the Brownian motion of a cylindrical particle with parameters:"
    !write(30,"(ES10.3,ES10.3)") "Radius = ",rayon," (m);   ","Temperature = ",temp," (K);   "
    !write(30,*) "Debye length = ",l_D," (m);","B =",B,"; L_typ = ",Lchrt
    !write(30,"(ES10.3,ES10.3)") "Min. potential height = ",q0,"Force k = ",kk
    write(30,*); write(30,*)
    !write(*,*) 'clight = ',clight

    
    open(unit=385,file="GBhdist.txt"); write(385,94) date(1), date(2), date(3), date(5), date(6); write(385,685)
    allocate(Upot(ratio_pdfz)); allocate(zseq(ratio_pdfz)); allocate(dltseq(ratio_pdfz)); allocate(Pzana(ratio_pdfz)) !! Add 20230228. Analytical prediction for G-B distribution
    do k=1,ratio_pdfz !! Add 20230228: compute GB distribution at eqlbm
        dltseq(k) = zmax*(k - 0.5d0) / ratio_pdfz ! Dimensionless variable
        zseq(k) = dltseq(k) !* reff !zmax*(k - 0.5d0)/ratio_pdf
        !Upot(k) = GBdist(zseq(k),0) / temp!* beta ! Dimensionless potential.
        Upot(k) = GBdist(zseq(k),0) * beta
        Pzana(k) = exp(-Upot(k))
    end do
    aint = 0.0d0; node=80 !! Add 20230228: Integration for normalized factor A of PDF
    j = node / 2
    do k=1,node
        tvl=Laguerre(node,1,k); wi=Laguerre(node,2,k)
        !if (i.le.j) then
        !    tvl = GLpoint(node,1,i); wi = GLpoint(node,2,i)
        !else
        !    tvl = -GLpoint(node,1,i-j); wi = GLpoint(node,2,i-j)
        !end if
        
        !!!aint = aint + wi * GBdist(tvl*reff,1) * exp(tvl) * reff !! \int P(z) dz = \int P(∆*reff) d∆ * reff = reff * \int P(∆*reff)*exp(±∆) d∆ 
        !!!aint = aint + wi * exp(-GBdist(tvl*reff,0) / temp) * exp(tvl) * reff
        
        aint = aint + wi * exp(-GBdist(tvl,0) * beta) * exp(tvl) * reff
        !tvl = (tvl + 1.0d0) * zmax / 2.0d0
        !aint = aint + wi * exp(-GBdist(tvl,0) * beta) * reff * zmax / 2.0d0
    end do
    write(385,*) " ## A = ",aint
    write(*,*) " ## A = ",aint
    do k=1,ratio_pdfz !! Output GB distribution prediction
        write(385,785) dltseq(k), Upot(k), Pzana(k), Pzana(k) / aint
        !Pzana(k) = Pzana(k) / aint
    end do
    close(385)
    bint = 0.0d0 !! Add 20230913: Integration for normalized factor A of PDF
    do k=1,node
        tvl=Laguerre(node,1,k); wi=Laguerre(node,2,k)
        !aint = aint + wi * GBdist(tvl*reff,1) * exp(tvl) * reff !! \int P(z) dz = \int P(∆*reff) d∆ * reff = reff * \int P(∆*reff)*exp(±∆) d∆ 
        !aint = aint + wi * exp(-GBdist(tvl*reff,0) / temp) * exp(tvl) * reff
        bint = bint + wi * exp(-GBdist(tvl,0) * beta) * exp(tvl) * reff / aint
    end do
    write(*,*) " ## B = ",bint
    
    !! Add 20230301 for velocity distribution
    open(unit=387,file="GBvdist.txt"); write(387,94) date(1), date(2), date(3), date(5), date(6); write(387,687)
    allocate(vseq(ratio_vz)); allocate(Pvana(ratio_vz))
    do k=1,ratio_vz
        vseq(k) = (k-0.5d0)*vzmax/ratio_vz
        Pvana(k) = GBvites(vseq(k))
    end do
    aint = 0.0d0; node=80
    do k=1,node
        tvl=Laguerre(node,1,k); wi=Laguerre(node,2,k)
        aint = aint + wi * GBvites(tvl) * exp(tvl) * clight !! \int P(v) dv = \int P(V*c) dV * c = c * \int P(V*reff)*exp(±V) dV 
        !write(*,*) tvl, GBdist(tvl*reff,1), wi, aint
    end do
    write(387,*) " ## A = ",aint
    do k=1,ratio_vz !! Output GB distribution prediction
        write(387,787) vseq(k), Pvana(k), Pvana(k) / aint
        !Pvana(k) = Pvana(k) / aint
    end do
    close(387)

    

    count_kappa = 0 !! Reset. Start from the minimum value
    1201 continue !! Loop label
    select case(loop_kappa) !! Determine if loop or not for kappa
    case(0) !! No loop, then take a constant
        kappa = min_kappa !1.6d-1 * 0.0d0
    case(1) !! Loop
        kappa = min_kappa + gap_kappa * count_kappa
    end select
    kxi = kappa*xi; kxe = kappa*xi*eps !! non-dimensional parameters w/o units.
    
        count_delta = 0 !! Reset. Start from the minimum value
        1202 continue !! Loop label
        select case(loop_delta) !! Determine if loop of not for ∆
        case(0) !! No loop, then take a constant
            ini_height = min_delta !0.15d1
        case(1) !! Loop
            ini_height = min_delta + gap_delta * count_delta
            !if ((ini_height.eq.0.2d0).or.(ini_height.eq.0.4e0).or.(ini_height.eq.0.5d0).or.(ini_height.eq.0.6d0)&
            !&.or.(ini_height.eq.0.7d0).or.(ini_height.eq.0.8d0).or.(ini_height.eq.0.9d0)) then
            !    count_delta = count_delta + 1
            !    goto 1202
        end select


            !rue_total  = "/Results/"//char(date(1))//char(date(2))//char(date(3))//"/"//&
            !&"kappa="//char(count_kappa)//"delta="//char(count_delta)
            write(no_kappa,*) count_kappa; no_kappa = ADJUSTL(no_kappa) !! Obtain the ordinal number of kappa as the string w/o spaces.
            write(no_delta,*) count_delta; no_delta = ADJUSTL(no_delta) !! Obtain the ordinal number of delta as the string w/o spaces.
            !write(*,*) "no_kappa = ", no_kappa
            !write(*,*) "no_delta = ", no_delta
            write(*,961) kappa, ini_height, temp; write(*,962) rayon, rhosty, rhosol
            write(*,963) kk, q0, num_particle; write(*,964) Nmax, dtmax
            rue_total = "k_"//trim(no_kappa)//"-"//"d_"//trim(no_delta)//"-"
            write(*,*) "rue_total = ", rue_total

            !! Add the below path to each file
            open(unit=31,file=TRIM(rue_total)//"debug.txt"); write(31,94) date(1), date(2), date(3), date(5), date(6)
            write(31,961) kappa, ini_height, temp; write(31,962) rayon, rhosty, rhosol; write(31,963) kk, q0, num_particle
            open(unit=32,file=TRIM(rue_total)//"positions.txt"); write(32,94) date(1), date(2), date(3), date(5), date(6)
            write(32,961) kappa, ini_height, temp; write(32,962) rayon, rhosty, rhosol; write(32,963) kk, q0, num_particle
            !open(unit=33,file=TRIM(rue_total)//"velocitys.txt"); write(33,94) date(1), date(2), date(3), date(5), date(6)
            !write(33,961) kappa, ini_height, temp
            open(unit=34,file=TRIM(rue_total)//"forces.txt"); write(34,94) date(1), date(2), date(3), date(5), date(6)
            write(34,961) kappa, ini_height, temp; write(34,962) rayon, rhosty, rhosol; write(34,963) kk, q0, num_particle
            !open(unit=35,file="inter_Euler_Maruyama.txt")
            !open(unit=36,file=TRIM(rue_total)//"Mass_Matrix_Inverse.txt"); write(36,961) kappa, ini_height, temp
            !open(unit=371,file=TRIM(rue_total)//"gamma0.txt"); write(371,961) kappa, ini_height, temp
            !open(unit=372,file=TRIM(rue_total)//"gamma1.txt"); write(372,961) kappa, ini_height, temp
            !open(unit=373,file=TRIM(rue_total)//"gamma1v.txt"); write(373,961) kappa, ini_height, temp
            open(unit=381,file=TRIM(rue_total)//"MSD.txt"); write(381,94) date(1), date(2), date(3), date(5), date(6)
            write(381,961) kappa, ini_height, temp; write(381,962) rayon, rhosty, rhosol; write(381,963) kk, q0, num_particle
            !open(unit=382,file=TRIM(rue_total)//"PDFz.txt"); write(382,94) date(1), date(2), date(3), date(5), date(6)
            !write(382,961) kappa, ini_height, temp
            open(unit=383,file=TRIM(rue_total)//"PDF_ztime.txt"); write(383,94) date(1), date(2), date(3), date(5), date(6)
            write(383,961) kappa, ini_height, temp; write(383,962) rayon, rhosty, rhosol; write(383,963) kk, q0, num_particle
            open(unit=384,file=TRIM(rue_total)//"PDF_xtime.txt"); write(384,94) date(1), date(2), date(3), date(5), date(6)
            write(384,961) kappa, ini_height, temp; write(384,962) rayon, rhosty, rhosol; write(384,963) kk, q0, num_particle
            open(unit=3861,file=TRIM(rue_total)//"PDF_vzt.txt"); write(3861,94) date(1), date(2), date(3), date(5), date(6)
            write(3861,961) kappa, ini_height, temp; write(3861,962) rayon, rhosty, rhosol; write(3861,963) kk, q0, num_particle
            open(unit=3862,file=TRIM(rue_total)//"PDF_vxt.txt"); write(3862,94) date(1), date(2), date(3), date(5), date(6)
            write(3862,961) kappa, ini_height, temp; write(3862,962) rayon, rhosty, rhosol; write(3862,963) kk, q0, num_particle
            !! Furnish Variable Description.
            write(31,61)
            write(32,62); write(34,64); write(381,6812)!; write(382,682); 
            write(383,683); write(384,684)
            write(3861,686); write(3862,686)
    

            !! Decalre this allocatable variables in advance, for ensemble average of many particles
            allocate(pdfzt(ratio_pdfz,pdft_max)) !! ratio_pdf defines the number of interval. 200 refers to the time section number.
            
            allocate(pdfxt(ratio_pdfx,pdft_max)) !! 20220822 added. For ratio_x = ceiling(2*xmax*50), 2 refers to ± two directions; xmax, the max absolute displacement; 50, interval number for unit length.
            allocate(pdfvz(ratio_vz,pdft_max)); allocate(pdfvx(ratio_vx,pdft_max)) !! 20230301 added for velocity distribution
            allocate(sumx(dtmax)); allocate(sumt(dtmax)); allocate(sumz(dtmax))

            
            
            write(30,961) kappa, ini_height, temp; write(30,962) rayon, rhosty, rhosol; write(30,963) kk, q0, num_particle
            write(30,*) ""
            !! Here we start LOOP from the 1st particle.
            count_particle = 1 
            pdfxt = 0; pdfzt = 0 !! Must put it here, outside the particle loop.
            pdfvx = 0; pdfvz = 0 !! 20230301 added.
            sumx = 0.0d0; sumt = 0.0d0; sumz = 0.0d0 !! 20220819 night add sumz
            1203 continue !! Loop label for different particles with the same initial conditions.

    
            call cpu_time(time_begin)
            !write(30,961) kappa, ini_height, temp
            !write(30,*) "     for the particle",count_particle

            1207 continue !! Restart if the simulation crash (the particle drops into the surface)
            count_sim = 0; count_msd = 0 !! Here are two variables to record the number of big loops for Euler method.
            crash = 0 !! Define the variable to determine whether the particle drops into the surface

            !! Initiation 
            allocate(position(3,Nmax+1)); allocate(velocity(3,Nmax+1)); allocate(force(Nmax+1)); allocate(momenta(3,Nmax+1))
            position=0.0d0; velocity=0.0d0; force=0.0d0; amplitude=0.0d0; momenta=0.0d0
            !allocate(pdfz(ratio_pdf)); pdfz = 0
            
            !allocate(Dcoefx(Nmax)); allocate(Dcoeft(Nmax))
            gmaeff=0.0d0; Minv=0.0d0; fext=0.0d0
            itm_q = 0.0d0; itm_v = 0.0d0; itm_f = 0.0d0; itm_v(1,2) = -0.0d0; itm_p = 0.0d0
            if (zfix.eq.2) then
                noise = normaldist(ini_height,sqrt(ini_height)/3.0d0,3)
                !itm_q(1,1) = noise(1,1)
                !write(*,*) "particle=",count_particle,itm_q(1,1)
                position(1,1) = abs(noise(1,1))
            else !! (zfix.eq.0) or (zfix.eq.1)
                position(1,1) = ini_height
            end if
            itm_q(1,2) = position(1,1)

            
            if (count_particle.le.outp_max) then
                write(31,*) ""; write(31,*) ""; write(32,*) ""; write(32,*) ""!; write(33,*) ""; write(33,*) ""
                write(32,60) count_particle; !write(33,60) count_particle
            end if
            !write(31,60); write(34,60); write(36,60); write(371,60); write(372,60); write(373,60)
            
    

            !1101 continue
            do i = 1 , ratio_sim * (Nmax - 1)
                noise = normaldist(0.0d0,dt,3)
                !! itm_q(:,1) refers to the position at present; itm_q(:,2) refers to the next step value. Same for itm_v
                itm_q(:,1) = itm_q(:,2); itm_v(:,1) = itm_v(:,2) !! Get the last values for recursion
                !itm_p(:,1) = itm_p(:,2)
                !! Note, we change variables: velocity(:,i) -> itm_v(:,1); position(:,i) -> itm_q(:,1)
                !! Note, we change variables: velocity(:,i+1) -> itm_v(:,2); position(:,i+1) -> itm_q(:,2)
                !! Then we can advance towards much longer time with limited storage space...
        

                !! Update all elements in the effective friction matrix, and mass matrix inverse
                do j=1,1!3
                    do k=1,1!3
                        do l=1,1!3
                            !gmaeff(j,k,l) = gammavalue(j,k,l,position(1,i),velocity(1,i),velocity(2,i),velocity(3,i))
                            gmaeff(j,k,l) = gammavalue(j,k,l,itm_q(1,1),itm_v(1,1))
                        end do
                        !Minv(j,k) = Minverse(position(1,i),j,k)
                        Minv(j,k) = Minverse(itm_q(1,1),j,k)
                    end do
                end do
                zsqrt = sqrt(itm_q(1,1)) !! 20230911 added
                !! 11/09/23 modified
                !amplitude(1,1)=sqrt(2.0d0/beta*(xi/(zsqrt*itm_q(1,1))+(63.0d0*kxi*itm_v(1,1))/(32.0d0*zsqrt*itm_q(1,1)**4)))*1.0d10
                amplitude(1,1)=sqrt((2.0d0*xi)/(beta*zsqrt*itm_q(1,1))*(1.0d0+(63.0d0*kappa*itm_v(1,1))/(32.0d0*itm_q(1,1)**3)))&
                &*sqrt(csquare/(reff))
                !!amplitude(1,1)=Minv(1,1)*sqrt(4.0d0/beta*mz*(gmaeff(1,1,1)) / (Lchrt * csquare * 1.3d0 * eps**2))
                !amplitude(2,2)=0.0d0
                !!amplitude(2,2)=Minv(2,2)*sqrt(4.0d0/beta*mx*(gmaeff(1,2,2)) / (Lchrt * csquare * eps))
                !amplitude(3,3)=0.0d0
                !amplitude(3,3)=Minv(3,3)*sqrt(4.0d0/beta*mt*(gmaeff(1,3,3)) / (Lchrt * csquare * eps))
                !if (fnoise.eq.1) then ! Cross correlation, DNF
                !    amplitude(1,2)=0.0d0
                !    amplitude(2,3)=0.0d0
                !    amplitude(3,1)=0.0d0
                !end if


                !force(i) = forcevalue(position(1,i))
                !itm_f(1) = repulsion(itm_q(1,1))
                !itm_f(2) = spurious_force(itm_q(1,1)) 
                !fext(1) = (2.0d0 - itm_q(1,1)) * 5.0d-2 !fext(2)=0.0d0; fext(3)=0.0d0
                fext(1) = repulsion(itm_q(1,1))


                !velocity(:,i+1) = velocity(:,i) + MATMUL(amplitude(:,:),noise(:,1)) + &
                !& dt * ( MATMUL(Minv(:,:),fext(:)) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),velocity(:,i)) )
                !itm_p(1,2) = itm_p(1,1) + amplitude(1,1) * noise(1,1) + dt * &
                !& 2.5d0*itm_q(1,1)*zsqrt*kxi*(84.0d0*itm_q(1,1)*(itm_p(1,1)/(8.0d0*itm_q(1,1)**3*zsqrt-15.0d0*kxi))**2) &
                !& + dt * fext(1) !+ external force?
                !itm_v(:,2) = MATMUL(Minv(:,:),itm_p(:,2)) !! p = M v
                
                !itm_v(:,2) = itm_v(:,1) + MATMUL(amplitude(:,:),noise(:,1)) + &
                !& dt * ( MATMUL(Minv(:,:),fext(:)) - MATMUL((gmaeff(1,:,:)+gmaeff(2,:,:)+gmaeff(3,:,:)),itm_v(:,1)) )
                itm_v(1,2) = itm_v(1,1) + 1.0d0 / Mass_eff(itm_q(1,1)) * ( amplitude(1,1)*noise(1,1) + dt * &
                & (fext(1) - (1.0d0+(21.0d0*kappa*itm_v(1,1))/(4.0d0 * itm_q(1,1)**3)) * (xi * itm_v(1,1))/(zsqrt * itm_q(1,1)) &
                & - (21.0d0*kxi)/(16.0d0 * beta * Mass_eff(itm_q(1,1)) * zsqrt * itm_q(1,1)**4) )) 
                !vext = spurious_velocity(itm_q(1,1)) * clight
                if (zfix.eq.0) itm_v(1,2) = 0.0d0 !! Fix v_∆ as 0 all the time if necessary.
                !if (xfix.eq.0) itm_v(2,2) = 0.0d0
                !if (tfix.eq.0) itm_v(3,2) = 0.0d0 !! Fix v_Θ as 0 all the time if necessary. Added on 20220823 to avoid affect of cross terms in gamma_1

                !position(:,i+1) = position(:,i) + velocity(:,i) * dt
                !itm_q(:,2) = itm_q(:,1) + itm_v(:,1) * dt
                itm_q(1,2) = itm_q(1,1) + itm_v(1,2) * dt
                !if (zfix.eq.0) itm_q(1,2) = itm_q(1,1) !! Fix ∆ all the time if necessary.
                !if (xfix.eq.0) itm_q(2,2) = itm_q(2,1)
                !if (tfix.eq.0) itm_q(3,2) = itm_q(3,1) !! Fix Θ all the time if necessary. Added on 20220823 same purpose above.
                if (itm_q(1,2).gt.0.0d0) then !! Extra treatment while the particle falls down.
                    continue !! 20220821凌晨重构，变判断为否定，避免随机到低海拔初值迅速坠落来不及小于零就NaN
                else
                    crash = 1
                    write(30,*) "   # Sim. failed",count_particle
                    write(*,*) "   # Sim. failed",count_particle
                    goto 1105 !! End the current Euler loop directly since it drops.
                end if

                !! Only take part of data based on « ratio_sim »
                !if (MOD(i,ratio_sim).eq.0) then
                continue
                    count_sim = count_sim + 1 !! count_sim = 0 in the beginning.
                    time_sim = i * dt
                    position(:,count_sim + 1) = itm_q(:,2)
                    velocity(:,count_sim + 1) = itm_v(:,2)
                    !momenta(:,count_sim + 1)  = itm_p(:,2)
                    force(count_sim) = fext(1) !itm_f
                    !write(*,*) itm_f

                    !! Calcul PDF_z
                    j = ceiling(ratio_pdfz * itm_q(1,2) / zmax) !! Determine the current height belongs to which interval
                    if (j.ge.ratio_pdfz) j = ratio_pdfz !! Avoid segmentation fault. j_max should not exceed ratio_pdf.
                    if (j.le.1) j = 1
                    !pdfz(j) = pdfz(j) + 1
                    !! Calcul PDF_t
                    if (count_sim.le.pdft_max) pdfzt(j,count_sim) = pdfzt(j,count_sim) + 1

                    !! Calculate PDF_x time-depending distribution
                    j = ceiling(ratio_pdfx * (itm_q(2,2) + xmax) * 0.5d0 / xmax) !! itm_q = - xmax, j = ratio_x
                    !if (j.le.1) j = 1
                    !if (j.ge.ratio_x) j = ratio_x
                    if ((j.ge.1).and.(j.le.ratio_pdfx)) then
                        if (count_sim.le.pdft_max) pdfxt(j,count_sim) = pdfxt(j,count_sim) + 1
                    end if

                    !! 20230301 added according to pdfzt. 
                    ! Calculate vz time-depending distribution
                    j = ceiling(ratio_vz * abs(itm_v(1,2)) / vzmax)
                    if (j.ge.ratio_vz) j = ratio_vz
                    if (j.le.1) j = 1
                    if (count_sim.le.pdft_max) pdfvz(j,count_sim) = pdfvz(j,count_sim) + 1
                    ! Calculate vx time-depending distribution
                    !j = ceiling(ratio_vx * abs(itm_v(2,2)) / vxmax)
                    !if (j.ge.ratio_vx) j = ratio_vx
                    !if (j.le.1) j = 1
                    !if (count_sim.le.pdft_max) pdfvx(j,count_sim) = pdfvx(j,count_sim) + 1


                    !! We only keep two sub-output files for shorter running time.
                    !write(31,71) time_sim,noise(1,1),noise(2,1),noise(3,1)
                    !write(32,72) time_sim,position(1,count_sim),position(2,count_sim),position(3,count_sim)
                    if (count_particle.le.outp_max) then
                        write(31,71) time_sim,position(1,count_sim),gmaeff(1,1,1),velocity(1,count_sim),force(count_sim),&
                        &itm_f(1),vext,amplitude(1,1)*noise(1,1)
                        write(32,72) time_sim,position(1,count_sim)!,position(2,count_sim),position(3,count_sim)
                        !write(33,73) time_sim,velocity(1,count_sim)!,velocity(2,count_sim),velocity(3,count_sim)
                        write(34,74) time_sim,mass,force(count_sim),amplitude(1,1)*noise(1,1),&
                        &amplitude(2,2)*noise(2,1),amplitude(3,3)*noise(3,1)
                    end if
                    !write(36,76) time_sim,Minv(1,1),Minv(2,2),Minv(1,2),Minv(2,1),Minv(2,2)
                continue
                !end if
            end do

            call cpu_time(time_end)
            !! Attention, the following expression would lead to some unexpected errors if no right MOD result.
            !! For example, if we set num_particle = 1, there would always be the error.
            !! Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.
            if (MOD(count_particle,ceiling(num_particle/6.0d0)).eq.0) then
                write(30,*) "   # Sim. done, particle",count_particle
                write(*,*) "   # Sim. done, particle",count_particle

                write(30,97) time_end-time_begin !! 20230301 move to here
                write(*,97) time_end-time_begin !! Print the time consumed for Euler method.
            end if
            !! Then we change to a better expression:
            !if
            !end if
            !1102 continue



            !! PDF_z calculs !! 20230301 add goto and abandon old pdfz
            goto 1109 !! Here we try to furnish theoretical expectation of Gibbs-Boltzmann distribution.
                call cpu_time(time_begin)
                write(382,*) ""; write(382,*) "" !! Insert two lines only with whitespace into the output file.
                !! 
                if (count_particle.le.0) then !! This sentence was written before inside the sim loop. To output all data, we comment that directly.
                    write(382,60) count_particle
                    do k=1,ratio_pdfz
                        itm_f = zmax*(k - 0.5d0)/(ratio_pdfz) * rayon * eps !! Just calculate the average height for k-th data. Do not want to define a new varaible.....
                        !! 20230301 comment the below setence
                        !write(382,782) dltseq(k), 1.0d0*pdfz(k)/Nmax, exp(-beta*mz*rhoa*rayon*grav*itm_f)/100.0d0
                        !& 4.0d-7/rayon*exp(-itm_f/l_D) + beta*mz*rayon*grav*itm_f
                        !&-beta*mz*rayon*grav*itm_f*1.0d0 - 4.0d-7*exp(-itm_f/l_D)!*0.0d0 + itm_f
                    end do
                end if
                write(30,*) "   # PDFz done for particle",count_particle
                !write(*,*) "   # PDFz done for particle",count_particle
                call cpu_time(time_end)
                write(30,982) time_end-time_begin !! Print the time consumed for PDF_z calcul/output.
            1109 continue
    

            !! MSD x/t calculs & MSD z added on 20220819 night
            continue !1103 continue
            !call cpu_time(time_begin)
            !write(381,*) ""
            !write(381,60) count_particle
            !goto 1104 !! MSD, Mean Square Displacement: <∆r**2(∆t)> = <∆r(t) ∆r(t+∆t)> = sum{∆r(t_i) ∆r(t_i + ∆t)}/N
                do i=1,dtmax !! Loop for ∆t
                    if (i.ge.Nmax) goto 1104 !! Avoid t_i + ∆t > t(Nmax)
                    if ((MOD(i,ratio_msd).eq.0).or.(i.le.25)) then !! Accelerate the MSD calculs.
                        !sumx(i) = 0.0d0; sumt(i) = 0.0d0
                        !time_sim = i * dt * ratio_sim
                        do j=1,Nmax-i !! Loop for t_i, t_i + ∆t < t(Nmax)
                            !msdx_num = (position(2,j) - position(2,i+j))**2
                            !sumx(i) = sumx(i) + msdx_num
                            !msdt_num = (position(3,j) - position(3,i+j)) * (position(2,j) - position(2,i+j))
                            !sumt(i) = sumt(i) + msdt_num
                            msdz_num = (position(1,j) - position(1,i+j))**2
                            sumz(i) = sumz(i) + msdz_num
                        end do
                        !write(381,7812) time_sim, sumz(l), sumz(l)/(num_particle*(Nmax-l))
                    end if
                end do
                !write(30,*) "   # MSD done for particle",count_particle
            1104 continue
            !call cpu_time(time_end)
            !write(30,980) time_end-time_begin !! Print the time consumed for MSD calculations


            1105 continue
            !! Deallocate most of arrays.
            deallocate(velocity); deallocate(position); deallocate(force); deallocate(momenta) !deallocate(pdfz)
            !! Restart the Euler loop until the particle always keeps a positive height.
            if (crash.eq.1) goto 1207 !! Return to the beginning of loop w/o adding count_particle
            !1106 continue 
    
            
            !! Here we continue LOOP toward the next particle.
            if (count_particle.lt.num_particle) then
                count_particle = count_particle + 1
                goto 1203
            end if
            

            !! MSD output for ensemble average
            call cpu_time(time_begin)
            !write(381,*) ""
            !goto 1107
                do l=1,dtmax
                    if (l.ge.Nmax) goto 1107 !! No output if t_i + ∆t > t(Nmax)
                    if ((MOD(l,ratio_msd).eq.0).or.(l.le.25)) then !! Only output non-zero results.
                        time_sim = l * dt * ratio_sim
                        !msdx_ana=exp(-time_sim*itm_g(1))/(mx*itm_g(1)**3)*(itm_g(1) - 2.0d0*itm_g(3) - time_sim*itm_g(1)*itm_g(3) &
                        !& + exp(time_sim*itm_g(1)) * (2.0d0*itm_g(3) + itm_g(1) * (time_sim * (itm_g(1) - itm_g(3)) - 1.0d0))) &
                        !& * 2.0d0*eps**3/(beta*rayon*clight**2) !rayon**2 * sqrt(2.0d0)

                        !write(381,781) time_sim, sumx(l), sumx(l)/(num_particle*(Nmax-l)), msdx_ana, &
                        !&sumt(l), sumt(l)/(Nmax-l), l*1.0d0, sumz(l), sumz(l)/(num_particle*(Nmax-l))
                        write(381,7812) time_sim, sumz(l), sumz(l)/(num_particle*(Nmax-l))
                    end if
                end do
            1107 continue
            call cpu_time(time_end)
            write(30,981) time_end-time_begin
            

            !! PDF_time output for ensemble average, "num_particle" in total
            call cpu_time(time_begin)
            write(383,*) ""
            3831 format("#",I6)
            3832 format(f12.6,",",4X,ES13.5,",")
            pdf_trapz = num_particle * zmax * reff / ratio_pdfz
            !pdf_trapx = num_particle * xmax * rayon * sqrt(2.0d0*eps) / ratio_pdfx ! + 20230301
            pdf_trapvz = num_particle * vzmax / ratio_vz * clight !* sqrt(eps/2.0d0) ! + 20230301
            !pdf_trapvx = num_particle * vxmax / ratio_vx * clight  ! + 20230301
            !write(*,*) "V norm = ", pdf_trapvz, pdf_trapvx
            do l=1,pdft_max
                write(383,*) ""; write(383,*) ""
                write(383,3831) l !! Furnish INDEX
                do k=1,ratio_pdfz
                    write(383,3832,advance='yes') dltseq(k), pdfzt(k,l)/pdf_trapz !1.0d0*(pdft(k,l)-0.0d0)/num_particle 
                    !! Must add "1.0d0*" else all results = 0.
                end do
                !! Below Added on 20220822 for PDF along x
                !write(384,*) ""; write(384,*) ""
                !write(384,3831) l !! Same, furnish INDEX
                !do k=1,ratio_pdfx
                !    write(384,3832,advance='yes') xmax*(k - 0.5d0)/ratio_pdfx*2.0d0-xmax, pdfxt(k,l)/pdf_trapx !! X position/parallel displacement; particle number;
                !end do
                !! Added on 20230301 for velocity dist
                write(3861,*) ""; write(3861,*) ""; write(3861,3831) l
                do k=1,ratio_vz
                    write(3861,3832,advance='yes') vzmax*(k-0.5d0)/ratio_vz, pdfvz(k,l)/pdf_trapvz
                end do
                !write(3862,*) ""; write(3862,*) ""; write(3862,3831) l
                !do k=1,ratio_vx
                !    write(3862,3832,advance='yes') vxmax*(k-0.5d0)/ratio_vx, pdfvx(k,l)/pdf_trapvx
                !end do
            end do
            call cpu_time(time_end)
            write(30,982) time_end-time_begin !! Print the time consumed for PDF output.
            !1110 continue
            

            !! Delete all intermediate data and close all sub- output files.
            deallocate(sumx); deallocate(sumt); deallocate(sumz); deallocate(pdfzt); deallocate(pdfxt) !deallocate(Dcoefx); deallocate(Dcoeft)
            deallocate(pdfvx); deallocate(pdfvz) !! Add 20230301
            !close(33); cclose(35); close(36); close(371); close(372); close(373); 
            close(31); close(32); close(34)
            close(381); close(382); close(383); close(384); close(3861); close(3862)

            write(30,*); write(30,*)
            write(30,*) "***** ***** ***** ***** ***** ***** ***** ***** ***** *****"; write(30,*)
    

        if ((loop_delta.eq.1).and.(ini_height.lt.max_delta)) then
            count_delta = count_delta + 1
            goto 1202
        end if

    if ((loop_kappa.eq.1).and.(kappa.lt.max_kappa)) then
        count_kappa = count_kappa + 1
        goto 1201
    end if
    
    deallocate(Upot); deallocate(zseq); deallocate(dltseq); deallocate(Pzana) !! Add 20230228
    deallocate(vseq); deallocate(Pvana) !! Add 20230301
    

    
    write(30,*); write(30,*); write(30,*) !! Blank rank for the output file.
    write(30,*) "  ## All simulation done with",count_particle,"paticles totally." !! Write the total particle number consumed into the output file.
    write(*,*) "  ## All simulation done with",count_particle,"particles totally." !! Write the total particle number consumed to the screen.
    call cpu_time(temps_fin)
    write(30,99) temps_fin-temps_debut !! Write the total time consumed into the output file.
    write(*,99) temps_fin-temps_debut !! Write the total time consumed to the screen.
    !write(*,*) "  ## Mass = ", mass
    !write(*,*) "  ## c = ", clight
    !write(*,*) "  ## fext = ", fext(1)
    !write(*,*) "beta' = ", beta
    !write(*,*) "pdf_trapz = ", pdf_trapz
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
        ruein="input_fpe.txt"
        open(unit=201,file=TRIM(ruein),form='formatted',status='old',action='read',iostat=status1,iomsg=msg)
        read(201,*) grav !! gravity acceleration
        read(201,*) temp !! ambient temperature
        read(201,*) rayon !! particle radius
        read(201,*) rhosty !! particle density
        read(201,*) rhosol !! solvent density
        read(201,*) kk
        read(201,*) q0 !! V(q) = kk/2*(q-q0)**2
        
        read(201,*) Nmax !! max step number for each simulation
        read(201,*) dtmax !! max ∆t number for MSD calculs
        read(201,*) num_particle !! particle number for one condition
        read(201,*) zfix !! Determine whether to fix height values or not. 0, fix; 1, no fix with constant initial height; 2, no fix with Gauss distributed initial height
        read(201,*) xfix !! Determine whether to fix horizontal value or note. 0, fix; 1, no fix.
        read(201,*) tfix !! Determine whether to fix rotation angles or not. 0, fix; 1 no fix.

        !! 03/02/23 added
        read(201,*) fspur !! Determine which type of spurious force used. 0, zero spurious force; 1, Gamma_z; 2, Gamma_z/x/t
        read(201,*) vspur !! Determine if add spurious velocity like Maxime; Added 170323. 0, zero; 1, ∂D(z)/∂z
        read(201,*) fnoise !! Determine which type of random noise used. 0, kT gamma0; 1, kT (gamma0 - gamma1); 2, gamma1v


        read(201,*) loop_kappa; read(201,*) min_kappa; read(201,*) max_kappa; read(201,*) gap_kappa
        read(201,*) loop_delta; read(201,*) min_delta; read(201,*) max_delta; read(201,*) gap_delta
        close(201)
    case(2)
        ruein="input_real.txt"
        !! Space left for possible parameters in reality such as viscosity.
    case default
        ruein="not_exist"
    end select
end subroutine



!*** Here is the function to calculate the gamma (effective friction matrix)
real*8 function gammavalue(no,i,j,z,vz)
    implicit none
    include "para.h"
    integer :: no,i,j !! no: gamma No.?; i/j: matrix index
    real*8 :: z,vz
    real*8 :: zroot
    real*8,external :: Mass_eff
    zroot = sqrt(z)

    select case(no)
    case(1) !! gamma_0
        if (i.eq.j) then
            if (i.eq.1) gammavalue = xi/(z*zroot)
            if (i.eq.2) gammavalue = 0.0d0 !(2.0d0*xi*eps)/(3.0d0*zroot)
            if (i.eq.3) gammavalue = 0.0d0 !(4.0d0*xi*eps)/(3.0d0*zroot)
        else
            gammavalue = 0.0d0
        end if
    case(2) !! gamma_1
        if (i.eq.1) then
            gammavalue = (21.0d0*kxi*vz)/(4.0d0*zroot*z**4)
        else
            gammavalue = 0.0d0
        end if
    case(3) !! gamma_1v
        if ((i.eq.1).and.(j.eq.1)) then
            gammavalue = (63.0d0 * kxi * vz) / (32.0d0 * zroot * z**4 * Mass_eff(z)) 
        else
            gammavalue = 0.0d0
        end if
    end select
    return
end function gammavalue




!*** Here is the function to calculate the force depending on z/∆.
real*8 function spurious_force(z)
    !This function would furnish the force on each direction
    !z: vertical position; vz: velocity vector
    implicit none
    real*8 :: z
    real*8 :: zroot,gabb!,gzzz,gzxx,gztt,
    include "para.h"
    include "loop.h"
    
    zroot = sqrt(z)
    !gzzz = (21.0d0*kxi)/(4.0d0*zroot*z**4)
    !gzxx = -(kxi)/(4.0d0*zroot*z**3)
    !!gztt = -(kxi)/(4.0d0*zroot*z**3)
    !gztt = gzxx

    ! 03/02/23 added
    if (fspur.eq.0) then 
        gabb = 0.0d0
    else !! fspur.eq.1 or 2
        gabb = kxi/(2.0d0*zroot*z**3) * (21.0d0/(2.0d0*z) - 1.0d0)
    end if

    
    !! Attention to the unit! As for the force, it should be [kg·m·s^{-2}]
    !! However, since m = rho * π * r**2 = [kg/m], the force unit would be [kg/s^2]
    !! Finally, we take all parameters dimensionless.
    !spurious = ( 1.0d0/beta * gabb/(rayon*eps) ) * (2.0d0*rayon)/(Lchrt*csquare)
    spurious_force = 0.0d0 !(2.0d0*rayon*gabb)/(beta*rayon*eps*Lchrt*csquare)
    !forcevalue =  - mz*grav*(1.0d0-rhosol/rhosty) * 2.0d0 * rayon / csquare + &
    !& 1.0d0/beta * (gzzz+gzxx+gztt)/(rayon*eps) * (2.0d0 * rayon) / (rayon * csquare) + &
    !& coulombmax/(rayon * beta * l_D) * exp(-z*rayon*eps/l_D) * 2.0d0 / csquare
    
    return
end function spurious_force
!*** Here is the function to calculate electrostatic force
real*8 function repulsion(Delta)
    !This function would furnish the force on each direction
    !z: vertical position; vz: velocity vector
    implicit none
    real*8 :: Delta!,B
    include "para.h"
    !! namely the "B" shown by Maxime, PHYSICAL REVIEW RESEARCH 3, L032011 (2021)
    !! Attention to the unit! As for the force, it should be [kg·m·s^{-2}]
    !! However, since m = rho * π * r**2 = [kg/m], the force unit would be [kg/s^2]
    !! Finally, we take all parameters dimensionless.
    repulsion = kk * (q0 - Delta) - 1.0d0
    !repulsion = exp(-Delta*rayon*eps) * (2.0d0*rayon)/(beta*csquare)
    return
end function repulsion


!*** Here is the function to calculate inverse mass matrix
real*8 function Minverse(z,i,j)
    implicit none
    include "para.h"
    real*8 :: z,zroot
    integer :: i,j !! matrix index
    zroot = sqrt(z)

    select case(i)
    case(1) !! ∆
        if (j.eq.1) then
            !Minverse = (1.0d0/mz) * (1.0d0 + (15.0d0*kxi)/(8.0d0*zroot*z**3)) !! Reviewed 26/07/22, ∆**{7/2} rather than 5/2
            Minverse = 1.0d0 / (1.0d0 - (15.0d0 * kxi)/(8.0d0 * zroot * z**3))
        else
            Minverse = 0.0d0
        end if
    case(2) !! X
        if (j.eq.1) Minverse = 0.0d0
        if (j.eq.2) Minverse = 0.0d0 !(1.0d0/mx) * (1.0d0 + kxe/(12.0d0*zroot*z**2))
        if (j.eq.3) Minverse = 0.0d0 !-kxe/(12.0d0*zroot*z**2*mt)
    case(3) !! THETA
        if (j.eq.1) Minverse = 0.0d0
        if (j.eq.2) Minverse = 0.0d0 !-kxe/(6.0d0*zroot*z**2*mx)
        if (j.eq.3) Minverse = 0.0d0 !(1.0d0/mt) * (1.0d0 + kxe/(6.0d0*zroot*z**2))
    end select
    return
end function Minverse
!*** Furnish effective Mass directly
real*8 function Mass_eff(z)
    implicit none
    include "para.h"
    real*8 :: z, zroot
    zroot = sqrt(z)
    Mass_eff = 1.0d0 - (15.0d0*kxi)/(8.0d0*zroot*z**3)
    return
end function Mass_eff



!*** Here is the function to compute Gibbs-Boltzmann distribution analytically. 20230228
real*8 function GBdist(z,xw) ! z: vertical height, unit (m)
    implicit none
    integer :: xw !! xw: mode selection. 0 for Udist & 1 for GBdist/Pdist
    real*8 :: z,Udist !! z = real vertical height; Udist = potential;; masstest = m or m*?
    real*8,external :: Minverse
    include "para.h"
    !masstest = 1.0/Minverse(z/(rayon*eps),1,1)
    !Udist = B/beta*exp(-z/l_D) + masstest*rhoa*Lchrt*grav*z!/Lchrt
    Udist = 0.5d0 * kk * (z-q0)**2 + z !(z-2.0d0)**2 * 5.0d-9
    !!Udist = B/beta*exp(-z/l_D) + rhoa*Lchrt*grav*z/Minverse(z/(rayon*eps),1,1)
    if (xw.eq.0) GBdist = Udist
    if (xw.eq.1) GBdist = exp(-Udist*beta)
    return
end function GBdist


!*** Here is the function to compute Gibbs-Boltzmann distribution analytically. 20230228
real*8 function GBvites(speed) ! z: vertical height; 
    use MATHS
    implicit none
    real*8 :: speed !! z = real vertical height; Udist = potential;
    real*8,external :: Minverse
    include "para.h"
    !GBvites = sqrt(2.0d0*beta*mass/pi) * exp(-0.5d0 * beta*mass*speed**2)
    GBvites = sqrt(2.0d0*beta/pi) * exp(-0.5d0 * beta*speed**2)
    return
end function GBvites


!***** Gauss-Laguerre Integration ***** 20230228 added from ntest.f90, for one-side infinite integral.
real*8 function Laguerre(Nnode,xw,j)
    implicit none
    integer :: Nnode,xw,j,i !! Nnode判断几阶多项式；xw判断零点/系数；j为第几个变量；i为函数内变量
    real*8,allocatable :: x(:),w(:)
    !Nnode=60; xw=1对应x & xw=2对应w；
    i=Nnode
    allocate(x(i)); allocate(w(i))

    !! 直接列出Legendre多项式零点；
    select case(Nnode)
        case(4)
        x(:)=(/3.22547689619392312d-1,1.74576110115834658d0,4.53602029692112798d0,9.39507091230113313d0/)
        w(:)=(/6.03154104341633602d-1,3.57418692437799687d-1,3.88879085150053843d-2,5.39294705561327450d-4/)
        case(8)
        x(:)=(/1.7027963230510100d-1,9.03701776799379912d-1,2.25108662986613069d0,4.26670017028765879d0,&
        &7.04590540239346570d0,1.07585160101809952d1,1.57406786412780046d1,2.28631317368892641d1/)
        w(:)=(/3.69188589341637530d-1,4.18780780814342956d-1,1.75794986637171806d-1,3.33434922612156515d-2,&
        &2.79453623522567252d-3,9.07650877335821310d-5,8.48574671627253154d-7,1.04800117487151038d-9/)
        case(12)
        x(:)=(/1.15722117358020675d-1,6.11757484515130065d-1,1.51261026977641879d0,2.83375133774350723d0,&
        &4.59922763941834848d0,6.84452545311517735d0,9.62131684245686704d0,1.30060549933063477d1,&
        &1.71168551874622557d1,2.21510903793970057d1,2.84879672509840003d1,3.70991210444069203d1/)
        w(:)=(/2.04731371055443190d-1,3.77759275873137982d-1,2.44082011319877564d-1,9.04492222116809307d-2,&
        &2.01023811540340965d-2,2.66397354186531588d-3,2.03231592662999392d-4,8.36505585681979875d-6,&
        &1.66849387654091026d-7,1.34239103051500415d-9,3.06160163503502078d-12,8.14807740742624168d-16/)
        case(16)
        x(:)=(/8.76494104789278403d-2,4.62696328915080832d-1,1.14105777483122686d0,2.12928364509838062d0,&
        &3.43708663389320665d0,5.07801861454976791d0,7.07033853504823413d0,9.43831433639193878d0,&
        &1.22142233688661587d1,1.54415273687816171d1,1.91801568567531349d1,2.35159056939919085d1,&
        &2.85787297428821404d1,3.45833987022866258d1,4.19404526476883326d1,5.17011603395433184d1/)
        w(:)=(/2.00151714957800994d-1,3.31057854950884166d-1,2.65795777644214153d-1,1.36296934296377540d-1,&
        &4.73289286941252190d-2,1.12999000803394532d-2,1.84907094352631086d-3,2.04271915308278460d-4,&
        &1.48445868739812988d-5,6.82831933087119956d-7,1.88102484107967321d-8,2.86235024297388162d-10,&
        &2.12707903322410297d-12,6.29796700251786779d-15,5.05047370003551282d-18,4.16146237037285519d-22/)
        case(20)
        x(:)=(/7.05398890919887534d-2,3.72126818001611444d-1,9.16582102483273565d-1,1.70730653102834388d0,&
        &2.74919925530943213d0,4.04892531385088692d0,5.61517497080161651d0,7.45901745367100331d0,&
        &9.59439286958109677d0,1.20388025409643163d1,1.48142934426307400d1,1.79488955205193700d1,&
        &2.14787882402850110d1,2.54517027931869055d1,2.99325546317006120d1,3.50134342404790000d1,&
        &4.08330570567285711d1,4.76199940473465021d1,5.58107957500638989d1,6.65244165256157538d1/)
        w(:)=(/1.68746801851113862d-1,2.91254362006008282d-1,2.66686102807001289d-1,1.66002453269506840d-1,&
        &7.48260646687923705d-2,2.49644173092832211d-2,6.20255084457223085d-3,1.14496238647690824d-3,&
        &1.55741773027811975d-4,1.54014408052249157d-5,1.08048036651798235d-6,5.33012090955671475d-8,&
        &1.75798117905058200d-9,3.72550240251232087d-11,4.76752925157819052d-13,3.37284424336243841d-15,&
        &1.15501433950039883d-17,1.53952214058234355d-20,5.28644272556915783d-24,1.65645661249902330d-28/)
        case(24)
        x(:)=(/5.90198521815079770d-2,3.11239146198483727d-1,7.66096905545936646d-1,1.42559759080361309d0,&
        2.29256205863219029d0,3.37077426420899772d0,4.66508370346717079d0,6.18153511873676541d0,&
        &7.92753924717215218d0,9.91209801507770602d0,1.21461027117297656d1,1.46427322895966743d1,&
        &1.74179926465089787d1,2.04914600826164247d1,2.38873298481697332d1,2.76359371743327174d1,&
        &3.17760413523747233d1,3.63581058016516217d1,4.14517204848707670d1,4.71531064451563230d1,&
        &5.36085745446950698d1,6.10585314472187616d1,6.99622400351050304d1,8.14982792339488854d1/)
        w(:)=(/1.42811973334781851d-1,2.58774107517423903d-1,2.58806707272869802d-1,1.83322688977778025d-1,&
        &9.81662726299188922d-2,4.07324781514086460d-2,1.32260194051201567d-2,3.36934905847830355d-3,&
        &6.72162564093547890d-4,1.04461214659275180d-4,1.25447219779933332d-5,1.15131581273727992d-6,&
        &7.96081295913363026d-8,4.07285898754999966d-9,1.50700822629258492d-10,3.91773651505845138d-12,&
        &6.89418105295808569d-14,7.81980038245944847d-16,5.35018881301003760d-18,2.01051746455550347d-20,&
        &3.60576586455295904d-23,2.45181884587840269d-26,4.08830159368065782d-30,5.57534578832835675d-34/)
        case(28)
        x(:)=(/5.07346248498738876d-2,2.67487268640741084d-1,6.58136628354791519d-1,1.22397180838490772d0,&
        &1.96676761247377770d0,2.88888332603032189d0,3.99331165925011414d0,5.28373606284344256d0,&
        &6.76460340424350515d0,8.44121632827132449d0,1.03198504629932601d1,1.24079034144606717d1,&
        &1.47140851641357488d1,1.72486634156080563d1,2.00237833299517127d1,2.30538901350302960d1,&
        &2.63562973744013176d1,2.99519668335961821d1,3.38666055165844592d1,3.81322544101946468d1,&
        &4.27896723707725763d1,4.78920716336227437d1,5.35112979596642942d1,5.97487960846412408d1,&
        &6.67569772839064696d1,7.47867781523391618d1,8.43178371072270431d1,9.65824206275273191d1/)
        w(:)=(/1.23778843954286428d-1,2.32279276900901161d-1,2.47511896036477212d-1,1.92307113132382827d-1,&
        &1.16405361721130006d-1,5.63459053644773065d-2,2.20663643262588079d-2,7.02588763558386773d-3,&
        &1.82060789269585487d-3,3.83344303857123177d-4,6.53508708069439831d-5,8.97136205341076834d-6,&
        &9.84701225624928887d-7,8.56407585267304245d-8,5.83683876313834429d-9,3.07563887784230228d-10,&
        &1.23259095272442282d-11,3.68217367410831200d-13,7.99879057596890965d-15,1.22492250032408341d-16,&
        &1.27112429503067374d-18,8.48859336768654320d-21,3.40245537942551185d-23,7.42015658886748513d-26,&
        &7.60041320580173769d-29,2.87391031794039581d-32,2.54182290388931800d-36,1.66137587802903396d-41/)
        case(32)
        x(:)=(/4.44893658332670184d-2,2.34526109519618537d-1,5.76884629301886426d-1,1.07244875381781763d0,&
        &1.72240877644464544d0,2.52833670642579488d0,3.49221327302199449d0,4.61645676974976739d0,&
        &5.90395850417424395d0,7.35812673318624111d0,8.98294092421259610d0,1.07830186325399721d1,&
        &1.27636979867427251d1,1.49311397555225573d1,1.72924543367153148d1,1.98558609403360547d1,&
        &2.26308890131967745d1, 2.56286360224592478d1,2.88621018163234747d1,3.23406291539647370d1,&
        &3.61004948057519738d1,4.01457197715394415d1,4.45092079957549380d1,4.92243949873086392d1,&
        &5.43337213333969073d1,5.98925091621340182d1,6.59753772879350528d1,7.26876280906627086d1,&
        &8.01874469779135231d1,8.87353404178923987d1,9.88295428682839726d1, 1.11751398097937695d2/)
        w(:)=(/1.09218341952384971d-1,2.10443107938813234d-1,2.35213229669848005d-1,1.95903335972881043d-1,&
        &1.29983786286071761d-1,7.05786238657174415d-2,3.17609125091750703d-2,1.19182148348385571d-2,&
        &3.73881629461152479d-3,9.80803306614955132d-4,2.14864918801364188d-4,3.92034196798794720d-5,&
        &5.93454161286863288d-6,7.41640457866755222d-7,7.60450787912078148d-8,0.35060222062580674d-9,&
        &4.28138297104092888d-10,2.30589949189133608d-11,9.79937928872709406d-13,3.23780165772926646d-14,&
        &8.17182344342071943d-16,1.54213383339382337d-17,2.11979229016361861d-19,2.05442967378804543d-21,&
        &1.34698258663739516d-23,5.66129413039735937d-26,1.41856054546303691d-28,1.91337549445422431d-31,&
        &1.19224876009822236d-34,2.67151121924013699d-38,1.33861694210025628d-42,4.51053019389897424d-48/)
        case(40)
            x(:)=(/0.035700394308888385122d0,0.18816228315869851600d0,0.46269428131457645357d0,0.85977296397293492226d0,&
            & 1.3800108205273371865d0,2.0242091359228267334d0,2.7933693535068164577d0,3.6887026779082702096d0,&
            & 4.7116411465549726936d0,5.8638508783437181143d0,7.1472479081022882507d0,8.5640170175861637627d0,&
            & 10.116634048451939407d0,11.807892294004584843d0,13.640933712537087228d0,15.619285893339073837d0,&
            & 17.746905950095663043d0,20.028232834574890530d0,22.468249983498418351d0,25.072560772426203794d0,&
            & 27.847480009168862721d0,30.800145739445462701d0,33.938657084913719609d0,37.272245880476004328d0,&
            & 40.811492823886920466d0,44.568603175334462707d0,48.557763533059992281d0,52.795611187216932969d0,&
            & 57.301863323393627495d0,62.100179072775111612d0,67.219370927126998799d0,72.695158847612462118d0,&
            & 78.572802911571309281d0,84.911231135704984543d0,91.789874671236376992d0,99.320808717446808250d0,&
            & 107.67244063938827252d0,117.12230951269068881d0,128.20184198825565119d0,142.28004446915999789d0/)
            w(:)=(/0.088412106190342440940d0,0.17681473909572229560d0,0.21136311701596243103d0,0.19408119531860179966d0,&
            & 0.14643428242412511441d0,0.093326798435770880507d0,0.050932204361044237026d0,0.023976193015684841840d0,&
            & 0.0097746252467144596189d0,0.0034579399930184868613d0,0.0010622468938968719350d0,0.00028327168532432471583d0,&
            & 0.000065509405003246292798d0,0.000013116069073267784125d0,2.2684528787793650545d-6, 3.3796264822006792108d-7,&
            & 4.3228213222820885689d-8, 4.7284937709907793279d-9, 4.4031741042328488129d-10, 3.4724414848038224856d-11,&
            & 2.3053815449168221616d-12, 1.2797725976766356072d-13, 5.8941771723511529447d-15, 2.2322175799045774184d-16,&
            & 6.8803364842843023409d-18, 1.7056037368180867485d-19, 3.3537119406661829355d-21, 5.1461995601366791408d-23,&
            & 6.0447625115876632890d-25, 5.3105847773213307528d-27, 3.3925280532805218961d-29, 1.5217354931814569975d-31,&
            & 4.5852916145026869176d-34, 8.7621586574862485610d-37, 9.8274157251479333061d-40, 5.8011520191697791085d-43,&
            & 1.5309086846066868536d-46, 1.3819863056493280997d-50, 2.5666336050123721838d-55, 2.7003609402170336406d-61/)
        case(60)
            x(:)=(/0.023897977262724994782d0,0.12593471888169076076d0,0.30957893432678988075d0,0.57499554209280526609d0,&
            & 0.92236948211666379157d0,1.3519383600081679397d0,1.8639963442992054802d0,2.4588958438224285850d0,&
            & 3.1370490097858959310d0,3.8989293872049917811d0,4.7450738001258887621d0,5.6760845082469168260d0,&
            & 6.6926316627865745904d0,7.7954560890310120110d0,8.9853724256576560911d0,10.263272655037909547d0,&
            & 11.630130063841871782d0,13.087003679350245100d0,14.635043234018346591d0,16.275494719209406879d0,&
            & 18.009706598857114263d0,19.839136765434035716d0,21.765360334373534885d0,23.790078389494180642d0,&
            & 25.915127811604900409d0,28.142492346079812161d0,30.474315093739509195d0,32.912912644080369003d0,&
            & 35.460791112322411216d0,38.120664393927127573d0,40.895475014812933693d0,43.788418035940638635d0,&
            & 46.802968571856479942d0,49.942913610317752243d0,53.212388982588307206d0,56.615922542696984717d0,&
            & 60.158484884500430204d0,63.845549279532241991d0,67.683162987059545208d0,71.678032714447406955d0,&
            & 75.837627854657063291d0,80.170306292607892439d0,84.685469194509277238d0,89.393753490252791880d0,&
            & 94.307274066118867443d0,99.439932542889869648d0,104.80781680747746147d0,110.42972668651628922d0,&
            & 116.32787889753132982d0,122.52887338413981661d0,129.06505218529826502d0,135.97646860411319811d0,&
            & 143.31384526024606358d0,151.14321669561511669d0,159.55362523885103398d0,168.67080654892222049d0,&
            & 178.68392501314637910d0,189.90524696213377268d0,202.93398795040067732d0,219.31811577379970553d0/)
            w(:)=(/0.059883611523733380661d0,0.12591096707540107193d0,0.16473078908210712807d0,0.17239118732674750403d0,&
            & 0.15442926800152203448d0,0.12180351302605122744d0,0.085807976879846703278d0,0.054435367226483746819d0,&
            & 0.031252789752223371974d0,0.016290259047635155781d0,0.0077247456605088567732d0,0.0033366894943076820164d0,&
            & 0.0013138658634496515358d0,0.00047178872610582552603d0,0.00015450073634605151289d0,0.000046133630291513734031d0,&
            & 0.000012555541586433208077d0,3.1126536116671386255d-6, 7.0238923169229563530d-7, 1.4413942338722032050d-7,&
            & 2.6871155382478738257d-8, 4.5453240465711692205d-9, 6.9667460870539940245d-10, 9.6611190516491905786d-11,&
            & 1.2101372902012551474d-11, 1.3666508971830879161d-12, 1.3887583568740821736d-13, 1.2670440927349749065d-14,&
            & 1.0354183850146325428d-15, 7.5590720583370492079d-17, 4.9160556832367860941d-18, 2.8393315957761979668d-19,&
            & 1.4514379644950183886d-20, 6.5427350200929260331d-22, 2.5902519096833059549d-23, 8.9664278435414917268d-25,&
            & 2.7006780092750217553d-26, 7.0398894915621409495d-28, 1.5787547853764471914d-29, 3.0258776894584834948d-31,&
            & 4.9201673552566388837d-33, 6.7316641160505118187d-35, 7.6780965388276276462d-37, 7.2246694240103543426d-39,&
            & 5.5415003983611333256d-41, 3.4176601279081430909d-43, 1.6681495220378452232d-45, 6.3255732716007591197d-48,&
            & 1.8231396385814367186d-50, 3.8906596692228002314d-53, 5.9551615457676981661d-56, 6.2854492261473107477d-59,&
            & 4.3523952400430152023d-62, 1.8533564849869102490d-65, 4.4482734830374030633d-69, 5.3225663149557690379d-73,&
            & 2.6412067805224611180d-77, 4.0165058425505468822d-82, 1.0516941039201472127d-87, 1.0909419486248200726d-94/)
        case(80)
            x(:)=(/0.017960423300698365554d0,0.094639912994353988811d0,0.23262286812586756921d0,0.43199254780238748026d0,&
            & 0.69282886135202183991d0,1.0152325561894714374d0,1.3993276878428727741d0,1.8452623038358451381d0,&
            & 2.3532088716092615245d0,2.9233646865554263248d0,3.5559523140461340594d0,4.2512200823098780832d0,&
            & 5.0094426336201647724d0,5.8309215386087190198d0,6.7159859778513171116d0,7.6649934948917730607d0,&
            & 8.6783308251677010954d0,9.7564148057429307132d0,10.899693371287855377d0,12.108646642365699901d0,&
            & 13.383788112778647370d0,14.725665943508585539d0,16.134864371662466579d0,17.612005243814437860d0,&
            & 19.157749684241247922d0,20.772799909792096092d0,22.457901204540458311d0,24.213844068958647377d0,&
            & 26.041466560165586693d0,27.941656841859465556d0,29.915355964900985501d0,31.963560902208920711d0,&
            & 34.087327864726189875d0,36.287775928781454459d0,38.566091009292210458d0,40.923530218031267200d0,&
            & 43.361426651731230296d0,45.881194661278886346d0,48.484335660833189136d0,51.172444544607010596d0,&
            & 53.947216789554447121d0,56.810456334636223134d0,59.764084342109954943d0,62.810148963926477204d0,&
            & 65.950836257456057343d0,69.188482420236277374d0,72.525587544263345359d0,75.964831127864174827d0,&
            & 79.509089629088836962d0,83.161456401053689663d0,86.925264419615623448d0,90.804112300940755952d0,&
            & 94.801894215947433207d0,98.922834446940579165d0,103.17152750803913023d0,107.55298497753990633d0,&
            & 112.07269048412833362d0,116.73666467350366632d0,121.55154249095262557d0,126.52466579651554034d0,&
            & 131.66419525212031087d0,136.97924668693697395d0,142.48005891216160193d0,148.17820245500444182d0,&
            & 154.08684228179869786d0,160.22107287009571594d0,166.59835193405391874d0,173.23907133424950383d0,&
            & 180.16732304903231798d0,187.41194967696377239d0,195.00802244153299145d0,202.99898419507493782d0,&
            & 211.43987049483646669d0,220.40236815173573965d0,229.98320607568000435d0,240.31908705584154042d0,&
            & 251.61587933049961117d0,264.21382388319910210d0,278.76673304600456365d0,296.96651199565134576d0/)
            w(:)=(/0.045272641464027452755d0,0.097622691129370068013d0,0.13365852798017004242d0,0.14937645827101128701d0,&
            & 0.14584706978660375008d0,0.12798046766667039213d0,0.10240364873547614395d0,0.075344057336652620107d0,&
            & 0.051240835821065727379d0,0.032323138479013890286d0,0.018956660110987408170d0,0.010353138932155252397d0,&
            & 0.0052716047497087076864d0,0.0025045007398191357512d0,0.0011108140479285224387d0,0.00046009904063076981317d0,&
            & 0.00017800427319655150323d0,0.000064328305181779756454d0,0.000021714152386661386056d0,6.8452422068302032728d-6,&
            & 2.0148387768652888554d-6, 5.5356096108309897758d-7, 1.4190662350227729006d-7, 3.3928319886067029298d-8,&
            & 7.5618521942413492988d-9, 1.5702103797210433426d-9, 3.0358719151260974984d-10, 5.4614960265651831348d-11,&
            & 9.1353117403076832507d-12, 1.4196250412935834672d-12, 2.0478152050918212355d-13, 2.7395267234969947724d-14,&
            & 3.3954734164377129687d-15, 3.8950030226439436306d-16, 4.1305652335626718890d-17, 4.0446805024741155105d-18,&
            & 3.6523826707208946996d-19, 3.0373390214298404880d-20, 2.3227574806659417900d-21, 1.6309293123767971508d-22,&
            & 1.0497107364507377730d-23, 6.1821869666743219930d-25, 3.3253390304932169409d-26, 1.6303421327356487754d-27,&
            & 7.2700624371570635928d-29, 2.9418210595485496140d-30, 1.0775663295033320128d-31, 3.5634814454894239535d-33,&
            & 1.0609024804179843587d-34, 2.8347893331015768975d-36, 6.7761153137539824132d-38, 1.4438133638965327344d-39,&
            & 2.7317419619723957224d-41, 4.5703808506247258942d-43, 6.7309359032141144705d-45, 8.6827148453882270607d-47,&
            & 9.7574467664989971107d-49, 9.4957661365042295333d-51, 7.9503372413637032776d-53, 5.6852252212220386817d-55,&
            & 3.4443689191523066463d-57, 1.7520838300614889115d-59, 7.4077279619841976241d-62, 2.5735398256614593841d-64,&
            & 7.2517056812244499503d-67, 1.6328052016335641917d-69, 2.8875122894699307129d-72, 3.9306497186688057597d-75,&
            & 4.0218963695440024828d-78, 3.0065656681535619953d-81, 1.5862629302776148181d-84, 5.6593984420409922697d-88,&
            & 1.2934555058728238527d-91, 1.7649977712530878877d-95, 1.3078621805548315574d-99, 4.6038641952371485236d-104,&
            & 6.2980356469014355246d-109, 2.4056926530320779390d-114, 1.3648444753359407889d-120, 2.2905062537183813302d-128/)
        case(100)
            x(:)=(/0.014386146995419669464d0,0.075803612023357124643d0,0.18631410205718717371d0,0.34596918099142909081d0,&
            & 0.55481093758091550960d0,0.81289128411566884504d0,1.1202738350075401486d0,1.4770343299238270697d0,&
            & 1.8832608263423947058d0,2.3390538496460341718d0,2.8445265427553590665d0,3.3998048274457119443d0,&
            & 4.0050275817586520175d0,4.6603468355689084595d0,5.3659279855851170149d0,6.1219500308040197891d0,&
            & 6.9286058293761730550d0,7.7861023778625174343d0,8.6946611139221679892d0,9.6545182435550807278d0,&
            & 10.665925094121675551d0,11.729148494472225072d0,12.844471183641030572d0,14.012192249694274649d0,&
            & 15.232627600466697843d0,16.506110468081989594d0,17.832991949326387430d0,19.213641584136066685d0,&
            & 20.648447974668349520d0,22.137819447656704407d0,23.682184763002376100d0,25.281993871834041539d0,&
            & 26.937718727574264944d0,28.649854153891292171d0,30.418918773790942920d0,32.245456004520665842d0,&
            & 34.130035123421516471d0,36.073252410379973226d0,38.075732373107093191d0,40.138129062115545763d0,&
            & 42.261127482984794178d0,44.445445114311806077d0,46.691833540651539419d0,49.001080210772436956d0,&
            & 51.374010332703994521d0,53.811488918355660438d0,56.314422991961708238d0,58.883763978282090294d0,&
            & 61.520510288396142647d0,64.225710123101560167d0,67.000464516419311597d0,69.845930644558374286d0,&
            & 72.763325428974586256d0,75.753929465939936193d0,78.819091319411472200d0,81.960232219060124004d0,&
            & 85.178851211218900200d0,88.476530817394611222d0,91.854943263049308931d0,95.315857348831720607d0,&
            & 98.861146047613527717d0,102.49279492391656096d0,106.21291148804681621d0,110.02373561603091696d0,&
            & 113.92765118897162371d0,117.92719913257322684d0,122.02509207044162108d0,126.22423084475038757d0,&
            & 130.52772320679941327d0,134.93890504022740758d0,139.46136455424015676d0,144.09896997721272416d0,&
            & 148.85590139775824946d0,153.73668754797303061d0,158.74624851171310443d0,163.88994558258723283d0,&
            & 169.17363981000302458d0,174.60376118237662671d0,180.18739094024569619d0,185.93236023966697143d0,&
            & 191.84736937224832912d0,197.94213310214325741d0,204.22755956703050772d0,210.71597286157694337d0,&
            & 217.42139327200148110d0,224.35989478887460815d0,231.55006802517248422d0,239.01362975131492447d0,&
            & 246.77624096724849046d0,254.86862925704743028d0,263.32816846915789311d0,272.20117002409253683d0,&
            & 281.54632828389738880d0,291.44013361637710726d0,301.98585525163915367d0,313.32953400407552434d0,&
            & 325.69126343702652000d0,339.43510192344961654d0,355.26131188853413247d0,374.98411283434267870d0/)
            w(:)=(/0.036392605883401356537d0,0.079676746212951398550d0,0.11211510334248694468d0,0.13035661297514618374d0,&
            & 0.13404333972846238040d0,0.12540709078066374997d0,0.10831411209726027355d0,0.087096638469959342035d0,&
            & 0.065551009312310614326d0,0.046340133582644259874d0,0.030846308627681445601d0,0.019367828113978789110d0,&
            & 0.011485442360179691617d0,0.0064389510016104295905d0,0.0034149799896926625930d0,0.0017143197401822081619d0,&
            & 0.00081487159158783785329d0,0.00036685483659948815078d0,0.00015645207417810679472d0,0.000063210870528885807631d0,&
            & 0.000024195752265189292644d0,8.7743097637554872419d-6,3.0142674860009475187d-6,9.8083358993452597680d-7,&
            & 3.0226387435322543080d-7,8.8200583952959386115d-8,2.4364258562006733674d-8,6.3697113739017570248d-9,&
            & 1.5756003204596680080d-9,3.6863292013461326689d-10,8.1547989242246112119d-11,1.7050625568260654077d-11,&
            & 3.3682141708667259378d-12,6.2835249553666573660d-13,1.1064980159833061707d-13,1.8383501754549013832d-14,&
            & 2.8801150572305857339d-15,4.2526128973372414067d-16,5.9144324672522439838d-17,7.7430891025001089332d-18,&
            & 9.5362494490747949743d-19,1.1040945169869171454d-19,1.2008469704304970569d-20,1.2260070074904774033d-21,&
            & 1.1740217146467280969d-22,1.0535940646087311407d-23,8.8532317684742975097d-25,6.9591987754834461444d-26,&
            & 5.1123698154658325241d-27,3.5062721488171387421d-28,2.2426462726939223396d-29,1.3362094234457496265d-30,&
            & 7.4074289575757510893d-32,3.8158601371363446901d-33,1.8242001976828543039d-34,8.0816323360641805778d-36,&
            & 3.3130638018623846737d-37,1.2548329976876779940d-38,4.3837909801617165881d-40,1.4101392292173474850d-41,&
            & 4.1688639996685860322d-43,1.1304819044771411726d-44,2.8060374967871239466d-46,6.3612705709821942695d-48,&
            & 1.3139811959781289076d-49,2.4668106921703939992d-51,4.1977470228717237504d-53,6.4562691049120031533d-55,&
            & 8.9473372486217303203d-57,1.1135669746997689184d-58,1.2402350422552966639d-60,1.2313765726083698151d-62,&
            & 1.0853674489051309034d-64,8.4549564909552701838d-67,5.7926307566097949606d-69,3.4718547763506601612d-71,&
            & 1.8098595335622393581d-73,8.1537526084108297035d-76,3.1524725487932880430d-78,1.0379005899152634183d-80,&
            & 2.8848755233073495343d-83,6.7048012259685912147d-86,1.2889637464360916588d-88,2.0248483451529085824d-91,&
            & 2.5633922199660732721d-94,2.5739615075159617883d-97,2.0126547874800776647d-100,1.1994844193823276437d-103,&
            & 5.3121327315398034324d-107,1.6959692594467450795d-110,3.7620972986344511541d-114,5.5396417544496093738d-118,&
            & 5.1106404770917605505d-122,2.7399654694003410753d-126,7.7136114926382004229d-131,9.8824946009588260463d-136,&
            & 4.6468630072942033152d-141,5.6260372950198530067d-147,8.9050314058891380744d-154,3.2465651634358090752d-162/)
    end select

    select case(xw)
    case(1)
        Laguerre=x(j)
    case(2)
        Laguerre=w(j)
    end select

    deallocate(x,w)
    return
end function Laguerre

!***** Gauss-Legendre Integration ***** 20230919 added from ntest.f90, for integral w/ cut-off.
real*8 function GLpoint(Nnode,xw,j)
    implicit none
    integer :: Nnode,xw,j,i !! Nnode判断几阶多项式；xw判断零点/系数；j为第几个变量；i为函数内变量
    real*8,allocatable :: x(:),w(:)
    !Nnode=60; xw=1对应x & xw=2对应w；
    i=Nnode/2
    allocate(x(i)); allocate(w(i))

    !! 直接列出Legendre多项式零点；
    select case(Nnode)
        case(2)
        x(:)=(/0.5773502691896258d0/)
        w(:)=(/1.d0/)
        case(4)
        x(:)=(/0.3399810435848563d0,0.8611363115940526d0/)
        w(:)=(/0.6521451548625462,0.3478548451374539/)
        case(6)
        x(:)=(/0.2386191860831969,0.6612093864662645,0.9324695142031520/)
        w(:)=(/0.467913934572691,0.3607615730481386,0.1713244923791704/)
        case(8)
        x(:)=(/0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975362/)
        w(:)=(/0.362683783378362,0.3137066458778873,0.2223810344533744,0.1012285362903763/)
        case(10)
        x(:)=(/0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717/)
        w(:)=(/0.2955242247147529,0.2692667193099964,0.219086362515982,0.1494513491505806,0.06667134430868817/)
        case(12)
        x(:)=(/0.1252334085114689,0.3678314989981802,0.5873179542866174,0.7699026741943047,0.9041172563704749,0.9815606342467193/)
        w(:)=(/0.2491470458134028,0.2334925365383548,0.2031674267230659,0.1600783285433462,0.1069393259953184,0.04717533638651188/)
        case(14)
        x(:)=(/0.1080549487073437,0.3191123689278898,0.5152486363581541,0.6872929048116855,&
        &0.8272013150697650,0.9284348836635735,0.9862838086968123/)
        w(:)=(/0.2152638534631578,0.2051984637212956,0.1855383974779378,0.1572031671581935,&
        &0.1215185706879032,0.0801580871597602,0.0351194603317518/)
        case(16)
        x(:)=(/0.09501250983763744,0.2816035507792589,0.4580167776572274,0.6178762444026437,0.7554044083550030,&
        &0.8656312023878317,0.9445750230732326,0.9894009349916499/)
        w(:)=(/0.1894506104550685,0.1826034150449236,0.1691565193950025,0.1495959888165768,0.1246289712555339,&
        &0.0951585116824928,0.06225352393864787,0.02715245941175404/)
        case(18)
        x(:)=(/0.08477501304173530,0.2518862256915055,0.4117511614628426,0.5597708310739475,0.6916870430603532,0.8037049589725231,&
        &0.8926024664975557,0.9558239495713978,0.9915651684209309/)
        w(:)=(/0.1691423829631436,0.1642764837458327,0.1546846751262652,0.1406429146706507,0.1225552067114785,0.1009420441062872,&
        &0.07642573025488906,0.04971454889496976,0.02161601352648329/)
        case(20)
        x(:)=(/0.07652652113349733,0.2277858511416451,0.3737060887154196,0.5108670019508271,0.6360536807265150,0.7463319064601508,&
        &0.8391169718222188,0.9122344282513259,0.9639719272779138,0.9931285991850949/)
        w(:)=(/0.1527533871307258,0.1491729864726037,0.1420961093183821,0.1316886384491766,0.1181945319615184,0.1019301198172404,&
        &0.0832767415767047,0.06267204833410905,0.04060142980038694,0.01761400713915213/)
        case(24)
        x(:)=(/0.06405689286260563,0.1911188674736163,0.3150426796961634,0.4337935076260451,0.5454214713888395,0.6480936519369756,&
        &0.7401241915785544,0.8200019859739029,0.8864155270044010,0.9382745520027328,0.9747285559713095,0.9951872199970214/)
        w(:)=(/0.1279381953467521,0.1258374563468283,0.1216704729278034,0.1155056680537256,0.1074442701159656,0.0976186521041139,&
        &0.0861901615319533,0.07334648141108032,0.05929858491543678,0.04427743881741979,0.02853138862893364,0.0123412297999872/)
        case(28)
        x(:)=(/0.05507928988403427,0.1645692821333808,0.2720616276351781,0.3762515160890787,0.4758742249551183,0.5697204718114017,&
        &0.6566510940388650,0.7356108780136318,0.8056413709171792,0.8658925225743950,0.9156330263921321,0.9542592806289382,&
        &0.9813031653708728,0.9964424975739544/)
        w(:)=(/0.1100470130164752,0.1087111922582941,0.1060557659228464,0.1021129675780608,0.0969306579979299,0.0905717443930328,&
        &0.0831134172289012,0.07464621423456879,0.06527292396699958,0.05510734567571676,0.04427293475900424,0.03290142778230436,&
        &0.02113211259277126,0.00912428259309449/)
        case(32)
        x(:)=(/0.04830766568773832,0.1444719615827965,0.2392873622521371,0.3318686022821276,0.4213512761306353,0.5068999089322294,&
        &0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,&
        &0.9349060759377397,0.9647622555875064,0.9856115115452683,0.9972638618494816/)
        w(:)=(/0.0965400885147278,0.0956387200792749,0.0938443990808046,0.0911738786957639,0.0876520930044038,0.0833119242269468,&
        &0.07819389578707031,0.07234579410884851,0.06582222277636184,0.05868409347853555,0.05099805926237617,0.04283589802222668,&
        &0.03427386291302144,0.02539206530926206,0.01627439473090566,0.007018610009470137/)
        case(36)
        x(:)=(/0.04301819847370861,0.1287361038093848,0.2135008923168656,0.2966849953440283,0.3776725471196892,0.4558639444334203,&
        &0.5306802859262452,0.6015676581359805,0.6680012365855211,0.7294891715935566,0.7855762301322065,0.8358471669924753,&
        &0.8799298008903971,0.9174977745156591,0.9482729843995075,0.9720276910496979,0.9885864789022122,0.9978304624840858/)
        w(:)=(/0.0859832756703948,0.0853466857393386,0.0840782189796619,0.0821872667043397,0.07968782891207159,0.07659841064587068,&
        &0.07294188500565308,0.06874532383573645,0.06403979735501551,0.05886014424532481,0.05324471397775992,0.04723508349026599,&
        &0.0408757509236449,0.03421381077030722,0.02729862149856879,0.02018151529773548,0.01291594728406556,0.005565719664245083/)
        case(40)
        x(:)=(/0.03877241750605082,0.1160840706752552,0.1926975807013711,0.2681521850072537,0.3419940908257585,&
        &0.4137792043716050,0.4830758016861787,0.5494671250951282,0.6125538896679802,0.6719566846141795,0.7273182551899271,&
        &0.7783056514265194,0.8246122308333117,0.8659595032122595,0.9020988069688743,0.9328128082786765,0.9579168192137917,&
        &0.9772599499837743,0.9907262386994570,0.9982377097105592/)
        w(:)=(/0.07750594797842481,0.07703981816424798,0.07611036190062625,0.07472316905796826,0.07288658239580405,&
        &0.07061164739128678,0.0679120458152339,0.06480401345660103,0.06130624249292894,0.05743976909939155,0.05322784698393682,&
        &0.04869580763507223,0.04387090818567327,0.038782167974472,0.03346019528254784,0.02793700698002339,0.02224584919416696,&
        &0.01642105838190788,0.0104982845311528,0.004521277098533155/)
        case(44)
        x(:)=(/0.03528923696413536,0.1056919017086532,0.1755680147755168,0.2445694569282013,0.3123524665027858,&
        &0.3785793520147071,0.4429201745254115,0.5050543913882023,0.5646724531854708,0.6214773459035758,0.6751860706661224,&
        &0.7255310536607170,0.7722614792487559,0.8151445396451350,0.8539665950047104,0.8885342382860432,0.9186752599841758,&
        &0.9442395091181941,0.9650996504224931,0.9811518330779140,0.9923163921385158,0.9985402006367742/)
        w(:)=(/0.07054915778935406,0.07019768547355821,0.06949649186157258,0.06844907026936666,0.06706063890629366,&
        &0.06533811487918141,0.06329007973320387,0.06092673670156198,0.0582598598775955,0.05530273556372805,0.05207009609170445,&
        &0.04857804644835204,0.04484398408197004,0.04088651231034623,0.03672534781380888,0.03238122281206983,0.027875782821281,&
        &0.0232314819020192,0.01847148173681476,0.01361958675557998,0.00870048136752482,0.003745404803112729/)
        case(48)
        x(:)=(/0.03238017096286936,0.09700469920946270,0.1612223560688917,0.2247637903946891,0.2873624873554556,0.3487558862921607,&
        &0.4086864819907167,0.4669029047509584,0.5231609747222330,0.5772247260839727,0.6288673967765136,0.6778723796326639,&
        &0.7240341309238147,0.7671590325157403,0.8070662040294426,0.8435882616243935,0.8765720202742479,0.9058791367155697,&
        &0.9313866907065543,0.9529877031604309,0.9705915925462473,0.9841245837228269,0.9935301722663508,0.9987710072524261/)
        w(:)=(/0.06473769681268394,0.06446616443595009,0.06392423858464817,0.06311419228625402,0.06203942315989267,&
        &0.06070443916589388,0.05911483969839564,0.05727729210040321,0.05519950369998416,0.05289018948519367,&
        &0.05035903555385447,0.04761665849249047,0.04467456085669429,0.04154508294346476,0.0382413510658307,&
        &0.03477722256477044,0.03116722783279808,0.02742650970835695,0.02357076083932438,0.01961616045735552,&
        &0.01557931572294385,0.01147723457923455,0.007327553901276245,0.003153346052305896/)
        case(52)
        x(:)=(/0.02991410979733877,0.08963524464890057,0.1490355086069492,0.2079022641563661,0.2660247836050018,0.3231950034348078,&
        &0.3792082691160937,0.4338640677187617,0.4869667456980961,0.5383262092858274,0.5877586049795791,0.6350869776952459,&
        &0.6801419042271677,0.7227620997499832,0.7627949951937450,0.8000972834304683,0.8345354323267345,0.8659861628460676,&
        &0.8943368905344953,0.9194861289164245,0.9413438536413591,0.9598318269330866,0.9748838842217445,0.9864461956515498,&
        &0.9944775909292160,0.9989511111039503/)
        w(:)=(/0.05981036574529185,0.05959626017124816,0.05916881546604297,0.05852956177181388,0.05768078745252684,&
        &0.05662553090236859,0.05536756966930265,0.05391140693275727,0.05226225538390699,0.05042601856634238,0.04840926974407489,&
        &0.04621922837278481,0.0438637342590004,0.04135121950056027,0.03869067831042397,0.03589163483509724,0.0329641090897188,&
        &0.02991858114714395,0.02676595374650401,0.02351751355398446,0.02018489150798079,0.01678002339630073,0.01331511498234096,&
        &0.00980263457946275,0.006255523962973276,0.002691316950047123/)
        case(56)
        x(:)=(/0.02779703528727544,0.08330518682243537,0.1385558468103762,0.1933782386352753,0.2476029094343372,0.3010622538672207,&
        &0.3535910321749545,0.4050268809270913,0.4552108148784596,0.5039877183843817,0.5512068248555346,0.5967221827706633,&
        &0.6403931068070069,0.6820846126944705,0.7216678344501881,0.7590204227051289,0.7940269228938665,0.8265791321428817,&
        &0.8565764337627486,0.8839261083278275,0.9085436204206555,0.9303528802474963,0.9492864795619626,0.9652859019054902,&
        &0.9783017091402564,0.9882937155401615,0.9952312260810697,0.9990943438014656/)
        w(:)=(/0.0555797463065144,0.05540795250324512,0.05506489590176242,0.05455163687088943,0.05386976186571449,&
        &0.05302137852401078,0.0520091091517414,0.05083608261779847,0.04950592468304757,0.04802274679360025,0.0463911333730019,&
        &0.04461612765269228,0.04270321608466707,0.04065831138474452,0.03848773425924767,0.03619819387231519,0.03379676711561176,&
        &0.03129087674731045,0.02868826847382273,0.02599698705839195,0.02322535156256532,0.02038192988240258,0.01747551291140094,&
        &0.01451508927802147,0.0115098243403834,0.0084690631633079,0.005402522246015325,0.002323855375773273/)
        case(60)
        x(:)=(/0.02595977230124780,0.07780933394953657,0.1294491353969450,0.1807399648734254,0.2315435513760293,0.2817229374232617,&
        &0.3311428482684482,0.3796700565767980,0.4271737415830784,0.4735258417617071,0.5186014000585697,0.5622789007539445,&
        &0.6044405970485104,0.6449728284894771,0.6837663273813554,0.7207165133557304,0.7557237753065857,0.7886937399322641,&
        &0.8195375261621458,0.8481719847859296,0.8745199226468983,0.8985103108100459,0.9200784761776276,0.9391662761164232,&
        &0.9557222558399961,0.9697017887650527,0.9810672017525982,0.9897878952222217,0.9958405251188382,0.9992101232274360/)
        w(:)=(/0.05190787763122064,0.05176794317491018,0.05148845150098093,0.05107015606985562,0.05051418453250937,&
        &0.04982203569055018,0.04899557545575683,0.04803703181997118,0.04694898884891221,0.04573437971611448,0.04439647879578711,&
        &0.04293889283593564,0.04136555123558475,0.03968069545238081,0.03788886756924344,0.03599489805108449,0.03400389272494643,&
        &0.03192121901929633,0.02975249150078894,0.02750355674992479,0.02518047762152125,0.02278951694399782,0.02033712072945729,&
        &0.01782990101420771,0.01527461859678479,0.01267816647681596,0.01004755718228797,0.007389931163345448,&
        &0.004712729926953579,0.00202681196887374/)
        case(66)
        x(:)=(/0.02361813338592,0.07080169886814,0.11782727807868,0.16458993856472,0.21098533454808,0.2569099397615,&
        &0.30226127845644,0.34693815406634,0.3908408750162,0.4338714771734,0.47593394244424,0.51693441302823,0.55678140085172,&
        &0.59538599171414,0.63266204369096,0.66852637935097,0.70289897135887,0.73570312104951,0.76686562957532,0.79631696124583,&
        &0.82399139869559,0.84982718953602,0.87376668416695,0.89575646444656,0.91574746294521,0.93369507254479,0.94955924620624,&
        &0.96330458685926,0.97490042774629,0.98432090493514,0.99154503038338,0.99655682108955,0.99934620987218/)
        w(:)=(/0.047227481263,0.0471220982876,0.0469115674876,0.0465963586396,0.04617717509792,0.04565495222527,0.04503085530544,&
        &0.044306276943153,0.0434828339567,0.04256236377006,0.0415469203132,0.040438769439,0.03924038386683,0.03795443766594,&
        &0.03658380028814,0.0351315301655,0.0336008678861,0.031995228964047,0.03031819621887,0.028573511782932,0.02676506875425,&
        &0.0248969025148,0.02297318173533,0.02099819909186,0.0189763617228,0.01691218147225,0.01481026500273,0.01267530398126,&
        &0.01051206598771,0.008325388765991,0.0061201920184,0.00390162564174,0.00167765374401/)
        case(72)
        x(:)=(/0.02166394603542,0.06495116631186,0.1081164475621,0.15107875148221,0.19375742083703,0.23607233088576,&
        &0.27794403980785,0.31929393784671,0.36004439489142,0.4001189062192,0.439442236125,0.4779405591689,0.515541598776,&
        &0.55217476292771,0.58777127669166,0.62226431133947,0.6555891098112,0.68768310829047,0.71848605366224,0.74794011663283,&
        &0.77599000029998,0.80258304396929,0.82766932202275,0.85120173765444,0.87313611129878,0.8934312635881,0.91204909268867,&
        &0.92895464588092,0.94411618527254,0.9575052475777,0.96909669799878,0.97886877855723,0.98680315237583,0.9928849510168,&
        &0.99710287164273,0.99944993445296/)
        w(:)=(/0.0433211121655,0.04323978130522,0.043077272274914,0.0428338901683,0.04251009191006,0.0421064853976,&
        &0.04162382836014,0.04106302693608,0.04042513397173,0.03971134704483,0.03892300621617,0.0380615915138,0.0371287201545,&
        &0.03612614350764,0.0350557438072,0.0339195306183,0.0327196370643,0.031458315822562,0.0301379348954,0.0287609731647,&
        &0.02733001573895,0.02584774910066,0.02431695606442,0.022740510555036,0.02112137221644,0.0194625808633,0.0177672507892,&
        &0.01603856495029,0.01427976905455,0.01249416561987,0.01068510816535,0.00885599607371,0.00701027232186,0.00515143601879,&
        &0.0032831697747,0.00141151639397/)
        case(80)
        x(:)=(/0.01951138325679,0.05850443715242,0.0974083984416,0.13616402280914,0.17471229183265,0.21299450285767,&
        &0.25095235839227,0.28852805488451,0.3256643707477,0.36230475349949,0.39839340588197,0.43387537083176,0.46869661517054,&
        &0.50280411188879,0.53614592089713,0.5686712681227,0.60033062282975,0.6310757730469,0.66085989898612,0.68963764434203,&
        &0.7173651853621,0.7440002975836,0.76950242013504,0.79383271750461,0.81695413868146,0.83883147358026,0.85943140666311,&
        &0.87872256767821,0.89667557943877,0.91326310257176,0.92845987717245,0.94224276130987,0.95459076634364,0.9654850890438,&
        &0.97490914058573,0.98284857273863,0.98929130249976,0.99422754096569,0.99764986439824,0.99955382265163/)
        w(:)=(/0.03901781365631,0.03895839596277,0.03883965105905,0.03866175977408,0.03842499300696,0.03812971131448,&
        &0.037776364362,0.0373654902387,0.03689771463828,0.036373749905836,0.0357943939534,0.03516052904475,0.0344731204518,&
        &0.0337332149846,0.0329419393976,0.0321004986735,0.0312101741881,0.0302723217596,0.02928836958327,0.02825981605728,&
        &0.02718822750049,0.02607523576757,0.0249225357641,0.02373188286593,0.02250509024633,0.0212440261158,0.01995061087814,&
        &0.0186268142083,0.0172746520563,0.0158961835837,0.01449350804051,0.0130687615924,0.0116241141208,0.010161766041103,&
        &0.00868394526926,0.00719290476812,0.005690922451403,0.004180313124695,0.00266353358951,0.0011449500032/)
        case(90)
        x(:)=(/0.0173557291463,0.052046275137207,0.08667410942073,0.12119750815392,0.15557487333053,0.18976478290338,&
        &0.2237260406947,0.2574177260344,0.29079924306617,0.3238303696623,0.3564713058886,0.3886827219595,0.4204258056282,&
        &0.45166230895187,0.4823545943777,0.51246568009303,0.54195928458591,0.5707998703612,0.59895268676074,0.62638381183505,&
        &0.65306019321684,0.6789496879466,0.7040211012024,0.72824422388739,0.75158986902964,0.77402990695033,0.79553729915825,&
        &0.81608613092948,0.83565164253338,0.85421025906707,0.8717396188629,0.88821860043475,0.9036273479313,0.91794729506659,&
        &0.93116118750043,0.94325310364536,0.9542084738815,0.96401409817151,0.97265816209019,0.98013025134515,0.98642136505783,&
        &0.99152392881106,0.99543181205835,0.99814037993857,0.99964697128664/)
        w(:)=(/0.034707972489,0.0346661520857,0.0345825616695,0.03445730196032,0.03429052388638,0.0340824284023,0.03383326624683,&
        &0.0335433376411,0.0332129919266,0.032842627144,0.0324326895543,0.03198367310022,0.031496118811819,0.03097061415408,&
        &0.0304077923193,0.029808331464,0.0291729538921,0.02850242518416,0.0277975532753,0.027059187481548,0.0262882174765,&
        &0.0254855722194,0.0246522188359,0.0237891614525,0.02289743998716,0.02197812889593,0.0210323358787,0.0200612005446,&
        &0.0190658930391,0.01804761263446,0.0170075862852,0.01594706715101,0.01486733308804,0.013769685112337,0.01265544583717,&
        &0.0115259578891,0.01038258230989,0.00922669695774,0.00805969494462,0.00688298320846,0.00569798156075,0.00450612361367,&
        &0.00330886724334,0.00210777877453,9.059323712E-04/)
        case(100)
        x(:)=(/0.015628984421543,0.04687168242159,0.07806858281344,0.10918920358006,0.1402031372361,0.1710800805386,&
        &0.20178986409574,0.23230248184497,0.2625881203715,0.29261718803847,0.32236034390053,0.35178852637242,0.38087298162463,&
        &0.4095852916783,0.43789740217203,0.46578164977336,0.4932107892082,0.52015801988176,0.54659701206509,0.5725019326214,&
        &0.59784747024718,0.62260886020371,0.64676190851413,0.6702830156031,0.6931491993558,0.71533811757306,0.736828089802,&
        &0.75759811851971,0.7776279096495,0.79689789239031,0.81538923833918,0.8330838798884,0.84996452787959,0.86601468849717,&
        &0.88121867938502,0.89556164497073,0.90902957098253,0.92160929814533,0.9332885350431,0.94405587013626,0.95390078292549,&
        &0.96281365425582,0.97078577576371,0.97780935848692,0.9838775407061,0.98898439524299,0.99312493703744,0.99629513473313,&
        &0.9984919506396,0.9997137267734/)
        w(:)=(/0.0312554234539,0.0312248842548,0.0311638356962,0.03107233742757,0.03095047885049,0.0307983790312,0.03061618658398,&
        &0.03040407952645,0.0301622651052,0.02989097959333,0.02959048805991,0.02926108411064,0.02890308960113,0.0285168543224,&
        &0.0281027556591,0.02766119822079,0.0271926134466,0.02669745918357,0.0261762192395,0.02562940291021,0.02505754448158,&
        &0.024461202708,0.02384096026597,0.0231974231853,0.02253122025634,0.02184300241625,0.02113344211253,&
        &0.0204032326462,0.01965308749444,0.0188837396134,0.0180959407221,0.01729046056832,0.0164680861761,&
        &0.0156296210775,0.01477588452744,0.0139077107037,0.013025947893,0.01213145766298,0.0112251140232,&
        &0.0103078025749,0.0093804196537,0.00844387146967,0.00749907325546,0.00654694845085,0.00558842800387,&
        &0.0046244500634,0.00365596120133,0.00268392537155,0.001709392653518,7.3463449051E-04/)
    end select

    select case(xw)
    case(1)
        GLpoint=x(j)
    case(2)
        GLpoint=w(j)
    end select

    deallocate(x,w)
    return
end function GLpoint





