!!!!! We define parameters used for equations
real*8 rayon,rhosty,rhosol,mass,mz,mx,mt,grav,kappa,eps,xi,kxi,kxe,temp,k_B,beta,clight,dt,l_D,csquare,Lchrt
common/para/rayon,rhosty,rhosol,mass,mz,mx,mt,grav,kappa,eps,xi,kxi,kxe,temp,k_B,beta,clight,dt,l_D,csquare,Lchrt

!! rayon: radius of particle
!! rhosty: density of styrene
!! rhosol: density of solvent (water)
!! mass: mass of particle
!! mz,mx,mt: mass along different direction
!! grav: gravity acceleration
!! kappa: compliance
!! eps: dimensionless parameter
!! xi: dimensionless viscosity
!! kxi = kappa*xi
!! kxe = kappa*xi*eps
!! temp: temperature
!! k_B: Boltzmann constant
!! beta = 1 / (temp*k_B)
!! clight: max falling speed of particle
!! dt: time gap for numerical simulations
!! l_D: Debye length for electrostatic repulsion
!! csquare = clight**2
!! Lchrt: Characteristic / Typical Length ~ \mathcal{L} ~ 10^{-6} m ?






!!!!! Then there leaves some ancient parameters useless.
!real*8 Minv,gamma,intma,intmb
!common/calcul/Minv,gamma,intma,intmb
real*8 noise(3,2),Minv(3,3),gmaeff(3,3,3)
common/calcul/noise,Minv,gmaeff
!! "noise" noise = random force; 3: z/x/Θ; 2: Box-Muller method
!! "Minv" effective mass matrix inverse; 3,3: matrix index
!! "gmaeff" effective friction matrix; 3: gamma0/gamma1/gamma1v; 3,3: matrix index
