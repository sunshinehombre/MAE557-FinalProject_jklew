program GEQN_SOLVER
  use CONSTS
  use AUX_SUBROUTINES
  implicit none

  ! gEqnSolver.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! GEQN_SOLVER solves the G-equation for the corrugated flamelets regime in a 2D
  ! incompressible channel flow (working fluid = air at STP).

  ! Spatial step size h=dx=dy
  real(WP), parameter :: h1 = 0.0625_WP ! [m], Grid is 16x16
  real(WP), parameter :: h2 = 0.03125_WP ! [m], Grid is 32x32
  real(WP), parameter :: h3 = 0.015625_WP ! [m], Grid is 64x64

  ! Length of domain (square) side
  real(WP), parameter :: lside1 = 1.0_WP ! [m]

  ! Time step size & total simulation time
  real(WP), parameter :: dt1 = 1E-6_WP ! [s]
  real(WP), parameter :: tfin1 = 2E-2_WP ! [s]

  ! Mach number based on bulk velocity
  real(WP), parameter :: ma1 = 0.05_WP

  ! Will need to change velocity ranges in plotting script if select Mach #s below
  ! real(WP), parameter :: ma2 = 0.1_WP
  ! real(WP), parameter :: ma3 = 0.3_WP

  ! Initial profile of flame
  integer, parameter :: gprofile1 = 1 ! Planar
  integer, parameter :: gprofile2 = 2 ! Sinusoidal

  ! For timing purposes
  integer :: t1,t2,clock_rate,clock_max

  ! ! Start timer
  ! call system_clock (t1,clock_rate,clock_max)

  ! ! Set number of OpenMP threads
  ! call OMP_SET_NUM_THREADS(NTHREADS)

  ! ! Solve for incompressible flow field with iso-surface from G-Equation
  ! call solveAll(h3,lside1,dt1,tfin1,ma1,gprofile2)

  ! ! Stop timer
  ! call system_clock (t2,clock_rate,clock_max)
  ! print *, 'Elapsed real time =', real(t2-t1)/real(clock_rate), 's'

  ! *** Below for temporal and spatial convergence ***
  call convergence_Error()

  ! ! *** Below for primary and secondary conservation ***
  ! call conservation_Error()
  
contains
  
  ! ********************************************************************************
  ! *** Incompressible Navier-Stokes solver with incorporated G-Equation ***********
  ! ********************************************************************************
  subroutine solveAll(h,lside,dt,tfin,ma,gprofile)
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), intent(in) :: h ! Spatial step size
    real(WP), intent(in) :: lside ! Length of square side
    real(WP), intent(in) :: dt ! Time step size
    real(WP), intent(in) :: tfin ! Total time of simulation
    real(WP), intent(in) :: ma ! Flow mach number
    integer, intent(in) :: gprofile ! Specify starting profile of flame

    real(WP) :: U ! Bulk velocity [m/s]
    real(WP) :: nu ! Kinematic viscosity [m^2/s]

    integer :: nxy ! Total spatial steps nx=ny
    integer :: nt ! Total time steps
    integer :: t_iter ! Loop iterator

    ! Solution matrix with indices (variable,timestep,xloc,yloc)
    ! where the first index can be the following values
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, parameter :: var = 3
    real(WP), dimension(:,:,:,:), allocatable :: soln
    
    integer :: AllocateStatus
    

    U = SOUND*ma
    nu = U*lside/RE
    nxy = int(lside/h)
    nt = int(tfin/dt)

    allocate(soln(var,2,nxy+1,nxy+1), STAT = AllocateStatus)
    if (AllocateStatus.ne.0) stop "Error: Not enough memory."  

    ! Initialize solution matrix at t=0 (initial zeta, psi, and G values)
    call soln_init(h,nxy,U,var,soln,gprofile)

    ! Write initial u, v, and G profiles to file
    soln(:,2,:,:) = soln(:,1,:,:)
    call write_solution(ma,h,nxy,dt,0,var,soln)
    soln(:,2,:,:) = 0.0_WP

    ! Iterate until final time
    do t_iter = 1,nt
       call vs_step(dt,t_iter,h,nxy,nu,var,soln)
       call gEqn_step(dt,t_iter,h,nxy,var,soln)

       ! Write u, v, and G every nt/10 time steps
       if (MOD(t_iter,nt/10).eq.0) then
          call write_solution(ma,h,nxy,dt,t_iter,var,soln)
       end if
    end do

    ! For testing purposes
    print *, 'Bulk velocity: ', U, ' m/s'

  end subroutine solveAll

  ! ********************************************************************************
  ! *** Convergence Analysis *******************************************************
  ! ********************************************************************************
  subroutine convergence_Error()
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), parameter :: h3 = 0.015625_WP ! Spatial step size
    real(WP) :: h ! Variable spatial step size
    real(WP), parameter :: lside = 1.0_WP ! Length of square side
    real(WP) :: dt ! Time step size
    real(WP), parameter :: ma = 0.05_WP! Flow mach number
    integer, parameter :: gprofile = 2 ! Specify starting profile of flame

    real(WP) :: U ! Bulk velocity [m/s]
    real(WP) :: nu ! Kinematic viscosity [m^2/s]

    integer :: nxy, nxy_whole, nxy_half ! Total spatial steps nx=ny
    integer :: t_iter, t_iter2, x_iter, x_iter2 ! Loop iterators
    integer, parameter :: numdt1 = 6 ! Number of time steps dt to use
    integer, parameter :: numdt2 = 4 ! Number of time steps dt to use
    integer, parameter :: numh1 = 4 ! Number of spatial steps h to use
    integer, parameter :: numh2 = 4 ! Number of spatial steps h to use

    ! Solution matrix with indices (variable,timestep,xloc,yloc)
    ! where the first index can be the following values
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, parameter :: var = 3
    real(WP), dimension(:,:,:,:), allocatable :: solnt1
    real(WP), dimension(:,:,:,:), allocatable :: solnt2
    real(WP), dimension(:,:,:,:), allocatable :: solnt3
    real(WP), dimension(:,:,:,:), allocatable :: solnt4
    real(WP), dimension(:,:,:,:), allocatable :: solns1
    real(WP), dimension(:,:,:,:), allocatable :: solns2
    real(WP), dimension(:,:,:,:), allocatable :: solns3
    real(WP), dimension(:,:,:,:), allocatable :: solns4

    ! Cumulative spatial error over all j,k
    real(WP) :: errsum1
    
    ! Store timestep and temporal error as terr(dt,error)
    real(WP), dimension(numdt1,2) :: terr1
    real(WP), dimension(numdt2,2) :: terr2

    ! Store spatial resolution and spatial error as serr(h,error)
    real(WP), dimension(numh1,2) :: serr1
    real(WP), dimension(numh2,2) :: serr2
    
    integer :: AS1,AS2,AS3,AS4,i,j,k

    U = SOUND*ma
    nu = U*lside/RE
    nxy = int(lside/h3)

    ! ! *** Verify temporal convergence **********************************************
    ! allocate(solnt1(var,2,nxy+1,nxy+1), STAT = AS1)
    ! allocate(solnt2(var,2,nxy+1,nxy+1), STAT = AS2)
    ! allocate(solnt3(var,2,nxy+1,nxy+1), STAT = AS3)
    ! allocate(solnt4(var,2,nxy+1,nxy+1), STAT = AS4)
    ! if (AS1.ne.0 .or. AS2.ne.0 .or. AS3.ne.0 .or. AS4.ne.0) stop "Error: Not enough memory."  

    ! ! Initialize solution matrix at t=0 (initial zeta, psi, and G values)
    ! call soln_init(h3,nxy,U,var,solnt1,gprofile)
    ! call soln_init(h3,nxy,U,var,solnt2,gprofile)
    ! call soln_init(h3,nxy,U,var,solnt3,gprofile)
    ! call soln_init(h3,nxy,U,var,solnt4,gprofile)


    ! do t_iter = 0,numdt1-1
    !    ! dt = 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001
    !    dt = 1E-8_WP*10.0_WP**t_iter

    !    do t_iter2 = 1,6
    !       call vs_step(dt,t_iter2,h3,nxy,nu,var,solnt1)
    !       call vs_step(dt/2.0_WP,t_iter2,h3,nxy,nu,var,solnt2)
    !    end do

    !    ! Pick arbitrary location (j,k)=(30,30) and sample temporal error
    !    terr1(t_iter+1,1) = dt
    !    terr1(t_iter+1,2) = abs(solnt2(2,2,30,30)-solnt1(2,2,30,30))
    ! end do

    ! do t_iter = 0,numdt2-1
    !    ! dt = 0.0000001, 0.000001, 0.00001, 0.0001
    !    dt = 1E-7_WP*10.0_WP**t_iter

    !    do t_iter2 = 1,6
    !       call gEqn_step(dt,t_iter2,h3,nxy,var,solnt3)
    !       call gEqn_step(dt/2.0_WP,t_iter2,h3,nxy,var,solnt4)
    !    end do

    !    ! Pick arbitrary location (j,k)=(30,30) and sample temporal error
    !    terr2(t_iter+1,1) = dt
    !    terr2(t_iter+1,2) = abs(solnt4(3,2,30,30)-solnt3(3,2,30,30))
    ! end do

    ! ! Write temporal errors to files
    ! open(unit=1,file="../data/timeErrVortSF.txt",status="replace",action="readwrite")
    ! do t_iter = 1,numdt1
    !    write(unit=1,fmt="(e10.3,1x,e10.3)") terr1(t_iter,1), terr1(t_iter,2)
    ! end do
    ! close(unit=1)

    ! open(unit=1,file="../data/timeErrGEqn.txt",status="replace",action="readwrite")
    ! do t_iter = 1,numdt2
    !    write(unit=1,fmt="(e10.3,1x,e10.3)") terr2(t_iter,1), terr2(t_iter,2)
    ! end do
    ! close(unit=1)


    ! *** Verify spatial convergence ***********************************************
    dt = 1E-6_WP
    serr1 = 0.0_WP
    serr2 = 0.0_WP
    
    do x_iter = 0,numh1-1
       if (MOD(x_iter,2).eq.0) then
          h = 0.005_WP*10.0_WP**(x_iter/2)
       else ! x_iter = 1, 3, 5
          h = 0.01_WP*10.0_WP**(x_iter/2)
       end if

       nxy_whole = int(lside/h)
       nxy_half = int(lside*2.0_WP/h)

       allocate(solns1(var,2,nxy_whole+1,nxy_whole+1), STAT = AS1)
       allocate(solns2(var,2,nxy_half+1,nxy_half+1), STAT = AS2)
       if (AS1.ne.0 .or. AS2.ne.0) stop "Error: Not enough memory."

       ! Initialize solution matrix at t=0 (initial zeta, psi, and G values)
       call soln_init(h,nxy_whole,U,var,solns1,gprofile)
       call soln_init(h/2.0_WP,nxy_half,U,var,solns2,gprofile)

       do t_iter = 1,150
          call vs_step(dt,t_iter,h,nxy_whole,nu,var,solns1)
          call vs_step(dt,t_iter,h/2.0_WP,nxy_half,nu,var,solns2)
       end do

       ! Select point in middle of domain & calculate spatial error
       serr1(x_iter+1,1) = h
       serr1(x_iter+1,2) = abs( solns2(2,1,(nxy_half+1)/2+1,(nxy_half+1)/2+1) - &
            solns1(2,1,(nxy_whole+1)/2+1,(nxy_whole+1)/2+1) )

       ! ! Calculate L1 Norm 
       ! serr1(x_iter+1,1) = h

       ! do j=1,nxy_whole+1
       !    do k=1,nxy_whole+1
       !       serr1(x_iter+1,2) = serr1(x_iter+1,2) + abs( solns2(2,2,2*j-1,2*k-1) - &
       !            solns1(2,2,j,k) )/(nxy_whole+1)**2
       !    end do
       ! end do

       deallocate(solns1)
       deallocate(solns2)
       
    end do
    
    ! do x_iter = 0,numh2-1
    !    ! h = 0.01, 0.05, 0.1, 0.5
    !    h = (5.0_WP**mod(x_iter,2))*0.001_WP*10.0_WP**(x_iter/2)

    !    nxy_whole = int(lside/h)
    !    nxy_half = int(lside*2.0_WP/h)

    !    allocate(solns3(var,2,nxy_whole+1,nxy_whole+1), STAT = AS3)
    !    allocate(solns4(var,2,nxy_half+1,nxy_half+1), STAT = AS4)
    !    if (AS3.ne.0 .or. AS4.ne.0) stop "Error: Not enough memory."

    !    ! Initialize solution matrix at t=0 (initial zeta, psi, and G values)
    !    call soln_init(h,nxy_whole,U,var,solns3,gprofile)
    !    call soln_init(h/2.0_WP,nxy_half,U,var,solns4,gprofile)

    !    do t_iter = 1,5
    !       call gEqn_step(dt,t_iter,h,nxy_whole,var,solns3)
    !       call gEqn_step(dt,t_iter,h/2.0_WP,nxy_half,var,solns4)
    !    end do

    !    ! Calculate L1 Norm
    !    serr2(x_iter+1,1) = h

    !    do j=1,nxy_whole+1
    !       do k=1,nxy_whole+1
    !          serr2(x_iter+1,2) = serr2(x_iter+1,2) + abs( solns4(3,2,2*j-1,2*k-1) - &
    !               solns3(3,2,j,k) )/(nxy_whole+1)**2
    !       end do
    !    end do
       
    !    deallocate(solns3)
    !    deallocate(solns4)
       
    ! end do

    ! Write spatial errors to files
    open(unit=1,file="../data/spceErrVortSF.txt",status="replace",action="readwrite")
    do x_iter = 1,numh1
       write(unit=1,fmt="(e10.3,1x,e10.3)") serr1(x_iter,1), serr1(x_iter,2)
    end do
    close(unit=1)

    ! open(unit=1,file="../data/spceErrGEqn.txt",status="replace",action="readwrite")
    ! do x_iter = 1,numh2
    !    write(unit=1,fmt="(e10.3,1x,e10.3)") serr2(x_iter,1), serr2(x_iter,2)
    ! end do
    ! close(unit=1)

  end subroutine convergence_Error
  

  ! ********************************************************************************
  ! *** Conservation Analysis ******************************************************
  ! ********************************************************************************
  subroutine conservation_Error()
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), parameter :: h = 0.015625_WP ! Spatial step size
    real(WP), parameter :: lside = 1.0_WP ! Length of square side
    real(WP), parameter :: dt = 1E-6_WP ! Time step size
    real(WP), parameter :: tfin = 3E-1_WP !1E-3_WP ! Total time of simulation
    real(WP), parameter :: ma = 0.05_WP! Flow mach number
    integer, parameter :: gprofile = 2 ! Specify starting profile of flame

    real(WP), parameter :: U = SOUND*ma ! Bulk velocity [m/s]
    real(WP), parameter :: nu = 0.0_WP ! No kinematic viscosity
    
    integer, parameter :: nxy = int(lside/h) ! Total spatial steps nx=ny
    integer, parameter :: nt = int(tfin/dt) ! Total time steps
    integer :: t_iter,j,k ! Loop iterators

    ! Solution matrix with indices (variable,timestep,xloc,yloc)
    ! where the first index can be the following values
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, parameter :: var = 3
    real(WP), dimension(var,2,nxy+1,nxy+1) :: soln
    real(WP), dimension(var,nxy+1,nxy+1) :: psoln ! Soln at previous t

    ! Variable index for velocity matrix
    !      var2 = 1: u
    !      var2 = 2: v
    integer, parameter :: var2 = 2
    real(WP), dimension(var2,nxy+1,nxy+1) :: uv
    real(WP), dimension(var2,nxy+1,nxy+1) :: puv ! u & v at previous t

    ! Arrays to store quantities for primary and secondary conservation
    real(WP), dimension(1,nt+1) :: primary
    real(WP), dimension(1,nt+1) :: secondary

    character(len=19) :: primaryoutfile
    character(len=21) :: secondaryoutfile
    
    
    ! Initialize solution matrix at t=0 (initial zeta, psi, and G values)
    call soln_init(h,nxy,U,var,soln,gprofile)
    psoln = soln(:,1,:,:)
    
    ! Initialize u & v at t=0
    do j=1,nxy+1
       do k=1,nxy+1
          if (k.eq.1 .or. k.eq.nxy+1) then
             puv(1,j,k) = 0.0_WP
             puv(2,j,k) = 0.0_WP
          else
             puv(1,j,k) = U
             puv(2,j,k) = 0.0_WP
          end if
       end do
    end do

    primary(1,1) = 0.0_WP
    secondary(1,1) = 0.0_WP

    ! Iterate until final time
    do t_iter = 2,nt+1
       call vs_step(dt,t_iter,h,nxy,nu,var,soln)
       call convert_streamfunc(h,nxy,var,soln,var2,uv)
       
       ! primary(1,t_iter) = SUM(sqrt(uv(1,:,:)**2+uv(2,:,:)**2)) - &
       ! SUM(sqrt(puv(1,:,:)**2+puv(2,:,:)**2))
       ! primary(1,t_iter) = SUM(uv(1,:,:)) - SUM(puv(1,:,:))
       primary(1,t_iter) = SUM(uv(1,:,:)) - SUM(puv(1,:,:))
       secondary(1,t_iter) = SUM(uv(1,:,:)**2+uv(2,:,:)**2) - &
            SUM(puv(1,:,:)**2+puv(2,:,:)**2)

       psoln = soln(:,2,:,:)
       puv = uv
    end do

    ! Write primary and secondary conservation values to files
    primaryoutfile = '../data/primary.txt'
    secondaryoutfile = '../data/secondary.txt'

    open(unit=1,file=primaryoutfile,status="replace",action="readwrite")
    do t_iter=1,nt+1
       write(unit=1,fmt="(2(e11.4,1x))") dt*(t_iter-1),primary(1,t_iter)
    end do
    close(unit=1)

    open(unit=1,file=secondaryoutfile,status="replace",action="readwrite")
    do t_iter=1,nt+1
       write(unit=1,fmt="(2(e11.4,1x))") dt*(t_iter-1),secondary(1,t_iter)
    end do
    close(unit=1)
    
  end subroutine conservation_Error
  
  
  ! ********************************************************************************
  ! *** Solve vorticity-stream function approach for a single timestep *************
  ! ********************************************************************************
  subroutine vs_step(dt,t_iter,h,nxy,nu,var,soln)
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: t_iter ! Current timestep index
    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), intent(in) :: nu ! Kinematic viscosity [m^2/s]

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln ! Solution matrix
    
    ! Variable index for solution matrix
    !      var2 = 1: u
    !      var2 = 2: v
    integer, parameter :: var2 = 2
    real(WP), dimension(var2,nxy+1,nxy+1) :: uv ! Values of u & v

    real(WP), dimension(nxy+1,nxy+1) :: psip ! Psi at t=t_iter
    
    integer :: i,j,k,subitercount ! Loop iterators
    real(WP), parameter :: tol = 0.001_WP !0.54_WP ! Tolerance
    real(WP) :: residual ! Residual
    real(WP), parameter :: beta = 1.5_WP !0.9_WP ! Coefficient in Successive Over Relaxation

    call OMP_SET_NUM_THREADS(NTHREADS)
    
    ! *** Substep 1: Iterate for new psi values ***
    ! Initialize subiteration counter and residual
    subitercount = 0
    residual = 1.0_WP
    do while (residual.gt.tol)
       residual = 0.0_WP
       
       ! Store psi values at t=t_iter
       psip = soln(2,1,:,:)
       
       ! Solve for psi at subiteration subitercount+1 with S.O.R.
       do j=1,nxy+1
          do k=2,nxy
             if (j.eq.1) then
                ! *** Periodic inflow B.C. ***
                soln(2,1,j,k) = 0.25_WP*beta*( soln(2,1,j+1,k)+soln(2,1,nxy+1,k) + &
                     soln(2,1,j,k+1) + soln(2,1,j,k-1) + h**2*soln(1,1,j,k) ) + &
                     (1.0_WP-beta)*soln(2,1,j,k)
                
             elseif (j.eq.nxy+1) then
                ! *** Periodic outflow B.C. ***
                soln(2,1,j,k) = 0.25_WP*beta*( soln(2,1,1,k) + soln(2,1,j-1,k) + &
                     soln(2,1,j,k+1) + soln(2,1,j,k-1) + h**2*soln(1,1,j,k) ) + &
                     (1.0_WP-beta)*soln(2,1,j,k)
                
             else
                soln(2,1,j,k) = 0.25_WP*beta*( soln(2,1,j+1,k) + soln(2,1,j-1,k) + &
                     soln(2,1,j,k+1) + soln(2,1,j,k-1) + h**2*soln(1,1,j,k) ) + &
                     (1.0_WP-beta)*soln(2,1,j,k)
                
             end if
          end do
       end do

       residual = sum(abs(soln(2,1,:,:) - psip))/sum(abs(psip))
       subitercount = subitercount + 1
       print *, "t_iter=", t_iter, "Subiteration=", subitercount, &
            "Residual=", residual
    end do
    
    ! Store converged psi (stream function) solution at t=t_iter+1
    soln(2,2,:,:) = soln(2,1,:,:)
    

    ! *** Substep 2: Determine zeta (vorticity) B.C.s using inner psi & zeta ***
    !$OMP PARALLEL
    
    !$OMP DO
    do j=1,nxy+1
       k = 1 ! Bottom wall
       soln(1,2,j,k) = 2.0_WP/h**2*(soln(2,2,j,k) - soln(2,2,j,k+1))
       
       k = nxy+1 ! Top wall
       soln(1,2,j,k) = 2.0_WP/h**2*(soln(2,2,j,k) - soln(2,2,j,k-1))
    end do
    !$OMP END DO   

    !$OMP DO
    do k=2,nxy
       j = 1 ! Left inflow

       ! *** Non-periodic inflow B.C. ***
       soln(1,2,j,k) = -(soln(2,2,j+2,k)-2.0_WP*soln(2,2,j+1,k)+soln(2,2,j,k))/ &
            h**2 - (soln(2,2,j,k+1)-2.0_WP*soln(2,2,j,k)+soln(2,2,j,k-1))/h**2
       
       j = nxy+1 ! Right outflow

       ! *** Non-periodic outflow B.C. ***
       soln(1,2,j,k) = -(soln(2,2,j,k)-2.0_WP*soln(2,2,j-1,k)+soln(2,2,j-2,k))/h**2 &
            - (soln(2,2,j,k+1)-2.0_WP*soln(2,2,j,k)+soln(2,2,j,k-1))/h**2
       
    end do
    !$OMP END DO

    
    ! *** Substep 3: Solve vorticity transport equation for inner cells ***
    !$OMP DO
    do j=2,nxy
       do k=2,nxy
          soln(1,2,j,k) = soln(1,1,j,k) + dt/(4.0_WP*h**2) * ( &
               -(soln(2,1,j,k+1)-soln(2,1,j,k-1))*(soln(1,1,j+1,k)-soln(1,1,j-1,k)) &
               +(soln(2,1,j+1,k)-soln(2,1,j-1,k))*(soln(1,1,j,k+1)-soln(1,1,j,k-1)) &
               +4.0_WP*nu*(soln(1,1,j+1,k) + soln(1,1,j-1,k) + soln(1,1,j,k+1) &
               +soln(1,1,j,k-1) - 4.0_WP*soln(1,1,j,k)) )
       end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

    ! Store solution at t=t_iter+1 for next timestep
    soln(1,1,:,:) = soln(1,2,:,:)
    
  end subroutine vs_step

  ! *********************************************************************************
  ! *** Solve the G-Equation over a single timestep *********************************
  ! *********************************************************************************
  subroutine gEqn_step(dt,t_iter,h,nxy,var,soln)
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: t_iter ! Current timestep index
    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln ! Solution matrix
    
    ! Variable index for solution matrix
    !      var2 = 1: u
    !      var2 = 2: v
    integer, parameter :: var2 = 2
    real(WP), dimension(var2,nxy+1,nxy+1) :: uv ! Values of u & v

    ! Store G values at timestep t_iter (the matrix soln(3,1,:,:) is duplicated so
    ! that I can type fewer characters in the expressions below
    real(WP), dimension(nxy+1,nxy+1) :: g

    real(WP) :: rhs1,rhs2,rhs3 ! Variables to represent terms on RHS of G-Equation

    ! Variables used in reinitialization process
    real(WP) :: tol, residual ! Tolerance & residual
    integer :: subiter ! Counter for reinitialization iterations
    real(WP), parameter :: epsilon = 0.01_WP ! For smoothing of sgn function approx.
    real(WP), dimension(nxy+1,nxy+1) :: gthalf ! G at t=t_iter+1/2
    
    integer :: i,j,k ! Loop iterators

    
    ! Convert stream function into velocities
    call convert_streamfunc(h,nxy,var,soln,var2,uv)

    ! Copy G values at timestep t_iter
    g(:,:) = soln(3,1,:,:)

    ! Set G values at top and bottom boundaries using grad(G) = <1,1>
    do j=1,nxy+1
       k = 1
       ! rhs1 = dt*SLNOT
       ! rhs2 = 0.0_WP
       ! rhs3 = -dt*ML/(2.0_WP*h)*( -3.0_WP*uv(1,j,k)+4.0_WP*uv(1,j,k+1) - &
       !            uv(1,j,k+2) - 3.0_WP*uv(2,j,k)+4.0_WP*uv(2,j,k+1)-uv(2,j,k+2) )
       ! soln(3,1,j,k) = g(j,k)-dt*(uv(1,j,k)+uv(2,j,k)) + rhs1 - rhs2 - rhs3

       soln(3,1,j,k) = soln(3,1,j,k+2) - 2.0_WP*h
       
       k = nxy+1
       ! rhs1 = dt*SLNOT
       ! rhs2 = 0.0_WP
       ! rhs3 = -dt*ML/(2.0_WP*h)*( 3.0_WP*uv(1,j,k)-4.0_WP*uv(1,j,k-1) + &
       !            uv(1,j,k-2) + 3.0_WP*uv(2,j,k)-4.0_WP*uv(2,j,k-1)+uv(2,j,k-2) )
       ! soln(3,1,j,k) = g(j,k)-dt*(uv(1,j,k)+uv(2,j,k)) + rhs1 - rhs2 - rhs3

       soln(3,1,j,k) = soln(3,1,j,k-2) + 2.0_WP*h
    end do

    ! Set G values in interior of domain
    !$OMP PARALLEL DO PRIVATE(j,k,rhs1,rhs2,rhs3)
    do j=1,nxy+1
       do k=2,nxy

          if (j.eq.1) then ! Left boundary with grad(G) = <1,1>
             
             ! SLNOT * |grad G|
             ! ! *** Centered difference version, x & y-dir ***
             ! rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(j+1,k)-g(nxy+1,k))**2 + &
             !      (g(j,k+1)-g(j,k-1))**2 )

             ! *** Centered difference version, x-dir only ***
             ! rhs1 = dt*SLNOT*(g(j+1,k)-g(nxy+1,k))/(2.0_WP*h)

             rhs1 = dt*SLNOT*sqrt(2.0_WP)

             ! ! *** Upwind/downwind version ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(nxy+1,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(nxy+1,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if
             ! endif

             ! SLNOT * ML * K * |grad G|
             ! rhs2 = -dt*SLNOT*ML/(4.0_WP*h**2)*(    &
             !      4.0_WP*(g(j+1,k)+g(nxy+1,k) + g(j,k+1)+g(j,k-1)-4.0_WP*g(j,k))- &
             !      1.0_WP/sqrt((g(j+1,k)-g(nxy+1,k))**2+(g(j,k+1)-g(j,k-1))**2)*(  &
             !      (g(j+1,k)-g(nxy+1,k))*( sqrt(16.0_WP*(g(j+1,k)-g(j,k))**2 + &
             !      (g(j,k+1)+g(j+1,k+1)-g(j,k-1)-g(j+1,k-1))**2) - &
             !      sqrt(16.0_WP*(g(j,k)-g(nxy+1,k))**2 + &
             !      (g(nxy+1,k+1)+g(j,k+1)-g(nxy+1,k-1)-g(j,k-1))**2) ) + &
             !      (g(j,k+1)-g(j,k-1)) * (  &
             !      sqrt( (g(j+1,k+1)+g(j+1,k)-g(nxy+1,k+1)-g(nxy+1,k))**2 + &
             !      16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
             !      sqrt( (g(j+1,k)+g(j+1,k-1)-g(nxy+1,k)-g(nxy+1,k-1))**2 + &
             !      16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             rhs2 = 0.0_WP

             ! ML * S * |grad G|
             ! rhs3 = -dt*ML/( 4.0_WP*h**2*sqrt((g(j+1,k)-g(nxy+1,k))**2+(g(j,k+1)- &
             !      g(j,k-1))**2) )*(    (g(j+1,k)-g(nxy+1,k))*(  4.0_WP*uv(1,j,k)* &
             !      (g(j+1,k)-2.0_WP*g(j,k)+g(nxy+1,k)) + (g(j+1,k)-g(nxy+1,k)) * &
             !      (uv(1,j+1,k)-uv(1,nxy+1,k)) + uv(2,j,k)*(g(j+1,k+1)-g(j+1,k-1)- &
             !      g(nxy+1,k+1)+g(nxy+1,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j+1,k) - &
             !      uv(2,nxy+1,k))-1.0_WP/sqrt((g(j+1,k)-g(nxy+1,k))**2+(g(j,k+1) - &
             !      g(j,k-1))**2)*(uv(1,j,k)*(g(j+1,k)-g(nxy+1,k)) + uv(2,j,k) * &
             !      (g(j,k+1)-g(j,k-1)))*( sqrt(16.0_WP*(g(j+1,k)-g(j,k))**2 + &
             !      (g(j,k+1)+g(j+1,k+1)-g(j,k-1)-g(j+1,k-1))**2) - &
             !      sqrt(16.0_WP*(g(j,k)-g(nxy+1,k))**2 + &
             !      (g(nxy+1,k+1)+g(j,k+1)-g(nxy+1,k-1)-g(j,k-1))**2) )  ) + &
             !      (g(j,k+1)-g(j,k-1))*(   uv(1,j,k) * &
             !      (g(j+1,k+1)-g(j+1,k-1)-g(nxy+1,k+1)+g(nxy+1,k-1)) + (g(j+1,k) - &
             !      g(nxy+1,k))*(uv(1,j,k+1)-uv(1,j,k-1))+4.0_WP*uv(2,j,k)*(g(j,k+1)- &
             !      2.0_WP*g(j,k)+g(j,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j,k+1) - &
             !      uv(2,j,k-1)) - 1.0_WP/sqrt((g(j+1,k)-g(nxy+1,k))**2 +(g(j,k+1)- &
             !      g(j,k-1))**2)*(uv(1,j,k)*(g(j+1,k)-g(nxy+1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1))) * (  &
             !      sqrt( (g(j+1,k+1)+g(j+1,k)-g(nxy+1,k+1)-g(nxy+1,k))**2 + &
             !      16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
             !      sqrt( (g(j+1,k)+g(j+1,k-1)-g(nxy+1,k)-g(nxy+1,k-1))**2 + &
             !      16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             rhs3 = -dt*ML/(2.0_WP*h)*( -3.0_WP*uv(1,j,k)+4.0_WP*uv(1,j+1,k) - &
                  uv(1,j+2,k) + uv(1,j,k+1) - uv(1,j,k-1) - 3.0_WP*uv(2,j,k) + &
                  4.0_WP*uv(2,j+1,k) - uv(2,j+2,k) + uv(2,j,k+1) - uv(2,j,k-1) )

             ! Substitute all terms in G-Eqn to obtain solution at t=t_iter+1/2
             ! ! *** Centered difference version, x & y-dir ***
             ! soln(3,1,j,k) = g(j,k)-dt/(2.0_WP*h)*( &
             !      uv(1,j,k)*(g(j+1,k)-g(nxy+1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3

             ! *** Centered difference version, x-dir only ***
             ! soln(3,1,j,k) =g(j,k)-dt*uv(1,j,k)*(g(j+1,k)-g(nxy+1,k))/(2.0_WP*h)+ &
             !      rhs1 !- rhs2 - rhs3

             soln(3,1,j,k) = g(j,k)-dt*(uv(1,j,k)+uv(2,j,k)) + rhs1 - rhs2 - rhs3
             
             ! ! *** Upwind/downwind version ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j+1,k)-g(j,k)) + &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j+1,k)-g(j,k)) + &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(nxy+1,k))+ &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(nxy+1,k))+ &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if
             ! endif
             
          elseif (j.eq.nxy+1) then ! Right boundary with grad(G) = <1,1>
             
             ! SLNOT * |grad G|
             ! ! *** Centered difference version, x & y-dir ***
             ! rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(1,k)-g(j-1,k))**2 + &
             !      (g(j,k+1)-g(j,k-1))**2 )

             ! *** Centered difference version, x-dir only ***
             ! rhs1 = dt*SLNOT*(g(1,k)-g(j-1,k))/(2.0_WP*h)

             rhs1 = dt*SLNOT*sqrt(2.0_WP)

             ! ! *** Upwind/downwind version ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(1,k)-g(j,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(1,k)-g(j,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if
             ! endif

             ! SLNOT * ML * K * |grad G|
             ! rhs2 = -dt*SLNOT*ML/(4.0_WP*h**2)*(    &
             !      4.0_WP*(g(1,k)+g(j-1,k) + g(j,k+1)+g(j,k-1)-4.0_WP*g(j,k)) - &
             !      1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1)-g(j,k-1))**2)*(  &
             !      (g(1,k)-g(j-1,k))*( sqrt(16.0_WP*(g(1,k)-g(j,k))**2 + &
             !      (g(j,k+1)+g(1,k+1)-g(j,k-1)-g(1,k-1))**2) - &
             !      sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
             !      (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) ) + &
             !      (g(j,k+1)-g(j,k-1)) * (  &
             !      sqrt( (g(1,k+1)+g(1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
             !      16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
             !      sqrt( (g(1,k)+g(1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
             !      16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             rhs2 = 0.0_WP

             ! ML * S * |grad G|
             ! rhs3 = -dt*ML/( 4.0_WP*h**2*sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
             !      g(j,k-1))**2) )*(    (g(1,k)-g(j-1,k))*(  4.0_WP*uv(1,j,k) * &
             !      (g(1,k)-2.0_WP*g(j,k)+g(j-1,k)) + (g(1,k)-g(j-1,k)) * &
             !      (uv(1,1,k)-uv(1,j-1,k)) + uv(2,j,k)*(g(1,k+1)-g(1,k-1) - &
             !      g(j-1,k+1)+g(j-1,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,1,k) - &
             !      uv(2,j-1,k)) - 1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
             !      g(j,k-1))**2)*(uv(1,j,k)*(g(1,k)-g(j-1,k)) + uv(2,j,k) * &
             !      (g(j,k+1)-g(j,k-1)))*( sqrt(16.0_WP*(g(1,k)-g(j,k))**2 + &
             !      (g(j,k+1)+g(1,k+1)-g(j,k-1)-g(1,k-1))**2) - &
             !      sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
             !      (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) )  ) + &
             !      (g(j,k+1)-g(j,k-1))*(   uv(1,j,k) * &
             !      (g(1,k+1)-g(1,k-1)-g(j-1,k+1)+g(j-1,k-1)) + (g(1,k) - &
             !      g(j-1,k))*(uv(1,j,k+1)-uv(1,j,k-1))+4.0_WP*uv(2,j,k)*(g(j,k+1)- &
             !      2.0_WP*g(j,k)+g(j,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j,k+1) - &
             !      uv(2,j,k-1)) - 1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
             !      g(j,k-1))**2)*(uv(1,j,k)*(g(1,k)-g(j-1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1))) * (  &
             !      sqrt( (g(1,k+1)+g(1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
             !      16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
             !      sqrt( (g(1,k)+g(1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
             !      16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             rhs3 = -dt*ML/(2.0_WP*h)*( 3.0_WP*uv(1,j,k) - 4.0_WP*uv(1,j-1,k) + &
                  uv(1,j-2,k) + uv(1,j,k+1) - uv(1,j,k-1) + 3.0_WP*uv(2,j,k) - &
                  4.0_WP*uv(2,j-1,k) + uv(2,j-2,k) + uv(2,j,k+1) - uv(2,j,k-1) )

             ! Substitute all terms in G-Eqn to obtain solution at t=t_iter+1/2
             ! ! *** Centered difference version, x & y-dir ***
             ! soln(3,1,j,k) = g(j,k)-dt/(2.0_WP*h)*( uv(1,j,k)*(g(1,k)-g(j-1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3

             ! *** Centered difference version, x-dir only ***
             ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k)*(g(1,k)-g(j-1,k))/(2.0_WP*h) + &
             !      rhs1 !- rhs2 - rhs3

             soln(3,1,j,k) = g(j,k)-dt*(uv(1,j,k)+uv(2,j,k)) + rhs1 - rhs2 - rhs3
             
             ! ! *** Upwind/downwind version ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(1,k)-g(j,k)) + &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(1,k)-g(j,k)) + &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if
             ! endif
             
          else ! Inner grid points
             
             ! SLNOT * |grad G|
             ! *** Centered difference version, x & y-dir ***
             rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(j+1,k)-g(j-1,k))**2 + &
                  (g(j,k+1)-g(j,k-1))**2 )

             ! ! ! *** Centered difference version, x-dir only ***
             ! ! rhs1 = dt*SLNOT*(g(j+1,k)-g(j-1,k))/(2.0_WP*h)

             ! ! *** Upwind/downwind version, x-dir only ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if
             ! endif

             ! ! *** Upwind/downwind version ***
             ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j+1,k)-g(j,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if          
             ! else ! u > 0
             !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k+1)-g(j,k))**2 )
             !    else ! v > 0
             !       rhs1 = dt*SLNOT/h * sqrt( (g(j,k)-g(j-1,k))**2 + &
             !            (g(j,k)-g(j,k-1))**2 )
             !    end if
             ! endif

             ! SLNOT * ML * K * |grad G|
             rhs2 = -dt*SLNOT*ML/(4.0_WP*h**2)*(    &
                  4.0_WP*(g(j+1,k)+g(j-1,k) + g(j,k+1)+g(j,k-1)-4.0_WP*g(j,k)) - &
                  1.0_WP/sqrt((g(j+1,k)-g(j-1,k))**2 + (g(j,k+1)-g(j,k-1))**2)*(  &
                  (g(j+1,k)-g(j-1,k))*( sqrt(16.0_WP*(g(j+1,k)-g(j,k))**2 + &
                  (g(j,k+1)+g(j+1,k+1)-g(j,k-1)-g(j+1,k-1))**2) - &
                  sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
                  (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) ) + &
                  (g(j,k+1)-g(j,k-1)) * (  &
                  sqrt( (g(j+1,k+1)+g(j+1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
                  16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
                  sqrt( (g(j+1,k)+g(j+1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
                  16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             ! ML * S * |grad G|
             rhs3 = -dt*ML/( 4.0_WP*h**2*sqrt((g(j+1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2) )*(    (g(j+1,k)-g(j-1,k))*(  4.0_WP*uv(1,j,k) * &
                  (g(j+1,k)-2.0_WP*g(j,k)+g(j-1,k)) + (g(j+1,k)-g(j-1,k)) * &
                  (uv(1,j+1,k)-uv(1,j-1,k)) + uv(2,j,k)*(g(j+1,k+1)-g(j+1,k-1) - &
                  g(j-1,k+1)+g(j-1,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j+1,k) - &
                  uv(2,j-1,k)) - 1.0_WP/sqrt((g(j+1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2)*(uv(1,j,k)*(g(j+1,k)-g(j-1,k)) + uv(2,j,k) * &
                  (g(j,k+1)-g(j,k-1)))*( sqrt(16.0_WP*(g(j+1,k)-g(j,k))**2 + &
                  (g(j,k+1)+g(j+1,k+1)-g(j,k-1)-g(j+1,k-1))**2) - &
                  sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
                  (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) )  ) + &
                  (g(j,k+1)-g(j,k-1))*(   uv(1,j,k) * &
                  (g(j+1,k+1)-g(j+1,k-1)-g(j-1,k+1)+g(j-1,k-1)) + (g(j+1,k) - &
                  g(j-1,k))*(uv(1,j,k+1)-uv(1,j,k-1)) + 4.0_WP*uv(2,j,k)*(g(j,k+1)- &
                  2.0_WP*g(j,k)+g(j,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j,k+1) - &
                  uv(2,j,k-1)) - 1.0_WP/sqrt((g(j+1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2)*(uv(1,j,k)*(g(j+1,k)-g(j-1,k)) + &
                  uv(2,j,k)*(g(j,k+1)-g(j,k-1))) * (  &
                  sqrt( (g(j+1,k+1)+g(j+1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
                  16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
                  sqrt( (g(j+1,k)+g(j+1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
                  16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             ! Substitute all terms in G-Eqn to obtain solution at t=t_iter+1/2
             ! ! *** Centered difference version, x & y-dir ***
             ! soln(3,1,j,k) = g(j,k)-dt/(2.0_WP*h)*( uv(1,j,k)*(g(j+1,k)-g(j-1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 - rhs2 - rhs3

             ! ! ! *** Centered difference version, x-dir only ***
             ! ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k)*(g(j+1,k)-g(j-1,k))/(2.0_WP*h) + &
             ! !      rhs1 !- rhs2 - rhs3

             ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k) + rhs1 !- rhs2 - rhs3

             ! *** Upwind/downwind version ***
             if(uv(1,j,k).le.0.0_WP) then ! u <= 0
                if(uv(2,j,k).le.0.0_WP) then ! v <= 0
                   soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j+1,k)-g(j,k)) + &
                        uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 - rhs2 - rhs3
                else ! v > 0
                   soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j+1,k)-g(j,k)) + &
                        uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 - rhs2 - rhs3
                end if          
             else ! u > 0
                if(uv(2,j,k).le.0.0_WP) then ! v <= 0
                   soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
                        uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 - rhs2 - rhs3
                else ! v > 0
                   soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
                        uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 - rhs2 - rhs3
                end if
             endif

          end if
          
       end do
    end do
    !$OMP END PARALLEL DO

    
    ! *** Reinitialize distances to flame *******************************************
    subiter = 0
    tol = 0.005_WP !0.5_WP ! Change for convergence stuff
    residual = 1.0_WP
    
    do while (residual.gt.tol)
       residual = 0.0_WP

       ! Store G values at t=t_iter+1/2
       gthalf = soln(3,1,:,:)

       ! Set g values at top and bottom boundaries using grad(g) = <1,1>
       do j=1,nxy+1
          k = 1
          soln(3,1,j,k) = gthalf(j,k+2) - 2.0_WP*h
       
          k = nxy+1
          soln(3,1,j,k) = gthalf(j,k-2) + 2.0_WP*h
       end do
       
       ! Solve for G at subiteration subitercount+1
       do j=1,nxy+1
          do k=2,nxy
             if (j.eq.1) then
                ! ! *** Periodic inflow B.C., x & y-dir ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - 0.5_WP/h*sqrt( (gthalf(j+1,k)-gthalf(nxy+1,k))**2 &
                !      + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! ! *** Periodic inflow B.C., x-dir only ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - (gthalf(j+1,k)-gthalf(nxy+1,k))/(2.0_WP*h) )

                ! ! *** Non-periodic inflow B.C. ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                     ! dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                     ! 1.0_WP - 0.5_WP/h*sqrt( (-3.0_WP*gthalf(j,k) + &
                     ! 4.0_WP*gthalf(j+1,k) - gthalf(j+2,k))**2 &
                ! + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! Left boundary where grad(g) = <1,1>
                soln(3,1,j,k) = gthalf(j,k) + &
                     dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(1.0_WP - &
                     sqrt(2.0_WP))
                
             elseif (j.eq.nxy+1) then
                ! ! *** Periodic outflow B.C., x & y-dir ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - 0.5_WP/h*sqrt( (gthalf(1,k)-gthalf(j-1,k))**2 &
                !      + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! ! *** Periodic inflow B.C., x-dir only ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - (gthalf(1,k)-gthalf(j-1,k))/(2.0_WP*h) )

                ! ! *** Non-periodic inflow B.C. ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                     ! dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                     ! 1.0_WP - 0.5_WP/h*sqrt( (3.0_WP*gthalf(j,k) - &
                     ! 4.0_WP*gthalf(j-1,k)+gthalf(j-2,k))**2 &
                ! + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! Right boundary where grad(g) = <1,1>
                soln(3,1,j,k) = gthalf(j,k) + &
                     dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(1.0_WP - &
                     sqrt(2.0_WP))
                
             else
                ! *** Inner cells, x & y-dir ***
                soln(3,1,j,k) = gthalf(j,k) + &
                     dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                     1.0_WP - 0.5_WP/h*sqrt( (gthalf(j+1,k)-gthalf(j-1,k))**2 &
                     + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! ! *** Inner cells, x-dir only ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - (gthalf(j+1,k)-gthalf(j-1,k))/(2.0_WP*h) )

                ! ! *** Upwind/downwind version ***
                ! if(uv(1,j,k).le.0.0_WP) then ! u <= 0
                !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
                !       soln(3,1,j,k) = gthalf(j,k) + &
                !            dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !            1.0_WP - 1.0_WP/h*sqrt( (gthalf(j+1,k)-gthalf(j,k))**2 &
                !            + (gthalf(j,k+1)-gthalf(j,k))**2 ) )
                !    else ! v > 0
                !       soln(3,1,j,k) = gthalf(j,k) + &
                !            dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !            1.0_WP - 1.0_WP/h*sqrt( (gthalf(j+1,k)-gthalf(j,k))**2 &
                !            + (gthalf(j,k)-gthalf(j,k-1))**2 ) )
                !    end if
                ! else ! u > 0
                !    if(uv(2,j,k).le.0.0_WP) then ! v <= 0
                !       soln(3,1,j,k) = gthalf(j,k) + &
                !            dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !            1.0_WP - 1.0_WP/h*sqrt( (gthalf(j,k)-gthalf(j-1,k))**2 &
                !            + (gthalf(j,k+1)-gthalf(j,k))**2 ) )
                !    else ! v > 0
                !       soln(3,1,j,k) = gthalf(j,k) + &
                !            dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !            1.0_WP - 1.0_WP/h*sqrt( (gthalf(j,k)-gthalf(j-1,k))**2 &
                !            + (gthalf(j,k)-gthalf(j,k-1))**2 ) )
                !    end if
                ! endif
   
             end if
          end do
       end do

       residual = sum(abs(soln(3,1,:,:) - gthalf))/sum(abs(gthalf))
       subiter = subiter + 1
       print *, "t_iter=", t_iter, "Subiteration=", subiter, &
            "Residual=", residual
    end do
    
    ! Store converged G values at t=t_iter+1
    soln(3,2,:,:) = soln(3,1,:,:)
    
  end subroutine gEqn_step
  

end program GEQN_SOLVER
