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
  real(WP), parameter :: tfin1 = 2E-2_WP !2E-3_WP !1E-1_WP !2.0_WP !1.0_WP ! [s]

  ! Mach number based on bulk velocity
  real(WP), parameter :: ma1 = 0.05_WP

  ! Will need to change velocity ranges in plotting script if select Mach #s below
  ! real(WP), parameter :: ma2 = 0.1_WP
  ! real(WP), parameter :: ma3 = 0.3_WP

  ! Initial profile of flame
  integer, parameter :: gprofile1 = 1 ! Planar
  integer, parameter :: gprofile2 = 2 ! Sinusoidal
  ! integer, parameter :: gprofile3 = 3 ! Not implemented yet

  ! For timing purposes
  integer :: t1,t2,clock_rate,clock_max

  ! Start timer
  call system_clock (t1,clock_rate,clock_max)

  ! Set number of OpenMP threads
  ! call OMP_SET_NUM_THREADS(NTHREADS)

  ! Solve for incompressible flow field with iso-surface from G-Equation
  call solveAll(h3,lside1,dt1,tfin1,ma1,gprofile2)

  ! Stop timer
  call system_clock (t2,clock_rate,clock_max)
  print *, 'Elapsed real time =', real(t2-t1)/real(clock_rate), 's'

  
contains
  

  ! *** Incompressible Navier-Stokes solver with incorporated G-Equation ************
  
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

  ! *********************************************************************************
  ! *** Solve vorticity-stream function approach for a single timestep **************
  ! *********************************************************************************
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
    ! subitercount = 0
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

       if (t_iter.le.2) then ! First and second timesteps where solution is 0
          residual = sum(abs(soln(2,1,:,:) - psip))
       else
          residual = sum(abs(soln(2,1,:,:) - psip))/sum(abs(psip))
       end if
       
       ! subitercount = subitercount + 1
       ! print *, "t_iter=", t_iter, "Subiteration=", subitercount, &
            ! "Residual=", residual
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

    !$OMP PARALLEL DO PRIVATE(j,k,rhs1,rhs2,rhs3)
    do j=1,nxy+1
       do k=2,nxy

          if (j.eq.1) then ! Left periodic boundary
             
             ! SLNOT * |grad G|
             ! ! *** Centered difference version, x & y-dir ***
             ! rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(j+1,k)-g(nxy+1,k))**2 + &
             !      (g(j,k+1)-g(j,k-1))**2 )

             ! *** Centered difference version, x-dir only ***
             ! rhs1 = dt*SLNOT*(g(j+1,k)-g(nxy+1,k))/(2.0_WP*h)

             rhs1 = dt*SLNOT

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
                  uv(1,j+2,k) + uv(1,j,k+1)-uv(1,j,k-1) )

             ! Substitute all terms in G-Eqn to obtain solution at t=t_iter+1/2
             ! ! *** Centered difference version, x & y-dir ***
             ! soln(3,1,j,k) = g(j,k)-dt/(2.0_WP*h)*( &
             !      uv(1,j,k)*(g(j+1,k)-g(nxy+1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3

             ! *** Centered difference version, x-dir only ***
             ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k)*(g(j+1,k)-g(nxy+1,k))/(2.0_WP*h) + &
             !      rhs1 !- rhs2 - rhs3

             soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k) + rhs1 !- rhs2 - rhs3

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
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(nxy+1,k)) + &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(nxy+1,k)) + &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if
             ! endif
             
          elseif (j.eq.nxy+1) then ! Right periodic boundary
             
             ! SLNOT * |grad G|
             ! ! *** Centered difference version, x & y-dir ***
             ! rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(1,k)-g(j-1,k))**2 + &
             !      (g(j,k+1)-g(j,k-1))**2 )

             ! *** Centered difference version, x-dir only ***
             ! rhs1 = dt*SLNOT*(g(1,k)-g(j-1,k))/(2.0_WP*h)

             rhs1 = dt*SLNOT

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
             rhs2 = -dt*SLNOT*ML/(4.0_WP*h**2)*(    &
                  4.0_WP*(g(1,k)+g(j-1,k) + g(j,k+1)+g(j,k-1)-4.0_WP*g(j,k)) - &
                  1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1)-g(j,k-1))**2)*(  &
                  (g(1,k)-g(j-1,k))*( sqrt(16.0_WP*(g(1,k)-g(j,k))**2 + &
                  (g(j,k+1)+g(1,k+1)-g(j,k-1)-g(1,k-1))**2) - &
                  sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
                  (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) ) + &
                  (g(j,k+1)-g(j,k-1)) * (  &
                  sqrt( (g(1,k+1)+g(1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
                  16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
                  sqrt( (g(1,k)+g(1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
                  16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             ! ML * S * |grad G|
             rhs3 = -dt*ML/( 4.0_WP*h**2*sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2) )*(    (g(1,k)-g(j-1,k))*(  4.0_WP*uv(1,j,k) * &
                  (g(1,k)-2.0_WP*g(j,k)+g(j-1,k)) + (g(1,k)-g(j-1,k)) * &
                  (uv(1,1,k)-uv(1,j-1,k)) + uv(2,j,k)*(g(1,k+1)-g(1,k-1) - &
                  g(j-1,k+1)+g(j-1,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,1,k) - &
                  uv(2,j-1,k)) - 1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2)*(uv(1,j,k)*(g(1,k)-g(j-1,k)) + uv(2,j,k) * &
                  (g(j,k+1)-g(j,k-1)))*( sqrt(16.0_WP*(g(1,k)-g(j,k))**2 + &
                  (g(j,k+1)+g(1,k+1)-g(j,k-1)-g(1,k-1))**2) - &
                  sqrt(16.0_WP*(g(j,k)-g(j-1,k))**2 + &
                  (g(j-1,k+1)+g(j,k+1)-g(j-1,k-1)-g(j,k-1))**2) )  ) + &
                  (g(j,k+1)-g(j,k-1))*(   uv(1,j,k) * &
                  (g(1,k+1)-g(1,k-1)-g(j-1,k+1)+g(j-1,k-1)) + (g(1,k) - &
                  g(j-1,k))*(uv(1,j,k+1)-uv(1,j,k-1)) + 4.0_WP*uv(2,j,k)*(g(j,k+1)- &
                  2.0_WP*g(j,k)+g(j,k-1)) + (g(j,k+1)-g(j,k-1))*(uv(2,j,k+1) - &
                  uv(2,j,k-1)) - 1.0_WP/sqrt((g(1,k)-g(j-1,k))**2 + (g(j,k+1) - &
                  g(j,k-1))**2)*(uv(1,j,k)*(g(1,k)-g(j-1,k)) + &
                  uv(2,j,k)*(g(j,k+1)-g(j,k-1))) * (  &
                  sqrt( (g(1,k+1)+g(1,k)-g(j-1,k+1)-g(j-1,k))**2 + &
                  16.0_WP*(g(j,k+1)-g(j,k))**2 ) - &
                  sqrt( (g(1,k)+g(1,k-1)-g(j-1,k)-g(j-1,k-1))**2 + &
                  16.0_WP*(g(j,k)-g(j,k-1))**2 )  )   )    )

             ! Substitute all terms in G-Eqn to obtain solution at t=t_iter+1/2
             ! ! *** Centered difference version, x & y-dir ***
             ! soln(3,1,j,k) = g(j,k)-dt/(2.0_WP*h)*( uv(1,j,k)*(g(1,k)-g(j-1,k)) + &
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3

             ! *** Centered difference version, x-dir only ***
             ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k)*(g(1,k)-g(j-1,k))/(2.0_WP*h) + &
             !      rhs1 !- rhs2 - rhs3

             soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k) + rhs1 !- rhs2 - rhs3
             
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
             ! ! *** Centered difference version, x & y-dir ***
             ! rhs1 = dt*SLNOT/(2.0_WP*h) * sqrt( (g(j+1,k)-g(j-1,k))**2 + &
             !      (g(j,k+1)-g(j,k-1))**2 )

             ! *** Centered difference version, x-dir only ***
             rhs1 = dt*SLNOT*(g(j+1,k)-g(j-1,k))/(2.0_WP*h)

             ! rhs1 = dt*SLNOT

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

             ! *** Upwind/downwind version ***
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
             !      uv(2,j,k)*(g(j,k+1)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3

             ! *** Centered difference version, x-dir only ***
             soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k)*(g(j+1,k)-g(j-1,k))/(2.0_WP*h) + &
                  rhs1 !- rhs2 - rhs3

             ! soln(3,1,j,k) = g(j,k)-dt*uv(1,j,k) + rhs1 !- rhs2 - rhs3

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
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
             !            uv(2,j,k)*(g(j,k+1)-g(j,k)) ) + rhs1 !- rhs2 - rhs3
             !    else ! v > 0
             !       soln(3,1,j,k) = g(j,k) - dt/h*( uv(1,j,k)*(g(j,k)-g(j-1,k)) + &
             !            uv(2,j,k)*(g(j,k)-g(j,k-1)) ) + rhs1 !- rhs2 - rhs3
             !    end if
             ! endif

          end if
          
       end do
    end do
    !$OMP END PARALLEL DO

    
    ! *** Reinitialize distances to flame *******************************************
    subiter = 0
    tol = 0.005_WP
    residual = 1.0_WP
    
    do while (residual.gt.tol)
       residual = 0.0_WP

       ! Store G values at t=t_iter+1/2
       gthalf = soln(3,1,:,:)
       
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
                
             else
                ! ! *** Inner cells, x & y-dir ***
                ! soln(3,1,j,k) = gthalf(j,k) + &
                !      dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                !      1.0_WP - 0.5_WP/h*sqrt( (gthalf(j+1,k)-gthalf(j-1,k))**2 &
                !      + (gthalf(j,k+1)-gthalf(j,k-1))**2 ) )

                ! *** Inner cells, x-dir only ***
                soln(3,1,j,k) = gthalf(j,k) + &
                     dt*gthalf(j,k)/sqrt(gthalf(j,k)**2 + epsilon**2)*(  &
                     1.0_WP - (gthalf(j+1,k)-gthalf(j-1,k))/(2.0_WP*h) )

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

       if (t_iter.le.2) then ! First and second timesteps where solution is 0
          residual = sum(abs(soln(3,1,:,:) - gthalf))
       else
          residual = sum(abs(soln(3,1,:,:) - gthalf))/sum(abs(gthalf))
       end if
       
       subiter = subiter + 1
       print *, "t_iter=", t_iter, "Subiteration=", subiter, &
            "Residual=", residual
    end do
    
    ! Store converged G values at t=t_iter+1
    soln(3,2,:,:) = soln(3,1,:,:)
    
  end subroutine gEqn_step
  

end program GEQN_SOLVER
