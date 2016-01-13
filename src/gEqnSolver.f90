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
  real(WP), parameter :: tfin1 = 1E-6_WP !0.01_WP ! [s]

  ! Mach number based on bulk velocity
  real(WP), parameter :: ma1 = 0.05_WP
  real(WP), parameter :: ma2 = 0.1_WP
  real(WP), parameter :: ma3 = 0.3_WP

  ! For timing purposes
  integer :: t1,t2,clock_rate,clock_max

  ! Set number of OpenMP threads
  CALL OMP_SET_NUM_THREADS(NTHREADS)

  ! Start timer
  call system_clock (t1,clock_rate,clock_max)

  ! Solve for incompressible flow field with iso-surface from G-Equation
  call solveGEqn(h3,lside1,dt1,tfin1,ma1)

  ! Stop timer
  call system_clock (t2,clock_rate,clock_max)
  print *, 'Elapsed real time =', real(t2-t1)/real(clock_rate), 's'

  
contains
  

  ! *** Incompressible Navier-Stokes solver with incorporated G-Equation ************
  
  subroutine solveGEqn(h,lside,dt,tfin,ma)
    use CONSTS
    use AUX_SUBROUTINES
    implicit none

    real(WP), intent(in) :: h ! Spatial step size
    real(WP), intent(in) :: lside ! Length of square side
    real(WP), intent(in) :: dt ! Time step size
    real(WP), intent(in) :: tfin ! Total time of simulation
    real(WP), intent(in) :: ma ! Flow mach number

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
    call soln_init(var,nxy,U,h,soln)

    ! Iterate until final time
    do t_iter = 1,nt
       call vs_step(dt,t_iter,U,nu,var,h,nxy,soln)
    end do

    ! Convert final stream function solution to velocities (var=1: u, var=2: v)
    call convert_streamfunc(dt,nt,U,var,h,nxy,soln)

    ! Write u, v, and G to file for plotting purposes
    call write_solution(ma,h,var,dt,nt,nxy,soln)

    ! For testing purposes
    print *, 'Bulk velocity: ', U, ' m/s'

  end subroutine solveGEqn


  ! *** Vorticity-stream function approach below ************************************
  subroutine vs_step(dt,t_iter,U,nu,var,h,nxy,soln)
    use CONSTS
    implicit none

    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: t_iter ! Current timestep index
    real(WP), intent(in) :: U ! Bulk velocity
    real(WP), intent(in) :: nu ! Kinematic viscosity [m^2/s]

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln ! Solution matrix

    real(WP), dimension(nxy+1,nxy+1) :: psip ! Psi at subiteration subitercount
    integer :: i,j,k,subitercount ! Loop iterators
    real(WP), parameter :: tol = 0.001_WP ! Tolerance
    real(WP) :: residual ! Residual
    real(WP), parameter :: beta = 1.5_WP ! Coefficient in Successive Over Relaxation

    ! *** Substep 1: Iterate for new psi values ***
    ! Initialize subiteration counter and residual
    subitercount = 0
    residual = 1.0_WP
    do while (residual.gt.tol)
       residual = 0.0_WP
       
       ! Initialize psi values at subiteration subitercount to psi values at t=t_iter
       psip = soln(2,1,:,:)
       
       ! Solve for psi at subiteration subitercount+1 with S.O.R.
       do j=1,nxy+1
          do k=2,nxy
             if (j.eq.1 .and. subitercount.gt.0) then
                soln(2,1,j,k) = 0.25_WP*beta*( soln(2,1,j+1,k) + soln(2,1,nxy+1,k)+ &
                     soln(2,1,j,k+1) + soln(2,1,j,k-1) + h**2*soln(1,1,j,k) ) + &
                     (1.0_WP-beta)*soln(2,1,j,k)
                
             elseif (j.eq.nxy+1) then
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
       
       subitercount = subitercount + 1
       print *, "t_iter=", t_iter, "Subiteration=", subitercount, &
            "Residual=", residual
    end do

    !$OMP PARALLEL
    
    ! Store converged psi (stream function) solution at t=t_iter+1
    !$OMP DO
    do j=1,nxy+1
       do k=1,nxy+1
          soln(2,2,j,k) = soln(2,1,j,k) 
       end do
    end do
    !$OMP END DO

    ! *** Substep 2: Determine zeta (vorticity) B.C.s using inner psi & zeta ***

    ! ! Left inflow periodic boundary conditions for top and bottom wall
    ! soln(1,2,1,1) = soln(1,2,nxy+1,1) ! Bottom wall
    ! soln(1,2,1,nxy+1) = soln(1,2,nxy+1,nxy+1) ! Top wall
    
    !$OMP DO
    do j=1,nxy+1
       k = 1 ! Bottom wall
       soln(1,2,j,k) = -2.0_WP/h**2*soln(2,2,j,k+1)
       
       k = nxy+1 ! Top wall
       soln(1,2,j,k) = -2.0_WP/h**2*soln(2,2,j,k-1)
    end do
    !$OMP END DO   

    !$OMP DO
    do k=2,nxy
       j = 1 ! Left inflow
       ! soln(1,2,j,k) = -(soln(2,2,nxy+1,k)-2.0_WP*soln(2,2,nxy,k) + &
            ! soln(2,2,nxy-1,k))/h**2 -(soln(2,2,nxy+1,k+1)-2.0_WP*soln(2,2,nxy+1,k)&
            ! + soln(2,2,nxy+1,k-1))/h**2

       soln(1,2,j,k) = -(soln(2,2,j+1,k)-2.0_WP*soln(2,2,j,k)+soln(2,2,nxy+1,k))/ &
            h**2 - (soln(2,2,j,k+1)-2.0_WP*soln(2,2,j,k)+soln(2,2,j,k-1))/h**2
       
       j = nxy+1 ! Right outflow
       ! soln(1,2,j,k) = -(soln(2,2,3,k)-2.0_WP*soln(2,2,2,k)+soln(2,2,1,k))/h**2 - &
            ! (soln(2,2,1,k+1)-2.0_WP*soln(2,2,1,k)+soln(2,2,1,k-1))/h**2

       soln(1,2,j,k) = -(soln(2,2,1,k)-2.0_WP*soln(2,2,j,k)+soln(2,2,j-1,k))/h**2 - &
            (soln(2,2,j,k+1)-2.0_WP*soln(2,2,j,k)+soln(2,2,j,k-1))/h**2
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
  
  
  ! *** Convert stream function to velocities u & v ***
  subroutine convert_streamfunc(dt,nt,U,var,h,nxy,soln)
    use CONSTS
    implicit none

    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: nt ! Total time steps
    real(WP), intent(in) :: U ! Bulk velocity

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln ! Solution matrix

    integer :: i,j,k

    ! Use space in soln containing old solution at t=dt*(nt-1) to store u & v fields
    ! This corresponds to u = soln(1,1,:,:) and v = soln(2,1,:,:)
    !$OMP PARALLEL

    ! Store u & v boundary conditions
    !$OMP DO
    do j=1,nxy+1
       k = 1 ! Bottom wall
       soln(1,1,j,k) = 0.0_WP ! u
       soln(2,1,j,k) = 0.0_WP ! v
       
       k = nxy+1 ! Top wall
       soln(1,1,j,k) = 0.0_WP ! u
       soln(2,1,j,k) = 0.0_WP ! v
    end do
    !$OMP END DO
    
    ! Store u & v at inner grid points
    !$OMP DO
    do j=1,nxy+1
       do k=2,nxy
          if (j.eq.1) then ! Left boundary
             soln(1,1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             ! soln(2,1,j,k) = -(-soln(2,2,j+2,k) + 4.0_WP*soln(2,2,j+1,k) - &
                  ! 3.0_WP*soln(2,2,j,k))/(2.0_WP*h) ! v

             soln(2,1,j,k) = -(soln(2,2,j+1,k) - soln(2,2,nxy+1,k))/(2.0_WP*h) ! v
             
          elseif (j.eq.nxy+1) then ! Right boundary
             soln(1,1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             ! soln(2,1,j,k) = -(3.0_WP*soln(2,2,j,k) - 4.0_WP*soln(2,2,j-1,k) + &
             ! soln(2,2,j-2,k))/(2.0_WP*h) ! v

             soln(2,1,j,k) = -(soln(2,2,1,k) - soln(2,2,j-1,k))/(2.0_WP*h) ! v
             
          elseif (j.gt.3 .and. j.lt.nxy-1) then
             soln(1,1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             soln(2,1,j,k) = -(soln(2,2,j+1,k) - soln(2,2,j-1,k))/(2.0_WP*h) ! v

          else
             
          end if
       end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL
    
  end subroutine convert_streamfunc
  

end program GEQN_SOLVER
