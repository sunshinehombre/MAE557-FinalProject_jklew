program GEQN_SOLVER
  use CONSTS
  use AUX_SUBROUTINES
  implicit none

  ! gEqnSolver.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! GEQN_SOLVER solves the G-equation for the corrugated flamelets regime in a 2D
  ! incompressible channel flow.
  

  ! Set number of OpenMP threads
  CALL OMP_SET_NUM_THREADS(NTHREADS)

  ! Spatial step size h=dx=dy
  real(WP), parameter :: h1 = 0.0625_WP ! [m], Grid is 16x16
  real(WP), parameter :: h2 = 0.03125_WP ! [m], Grid is 32x32
  real(WP), parameter :: h3 = 0.015625_WP ! [m], Grid is 64x64

  ! Length of domain (square) side
  real(WP), parameter :: lside1 = 1.0_WP ! [m]

  ! Time step size & total simulation time
  real(WP), parameter :: dt1 = 1E-6_WP ! [s]
  real(WP), parameter :: tott1 = 5.0_WP ! [s]

  ! Flow Mach number
  real(WP), parameter :: ma1 = 0.05_WP
  real(WP), parameter :: ma2 = 0.1_WP
  real(WP), parameter :: ma3 = 0.3_WP

  ! For timing purposes
  integer :: t1,t2,clock_rate,clock_max

  ! Start timer
  call system_clock (t1,clock_rate,clock_max)

  ! Solve for incompressible flow field with iso-surface from G-Equation
  call solveGEqn(h1,lside1,dt1,tott1,ma1)

  ! Stop timer
  call system_clock (t2,clock_rate,clock_max)
  print *, 'Elapsed real time =', real(t2-t1)/real(clock_rate), 's'
  
contains

  ! *** Incompressible Navier-Stokes solver with incorporated G-Equation ************
  
  subroutine solveGEqn(h,lside,dt,tott,ma)
    use CONSTS
    use AUX_SUBROUTINES
    implicit none


  end subroutine solveGEqn

end program GEQN_SOLVER
