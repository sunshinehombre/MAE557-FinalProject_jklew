module CONSTS
  use omp_lib
  implicit none

  ! consts.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! module CONSTS provides constants that are used in several files

  ! Define working precision, pi, and number of OpenMP threads to use
  integer, parameter :: WP = kind(1.0d0)
  real(WP), parameter :: PI = 3.14159265358979323846_WP
  integer, parameter :: NTHREADS = 1

  ! Burning velocity of unstretched planar flame [m/s]
  real(WP), parameter :: SLNOT = 25.0_WP !5.0_WP

  ! Markstein length [m]
  real(WP), parameter :: ML = 0.0001_WP !0.1_WP !0.01_WP !-0.0001_WP !-0.01_WP 

  ! Reynolds number constant
  real(WP), parameter :: RE = 100.0_WP

  ! Specific heat constants & specific gas constant
  real(WP), parameter :: GAMMA = 1.4_WP
  real(WP), parameter :: RGAS = 287.0_WP ! J/(kg K)
  real(WP), parameter :: CP = GAMMA/(GAMMA-1.0_WP)*RGAS ! J/(kg K)
  
  ! Reference pressure, temperature, density
  real(WP), parameter :: PREF = 100000.0_WP ! Pa
  real(WP), parameter :: TREF = 300.0_WP ! K
  real(WP), parameter :: DREF = PREF/(RGAS*TREF) ! kg/m^3

  ! Reference speed of sound
  real(WP), parameter :: SOUND = sqrt(GAMMA*RGAS*TREF)

end module CONSTS
