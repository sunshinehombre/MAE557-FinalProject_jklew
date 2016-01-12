module CONSTS
  use omp_lib
  implicit none

  ! consts.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! module CONSTS provides constants that are used in several files

  ! Define working precision, pi, and number of OpenMP threads to use
  integer, parameter :: WP = kind(1.0d0)
  real(WP), parameter :: PI = 3.14159265358979323846_WP
  integer, parameter :: NTHREADS = 8

  ! Burning velocity of unstretched planar flame [m/s]
  real(WP), parameter :: SLNOT = 0.1_WP

  ! Markstein length [m]
  real(WP), parameter :: ML = 0.0001_WP

end module CONSTS
