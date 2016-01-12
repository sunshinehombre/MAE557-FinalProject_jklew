module AUX_SUBROUTINES
  use CONSTS
  implicit none

  ! aux_subroutines.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! module AUX_SUBROUTINES provides subroutines for initializing a solution matrix
  ! as well as writing data to files

contains

  ! *** Initialize grid and solution matrix *****************************************
  subroutine soln_init(var,nxy,soln)
    use CONSTS
    implicit none

    integer, intent(in) :: var
    integer, intent(in) :: nxy
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln

    integer :: i,j,k


  end subroutine soln_init


  ! *** Write solution to file ******************************************************
  subroutine write_solution()
    use CONSTS
    implicit none


  end subroutine write_solution


end module AUX_SUBROUTINES


    
