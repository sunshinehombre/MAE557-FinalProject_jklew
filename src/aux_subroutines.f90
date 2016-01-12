module AUX_SUBROUTINES
  use CONSTS
  implicit none

  ! aux_subroutines.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! module AUX_SUBROUTINES provides subroutines for initializing a solution matrix
  ! as well as writing data to files

contains

  ! *** Initialize grid and solution matrix *****************************************
  subroutine soln_init(var,nxy,U,h,soln)
    use CONSTS
    implicit none

    integer, intent(in) :: var
    integer, intent(in) :: nxy
    real(WP), intent(in) :: U
    real(WP), intent(in) :: h
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln

    integer :: i,j,k

    do i=1,var 
       select case(i)
       case(1) ! zeta (vorticity)
          do j=1,nxy+1
             do k=1,nxy+1
                soln(i,1,j,k) = 0.0_WP
             end do
          end do
       case(2) ! psi (stream function)
          do j=1,nxy+1
             do k=1,nxy+1
                if (j.eq.1 .and. k.ge.2 .and. k.le.nxy) then ! Uniform vel. profile
                   soln(i,1,j,k) = U*(k-1)*h
                else
                   soln(i,1,j,k) = 0.0_WP
                end if
             end do
          end do
       case(3) ! G (iso-surface)
          do j=1,nxy+1
             do k=1,nxy+1
                if (j.eq.int((nxy+1)/2.0_WP)) then ! Halfway along x-dir in channel
                   soln(i,1,j,k) = 1.0_WP
                else
                   soln(i,1,j,k) = 0.0_WP
                end if
             end do
          end do
       case default
          print *, 'Attempted to initialize unknown variable.'
       end select
    end do
    
  end subroutine soln_init


  ! *** Write solution to file ******************************************************
  subroutine write_solution(ma,h,var,dt,nt,nxy,soln)
    use CONSTS
    implicit none

    real(WP), intent(in) :: ma ! Mach number
    real(WP), intent(in) :: h ! Spatial step size

    ! Variable index for solution matrix
    !      var = 1: u
    !      var = 2: v
    !      var = 3: G
    integer, intent(in) :: var

    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: nt ! Total time steps
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(in) :: soln

    integer :: i,j,k ! Loop iterators

    ! Strings used for output filenames
    character(len=8), parameter :: wd = '../data/'
    character(len=54) :: outfile ! Output filenames
    
    character(len=3), parameter :: tstep = 'dt_'
    character(len=5) :: mach
    character(len=2), parameter :: time = '_t', sstep = '_h'
    character(len=4), parameter :: varnum = '_var', endname = '.txt'

    ! Determine the mach number under consideration 
    select case(int(ma*100.0_WP))
       case(5)
          mach = 'ma005'
       case(10)
          mach = 'ma010'
       case(30)
          mach = 'ma030'
       case default
          print *, 'Unknown Mach number.'
       end select
    
       ! Write data to file
       do i=1,var
          write(outfile,"(a8,a5,a3,i10.10,a2,i9.9,a2,i6.6,a4,i1,a4)") wd,mach, &
               tstep,int(dt*1E9_WP),time,nt+1,sstep,int(h*1E5_WP),varnum,i,endname

          open(unit=1,file=outfile,status="replace",action="readwrite")
          ! Write u, v, and G
          do k=1,nxy+1
             write(unit=1,fmt="(1000000(e11.4,1x))") &
                  (soln(i,1,j,k), j=1,nxy+1)
          end do
          close(unit=1)

       end do

  end subroutine write_solution


end module AUX_SUBROUTINES


    
