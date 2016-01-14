module AUX_SUBROUTINES
  use CONSTS
  implicit none

  ! aux_subroutines.f90 written by Jeffry Lew, jklew@princeton.edu
  !
  ! module AUX_SUBROUTINES provides subroutines for initializing a solution matrix
  ! as well as writing data to files

contains

  ! *** Initialize grid and solution matrix *****************************************
  subroutine soln_init(h,nxy,U,var,soln)
    use CONSTS
    implicit none

    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), intent(in) :: U ! Bulk velocity [m/s]

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(inout) :: soln ! Solution matrix

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
                ! if (j.eq.1 .or. j.eq.2 .or. j.eq.nxy .or. j.eq.nxy+1) then
                ! if (j.eq.1 .and. (k.ge.2 .and. k.le.nxy)) then
                if (j.eq.1) then
                   soln(i,2,j,k) = U*(k-1)*h ! Uniform vel.
                   
                   ! soln(i,1,j,k) = -2.0_WP*U*((k-1)*h)**3/(nxy*h)**2 + &
                        ! 1.5_WP*U*(k-1)*h ! Parabolic vel. (fully developed)
                   
                else
                   soln(i,2,j,k) = 0.0_WP
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


  ! *** Convert stream function to velocities u & v ***
  subroutine convert_streamfunc(h,nxy,var,soln,var2,uv)
    use CONSTS
    implicit none

    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(in) :: soln ! Solution matrix

    ! Variable index for velocity matrix
    !      var2 = 1: u
    !      var2 = 2: v
    integer, intent(in) :: var2
    real(WP), dimension(var2,nxy+1,nxy+1), intent(inout) :: uv ! Values of u & v

    integer :: i,j,k ! Loop iterators
    

    !$OMP PARALLEL

    ! Store u & v boundary conditions
    !$OMP DO
    do j=1,nxy+1
       k = 1 ! Bottom wall
       uv(1,j,k) = 0.0_WP ! u
       uv(2,j,k) = 0.0_WP ! v

       k = nxy+1 ! Top wall
       uv(1,j,k) = 0.0_WP ! u
       uv(2,j,k) = 0.0_WP ! v
    end do
    !$OMP END DO

    ! Store u & v at inner grid points
    !$OMP DO
    do j=1,nxy+1
       do k=2,nxy
          if (j.eq.1) then ! Left boundary
             uv(1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             ! uv(2,j,k) = -(-soln(2,2,j+2,k) + 4.0_WP*soln(2,2,j+1,k) - &
             ! 3.0_WP*soln(2,2,j,k))/(2.0_WP*h) ! v

             uv(2,j,k) = -(soln(2,2,j+1,k) - soln(2,2,nxy+1,k))/(2.0_WP*h) ! v

          elseif (j.eq.nxy+1) then ! Right boundary
             uv(1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             uv(2,j,k) = -(3.0_WP*soln(2,2,j,k) - 4.0_WP*soln(2,2,j-1,k) + &
                  soln(2,2,j-2,k))/(2.0_WP*h) ! v

             ! uv(2,j,k) = -(soln(2,2,1,k) - soln(2,2,j-1,k))/(2.0_WP*h) ! v

          else
             uv(1,j,k) = (soln(2,2,j,k+1) - soln(2,2,j,k-1))/(2.0_WP*h) ! u
             uv(2,j,k) = -(soln(2,2,j+1,k) - soln(2,2,j-1,k))/(2.0_WP*h) ! v

          end if
       end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine convert_streamfunc


  ! *** Write solution to file ******************************************************
  subroutine write_solution(ma,h,nxy,dt,nt,var,soln)
    use CONSTS
    implicit none

    real(WP), intent(in) :: ma ! Mach number
    real(WP), intent(in) :: h ! Spatial step size
    integer, intent(in) :: nxy ! Total spatial steps nx=ny
    real(WP), intent(in) :: dt ! Time step size
    integer, intent(in) :: nt ! Total time steps

    ! Variable index for solution matrix
    !      var = 1: zeta (vorticity)
    !      var = 2: psi (stream function)
    !      var = 3: G (iso-surface)
    integer, intent(in) :: var
    real(WP), dimension(var,2,nxy+1,nxy+1), intent(in) :: soln
    
    ! Variable index for velocity matrix
    !      var2 = 1: u
    !      var2 = 2: v
    integer, parameter :: var2 = 2
    real(WP), dimension(var2,nxy+1,nxy+1) :: uv

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

    ! Convert final stream function solution to velocities
    call convert_streamfunc(h,nxy,var,soln,var2,uv)

    ! Write data to file
    do i=1,var
       write(outfile,"(a8,a5,a3,i10.10,a2,i9.9,a2,i6.6,a4,i1,a4)") wd,mach,tstep, &
            int(dt*1E9_WP),time,nt+1,sstep,int(h*1E5_WP),varnum,i,endname

       open(unit=1,file=outfile,status="replace",action="readwrite")

       ! Write u, v, and G values
       select case(i)
       case(1:2) ! Write u and v values
          do k=1,nxy+1
             write(unit=1,fmt="(1000000(e11.4,1x))") (uv(i,j,k), j=1,nxy+1)
          end do
       case(3) ! Write G values
          do k=1,nxy+1
             write(unit=1,fmt="(1000000(e11.4,1x))") (soln(i,2,j,k), j=1,nxy+1)
          end do
       case default
          print *, 'Attempted to write unknown variable.'
       end select

       close(unit=1)
    end do

  end subroutine write_solution

  
end module AUX_SUBROUTINES


    
