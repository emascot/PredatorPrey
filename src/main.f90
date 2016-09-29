! Common constants
module constants
  implicit none
  ! Double precision
  integer, parameter :: dp = kind(0.d0)
  ! Total number of test velocities
  integer, parameter :: nvel = 100
  ! Total number of time steps
  integer, parameter :: nsteps = 1000
  ! Total amount of time
  real(dp), parameter :: totaltime = 30.d0
  ! Time step
  real(dp), parameter :: tstep = totaltime / nsteps
  ! Velocity of prey
  real(dp), parameter :: prey_vel = 1.d0
  ! Start Distance
  real(dp), parameter :: start_dist = 10.d0 * prey_vel
  ! Minimum velocity of predator
  real(dp), parameter :: min_vel = 1.d0
  ! Maximum velocity of predator
  real(dp), parameter :: max_vel = 2.d0
  ! Distance to catch prey
  real(dp), parameter :: range = 0.10
end module constants

! Record trajectories for different velocities
program main
  use mpi
  use constants
  use parallel_tasks
  implicit none
  ! Command line argument
  character(len=32) :: cmd
  ! Granularity
  integer :: nvec = 1
  ! Start, wall, and cpu time
  real(dp) :: start, elapsed, total
  ! Find minimum velocity
  real(dp) :: find_min_vel
  ! Position as time vs velocity
  complex(dp) :: position(nsteps+1,nvel)
  ! Time and velocity iterator
  integer :: it, ivel
  ! Velocity array
  real(dp) :: velocities(nvel)
  ! Calculate position at each iteration
  external :: chase
  ! MPI ID and error status
  integer :: rank, size, ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  call mpi_comm_size(mpi_comm_world,size,ierr)
  start = mpi_wtime()

  ! Read command line argument
  call get_command_argument(1,cmd)
  ! Convert to integer
  if (len_trim(cmd) > 0) read(cmd, *) nvec

  ! Find minimum velocity
  if (rank.eq.0) write(*,*) 'Minimum Velocity:', find_min_vel()

  ! Array of test velocities
  do ivel = 0,nvel-1
    ! Velocities ranging from min to max
    velocities(ivel+1) = min_vel + (max_vel-min_vel) * ivel / real(nvel-1, dp)
  enddo

#ifdef VERBOSE
  ! Show progress bar
  prog_bar = .true.
#endif

  call task_manager(nvel,1,nsteps+1,nvec,velocities,chase,position,ierr)

#ifdef EXPORT
  ! Export to file
  open(10, file='t.txt', status='replace')
  open(11, file='x.txt', status='replace')
  open(12, file='y.txt', status='replace')
  open(13, file='data.txt', status='replace')
  write(13,'(3X,A,18X,A,22X,A,25X,A)') 'Velocity', 'Time', 'x', 'y'
  do it=1,nsteps+1
    write(10,*) (it-1) * tstep
    write(11,*) real(position(it,:))
    write(12,*) aimag(position(it,:))
    do ivel=1,nvel
      write(13,*) velocities(ivel), (it-1) * tstep, real(position(it,ivel)), aimag(position(it,ivel))
    enddo
  enddo
  close(10)
  close(11)
  close(12)
  close(13)
#endif

  elapsed = MPI_Wtime() - start
  call MPI_Reduce(elapsed, total, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (rank.eq.0) then
    write(*,*) size, nvec, total, elapsed
  endif
  call mpi_finalize(ierr)
end program main

subroutine chase(ndim,velocities,nfun,positions,nvec)
  use constants
  implicit none
  integer, intent(in) :: ndim, nfun, nvec
  real(dp), intent(in) :: velocities(ndim,nvec)
  complex(dp), intent(out) :: positions(nfun,nvec)
  ! Time and velocity iterators
  integer :: it, ivel
  ! Time
  real(dp) :: time
  ! Target position
  complex(dp) :: prey_position
  ! Difference in positions
  complex(dp) :: diff
  
  positions(1,:) = cmplx(0.d0, start_dist, dp)
  do ivel = 1,ndim
    do it = 1,nsteps
      ! Find prey
      time = it * tstep
      prey_position = cmplx(prey_vel * time, 0.d0, dp)
      ! Aim towards prey
      diff = prey_position - positions(it,ivel)
      ! Check if within reach
      if (abs(diff) < velocities(1,ivel) * tstep + range) then
        positions(it:nsteps+1,ivel) = prey_position
        exit
      else
        ! Move towards prey
        positions(it+1,ivel) = positions(it,ivel) + velocities(1,ivel) * tstep * diff / abs(diff)
      endif
    enddo
  enddo
end subroutine chase

! Converge onto minimum velocity
real(dp) function find_min_vel()
  use constants
  implicit none
  ! Number of tests
  integer, parameter :: ntest = 1000
  ! Test iterator
  integer :: itest
  ! Current minimum velocity
  real(dp) :: min = min_vel
  ! Current maximum velocity
  real(dp) :: max = max_vel
  ! Current velocity
  real(dp) :: velocity = max_vel
  ! Predator catches prey test function
  logical :: will_catch

  ! Repeat ntest times
  do itest=1,ntest
    ! Bisect range
    velocity = (min + max) * 0.5d0
    if (will_catch(velocity)) then
      ! Lower maximum
      max = velocity
    else
      ! Increase minimum
      min = velocity
    endif
    if (abs(min - max) < 5.0e-16) exit
  enddo
  find_min_vel = velocity
end function find_min_vel

! Return true if fast enough to catch prey
logical function will_catch(velocity)
  use constants
  implicit none
  real(dp), intent(in) :: velocity
  ! Time iterator
  integer :: it
  ! Time
  real(dp) :: time
  ! Current position
  complex(dp) :: position
  ! Target position
  complex(dp) :: prey_position
  ! Difference in positions
  complex(dp) :: diff

  will_catch = .false.
  position = cmplx(0.d0, 10.d0, dp)
  do it = 1,nsteps
    ! Find prey
    time = it * tstep
    prey_position = cmplx(prey_vel * time, 0.d0, dp)
    ! Aim towards prey
    diff = prey_position - position
    ! Check if within reach
    if (abs(diff) < velocity * tstep + range) then
      will_catch = .true.
      return
    else
      ! Move towards prey
      position = position + velocity * tstep * diff / abs(diff)
    endif
  enddo
end function will_catch
