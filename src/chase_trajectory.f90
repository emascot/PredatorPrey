! Record trajectories for different velocities
program main
  use mpi
  use chase
  use parallel_tasks
  implicit none
  ! MPI ID, size, and error status
  integer :: rank, size, ierr
  ! Start, wall, and cpu time
  real(dp) :: start, elapsed, total
  ! Command line argument
  character(len=32) :: cmd
  ! Granularity
  integer :: nvec = 1
  ! Time and velocity iterator
  integer :: it, ivel
  ! Velocity array
  real(dp) :: velocities(nvel)
  ! Position at each time for each velocity
  complex(dp) :: positions(nsteps+1,nvel)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  call mpi_comm_size(mpi_comm_world,size,ierr)
  start = mpi_wtime()

  ! Read command line argument
  call get_command_argument(1,cmd)
  ! Convert to integer
  if (len_trim(cmd) > 0) read(cmd, *) nvec

  ! Array of test velocities
  do ivel = 1,nvel
    ! Velocities ranging from min to max
    velocities(ivel) = min_vel + (max_vel-min_vel) * (ivel-1) / real(nvel-1, dp)
  enddo

#ifdef VERBOSE
  ! Find minimum velocity
  if (rank.eq.0) write(*,*) 'Minimum Velocity:', find_min_vel()
  ! Show progress bar
  prog_bar = .true.
#endif

  call task_manager(nvel,1,nsteps+1,nvec,velocities,chase_trajectory,positions,ierr)

#ifdef EXPORT
  ! Export to file
  open(10, file='t.txt', status='replace')
  open(11, file='x.txt', status='replace')
  open(12, file='y.txt', status='replace')
  do it=1,nsteps+1
    write(10,*) (it-1) * tstep
    write(11,*) real(positions(it,:))
    write(12,*) aimag(positions(it,:))
  enddo
  close(10)
  close(11)
  close(12)

  open(13, file='data.txt', status='replace')
  write(13,'(3X,A,18X,A,22X,A,25X,A)') 'Velocity', 'Time', 'x', 'y'
  do ivel=1,nvel
    do it=1,nsteps+1
      write(13,*) velocities(ivel), (it-1) * tstep, &
        real(positions(it,ivel)), aimag(positions(it,ivel))
    enddo
  enddo
  close(13)
#endif

  elapsed = MPI_Wtime() - start
  call MPI_Reduce(elapsed, total, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (rank.eq.0) then
    write(*,*) size, nvec, total, elapsed
  endif
  call mpi_finalize(ierr)
end program main