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
  ! Velocity iterator
  integer :: ivel
  ! Velocity array
  real(dp) :: velocities(nvel)
  ! Catch time for each velocity
  complex(dp) :: times(nvel)

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

  call task_manager(nvel,1,1,nvec,velocities,chase_time,times,ierr)

#ifdef EXPORT
  open(14, file='tf.txt', status='replace')
  do ivel=1,nvel
    write(14,*) velocities(ivel), real(times(ivel))
  enddo
  close(14)
#endif

  elapsed = MPI_Wtime() - start
  call MPI_Reduce(elapsed, total, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (rank.eq.0) then
    write(*,*) size, nvec, total, elapsed
  endif
  call mpi_finalize(ierr)
end program main