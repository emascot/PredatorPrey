#define FUNTYPE real(dp)
#define RESTYPE complex(dp)
#define MPIFUNTYPE MPI_REAL8
#define MPIRESTYPE MPI_DOUBLE_COMPLEX

!==============================================================
!> Distribute tasks to processes.
!> Tasks are an array of FUNTYPE sent to a subroutine.
!> |option    |description                                |
!> |----------|-------------------------------------------|
!> |stream    |Write output as byte stream instead of text|
!> |sequential|Write after each task is completed         |
!> |prog_bar  |Show a progress bar with estimated time    |
!> |bcast     |Send results to all processes              |
!> \date September 2016
!> \author Eric Mascot
!==============================================================
!> \par Example code:
!> \code
!> program example
!>   use mpi
!>   use parallel_tasks
!>   implicit none
!>   integer, parameter :: Ntasks=10, Nfun=2, Nres=2, Nvec=1, dp=kind(0.d0)
!>   real(dp) :: tasks(Nfun,Ntasks), results(Nres,Ntasks)
!>   integer :: ierr
!> 
!>   ! Initialize MPI
!>   call MPI_INIT(ierr)
!> 
!>   ! Make task list
!>   call random_number(tasks)
!>
!>   ! Turn on progress bar
!>   prog_bar = .true.
!>
!>   ! Assign tasks and write results to "example.dat"
!>   call task_manager(Ntasks,Nfun,Nres,Nvec,tasks,sum_prod,results,ierr,"example.dat")
!> 
!>   ! Finalize MPI
!>   call MPI_FINALIZE(ierr)
!> 
!> contains
!>
!>   ! Define task
!>   subroutine sum_prod(Nfun,fun,Nres,res,Nvec)
!>     integer, intent(in) :: Nfun, Nres, Nvec
!>     real(dp), intent(in) :: fun(Nfun,Nvec)
!>     real(dp), intent(out) :: res(Nres,Nvec)
!>     integer :: iv, if
!>     do iv=1,Nvec
!>       do if=1,Nfun
!>         res(1,iv) = fun(if,iv) + fun(if,iv)
!>         res(2,iv) = fun(if,iv) * fun(if,iv)
!>       enddo
!>     enddo
!>   end subroutine sum_prod
!> end program example
!> \endcode
!==============================================================

module parallel_tasks
  use mpi
  implicit none
  integer, parameter :: dp=kind(0.d0)
  ! true  - Write as byte stream (fast)
  ! false - Write in human readable format
  logical :: stream=.false.
  ! true  - Sequentially write after each task
  ! false - Do not do sequential write
  logical :: sequential=.false.
  ! true  - Show progress bar
  ! false - Do not show progress
  logical :: prog_bar=.false.
  ! true  - Broadcast results
  ! false - Only root has results
  logical :: bcast=.false.
  save
  private
  public :: task_manager, task_static, task_dynamic,&
            stream, sequential, prog_bar, bcast
contains


!==============================================================
!  task_manager
!> Assigns a new task to a process when process finishes task.
!> Takes an array of inputs (tasks) and passes it through
!> an externally declared subroutine (func).
!
!> \param[in]  Ntasks  Number of tasks \n
!>                     Size of tasks array
!> \param[in]  Nfun    Number of function arguments \n
!>                     Size of arguments array for subroutine "func"
!> \param[in]  Nres    Number of results \n
!>                     Size of results array from subroutine "func"
!> \param[in]  Nvec    Size of vector \n
!>                     Number of tasks to do at a time
!> \param[in]  tasks   Array of function inputs \n
!>                     Arguments to use for subroutine "func"
!> \param[in]  func    externally declared subroutine to calculate tasks \n
!>                     Must have parameters (Nfun,fun,Nres,res)
!>                     - Inputs:
!>                       - Nfun - Size of fun array
!>                       - fun - Array of inputs
!>                       - Nres - Size of res array
!>                     - Output:
!>                       - res - Array of results
!> \param[out] results Array of results
!>                     Output of subroutine "func"
!> \param[out] ierr    MPI error return value
!>   - MPI_SUCCESS
!>   - No error; MPI routine completed successfully.
!>   - MPI_ERR_COMM
!>   - Invalid communicator. A common error is to use a null communicator in a call
!>     (not even allowed in MPI_Comm_rank).
!>   - MPI_ERR_COUNT
!>   - Invalid count argument. Count arguments must be non-negative; a count of zero is often valid.
!>   - MPI_ERR_TYPE
!>   - Invalid datatype argument. May be an uncommitted MPI_Datatype (see MPI_Type_commit).
!>   - MPI_ERR_TAG
!>   - Invalid tag argument. Tags must be non-negative; tags in a receive
!>     (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_TAG.
!>     The largest tag value is available through the the attribute MPI_TAG_UB.
!>   - MPI_ERR_RANK
!>   - Invalid source or destination rank. Ranks must be between zero and the size of the communicator minus one;
!>     ranks in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_SOURCE.
!> \param[in] fname    Filename of output \n
!>                     Skips if fname is not provided \n
!>                     File has first columns as tasks and last columns as results
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine task_manager(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
  implicit none
  integer, intent(in) :: Ntasks, Nfun, Nres, Nvec
  FUNTYPE, intent(in) :: tasks(Nfun,Ntasks)
  external :: func
  RESTYPE, intent(out) :: results(Nres,Ntasks)
  integer, intent(out) :: ierr
  character(len=*), optional, intent(in) :: fname
  integer :: size

  ! Get number of processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

  ! For small number of processes, use static load balancing
  if (size.lt.3) then
    call task_static(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
  else
#ifdef STATIC
    call task_static(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
#else
    call task_dynamic(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
#endif
  end if
end subroutine task_manager


!==============================================================
!  task_dynamic
!> Assigns a new task to a process when process finishes task.
!> Takes an array of inputs (tasks) and passes it through
!> an externally declared subroutine (func).
!
!> \param[in]  Ntasks  Number of tasks \n
!>                     Size of tasks array
!> \param[in]  Nfun    Number of function Arguments \n
!>                     Size of arguments array for subroutine "func"
!> \param[in]  Nres    Number of results \n
!>                     Size of results array from subroutine "func"
!> \param[in]  Nvec    Size of vector \n
!>                     Number of tasks to do at a time
!> \param[in]  tasks   Array of function inputs \n
!>                     Arguments to use for subroutine "func"
!> \param[in]  func    externally declared subroutine to calculate tasks \n
!>                     Must have parameters (Nfun,fun,Nres,res)
!>                     - Inputs:
!>                       - Nfun - Size of fun array
!>                       - fun - Array of inputs
!>                       - Nres - Size of res array
!>                     - Output:
!>                       - res - Array of results
!> \param[out] results Array of results
!>                     Output of subroutine "func"
!> \param[out] ierr    MPI error return value
!>   - MPI_SUCCESS
!>   - No error; MPI routine completed successfully.
!>   - MPI_ERR_COMM
!>   - Invalid communicator. A common error is to use a null communicator in a call
!>     (not even allowed in MPI_Comm_rank).
!>   - MPI_ERR_COUNT
!>   - Invalid count argument. Count arguments must be non-negative; a count of zero is often valid.
!>   - MPI_ERR_TYPE
!>   - Invalid datatype argument. May be an uncommitted MPI_Datatype (see MPI_Type_commit).
!>   - MPI_ERR_TAG
!>   - Invalid tag argument. Tags must be non-negative; tags in a receive
!>     (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_TAG.
!>     The largest tag value is available through the the attribute MPI_TAG_UB.
!>   - MPI_ERR_RANK
!>   - Invalid source or destination rank. Ranks must be between zero and the size of the communicator minus one;
!>     ranks in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_SOURCE.
!> \param[in] fname    Filename of output \n
!>                     Skips if fname is not provided \n
!>                     File has first columns as tasks and last columns as results
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine task_dynamic(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
  implicit none
  integer, intent(in) :: Ntasks, Nfun, Nres, Nvec
  FUNTYPE, intent(in) :: tasks(Nfun,Ntasks)
  external :: func
  RESTYPE, intent(out) :: results(Nres,Ntasks)
  integer, intent(out) :: ierr
  character(len=*), optional, intent(in) :: fname
  integer :: size, rank
  real(dp) :: start

  ! Start timer
  start = MPI_Wtime()

  ! Get number of processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  ! Get number of this process
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ! Otherwise use root as master
  if ( rank.eq.0 ) then
    ! Root process assigns task
    call assign_tasks
  else
    ! Worker processes receive tasks
    call receive_tasks
  end if

  ! Broadcast results so all processes have same results
  if (bcast) call MPI_BCAST(results,Ntasks*Nres,MPIRESTYPE,0,MPI_COMM_WORLD,ierr)

  ! Write to file if present
  if (present(fname)) call export(Ntasks,Nfun,Nres,tasks,results,fname)
  
contains

  subroutine assign_tasks
    implicit none
    integer :: i, i0, i1, tag, source, status(MPI_STATUS_SIZE)
    RESTYPE :: buffer(Nres,Nvec)

    ! Assign processes first task
    do i=1,size-1
      i0 = min(1+(i-1)*Nvec, Ntasks)
      i1 = min(i*Nvec, Ntasks)
      call MPI_SEND(tasks(:,i0:i1), Nfun, MPIFUNTYPE, i, i, MPI_COMM_WORLD, ierr)
    end do

    ! Assign rest of tasks
    do while (1+(i-1)*Nvec .le. Ntasks)
      ! Wait for process to finish
      call MPI_RECV(buffer, Nres*Nvec, MPIRESTYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      call check_error(ierr)
      ! Get number of process
      source = status(MPI_SOURCE)
      ! Get number of task
      tag = status(MPI_TAG)
      ! Save result
      i0 = min(1+(tag-1)*Nvec, Ntasks)
      i1 = min(tag*Nvec, Ntasks)
      results(:,i0:i1) = buffer
      ! Write result
      if (sequential) call seq_write(Nfun,Nres,Nvec,tasks(:,i0:i1),buffer,fname)
      ! Show progress bar
      if (prog_bar) call progress(i+1-size, 1+(Ntasks-1)/Nvec, start)
      ! Send process next task
      i0 = min(1+(i-1)*Nvec, Ntasks)
      i1 = min(i*Nvec, Ntasks)
      call MPI_SEND(tasks(:,i0:i1), Nfun, MPIFUNTYPE, source, i, MPI_COMM_WORLD, ierr)
      call check_error(ierr)
      i = i + 1
    end do

    ! Receive last results and send finish signal
    do i=1,size-1
      ! Wait for process to finish
      call MPI_RECV(buffer, Nres*Nvec, MPIRESTYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      call check_error(ierr)
      ! Get number of process
      source = status(MPI_SOURCE)
      ! Get number of task
      tag = status(MPI_TAG)
      ! Save result
      i0 = min(1+(tag-1)*Nvec, Ntasks)
      i1 = min(tag*Nvec, Ntasks)
      results(:,i0:i1) = buffer(:,1:i1-i0+1)
      ! Write result
      if (sequential) call seq_write(Nfun,Nres,Nvec,tasks(:,i0:i1),buffer,fname)
      ! Show progress bar
      if (prog_bar) call progress(i+2-size+(Ntasks-1)/Nvec, 1+(Ntasks-1)/Nvec, start)
      ! Send finish signal to process
      call MPI_SEND(tasks(:,1:Nvec), Nfun, MPIFUNTYPE, source, 0, MPI_COMM_WORLD, ierr)
      call check_error(ierr)
    end do
  end subroutine assign_tasks

  subroutine receive_tasks
    implicit none
    integer :: tag, status(MPI_STATUS_SIZE)
    FUNTYPE :: task(Nfun,Nvec)
    RESTYPE :: result(Nres,Nvec)

    do
      ! Receive task
      call MPI_RECV(task, Nfun*Nvec, MPIFUNTYPE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      call check_error(ierr)
      ! Get number of task
      tag = status(MPI_TAG)
      ! Exit if all tasks are finished
      if (tag.eq.0) exit
      ! Do task
      call func(Nfun,task,Nres,result,Nvec)
      ! Send results
      call MPI_SEND(result, Nres*Nvec, MPIRESTYPE, 0, tag, MPI_COMM_WORLD, ierr)
      call check_error(ierr)
    end do
  end subroutine receive_tasks
end subroutine task_dynamic

!==============================================================
!  task_static
!> Assigns each process an equal number of tasks.
!> Takes an array of inputs (tasks) and passes it through
!> an externally declared subroutine (func).
!
!> \param[in]  Ntasks  Number of tasks \n
!>                     Size of tasks array
!> \param[in]  Nfun    Number of function Arguments \n
!>                     Size of arguments array for subroutine "func"
!> \param[in]  Nres    Number of results \n
!>                     Size of results array from subroutine "func"
!> \param[in]  Nvec    Size of vector \n
!>                     Number of tasks to do at a time
!> \param[in]  tasks   Array of function inputs \n
!>                     Arguments to use for subroutine "func"
!> \param[in]  func    externally declared subroutine to calculate tasks \n
!>                     Must have parameters (Nfun,fun,Nres,res)
!>                     - Inputs:
!>                       - Nfun - Size of fun array
!>                       - fun - Array of inputs
!>                       - Nres - Size of res array
!>                     - Output:
!>                       - res - Array of results
!> \param[out] results Array of results
!>                     Output of subroutine "func"
!> \param[out] ierr    MPI error return value
!>   - MPI_SUCCESS
!>   - No error; MPI routine completed successfully.
!>   - MPI_ERR_COMM
!>   - Invalid communicator. A common error is to use a null communicator in a call
!>     (not even allowed in MPI_Comm_rank).
!>   - MPI_ERR_COUNT
!>   - Invalid count argument. Count arguments must be non-negative; a count of zero is often valid.
!>   - MPI_ERR_TYPE
!>   - Invalid datatype argument. May be an uncommitted MPI_Datatype (see MPI_Type_commit).
!>   - MPI_ERR_TAG
!>   - Invalid tag argument. Tags must be non-negative; tags in a receive
!>     (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_TAG.
!>     The largest tag value is available through the the attribute MPI_TAG_UB.
!>   - MPI_ERR_RANK
!>   - Invalid source or destination rank. Ranks must be between zero and the size of the communicator minus one;
!>     ranks in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_SOURCE.
!> \param[in] fname    Filename of output \n
!>                     Skips if fname is not provided \n
!>                     File has first columns as tasks and last columns as results
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine task_static(Ntasks,Nfun,Nres,Nvec,tasks,func,results,ierr,fname)
  implicit none
  integer, intent(in) :: Ntasks, Nfun, Nres, Nvec
  FUNTYPE, intent(in) :: tasks(Nfun,Ntasks)
  external :: func
  RESTYPE, intent(out) :: results(Nres,Ntasks)
  integer, intent(out) :: ierr
  character(len=*), optional, intent(in) :: fname
  integer :: i, i1, j, size, rank, source, status(MPI_STATUS_SIZE)
  real(dp) :: start
  FUNTYPE :: task(Nfun,Nvec)
  RESTYPE :: result(Nres,Nvec), buffer(Nres,Ntasks)

  ! Start timer
  start = MPI_Wtime()

  ! Get number of processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  ! Get number of this process
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ! Divide tasks
  do i=1+rank*Nvec,Ntasks,size*Nvec
    i1 = min(i+Nvec-1, Ntasks)
    task = tasks(:,i:i1)
    ! Do task
    call func(Nfun,task,Nres,result,i1-i+1)
    ! Save result
    results(:,i:i1) = result(:,1:i1-i+1)
    ! Write result
    if (sequential) call seq_write(Nfun,Nres,Nvec,task,result,fname)
    ! Show progress bar
    if (prog_bar) call progress(i, Ntasks, start)
  end do

  ! Merge results
  if (rank.eq.0) then
    do j=1,size-1
      ! Receive results
      call MPI_RECV(buffer, Nres*Ntasks, MPIRESTYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      call check_error(ierr)
      ! Get number of process
      source = status(MPI_SOURCE)
      ! Save results
      do i=1+source*Nvec,Ntasks,size*Nvec
        i1 = min(i+Nvec-1, Ntasks)
        results(:,i:i1) = buffer(:,i:i1)
      end do
    end do
  else
    ! Send results to root
    call MPI_SEND(results, Nres*Ntasks, MPIRESTYPE, 0, 0, MPI_COMM_WORLD, ierr)
    call check_error(ierr)
  end if

  ! Broadcast results so all processes have same results
  if (bcast) call MPI_BCAST(results,Ntasks*Nres,MPIRESTYPE,0,MPI_COMM_WORLD,ierr)

  ! Write to file if present
  if (present(fname)) call export(Ntasks,Nfun,Nres,tasks,results,fname)
end subroutine task_static


!==============================================================
!  seq_write
!> Sequential write of tasks and results.
!> Write when tasks finish.
!
!> \param[in]  Nfun    Number of function Arguments \n
!>                     Size of arguments array for subroutine "func"
!> \param[in]  Nres    Number of results \n
!>                     Size of results array from subroutine "func"
!> \param[in]  task    Array of function inputs \n
!>                     Arguments to use for subroutine "func"
!> \param[out] result  Array of results
!>                     Output of subroutine "func"
!> \param[in]  fname   Filename of output \n
!>                     Skips if fname is not provided \n
!>                     File has first columns as tasks and last columns as results
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine seq_write(Nfun,Nres,Nvec,task,result,fname)
  integer, intent(in) :: Nfun, Nres, Nvec
  FUNTYPE, intent(in) :: task(Nfun,Nvec)
  RESTYPE, intent(in) :: result(Nres,Nvec)
  character(len=*), optional, intent(in) :: fname
  integer :: i, rank, ierr
  character(len=64) :: fname_MPI

  ! Stop if sequential writing turned off
  if (.not. present(fname)) return

  ! Get number of this process
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ! Add rank extension to filename
  write(fname_MPI,'(A,".",I0.3)') trim(fname), rank

  ! Export data to file
  if (stream) then
    ! Write as unformatted data
    open(10, file=fname_MPI, position='append', form='unformatted')
      do i=1,Nvec
        write(10) task(:,i), result(:,i)
      end do
    close(10)
  else
    ! Write as human readable data
    open(10, file=fname_MPI, position='append')
      do i=1,Nvec
        write(10,*) task(:,i), result(:,i)
      end do
    close(10)
  end if
end subroutine seq_write


!==============================================================
!  export
!> Write tasks and results to file
!
!> \param[in]  Ntasks  Number of tasks \n
!>                     Size of tasks array
!> \param[in]  Nfun    Number of function Arguments \n
!>                     Size of arguments array for subroutine "func"
!> \param[in]  Nres    Number of results \n
!>                     Size of results array from subroutine "func"
!> \param[in]  tasks   Array of function inputs \n
!>                     Arguments to use for subroutine "func"
!> \param[out] results Array of results
!>                     Output of subroutine "func"
!> \param[in]  fname   Filename of output \n
!>                     Skips if fname is not provided \n
!>                     File has first columns as tasks and last columns as results
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine export(Ntasks,Nfun,Nres,tasks,results,fname)
  integer, intent(in) :: Ntasks, Nfun, Nres
  FUNTYPE, intent(in) :: tasks(Nfun,Ntasks)
  RESTYPE, intent(in) :: results(Nres,Ntasks)
  character(len=*), optional, intent(in) :: fname
  integer :: i, rank, ierr

  ! Stop if no filename is present
  if (.not. present(fname)) return

  ! Get number of this process
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ! Export data to file
  if (rank.eq.0) then
    if (stream) then
      ! Write as unformatted data
      open(10, file=fname, status='replace', access='stream')
        do i=1,Ntasks
          write(10) tasks(:,i), results(:,i)
        end do
      close(10)
    else
      ! Write as human readable data
      open(10, file=fname, status='replace')
        do i=1,Ntasks
          write(10,*) tasks(:,i), results(:,i)
        end do
      close(10)
    end if
  end if
end subroutine export

!==============================================================
!  progress
!> Show progress bar
!
!> \param[in]  percent Percent completed
!> \param[in]  start   Time of start from system_clock
!> \date September 2016
!> \author Eric Mascot
!==============================================================

subroutine progress(i,imax,start)
  use iso_fortran_env
  implicit none
  integer, parameter :: w=30
  integer, intent(in) :: i, imax
  real(dp), intent(in) :: start
  integer :: ticks
  real(dp) :: percent, elapsed, remaining
  character(len=w+2) :: bar

  ticks = w*i/imax
  if (ticks > w) ticks = w
  if (ticks < 0) ticks = 0

  percent = real(i,dp)/real(imax,dp)
  elapsed = MPI_Wtime()-start
  remaining = int(elapsed*(1.0/percent-1.0))
  bar = "["//repeat("=",ticks)//repeat(" ",w-ticks)//"]"
  
  write(OUTPUT_UNIT,"(A,I3,'% ',I4,':',I2.2,' elapsed',I4,':',I2.2,' remaining')") &
    bar, int(percent*100), &
    int(elapsed)/3600, mod(int(elapsed)/60,60), &
    int(remaining)/3600, mod(int(remaining)/60,60)
  flush(OUTPUT_UNIT)
end subroutine progress

subroutine check_error(ierr)
  use iso_fortran_env
  implicit none
  integer :: ierr
  select case (ierr)
  case (MPI_ERR_COMM)
    write(ERROR_UNIT,*) 'MPI_ERR_COMM'
  case (MPI_ERR_TYPE)
    write(ERROR_UNIT,*) 'MPI_ERR_TYPE'
  case (MPI_ERR_COUNT)
    write(ERROR_UNIT,*) 'MPI_ERR_COUNT'
  case (MPI_ERR_TAG)
    write(ERROR_UNIT,*) 'MPI_ERR_TAG'
  case (MPI_ERR_RANK)
    write(ERROR_UNIT,*) 'MPI_ERR_RANK'
  end select
  flush(ERROR_UNIT)
end subroutine check_error
end module parallel_tasks
