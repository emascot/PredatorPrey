module chase
  implicit none
  ! Double precision
  integer, parameter :: dp = kind(0.d0)
  ! Total number of test velocities
  integer, parameter :: nvel = 101
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

contains

  ! Calculate position at each time for given velocities
  subroutine chase_trajectory(ndim,velocities,nfun,positions,nvec)
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
    do ivel = 1,nvec
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
  end subroutine chase_trajectory

  ! Calculate catch time for given velocities
  subroutine chase_time(ndim,velocity,nfun,time,nvec)
    use constants, only: dp
    implicit none
    integer, intent(in) :: ndim, nfun, nvec
    real(dp), intent(in) :: velocity(ndim,nvec)
    complex(dp), intent(out) :: time(nfun,nvec)
    integer :: ivec

    do ivec=1,nvec
      time(1,ivec) = catch_time(velocity(1,ivec))
    enddo
  end subroutine chase_time

  ! Return time when prey is caught
  real(dp) function catch_time(velocity)
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
    ! Distance between position and target
    real(dp) :: dist

    catch_time = 0._dp
    position = cmplx(0.d0, 10.d0, dp)
    do it = 1,nsteps
      ! Find prey
      time = it * tstep
      prey_position = cmplx(prey_vel * time, 0.d0, dp)
      ! Aim towards prey
      diff = prey_position - position
      dist = abs(diff)
      ! Check if within reach
      if (dist < velocity * tstep + range) then
        catch_time = time
        return
      else
        ! Move towards prey
        position = position + velocity * tstep * diff / dist
      endif
    enddo
  end function catch_time

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

    ! Repeat ntest times
    do itest=1,ntest
      ! Bisect range
      velocity = (min + max) * 0.5d0
      if (catch_time(velocity) > 0.d0) then
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
end module chase
