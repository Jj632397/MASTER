module mdl_steps
  use mdl_input_manager
  implicit none

  !Physical time step
  integer, private, save :: nstep  = 0
  integer, private, save :: nstart = 0
  integer, private, save :: nfinal = 0

  !Inner time step
  integer, private, save :: nstep_it  = 0
  integer, private, save :: nstart_it = 0
  integer, private, save :: nfinal_it = 0

contains

  subroutine steps_init
    integer :: mode_restart
    call input_restart(nstart,mode_restart)
    nstep  = nstart
    call input_steps(nfinal)
    nfinal = nfinal+nstart
  end subroutine steps_init

  subroutine steps_final
    nstep = nstart
  end subroutine steps_final

  subroutine steps_get_step_now(ns)
    integer, intent(out) :: ns
    ns = nstep
  end subroutine steps_get_step_now

  subroutine steps_get_step_start(ns)
    integer, intent(out) :: ns
    ns = nstart
  end subroutine steps_get_step_start

  subroutine steps_get_step_final(ns)
    integer, intent(out) :: ns
    ns = nfinal
  end subroutine steps_get_step_final

  subroutine steps_increment
    nstep = nstep+1
  end subroutine steps_increment

  subroutine steps_isend(is)
    integer, intent(out) :: is
    if (nstep.ge.nfinal) then
      is = 1
    else
      is = 0
    endif
  end subroutine steps_isend
!-------------------------------------------------------------------------------
  subroutine steps_it_init
    integer :: mode,inner_iteration
    real(8) :: relaxation_factor
    call input_implicit(mode,inner_iteration,relaxation_factor)
    nstart_it = 0
    nstep_it  = nstart_it
    nfinal_it = inner_iteration
  end subroutine steps_it_init

  subroutine steps_it_final
    nstep_it  = nstart_it
  end subroutine steps_it_final

  subroutine steps_it_get_step_now(ns)
    integer, intent(out) :: ns
    ns = nstep_it
  end subroutine steps_it_get_step_now

  subroutine steps_it_get_step_start(ns)
    integer, intent(out) :: ns
    ns = nstart_it
  end subroutine steps_it_get_step_start

  subroutine steps_it_get_step_final(ns)
    integer, intent(out) :: ns
    ns = nfinal_it
  end subroutine steps_it_get_step_final

  subroutine steps_it_increment
    nstep_it = nstep_it+1
  end subroutine steps_it_increment

  subroutine steps_it_isend(is)
    integer, intent(out) :: is
    if (nstep_it.ge.nfinal_it) then
      is = 1
    else
      is = 0
    endif
  end subroutine steps_it_isend

end module mdl_steps
