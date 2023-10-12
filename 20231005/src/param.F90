module mdl_param
  implicit none

  real(8), parameter, public :: physical_time_step = 1.0d-3 !Nondimensional value.
                                                            !For steady state calculation,
                                                            !set an arbitrary large value such as 10^16.

  real(8), parameter, public :: precon_mach_cutoff = 1.0d-6 !For disabling preconditioning,
                                                            !set an arbitrary large value such as 10^16.

  !<positive>:2-dimension <0 or negative>:3-dimension
  integer, parameter, public :: i2dimension = 0

  !<0>:inactive <1>:active
  integer, parameter, public :: iact_lp_filter = 0

  integer, parameter, public :: neqns = 7  ! Furusawa

  integer, parameter, public :: halo_width = 3

  integer, parameter, public :: iact_tbl_sa = 0
  integer, parameter, public :: neqn_tbl_sa = 6
  integer, parameter, public :: step_tbl_sa = 50
  real(8), parameter, public :: tkin_tbl_sa = 5.0d0

  integer, parameter, public :: iact_passive_scalar = 1
  integer, parameter, public :: neqn_passive_scalar = 6    ! Furusawa 6, 7 are passive scalar 

  !parameters for rotating reference frame
  real(8), public, parameter               :: frame_rpm    = 0.0d0           ! rotating speed [rpm]
  integer, public, parameter               :: frame_axis   = 1               ! <1>:x <2>:y <3>:z
  real(8), public, parameter, dimension(3) :: frame_origin = (/0.0,0.0,0.0/) ! Non-dimensional value

  !Specific heat at constant pressure for ideal gas
  real(8), parameter, public :: cp_i = 1.004d3

  !Gas constant for ideal gas
  real(8), parameter, public :: cr_i = 8.3144d0/28.96d-3

  !Prandtl number
  real(8), parameter, public :: prntl  = 0.72d0
  real(8), parameter, public :: tprntl = 0.98d0

  !Schmidt number
!!  real(8), parameter, public :: schmt  = 1.0d0  !! FURUSAWA
  real(8), parameter, public :: tschmt = 0.7d0

end module mdl_param
