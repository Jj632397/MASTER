module mdl_refs
  implicit none

  real(8), public, save :: XYREF
  real(8), public, save :: RHREF
  real(8), public, save :: UUREF
  real(8), public, save :: TTREF
  real(8), public, save :: CPREF
  real(8), public, save :: CKREF
  real(8), public, save :: VSREF
  real(8), public, save :: DFREF
  real(8), public, save :: DTREF
  real(8), public, save :: PPREF
  !-----------------------------
  real(8), public, save :: MACHN
  real(8), public, save :: RYNLD
  !-----------------------------
  real(8), public, save :: MACHP !Reference Mach number for precondition
  real(8), public, save :: UUPRE !Reference velocity for precondition

end module mdl_refs
