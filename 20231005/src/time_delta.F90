module mdl_time_delta
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  implicit none

  !input------------------------------------------------------------------------
  integer, private, save :: mode !<0>:global time step  <1>:local time step
  real(8), private, save :: cfl
  real(8), private, save :: dt_minlimit
  real(8), private, save :: dt_maxlimit
  !-----------------------------------------------------------------------------

contains

  subroutine preconditioning_velocity_parameter(BLK)
    type(block), intent(inout) :: BLK

    integer :: i,j,k
    real(8) :: u,v,w,c,vu,vv

    real(8) :: pi
    real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
    real(8) :: urot,vrot,wrot

    pi = 4.0d0*datan(1.0d0)

    nx = 0.0d0
    ny = 0.0d0
    nz = 0.0d0
    if (frame_axis==1) nx = 1.0d0
    if (frame_axis==2) ny = 1.0d0
    if (frame_axis==3) nz = 1.0d0

    ox = 2.0d0*pi*frame_rpm/60.0d0*nx*DTREF
    oy = 2.0d0*pi*frame_rpm/60.0d0*ny*DTREF
    oz = 2.0d0*pi*frame_rpm/60.0d0*nz*DTREF

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          rx = BLK%XCRD(i,j,k)-frame_origin(1)
          ry = BLK%YCRD(i,j,k)-frame_origin(2)
          rz = BLK%ZCRD(i,j,k)-frame_origin(3)

          urot = oy*rz-oz*ry
          vrot = oz*rx-ox*rz
          wrot = ox*ry-oy*rx

          !Relative velocity
          u = BLK%ULCT(i,j,k) - urot
          v = BLK%VLCT(i,j,k) - vrot
          w = BLK%WLCT(i,j,k) - wrot

          c = BLK%CSOS(i,j,k)

          vu = UUPRE/UUREF
          vv = dsqrt(u*u+v*v+w*w)

!          BLK%VPRE(i,j,k) = dmax1(vu,dsqrt(3.0d0)*vv,c*precon_mach_cutoff)
          BLK%VPRE(i,j,k) = dmax1(3.0d0*vu,c*precon_mach_cutoff)
        enddo
      enddo
    enddo

  end subroutine preconditioning_velocity_parameter
!-------------------------------------------------------------------------------
  subroutine time_delta(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: i,j,k
    real(8) :: dx,dy,dz
    real(8) :: u,v,w,qq
    real(8) :: dxsi,deta,dzta,dl,ul,dt,dt_min,c,c2,vr,vr2

    real(8) :: pi
    real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
    real(8) :: urot,vrot,wrot

    call input_time_delta(mode,cfl,dt_minlimit,dt_maxlimit)

    call mpisub_sbsp_get_myrank_world(myrank)

    pi = 4.0d0*datan(1.0d0)

    nx = 0.0d0
    ny = 0.0d0
    nz = 0.0d0
    if (frame_axis==1) nx = 1.0d0
    if (frame_axis==2) ny = 1.0d0
    if (frame_axis==3) nz = 1.0d0

    ox = 2.0d0*pi*frame_rpm/60.0d0*nx*DTREF
    oy = 2.0d0*pi*frame_rpm/60.0d0*ny*DTREF
    oz = 2.0d0*pi*frame_rpm/60.0d0*nz*DTREF

    dt_min = 1.0d16

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          rx = BLK%XCRD(i,j,k)-frame_origin(1)
          ry = BLK%YCRD(i,j,k)-frame_origin(2)
          rz = BLK%ZCRD(i,j,k)-frame_origin(3)

          urot = oy*rz-oz*ry
          vrot = oz*rx-ox*rz
          wrot = ox*ry-oy*rx

          dx = BLK%XCRD_FACE_ETAZTA(i+1,j,k)-BLK%XCRD_FACE_ETAZTA(i,j,k)
          dy = BLK%YCRD_FACE_ETAZTA(i+1,j,k)-BLK%YCRD_FACE_ETAZTA(i,j,k)
          dZ = BLK%ZCRD_FACE_ETAZTA(i+1,j,k)-BLK%ZCRD_FACE_ETAZTA(i,j,k)
          if (i2dimension>0) then
            dxsi = dsqrt(dx*dx+dy*dy)
          else
            dxsi = dsqrt(dx*dx+dy*dy+dz*dz)
          endif

          dx = BLK%XCRD_FACE_ZTAXSI(i,j+1,k)-BLK%XCRD_FACE_ZTAXSI(i,j,k)
          dy = BLK%YCRD_FACE_ZTAXSI(i,j+1,k)-BLK%YCRD_FACE_ZTAXSI(i,j,k)
          dZ = BLK%ZCRD_FACE_ZTAXSI(i,j+1,k)-BLK%ZCRD_FACE_ZTAXSI(i,j,k)
          if (i2dimension>0) then
            deta = dsqrt(dx*dx+dy*dy)
          else
            deta = dsqrt(dx*dx+dy*dy+dz*dz)
          endif

          if (i2dimension>0) then
            dzta = 0.0d0
          else
            dx = BLK%XCRD_FACE_XSIETA(i,j,k+1)-BLK%XCRD_FACE_XSIETA(i,j,k)
            dy = BLK%YCRD_FACE_XSIETA(i,j,k+1)-BLK%YCRD_FACE_XSIETA(i,j,k)
            dZ = BLK%ZCRD_FACE_XSIETA(i,j,k+1)-BLK%ZCRD_FACE_XSIETA(i,j,k)
            dzta = dsqrt(dx*dx+dy*dy+dz*dz)
          endif

          if (i2dimension>0) then
            dl = dmin1(dxsi,deta)
          else
            dl = dmin1(dxsi,deta,dzta)
          endif

          !Relative velocity
          u = BLK%ULCT(i,j,k) - urot
          v = BLK%VLCT(i,j,k) - vrot
          w = BLK%WLCT(i,j,k) - wrot

          c = BLK%CSOS(i,j,k)

          vr = dmin1(BLK%VPRE(i,j,k),c)

          c2 = c*c
          vr2= vr*vr

          qq = u*u+v*v+w*w

          ul = 0.5d0*(dsqrt(qq)*(1.0d0+vr2/c2)+dsqrt(qq*(1.0d0-vr2/c2)*(1.0d0-vr2/c2)+4.0d0*vr2))

          dt = cfl*dl/ul

          dt = dmax1(dt,dt_minlimit)
          dt = dmin1(dt,dt_maxlimit)

          BLK%DTLC(i,j,k) = dmin1(dt,physical_time_step)

          dt_min = dmin1(dt_min,dt)
        enddo
      enddo
    enddo


    if (mode/=1) then
      call mpisub_sbsp_allreduce_min(dt_min)

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            BLK%DTLC(i,j,k) = dt_min
          enddo
        enddo
      enddo

    endif

  end subroutine time_delta

end module mdl_time_delta
