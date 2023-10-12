module mdl_output
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_steps
  use mdl_exchange_domain_halo
  use mdl_halo
  implicit none

contains

  subroutine output_grid(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    fname = 'data/grid.dat'

    lmax = 3

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QB(ii,jj,kk,1) = BLK%XCRD(i,j,k)
          QB(ii,jj,kk,2) = BLK%YCRD(i,j,k)
          QB(ii,jj,kk,3) = BLK%ZCRD(i,j,k)
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d(fname,0,QB,imax,jmax,kmax,lmax)

    deallocate(QB)

  end subroutine output_grid
!-------------------------------------------------------------------------------
  subroutine output_grid_planes(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    integer :: ndomain,id_dm
    integer, allocatable, dimension(:) :: id_dms
    integer, allocatable, dimension(:) :: id_dirs
    integer, allocatable, dimension(:) :: idxs

    fname = 'data/ani/grid_2d.dat'

    lmax = 3

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          QB(ii,jj,kk,1) = BLK%XCRD(i,j,k)
          QB(ii,jj,kk,2) = BLK%YCRD(i,j,k)
          QB(ii,jj,kk,3) = BLK%ZCRD(i,j,k)
        enddo
      enddo
    enddo

    call Get_DomainNum(ndomain)

    allocate(id_dms(ndomain))
    allocate(id_dirs(ndomain))
    allocate(idxs(ndomain))

    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      id_dms(id_dm) = id_dm
      id_dirs(id_dm)= 3
      idxs(id_dm) = kmax/2+mod(kmax,2)
    enddo

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call mpisub_sbsp_out_plot3d_planes(fname,0,QB,imax,jmax,kmax,lmax, &
         id_dms,id_dirs,idxs,ndomain)

    deallocate(id_dms)
    deallocate(id_dirs)
    deallocate(idxs)

    deallocate(QB)

  end subroutine output_grid_planes
!-------------------------------------------------------------------------------
  subroutine output_grid_block_with_halo(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: lmax
    integer :: i,j,k

    real(8), allocatable, dimension(:,:,:,:) :: QB_WH

800 format('data/block/grid_b',i5.5,'.dat')

    lmax = 3

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    write(fname,800) id_block

    allocate(QB_WH(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax))

    do k=1,BLK%kmax_wh
      do j=1,BLK%jmax_wh
        do i=1,BLK%imax_wh
          QB_WH(i,j,k,1) = BLK%XCRD(i,j,k)
          QB_WH(i,j,k,2) = BLK%YCRD(i,j,k)
          QB_WH(i,j,k,3) = BLK%ZCRD(i,j,k)
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d_block_with_halo(fname,0,halo_width,IFLAG_MTRC, &
         QB_WH,BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax)

    deallocate(QB_WH)

  end subroutine output_grid_block_with_halo
!-------------------------------------------------------------------------------
  subroutine output_phys(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk
    integer :: nstep,nstep_it

    real(8) :: xsix,xsiy,xsiz
    real(8) :: etax,etay,etaz
    real(8) :: ztax,ztay,ztaz
    real(8) :: gxsi,geta,gzta

    real(8) :: dt
    real(8) :: rh,u,v,w,p,t,c,ys,knet,di
    real(8) :: c2,vp,vp2,vpd,ggdxsi,ggdeta,ggdzta
    real(8) :: rhx,rhy,rhz,ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(8) :: vrtxx,vrtxy,vrtxz
    real(8) :: smgs,smga
    real(8) :: uc,vc,wc,ul,cfl

    real(8) :: pi
    real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
    real(8) :: urot,vrot,wrot
    real(8) :: ac,bc,gc

    real(8), allocatable, dimension(:,:,:,:) :: QB

800 format('data/phys_cell_n',i7.7,'.dat')

900 format('data/phys_cell_n',i7.7,'_p',i7.7,'.dat')

    call steps_get_step_now(nstep)
    call steps_it_get_step_now(nstep_it)

    if (nstep_it==0) then
      write(fname,800) nstep
    else
      write(fname,900) nstep,nstep_it
    endif

    lmax = 14

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    !call schs_source(BLK,RHS)

    allocate(QB(imax,jmax,kmax,lmax))

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
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          rx = BLK%XCRD(i,j,k)-frame_origin(1)
          ry = BLK%YCRD(i,j,k)-frame_origin(2)
          rz = BLK%ZCRD(i,j,k)-frame_origin(3)

          urot = oy*rz-oz*ry
          vrot = oz*rx-ox*rz
          wrot = ox*ry-oy*rx

          xsix = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          xsiy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
          xsiz = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)
          etax = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          etay = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
          etaz = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)
          ztax = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
          ztay = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
          ztaz = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)

          rh = BLK%DNST(i,j,k)
          u  = BLK%ULCT(i,j,k)
          v  = BLK%VLCT(i,j,k)
          w  = BLK%WLCT(i,j,k)

          p  = BLK%PRSS(i,j,k)
          t  = BLK%TMPR(i,j,k)
          c  = BLK%CSOS(i,j,k)
          knet=BLK%KNET(i,j,k)
          di = BLK%DIEL(i,j,k)

          rhx = BLK%DDNX(i,j,k)
          rhy = BLK%DDNY(i,j,k)
          rhz = BLK%DDNZ(i,j,k)
          ux  = BLK%DULX(i,j,k)
          uy  = BLK%DULY(i,j,k)
          uz  = BLK%DULZ(i,j,k)
          vx  = BLK%DVLX(i,j,k)
          vy  = BLK%DVLY(i,j,k)
          vz  = BLK%DVLZ(i,j,k)
          wx  = BLK%DWLX(i,j,k)
          wy  = BLK%DWLY(i,j,k)
          wz  = BLK%DWLZ(i,j,k)

          dt = BLK%DTLC(i,j,k)

          vrtxx = wy-vz
          vrtxy = uz-wx
          vrtxz = vx-uy

          smgs = dsqrt(ux*ux+vy*vy+wz*wz &
                      +0.5d0*(uy+vx)**2+0.5d0*(wx+uz)**2+0.5d0*(vz+wy)**2)
          smga = dsqrt(0.25d0*((uy-vx)**2+(uz-wx)**2+(vx-uy)**2 &
                              +(vz-wy)**2+(wx-uz)**2+(wy-vz)**2))

          c2 = c*c

          vp = BLK%VPRE(i,j,k)
          vp = dmin1(vp,c)
          vp2= vp*vp

          vpd= vp2/c2

          uc = u*xsix+v*xsiy+w*xsiz
          vc = u*etax+v*etay+w*etaz
          wc = u*ztax+v*ztay+w*ztaz

          ac = urot*xsix+vrot*xsiy+wrot*xsiz
          bc = urot*etax+vrot*etay+wrot*etaz
          gc = urot*ztax+vrot*ztay+wrot*ztaz

          !Relative velocity
          uc = dabs(uc-ac)
          vc = dabs(vc-bc)
          wc = dabs(wc-gc)

          gxsi = dsqrt(xsix*xsix+xsiy*xsiy+xsiz*xsiz)
          geta = dsqrt(etax*etax+etay*etay+etaz*etaz)
          gzta = dsqrt(ztax*ztax+ztay*ztay+ztaz*ztaz)

          ggdxsi = dsqrt(uc*uc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gxsi*gxsi)
          ggdeta = dsqrt(vc*vc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*geta*geta)
          ggdzta = dsqrt(wc*wc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gzta*gzta)

          ul= dmax1( 0.5d0*(uc*(1.0d0+vpd)+ggdxsi) &
                    ,0.5d0*(vc*(1.0d0+vpd)+ggdeta) &
                    ,0.5d0*(wc*(1.0d0+vpd)+ggdzta) )

          cfl = ul*dt

          QB(ii,jj,kk,1 ) = rh*RHREF
          QB(ii,jj,kk,2 ) = u *UUREF
          QB(ii,jj,kk,3 ) = v *UUREF
          QB(ii,jj,kk,4 ) = w *UUREF
          QB(ii,jj,kk,5 ) = p *PPREF
          QB(ii,jj,kk,6 ) = t *TTREF
          QB(ii,jj,kk,7 ) = dsqrt(u*u+v*v+w*w)/c
          QB(ii,jj,kk,8 ) = vrtxz
          QB(ii,jj,kk,9 ) = 0.5d0*(smga*smga-smgs*smgs)
          QB(ii,jj,kk,10) = BLK%VSST(i,j,k)
          QB(ii,jj,kk,11) = dt
!!F          QB(ii,jj,kk,12) = cfl
          if (iact_passive_scalar>0) then
            QB(ii,jj,kk,11) = BLK%QCS1(i,j,k,neqn_passive_scalar)
            QB(ii,jj,kk,12) = BLK%QCS1(i,j,k,neqn_passive_scalar+1)
          endif
           QB(ii,jj,kk,13) = knet
           QB(ii,jj,kk,14) = di
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d(fname,1,QB,imax,jmax,kmax,lmax)

    deallocate(QB)

  end subroutine output_phys
!-------------------------------------------------------------------------------
  subroutine output_phys_planes(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk
    integer :: nstep,nstep_it

    real(8) :: xsix,xsiy,xsiz
    real(8) :: etax,etay,etaz
    real(8) :: ztax,ztay,ztaz
    real(8) :: gxsi,geta,gzta

    real(8) :: dt
    real(8) :: rh,u,v,w,p,t,c
    real(8) :: c2,vp,vp2,vpd,ggdxsi,ggdeta,ggdzta
    real(8) :: rhx,rhy,rhz,ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(8) :: vrtxx,vrtxy,vrtxz
    real(8) :: smgs,smga
    real(8) :: uc,vc,wc,ul,cfl

    real(8) :: pi
    real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
    real(8) :: urot,vrot,wrot
    real(8) :: ac,bc,gc

    real(8), allocatable, dimension(:,:,:,:) :: QB

    integer :: ndomain,id_dm
    integer, allocatable, dimension(:) :: id_dms
    integer, allocatable, dimension(:) :: id_dirs
    integer, allocatable, dimension(:) :: idxs

800 format('data/ani/phys_2d_n',i7.7,'.dat')

900 format('data/ani/phys_2d_n',i7.7,'_p',i7.7,'.dat')

    call steps_get_step_now(nstep)
    call steps_it_get_step_now(nstep_it)

    if (nstep_it==0) then
      write(fname,800) nstep
    else
      write(fname,900) nstep,nstep_it
    endif

    lmax = 12

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

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
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          rx = BLK%XCRD(i,j,k)-frame_origin(1)
          ry = BLK%YCRD(i,j,k)-frame_origin(2)
          rz = BLK%ZCRD(i,j,k)-frame_origin(3)

          urot = oy*rz-oz*ry
          vrot = oz*rx-ox*rz
          wrot = ox*ry-oy*rx

          xsix = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          xsiy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
          xsiz = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)
          etax = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          etay = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
          etaz = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)
          ztax = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
          ztay = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
          ztaz = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)

          rh = BLK%DNST(i,j,k)
          u  = BLK%ULCT(i,j,k)
          v  = BLK%VLCT(i,j,k)
          w  = BLK%WLCT(i,j,k)

          p  = BLK%PRSS(i,j,k)
          t  = BLK%TMPR(i,j,k)
          c  = BLK%CSOS(i,j,k)

          rhx = BLK%DDNX(i,j,k)
          rhy = BLK%DDNY(i,j,k)
          rhz = BLK%DDNZ(i,j,k)
          ux  = BLK%DULX(i,j,k)
          uy  = BLK%DULY(i,j,k)
          uz  = BLK%DULZ(i,j,k)
          vx  = BLK%DVLX(i,j,k)
          vy  = BLK%DVLY(i,j,k)
          vz  = BLK%DVLZ(i,j,k)
          wx  = BLK%DWLX(i,j,k)
          wy  = BLK%DWLY(i,j,k)
          wz  = BLK%DWLZ(i,j,k)

          dt = BLK%DTLC(i,j,k)

          vrtxx = wy-vz
          vrtxy = uz-wx
          vrtxz = vx-uy

          smgs = dsqrt(ux*ux+vy*vy+wz*wz &
                      +0.5d0*(uy+vx)**2+0.5d0*(wx+uz)**2+0.5d0*(vz+wy)**2)
          smga = dsqrt(0.25d0*((uy-vx)**2+(uz-wx)**2+(vx-uy)**2 &
                              +(vz-wy)**2+(wx-uz)**2+(wy-vz)**2))

          c2 = c*c

          vp = BLK%VPRE(i,j,k)
          vp = dmin1(vp,c)
          vp2= vp*vp

          vpd= vp2/c2

          uc = u*xsix+v*xsiy+w*xsiz
          vc = u*etax+v*etay+w*etaz
          wc = u*ztax+v*ztay+w*ztaz

          ac = urot*xsix+vrot*xsiy+wrot*xsiz
          bc = urot*etax+vrot*etay+wrot*etaz
          gc = urot*ztax+vrot*ztay+wrot*ztaz

          !Relative velocity
          uc = dabs(uc-ac)
          vc = dabs(vc-bc)
          wc = dabs(wc-gc)

          gxsi = dsqrt(xsix*xsix+xsiy*xsiy+xsiz*xsiz)
          geta = dsqrt(etax*etax+etay*etay+etaz*etaz)
          gzta = dsqrt(ztax*ztax+ztay*ztay+ztaz*ztaz)

          ggdxsi = dsqrt(uc*uc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gxsi*gxsi)
          ggdeta = dsqrt(vc*vc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*geta*geta)
          ggdzta = dsqrt(wc*wc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gzta*gzta)

          ul= dmax1( 0.5d0*(uc*(1.0d0+vpd)+ggdxsi) &
                    ,0.5d0*(vc*(1.0d0+vpd)+ggdeta) &
                    ,0.5d0*(wc*(1.0d0+vpd)+ggdzta) )

          cfl = ul*dt

          QB(ii,jj,kk,1 ) = rh*RHREF
          QB(ii,jj,kk,2 ) = u *UUREF
          QB(ii,jj,kk,3 ) = v *UUREF
          QB(ii,jj,kk,4 ) = w *UUREF
          QB(ii,jj,kk,5 ) = p *PPREF
          QB(ii,jj,kk,6 ) = t *TTREF
          QB(ii,jj,kk,7 ) = dsqrt(u*u+v*v+w*w)/c
          QB(ii,jj,kk,8 ) = vrtxz
          QB(ii,jj,kk,9 ) = 0.5d0*(smga*smga-smgs*smgs)
          QB(ii,jj,kk,10) = BLK%VSST(i,j,k)
          QB(ii,jj,kk,11) = dt
          if (iact_passive_scalar>0) then
            QB(ii,jj,kk,11) = BLK%QCS1(i,j,k,neqn_passive_scalar)
          endif
          QB(ii,jj,kk,12) = cfl
        enddo
      enddo
    enddo

    call Get_DomainNum(ndomain)

    allocate(id_dms(ndomain))
    allocate(id_dirs(ndomain))
    allocate(idxs(ndomain))

    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      id_dms(id_dm) = id_dm
      id_dirs(id_dm)= 3
      idxs(id_dm) = kmax/2+mod(kmax,2)
    enddo

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    call mpisub_sbsp_out_plot3d_planes(fname,1,QB,imax,jmax,kmax,lmax, &
         id_dms,id_dirs,idxs,ndomain)

    deallocate(id_dms)
    deallocate(id_dirs)
    deallocate(idxs)

    deallocate(QB)

  end subroutine output_phys_planes
!-------------------------------------------------------------------------------
  subroutine output_phys_block_with_halo(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: lmax
    integer :: i,j,k
    integer :: nstep,nstep_it

    real(8) :: xsix,xsiy,xsiz
    real(8) :: etax,etay,etaz
    real(8) :: ztax,ztay,ztaz
    real(8) :: gxsi,geta,gzta

    real(8) :: dt
    real(8) :: rh,u,v,w,p,t,c,ys
    real(8) :: c2,vp,vp2,vpd,ggdxsi,ggdeta,ggdzta
    real(8) :: rhx,rhy,rhz,ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(8) :: vrtxx,vrtxy,vrtxz
    real(8) :: smgs,smga
    real(8) :: uc,vc,wc,ul,cfl

    real(8) :: pi
    real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
    real(8) :: urot,vrot,wrot
    real(8) :: ac,bc,gc

    real(8), allocatable, dimension(:,:,:,:) :: QB_WH

    integer, allocatable, dimension(:,:,:) :: MASK

800 format('data/block/phys_cell_b',i5.5,'_n',i7.7,'.dat')

900 format('data/block/phys_cell_b',i5.5,'_n',i7.7,'_p',i7.7,'.dat')

    allocate(MASK(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh))

    MASK(:,:,:) = 0

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          MASK(i,j,k) = 1
        enddo
      enddo
    enddo

    if (IFLAG_PRIM(1)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(2)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(3)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(4)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_PRIM(5)==1)) then
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_PRIM(6)==1)) then
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call steps_get_step_now(nstep)

    call steps_it_get_step_now(nstep_it)

    if (nstep_it==0) then
      write(fname,800) id_block,nstep
    else
      write(fname,900) id_block,nstep,nstep_it
    endif

    lmax = 12

    allocate(QB_WH(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax))
    QB_WH(:,:,:,:) = 0.0d0

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

    do k=1,BLK%kmax_wh
      do j=1,BLK%jmax_wh
        do i=1,BLK%imax_wh
          if (MASK(i,j,k)/=1) cycle

          rx = BLK%XCRD(i,j,k)-frame_origin(1)
          ry = BLK%YCRD(i,j,k)-frame_origin(2)
          rz = BLK%ZCRD(i,j,k)-frame_origin(3)

          urot = oy*rz-oz*ry
          vrot = oz*rx-ox*rz
          wrot = ox*ry-oy*rx

          xsix = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          xsiy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
          xsiz = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)
          etax = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          etay = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
          etaz = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)
          ztax = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
          ztay = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
          ztaz = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)

          rh = BLK%DNST(i,j,k)
          u  = BLK%ULCT(i,j,k)
          v  = BLK%VLCT(i,j,k)
          w  = BLK%WLCT(i,j,k)

          p  = BLK%PRSS(i,j,k)
          t  = BLK%TMPR(i,j,k)
          c  = BLK%CSOS(i,j,k)

          rhx = BLK%DDNX(i,j,k)
          rhy = BLK%DDNY(i,j,k)
          rhz = BLK%DDNZ(i,j,k)
          ux  = BLK%DULX(i,j,k)
          uy  = BLK%DULY(i,j,k)
          uz  = BLK%DULZ(i,j,k)
          vx  = BLK%DVLX(i,j,k)
          vy  = BLK%DVLY(i,j,k)
          vz  = BLK%DVLZ(i,j,k)
          wx  = BLK%DWLX(i,j,k)
          wy  = BLK%DWLY(i,j,k)
          wz  = BLK%DWLZ(i,j,k)

          dt = BLK%DTLC(i,j,k)

          vrtxx = wy-vz
          vrtxy = uz-wx
          vrtxz = vx-uy

          smgs = dsqrt(ux*ux+vy*vy+wz*wz &
                      +0.5d0*(uy+vx)**2+0.5d0*(wx+uz)**2+0.5d0*(vz+wy)**2)
          smga = dsqrt(0.25d0*((uy-vx)**2+(uz-wx)**2+(vx-uy)**2 &
                              +(vz-wy)**2+(wx-uz)**2+(wy-vz)**2))

          c2 = c*c

          vp = BLK%VPRE(i,j,k)
          vp = dmin1(vp,c)
          vp2= vp*vp

          vpd= vp2/c2

          uc = u*xsix+v*xsiy+w*xsiz
          vc = u*etax+v*etay+w*etaz
          wc = u*ztax+v*ztay+w*ztaz

          ac = urot*xsix+vrot*xsiy+wrot*xsiz
          bc = urot*etax+vrot*etay+wrot*etaz
          gc = urot*ztax+vrot*ztay+wrot*ztaz

          !Relative velocity
          uc = dabs(uc-ac)
          vc = dabs(vc-bc)
          wc = dabs(wc-gc)

          gxsi = dsqrt(xsix*xsix+xsiy*xsiy+xsiz*xsiz)
          geta = dsqrt(etax*etax+etay*etay+etaz*etaz)
          gzta = dsqrt(ztax*ztax+ztay*ztay+ztaz*ztaz)

          ggdxsi = dsqrt(uc*uc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gxsi*gxsi)
          ggdeta = dsqrt(vc*vc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*geta*geta)
          ggdzta = dsqrt(wc*wc*(1.0d0-vpd)*(1.0d0-vpd)+4.0d0*vp2*gzta*gzta)

          ul= dmax1( 0.5d0*(uc*(1.0d0+vpd)+ggdxsi) &
                    ,0.5d0*(vc*(1.0d0+vpd)+ggdeta) &
                    ,0.5d0*(wc*(1.0d0+vpd)+ggdzta) )

          cfl = ul*dt

          QB_WH(i,j,k,1 ) = rh*RHREF
          QB_WH(i,j,k,2 ) = u *UUREF
          QB_WH(i,j,k,3 ) = v *UUREF
          QB_WH(i,j,k,4 ) = w *UUREF
          QB_WH(i,j,k,5 ) = p *PPREF
          QB_WH(i,j,k,6 ) = t *TTREF
          QB_WH(i,j,k,7 ) = dsqrt(u*u+v*v+w*w)/c
          QB_WH(i,j,k,8 ) = vrtxz
          QB_WH(i,j,k,9 ) = 0.5d0*(smga*smga-smgs*smgs)
          QB_WH(i,j,k,10) = BLK%VSST(i,j,k)
          QB_WH(i,j,k,11) = dt
          if (iact_passive_scalar>0) then
            QB_WH(i,j,k,11) = BLK%QCS1(i,j,k,neqn_passive_scalar)
          endif
          QB_WH(i,j,k,12) = cfl
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d_block_with_halo(fname,1,halo_width,IFLAG_PRIM, &
         QB_WH,BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax)

    deallocate(QB_WH)
    deallocate(MASK)

  end subroutine output_phys_block_with_halo
!-------------------------------------------------------------------------------
  subroutine output_qout(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk,l
    integer :: nstep

    real(8), allocatable, dimension(:,:,:,:) :: QB

908 format('data/qout_n',i7.7,'.dat')

    call steps_get_step_now(nstep)

    write(fname,908) nstep

    lmax = neqns*3

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QB(ii,jj,kk,l) = BLK%QOUT(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d(fname,1,QB,imax,jmax,kmax,lmax)

    deallocate(QB)

  end subroutine output_qout
!-------------------------------------------------------------------------------
  subroutine output_qout_block_with_halo(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: lmax
    integer :: i,j,k,l
    integer :: nstep

    real(8), allocatable, dimension(:,:,:,:) :: QB_WH

    integer, allocatable, dimension(:,:,:) :: MASK

908 format('data/block/qout_b',i5.5,'_n',i7.7,'.dat')

    allocate(MASK(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh))

    MASK(:,:,:) = 0

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          MASK(i,j,k) = 1
        enddo
      enddo
    enddo

    if (IFLAG_PRIM(1)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(2)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(3)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_PRIM(4)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_PRIM(5)==1)) then
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_PRIM(6)==1)) then
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    lmax = neqns*3

    call steps_get_step_now(nstep)

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    allocate(QB_WH(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax))

    QB_WH(:,:,:,:) = 0.0d0

    do l=1,lmax
      do k=1,BLK%kmax_wh
        do j=1,BLK%jmax_wh
          do i=1,BLK%imax_wh
            if (MASK(i,j,k)/=1) cycle

            QB_WH(i,j,k,l) = BLK%QOUT(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    write(fname,908) id_block,nstep

    call mpisub_sbsp_out_plot3d_block_with_halo(fname,1,halo_width,IFLAG_PRIM, &
         QB_WH,BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax)

    deallocate(QB_WH)
    deallocate(MASK)

  end subroutine output_qout_block_with_halo
!-------------------------------------------------------------------------------
  subroutine output_xsix(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    fname = 'data/xsix.dat'

    lmax = 10

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QB(ii,jj,kk,1 ) = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,2 ) = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,3 ) = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,4 ) = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,5 ) = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,6 ) = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,7 ) = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,8 ) = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,9 ) = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)
          QB(ii,jj,kk,10) = BLK%RJCB(i,j,k)
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d(fname,1,QB,imax,jmax,kmax,lmax)

    deallocate(QB)

  end subroutine output_xsix
!-------------------------------------------------------------------------------
  subroutine output_xsix_block_with_halo(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: lmax
    integer :: i,j,k

    real(8), allocatable, dimension(:,:,:,:) :: QB_WH

    integer, allocatable, dimension(:,:,:) :: MASK

800 format('data/block/xsix_b',i5.5,'.dat')

    allocate(MASK(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh))

    MASK(:,:,:) = 0

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          MASK(i,j,k) = 1
        enddo
      enddo
    enddo

    if (IFLAG_MTRC(1)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista_wh,BLK%ista_wh+halo_width-1
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_MTRC(2)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%iend_wh-halo_width+1,BLK%iend_wh
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_MTRC(3)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta_wh,BLK%jsta_wh+halo_width-1
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if (IFLAG_MTRC(4)==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jend_wh-halo_width+1,BLK%jend_wh
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_MTRC(5)==1)) then
      do k=BLK%ksta_wh,BLK%ksta_wh+halo_width-1
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    if ((.not.(i2dimension>0)).and.(IFLAG_MTRC(6)==1)) then
      do k=BLK%kend_wh-halo_width+1,BLK%kend_wh
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            MASK(i,j,k) = 1
          enddo
        enddo
      enddo
    endif

    lmax = 10

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    write(fname,800) id_block

    allocate(QB_WH(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax))

    QB_WH(:,:,:,:) = 0.0d0

    do k=1,BLK%kmax_wh
      do j=1,BLK%jmax_wh
        do i=1,BLK%imax_wh
          if (MASK(i,j,k)/=1) cycle

          QB_WH(i,j,k,1 ) = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,2 ) = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,3 ) = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,4 ) = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,5 ) = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,6 ) = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,7 ) = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,8 ) = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,9 ) = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)
          QB_WH(i,j,k,10) = BLK%RJCB(i,j,k)
        enddo
      enddo
    enddo

    call mpisub_sbsp_out_plot3d_block_with_halo(fname,1,halo_width,IFLAG_MTRC, &
         QB_WH,BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,lmax)

    deallocate(QB_WH)
    deallocate(MASK)

  end subroutine output_xsix_block_with_halo
!-------------------------------------------------------------------------------
  subroutine output_l2norm(BLK)
    type(block), intent(in) :: BLK

    character(60) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: time_integral_mode
    integer :: i,j,k,l,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    integer :: nstep,step_start
    real(8), allocatable, dimension(:) :: l2norm

    805 format(6A13)
    806 format(7A13)
    807 format(8A13)

    905 format(I13,5E13.5)
    906 format(I13,6E13.5)
    907 format(I13,7E13.5)

    fname = 'data/l2norm_of_physical_time_step.dat'

    call input_time_integral(time_integral_mode)

    call mpisub_sbsp_get_myrank_world(myrank)

    call steps_get_step_now(nstep)

    call steps_get_step_start(step_start)

    if (nstep==step_start+1) then
      if (myrank==0) then
        open(65,file=fname,form='formatted',status='replace')
        if (neqns==5) write(65,805) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5'
        if (neqns==6) write(65,806) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5','EQ6'
        if (neqns==7) write(65,807) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5','EQ6','EQ7'
        close(65)
      endif
    endif

    if (.not.(nstep>step_start.and.mod(nstep,100)==0)) return

    lmax = neqns

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    if (time_integral_mode==0) then
      do l=1,lmax
        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              ii = i-BLK%ista+1
              jj = j-BLK%jsta+1
              kk = k-BLK%ksta+1
              QB(ii,jj,kk,l) = BLK%QCS3(i,j,k,l)-BLK%QCS1(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
    else if (time_integral_mode==1) then
      do l=1,lmax
        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              ii = i-BLK%ista+1
              jj = j-BLK%jsta+1
              kk = k-BLK%ksta+1
              QB(ii,jj,kk,l) = BLK%QCS3(i,j,k,l)-BLK%QCS2(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
    endif

    allocate(l2norm(lmax))
    l2norm = 0.0d0

    call mpisub_sbsp_calc_l2norm(QB,imax,jmax,kmax,lmax,l2norm)

    if (myrank==0) then
      open(66,file=fname,form='formatted',position='append')
      if (neqns==5) write(66,905) nstep,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5)
      if (neqns==6) write(66,906) nstep,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5),l2norm(6)
      if (neqns==7) write(66,907) nstep,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5),l2norm(6),l2norm(7)
      close(66)
    endif

    deallocate(l2norm)

    deallocate(QB)

  end subroutine output_l2norm
!-------------------------------------------------------------------------------
  subroutine output_l2norm_inner_iteration(BLK)
    type(block), intent(in) :: BLK

    character(100) :: fname
    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    integer :: nstep,step_start
    integer :: nstep_it,step_start_it,step_final_it
    integer :: nperiod
    real(8), allocatable, dimension(:) :: l2norm

    700 format('./data/l2norm_of_inner_iteration_in_step_',i7.7,'.txt')

    805 format(6A13)
    806 format(7A13)
    807 format(8A13)

    905 format(I13,5E13.5)
    906 format(I13,6E13.5)
    907 format(I13,7E13.5)

    call mpisub_sbsp_get_myrank_world(myrank)

    call steps_get_step_now(nstep)

    call steps_get_step_start(step_start)

    if (.not.(nstep==step_start.or.mod(nstep,1000)==0)) return

    call steps_it_get_step_now(nstep_it)

    call steps_it_get_step_start(step_start_it)

    call steps_it_get_step_final(step_final_it)

    if (nstep_it==step_start_it+1) then
      if (myrank==0) then
        write(fname,700) nstep
        open(65,file=fname,form='formatted',status='replace')
        if (neqns==5) write(65,805) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5'
        if (neqns==6) write(65,806) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5','EQ6'
        if (neqns==7) write(65,807) 'STEP','EQ1','EQ2','EQ3','EQ4','EQ5','EQ6','EQ7'
        close(65)
      endif
    endif

    nperiod = 1
    if (step_final_it>100 ) nperiod = 10
    if (step_final_it>1000) nperiod = 100

    if (.not.(nstep_it>step_start_it.and.mod(nstep_it,nperiod)==0)) return

    lmax = neqns

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do l=1,lmax
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QB(ii,jj,kk,l) = BLK%RHSS(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    allocate(l2norm(lmax))
    l2norm = 0.0d0

    call mpisub_sbsp_calc_l2norm(QB,imax,jmax,kmax,lmax,l2norm)

    if (myrank==0) then
      write(fname,700) nstep
      open(66,file=fname,form='formatted',position='append')
      if (neqns==5) write(66,905) nstep_it,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5)
      if (neqns==6) write(66,906) nstep_it,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5),l2norm(6)
      if (neqns==7) write(66,907) nstep_it,l2norm(1),l2norm(2),l2norm(3),l2norm(4),l2norm(5),l2norm(6),l2norm(7)
      close(66)
    endif

    deallocate(l2norm)

    deallocate(QB)

  end subroutine output_l2norm_inner_iteration
!-------------------------------------------------------------------------------
end module
