module mdl_implicit_sgs
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_halo
  implicit none

contains

  subroutine lu_sgs(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block
    integer :: imax,jmax,kmax
    integer :: i,j,k
    integer, allocatable, dimension(:,:,:) :: MASK

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

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

!    if (i2dimension>0) then
!      call lu_sgs_2d
!    else
      call lu_sgs_3d
!    endif

    deallocate(MASK)

  contains
  !-----------------------------------------------------------------------------
    subroutine lu_sgs_2d
      integer :: i,j,k,l,m,mm,mmin,mmax,idiv1,idiv2,ista,iend
      integer :: im,jm,ip,jp,ii,jj,kk
      integer :: iista,iiend,jjsta,jjend
      real(8) :: sgn,dt,bb,h
      real(8) :: rh,u,v,hst,hsp,rrh
      real(8) :: sx,sy,sxy,sxy2,c,c2,uc,vp,vp2,beta,ggd
      real(8) :: eig1,eig2,eig3,rada,radb
      real(8) :: qm01,qm02,qm03,qm04,qm05,qm06,qm07,qm08,qm09,qm10,qm11,qm12,dinv
      real(8) :: duna01,duna02,duna03,duna04,duna05,duna06,duna07,duna08,duna09,duna10,duna11,duna12
      real(8) :: dunb01,dunb02,dunb03,dunb04,dunb05,dunb06,dunb07,dunb08,dunb09,dunb10,dunb11,dunb12
      real(8), allocatable, dimension(:,:,:,:) :: DQSTAR

      allocate(DQSTAR(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,neqns))

      DQSTAR(:,:,:,:) = 0.0d0

      do l=1,neqns
        do k=BLK%ksta_wh,BLK%kend_wh
          do j=BLK%jsta_wh,BLK%jend_wh
            do i=BLK%ista_wh,BLK%iend_wh
              BLK%RHSS(i,j,k,l) = BLK%RHSS(i,j,k,l)*dble(MASK(i,j,k))
            enddo
          enddo
        enddo
      enddo

      DQSTAR(:,:,:,:) = BLK%RHSS(:,:,:,:)

      iista = 2
      iiend = BLK%imax_wh-1
      jjsta = 2
      jjend = BLK%jmax_wh-1

      mmin = iista+jjsta
      mmax = iiend+jjend

      !Forward sweep
      do m=mmin,mmax
        idiv1= m-jjend
        idiv2= m-jjsta
        ista = max0(idiv1,iista)
        iend = min0(idiv2,iiend)
        do i=ista,iend
          j = m-i
          k = 1
          im= i-1
          jm= j-1

          if (MASK(i,j,k)/=1) cycle

          sgn = 1.0d0

          dt = BLK%DTLC(i,j,k)

          bb = 1.0d0+1.5d0*dt/physical_time_step
          bb = 1.0d0/bb

          h = bb*dt

          rh= BLK%DNST(i,j,k)
          u = BLK%ULCT(i,j,k)
          v = BLK%VLCT(i,j,k)

          hst = BLK%SHCP(i,j,k)*CPREF
          hsp = BLK%DHPR(i,j,k)

          rrh = 1.0d0/rh

          c  = BLK%CSOS(i,j,k)
          vp = dmin1(BLK%VPRE(i,j,k),c)

          c2 = c*c
          vp2= vp*vp

          beta = (1.0d0-rh*hsp)/(rh*hst)

          !JACOBI A+ (i-1,j,k)
          ii = im
          jj = j
          kk = k

          sx = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          sy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)

          uc  = sx*u+sy*v
          sxy2= sx*sx+sy*sy
          sxy = dsqrt(sxy2)

          ggd = dsqrt(uc*uc*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxy2)

          eig1 = uc
          eig2 = 0.5d0*(uc*(1.0d0+vp2/c2)+ggd)
          eig3 = 0.5d0*(uc*(1.0d0+vp2/c2)-ggd)

          eig1 = dabs(eig1)
          eig2 = dabs(eig2)
          eig3 = dabs(eig3)

          rada = 1.01D0*dmax1(eig1,eig2,eig3)

          qm01 = DQSTAR(ii,jj,kk,1)
          qm02 = DQSTAR(ii,jj,kk,2)
          qm03 = DQSTAR(ii,jj,kk,3)
          qm04 = DQSTAR(ii,jj,kk,4)
          qm05 = DQSTAR(ii,jj,kk,5)
          if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
          if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
          if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
          if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
          if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
          if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
          if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

          duna01 = (sgn*rada+uc           )*qm01 &
                  +(beta*vp2*sx*rh        )*qm02 &
                  +(beta*vp2*sy*rh        )*qm03 &
                  +(0.0d0                 )*qm04 &
                  -(beta*(1.0d0-vp2/c2)*uc)*qm05
          duna02 = (0.0d0                 )*qm01 &
                  +(sgn*rada+uc           )*qm02 &
                  +(0.0d0                 )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sx*rrh                )*qm05
          duna03 = (0.0d0                 )*qm01 &
                  +(0.0d0                 )*qm02 &
                  +(sgn*rada+uc           )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sy*rrh                )*qm05
          duna04 = 0.0d0
          duna05 = (0.0d0                 )*qm01 &
                  +(vp2*sx*rh             )*qm02 &
                  +(vp2*sy*rh             )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sgn*rada+vp2/c2*uc    )*qm05
          if (neqns>=6) duna06 = (sgn*rada+uc)*qm06
          if (neqns>=7) duna07 = (sgn*rada+uc)*qm07
          if (neqns>=8) duna08 = (sgn*rada+uc)*qm08
          if (neqns>=9) duna09 = (sgn*rada+uc)*qm09
          if (neqns>=10)duna10 = (sgn*rada+uc)*qm10
          if (neqns>=11)duna11 = (sgn*rada+uc)*qm11
          if (neqns>=12)duna12 = (sgn*rada+uc)*qm12

          duna01 = 0.5D0*duna01
          duna02 = 0.5D0*duna02
          duna03 = 0.5D0*duna03
          duna04 = 0.5D0*duna04
          duna05 = 0.5D0*duna05
          if (neqns>=6) duna06 = 0.5D0*duna06
          if (neqns>=7) duna07 = 0.5D0*duna07
          if (neqns>=8) duna08 = 0.5D0*duna08
          if (neqns>=9) duna09 = 0.5D0*duna09
          if (neqns>=10)duna10 = 0.5D0*duna10
          if (neqns>=11)duna11 = 0.5D0*duna11
          if (neqns>=12)duna12 = 0.5D0*duna12

          !JACOBI B+ (i,j-1,k)
          ii = i
          jj = jm
          kk = k

          sx = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          sy = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)

          uc  = sx*u+sy*v
          sxy2= sx*sx+sy*sy
          sxy = dsqrt(sxy2)

          ggd = dsqrt(uc*uc*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxy2)

          eig1 = uc
          eig2 = 0.5d0*(uc*(1.0d0+vp2/c2)+ggd)
          eig3 = 0.5d0*(uc*(1.0d0+vp2/c2)-ggd)

          eig1 = dabs(eig1)
          eig2 = dabs(eig2)
          eig3 = dabs(eig3)

          radb = 1.01D0*dmax1(eig1,eig2,eig3)

          qm01 = DQSTAR(ii,jj,kk,1)
          qm02 = DQSTAR(ii,jj,kk,2)
          qm03 = DQSTAR(ii,jj,kk,3)
          qm04 = DQSTAR(ii,jj,kk,4)
          qm05 = DQSTAR(ii,jj,kk,5)
          if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
          if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
          if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
          if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
          if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
          if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
          if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

          dunb01 = (sgn*radb+uc           )*qm01 &
                  +(beta*vp2*sx*rh        )*qm02 &
                  +(beta*vp2*sy*rh        )*qm03 &
                  +(0.0d0                 )*qm04 &
                  -(beta*(1.0d0-vp2/c2)*uc)*qm05
          dunb02 = (0.0d0                 )*qm01 &
                  +(sgn*radb+uc           )*qm02 &
                  +(0.0d0                 )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sx*rrh                )*qm05
          dunb03 = (0.0d0                 )*qm01 &
                  +(0.0d0                 )*qm02 &
                  +(sgn*radb+uc           )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sy*rrh                )*qm05
          dunb04 = 0.0d0
          dunb05 = (0.0d0                 )*qm01 &
                  +(vp2*sx*rh             )*qm02 &
                  +(vp2*sy*rh             )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sgn*radb+vp2*c2*uc    )*qm05
          if (neqns>=6) dunb06 = (sgn*radb+uc)*qm06
          if (neqns>=7) dunb07 = (sgn*radb+uc)*qm07
          if (neqns>=8) dunb08 = (sgn*radb+uc)*qm08
          if (neqns>=9) dunb09 = (sgn*radb+uc)*qm09
          if (neqns>=10)dunb10 = (sgn*radb+uc)*qm10
          if (neqns>=11)dunb11 = (sgn*radb+uc)*qm11
          if (neqns>=12)dunb12 = (sgn*radb+uc)*qm12

          dunb01 = 0.5D0*dunb01
          dunb02 = 0.5D0*dunb02
          dunb03 = 0.5D0*dunb03
          dunb04 = 0.5D0*dunb04
          dunb05 = 0.5D0*dunb05
          if (neqns>=6) dunb06 = 0.5D0*dunb06
          if (neqns>=7) dunb07 = 0.5D0*dunb07
          if (neqns>=8) dunb08 = 0.5D0*dunb08
          if (neqns>=9) dunb09 = 0.5D0*dunb09
          if (neqns>=10)dunb10 = 0.5D0*dunb10
          if (neqns>=11)dunb11 = 0.5D0*dunb11
          if (neqns>=12)dunb12 = 0.5D0*dunb12

          dinv = 1.0D0/(1.0D0+h*(rada+radb))

          DQSTAR(i,j,k,1) = dinv*(BLK%RHSS(i,j,k,1)+h*(duna01+dunb01))
          DQSTAR(i,j,k,2) = dinv*(BLK%RHSS(i,j,k,2)+h*(duna02+dunb02))
          DQSTAR(i,j,k,3) = dinv*(BLK%RHSS(i,j,k,3)+h*(duna03+dunb03))
          DQSTAR(i,j,k,4) = dinv*(BLK%RHSS(i,j,k,4)+h*(duna04+dunb04))
          DQSTAR(i,j,k,5) = dinv*(BLK%RHSS(i,j,k,5)+h*(duna05+dunb05))
          if (neqns>=6)  DQSTAR(i,j,k,6) = dinv*(BLK%RHSS(i,j,k, 6)+h*(duna06+dunb06))
          if (neqns>=7)  DQSTAR(i,j,k,7) = dinv*(BLK%RHSS(i,j,k, 7)+h*(duna07+dunb07))
          if (neqns>=8)  DQSTAR(i,j,k,8) = dinv*(BLK%RHSS(i,j,k, 8)+h*(duna08+dunb08))
          if (neqns>=9)  DQSTAR(i,j,k,9) = dinv*(BLK%RHSS(i,j,k, 9)+h*(duna09+dunb09))
          if (neqns>=10) DQSTAR(i,j,k,10)= dinv*(BLK%RHSS(i,j,k,10)+h*(duna10+dunb10))
          if (neqns>=11) DQSTAR(i,j,k,11)= dinv*(BLK%RHSS(i,j,k,11)+h*(duna11+dunb11))
          if (neqns>=12) DQSTAR(i,j,k,12)= dinv*(BLK%RHSS(i,j,k,12)+h*(duna12+dunb12))
        enddo
      enddo


      !Backward sweep
      do mm=mmin,mmax
        m = mmax-mm+mmin
        idiv1= m-jjend
        idiv2= m-jjsta
        ista = max0(idiv1,iista)
        iend = min0(idiv2,iiend)
        do i=ista,iend
          j = m-i
          k = 1
          ip= i+1
          jp= j+1

          if (MASK(i,j,k)/=1) cycle

          sgn = -1.0d0

          dt = BLK%DTLC(i,j,k)

          bb = 1.0d0+1.5d0*dt/physical_time_step
          bb = 1.0d0/bb

          h = bb*dt

          rh= BLK%DNST(i,j,k)
          u = BLK%ULCT(i,j,k)
          v = BLK%VLCT(i,j,k)

          hst = BLK%SHCP(i,j,k)*CPREF
          hsp = BLK%DHPR(i,j,k)

          rrh = 1.0d0/rh

          c  = BLK%CSOS(i,j,k)
          vp = dmin1(BLK%VPRE(i,j,k),c)

          c2 = c*c
          vp2= vp*vp

          beta = (1.0d0-rh*hsp)/(rh*hst)

          !JACOBI A- (i+1,j,k)
          ii = ip
          jj = j
          kk = k

          sx = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
          sy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)

          uc = sx*u+sy*v
          sxy2= sx*sx+sy*sy
          sxy = dsqrt(sxy2)

          ggd = dsqrt(uc*uc*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxy2)

          eig1 = uc
          eig2 = 0.5d0*(uc*(1.0d0+vp2/c2)+ggd)
          eig3 = 0.5d0*(uc*(1.0d0+vp2/c2)-ggd)

          eig1 = dabs(eig1)
          eig2 = dabs(eig2)
          eig3 = dabs(eig3)

          rada = 1.01D0*dmax1(eig1,eig2,eig3)

          qm01 = BLK%RHSS(ii,jj,kk,1)
          qm02 = BLK%RHSS(ii,jj,kk,2)
          qm03 = BLK%RHSS(ii,jj,kk,3)
          qm04 = BLK%RHSS(ii,jj,kk,4)
          qm05 = BLK%RHSS(ii,jj,kk,5)
          if (neqns>=6) qm06 = BLK%RHSS(ii,jj,kk,6)
          if (neqns>=7) qm07 = BLK%RHSS(ii,jj,kk,7)
          if (neqns>=8) qm08 = BLK%RHSS(ii,jj,kk,8)
          if (neqns>=9) qm09 = BLK%RHSS(ii,jj,kk,9)
          if (neqns>=10)qm10 = BLK%RHSS(ii,jj,kk,10)
          if (neqns>=11)qm11 = BLK%RHSS(ii,jj,kk,11)
          if (neqns>=12)qm12 = BLK%RHSS(ii,jj,kk,12)

          duna01 = (sgn*rada+uc           )*qm01 &
                  +(beta*vp2*sx*rh        )*qm02 &
                  +(beta*vp2*sy*rh        )*qm03 &
                  +(0.0d0                 )*qm04 &
                  -(beta*(1.0d0-vp2/c2)*uc)*qm05
          duna02 = (0.0d0                 )*qm01 &
                  +(sgn*rada+uc           )*qm02 &
                  +(0.0d0                 )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sx*rrh                )*qm05
          duna03 = (0.0d0                 )*qm01 &
                  +(0.0d0                 )*qm02 &
                  +(sgn*rada+uc           )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sy*rrh                )*qm05
          duna04 = 0.0d0
          duna05 = (0.0d0                 )*qm01 &
                  +(vp2*sx*rh             )*qm02 &
                  +(vp2*sy*rh             )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sgn*rada+vp2/c2*uc    )*qm05
          if (neqns>=6) duna06 = (sgn*rada+uc)*qm06
          if (neqns>=7) duna07 = (sgn*rada+uc)*qm07
          if (neqns>=8) duna08 = (sgn*rada+uc)*qm08
          if (neqns>=9) duna09 = (sgn*rada+uc)*qm09
          if (neqns>=10)duna10 = (sgn*rada+uc)*qm10
          if (neqns>=11)duna11 = (sgn*rada+uc)*qm11
          if (neqns>=12)duna12 = (sgn*rada+uc)*qm12

          duna01 = 0.5D0*duna01
          duna02 = 0.5D0*duna02
          duna03 = 0.5D0*duna03
          duna04 = 0.5D0*duna04
          duna05 = 0.5D0*duna05
          if (neqns>=6) duna06 = 0.5D0*duna06
          if (neqns>=7) duna07 = 0.5D0*duna07
          if (neqns>=8) duna08 = 0.5D0*duna08
          if (neqns>=9) duna09 = 0.5D0*duna09
          if (neqns>=10)duna10 = 0.5D0*duna10
          if (neqns>=11)duna11 = 0.5D0*duna11
          if (neqns>=12)duna12 = 0.5D0*duna12

          !JACOBI B- (i,j+1,k)
          ii = i
          jj = jp
          kk = k

          sx = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
          sy = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)

          uc  = sx*u+sy*v
          sxy2= sx*sx+sy*sy
          sxy = dsqrt(sxy2)

          ggd = dsqrt(uc*uc*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxy2)

          eig1 = uc
          eig2 = 0.5d0*(uc*(1.0d0+vp2/c2)+ggd)
          eig3 = 0.5d0*(uc*(1.0d0+vp2/c2)-ggd)

          eig1 = dabs(eig1)
          eig2 = dabs(eig2)
          eig3 = dabs(eig3)

          radb = 1.01D0*dmax1(eig1,eig2,eig3)

          qm01 = BLK%RHSS(ii,jj,kk,1)
          qm02 = BLK%RHSS(ii,jj,kk,2)
          qm03 = BLK%RHSS(ii,jj,kk,3)
          qm04 = BLK%RHSS(ii,jj,kk,4)
          qm05 = BLK%RHSS(ii,jj,kk,5)
          if (neqns>=6) qm06 = BLK%RHSS(ii,jj,kk,6)
          if (neqns>=7) qm07 = BLK%RHSS(ii,jj,kk,7)
          if (neqns>=8) qm08 = BLK%RHSS(ii,jj,kk,8)
          if (neqns>=9) qm09 = BLK%RHSS(ii,jj,kk,9)
          if (neqns>=10)qm10 = BLK%RHSS(ii,jj,kk,10)
          if (neqns>=11)qm11 = BLK%RHSS(ii,jj,kk,11)
          if (neqns>=12)qm12 = BLK%RHSS(ii,jj,kk,12)

          dunb01 = (sgn*radb+uc           )*qm01 &
                  +(beta*vp2*sx*rh        )*qm02 &
                  +(beta*vp2*sy*rh        )*qm03 &
                  +(0.0d0                 )*qm04 &
                  -(beta*(1.0d0-vp2/c2)*uc)*qm05
          dunb02 = (0.0d0                 )*qm01 &
                  +(sgn*radb+uc           )*qm02 &
                  +(0.0d0                 )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sx*rrh                )*qm05
          dunb03 = (0.0d0                 )*qm01 &
                  +(0.0d0                 )*qm02 &
                  +(sgn*radb+uc           )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sy*rrh                )*qm05
          dunb04 = 0.0d0
          dunb05 = (0.0d0                 )*qm01 &
                  +(vp2*sx*rh             )*qm02 &
                  +(vp2*sy*rh             )*qm03 &
                  +(0.0d0                 )*qm04 &
                  +(sgn*radb+vp2/c2*uc    )*qm05
          if (neqns>=6) dunb06 = (sgn*radb+uc)*qm06
          if (neqns>=7) dunb07 = (sgn*radb+uc)*qm07
          if (neqns>=8) dunb08 = (sgn*radb+uc)*qm08
          if (neqns>=9) dunb09 = (sgn*radb+uc)*qm09
          if (neqns>=10)dunb10 = (sgn*radb+uc)*qm10
          if (neqns>=11)dunb11 = (sgn*radb+uc)*qm11
          if (neqns>=12)dunb12 = (sgn*radb+uc)*qm12

          dunb01 = 0.5D0*dunb01
          dunb02 = 0.5D0*dunb02
          dunb03 = 0.5D0*dunb03
          dunb04 = 0.5D0*dunb04
          dunb05 = 0.5D0*dunb05
          if (neqns>=6) dunb06 = 0.5D0*dunb06
          if (neqns>=7) dunb07 = 0.5D0*dunb07
          if (neqns>=8) dunb08 = 0.5D0*dunb08
          if (neqns>=9) dunb09 = 0.5D0*dunb09
          if (neqns>=10)dunb10 = 0.5D0*dunb10
          if (neqns>=11)dunb11 = 0.5D0*dunb11
          if (neqns>=12)dunb12 = 0.5D0*dunb12

          dinv = 1.0D0/(1.0D0+h*(rada+radb))

          BLK%RHSS(i,j,k,1) = DQSTAR(i,j,k,1)-dinv*h*(duna01+dunb01)
          BLK%RHSS(i,j,k,2) = DQSTAR(i,j,k,2)-dinv*h*(duna02+dunb02)
          BLK%RHSS(i,j,k,3) = DQSTAR(i,j,k,3)-dinv*h*(duna03+dunb03)
          BLK%RHSS(i,j,k,4) = DQSTAR(i,j,k,4)-dinv*h*(duna04+dunb04)
          BLK%RHSS(i,j,k,5) = DQSTAR(i,j,k,5)-dinv*h*(duna05+dunb05)
          if (neqns>=6) BLK%RHSS(i,j,k,6) = DQSTAR(i,j,k, 6)-dinv*h*(duna06+dunb06)
          if (neqns>=7) BLK%RHSS(i,j,k,7) = DQSTAR(i,j,k, 7)-dinv*h*(duna07+dunb07)
          if (neqns>=8) BLK%RHSS(i,j,k,8) = DQSTAR(i,j,k, 8)-dinv*h*(duna08+dunb08)
          if (neqns>=9) BLK%RHSS(i,j,k,9) = DQSTAR(i,j,k, 9)-dinv*h*(duna09+dunb09)
          if (neqns>=10)BLK%RHSS(i,j,k,10)= DQSTAR(i,j,k,10)-dinv*h*(duna10+dunb10)
          if (neqns>=11)BLK%RHSS(i,j,k,11)= DQSTAR(i,j,k,11)-dinv*h*(duna11+dunb11)
          if (neqns>=12)BLK%RHSS(i,j,k,12)= DQSTAR(i,j,k,12)-dinv*h*(duna12+dunb12)
        enddo
      enddo

      deallocate(DQSTAR)

    end subroutine lu_sgs_2d
  !-----------------------------------------------------------------------------
    subroutine lu_sgs_3d
      integer :: i,j,k,l,m,mm,mmin,mmax,idiv1,idiv2,ista,iend
      integer :: im,jm,km,ip,jp,kp,ii,jj,kk
      integer :: iista,iiend,jjsta,jjend,kksta,kkend
      real(8) :: sgn,dt,bb,h
      real(8) :: rh,u,v,w,hst,hsp,rrh,vp,vp2,beta,ggd
      real(8) :: sx,sy,sz,sxyz,sxyz2,c,c2,uc
      real(8) :: eig1,eig2,eig3,rada,radb,radc
      real(8) :: qm01,qm02,qm03,qm04,qm05,qm06,qm07,qm08,qm09,qm10,qm11,qm12,dinv
      real(8) :: duna01,duna02,duna03,duna04,duna05,duna06,duna07,duna08,duna09,duna10,duna11,duna12
      real(8) :: dunb01,dunb02,dunb03,dunb04,dunb05,dunb06,dunb07,dunb08,dunb09,dunb10,dunb11,dunb12
      real(8) :: dunc01,dunc02,dunc03,dunc04,dunc05,dunc06,dunc07,dunc08,dunc09,dunc10,dunc11,dunc12
      real(8), allocatable, dimension(:,:,:,:) :: DQSTAR

      real(8) :: pi
      real(8) :: nx,ny,nz,ox,oy,oz,rx,ry,rz
      real(8) :: urot,vrot,wrot,ac

      allocate(DQSTAR(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,neqns))

      DQSTAR(:,:,:,:) = 0.0d0

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

      if (i2dimension>0) then
        do k=BLK%ksta_wh,BLK%kend_wh
          do j=BLK%jsta_wh,BLK%jend_wh
            do i=BLK%ista_wh,BLK%iend_wh
              BLK%RHSS(i,j,k,4) = 0.0d0
            enddo
          enddo
        enddo
      endif

      do l=1,neqns
        do k=BLK%ksta_wh,BLK%kend_wh
          do j=BLK%jsta_wh,BLK%jend_wh
            do i=BLK%ista_wh,BLK%iend_wh
              BLK%RHSS(i,j,k,l) = BLK%RHSS(i,j,k,l)*dble(MASK(i,j,k))
            enddo
          enddo
        enddo
      enddo

      DQSTAR(:,:,:,:) = BLK%RHSS(:,:,:,:)

      iista = 2
      iiend = BLK%imax_wh-1
      jjsta = 2
      jjend = BLK%jmax_wh-1
      kksta = 2
      kkend = BLK%kmax_wh-1
      if (i2dimension>0) then
        kksta = 1
        kkend = 1
      endif

      mmin = iista+jjsta
      mmax = iiend+jjend

      !Forward sweep
      do k=kksta,kkend
        do m=mmin,mmax
          idiv1= m-jjend
          idiv2= m-jjsta
          ista = max0(idiv1,iista)
          iend = min0(idiv2,iiend)
          do i=ista,iend
            j = m-i
            im= i-1
            jm= j-1
            km= k-1

            if (MASK(i,j,k)/=1) cycle

            sgn = 1.0d0

            dt = BLK%DTLC(i,j,k)

            bb = 1.0d0+1.5d0*dt/physical_time_step
            bb = 1.0d0/bb

            h = bb*dt

            rh= BLK%DNST(i,j,k)
            u = BLK%ULCT(i,j,k)
            v = BLK%VLCT(i,j,k)
            w = BLK%WLCT(i,j,k)

            hst = BLK%SHCP(i,j,k)*CPREF
            hsp = BLK%DHPR(i,j,k)

            rrh = 1.0d0/rh

            c  = BLK%CSOS(i,j,k)
            vp = dmin1(BLK%VPRE(i,j,k),c)

            c2 = c*c
            vp2= vp*vp

            beta = (1.0d0-rh*hsp)/(rh*hst)

            rx = BLK%XCRD(i,j,k)-frame_origin(1)
            ry = BLK%YCRD(i,j,k)-frame_origin(2)
            rz = BLK%ZCRD(i,j,k)-frame_origin(3)

            urot = oy*rz-oz*ry
            vrot = oz*rx-ox*rz
            wrot = ox*ry-oy*rx

            !JACOBI A+ (i-1,j,k)
            ii = im
            jj = j
            kk = k

            sx = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
            sy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
            sz = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)

            uc = sx*u   +sy*v   +sz*w
            ac = sx*urot+sy*vrot+sz*wrot

            sxyz2= sx*sx+sy*sy+sz*sz
            sxyz = dsqrt(sxyz2)

            ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

            eig1 = uc-ac
            eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
            eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

            eig1 = dabs(eig1)
            eig2 = dabs(eig2)
            eig3 = dabs(eig3)

            rada = 1.01D0*dmax1(eig1,eig2,eig3)

            qm01 = DQSTAR(ii,jj,kk,1)
            qm02 = DQSTAR(ii,jj,kk,2)
            qm03 = DQSTAR(ii,jj,kk,3)
            qm04 = DQSTAR(ii,jj,kk,4)
            qm05 = DQSTAR(ii,jj,kk,5)
            if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
            if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
            if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
            if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
            if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
            if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
            if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

            duna01 = (sgn*rada+uc-ac             )*qm01 &
                    +(beta*vp2*sx*rh             )*qm02 &
                    +(beta*vp2*sy*rh             )*qm03 &
                    +(beta*vp2*sz*rh             )*qm04 &
                    -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
            duna02 = (0.0d0                      )*qm01 &
                    +(sgn*rada+uc-ac             )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sx*rrh                     )*qm05
            duna03 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(sgn*rada+uc-ac             )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sy*rrh                     )*qm05
            duna04 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(sgn*rada+uc-ac             )*qm04 &
                    +(sz*rrh                     )*qm05
            duna05 = (0.0d0                      )*qm01 &
                    +(vp2*sx*rh                  )*qm02 &
                    +(vp2*sy*rh                  )*qm03 &
                    +(vp2*sz*rh                  )*qm04 &
                    +(sgn*rada+vp2/c2*(uc-ac)    )*qm05
            if (neqns>=6) duna06 = (sgn*rada+uc-ac)*qm06
            if (neqns>=7) duna07 = (sgn*rada+uc-ac)*qm07
            if (neqns>=8) duna08 = (sgn*rada+uc-ac)*qm08
            if (neqns>=9) duna09 = (sgn*rada+uc-ac)*qm09
            if (neqns>=10)duna10 = (sgn*rada+uc-ac)*qm10
            if (neqns>=11)duna11 = (sgn*rada+uc-ac)*qm11
            if (neqns>=12)duna12 = (sgn*rada+uc-ac)*qm12

            duna01 = 0.5D0*duna01
            duna02 = 0.5D0*duna02
            duna03 = 0.5D0*duna03
            duna04 = 0.5D0*duna04
            duna05 = 0.5D0*duna05
            if (neqns>=6) duna06 = 0.5D0*duna06
            if (neqns>=7) duna07 = 0.5D0*duna07
            if (neqns>=8) duna08 = 0.5D0*duna08
            if (neqns>=9) duna09 = 0.5D0*duna09
            if (neqns>=10)duna10 = 0.5D0*duna10
            if (neqns>=11)duna11 = 0.5D0*duna11
            if (neqns>=12)duna12 = 0.5D0*duna12

            !JACOBI B+ (i,j-1,k)
            ii = i
            jj = jm
            kk = k

            sx = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
            sy = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
            sz = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)

            uc = sx*u   +sy*v   +sz*w
            ac = sx*urot+sy*vrot+sz*wrot

            sxyz2= sx*sx+sy*sy+sz*sz
            sxyz = dsqrt(sxyz2)

            ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

            eig1 = uc-ac
            eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
            eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

            eig1 = dabs(eig1)
            eig2 = dabs(eig2)
            eig3 = dabs(eig3)

            radb = 1.01D0*dmax1(eig1,eig2,eig3)

            qm01 = DQSTAR(ii,jj,kk,1)
            qm02 = DQSTAR(ii,jj,kk,2)
            qm03 = DQSTAR(ii,jj,kk,3)
            qm04 = DQSTAR(ii,jj,kk,4)
            qm05 = DQSTAR(ii,jj,kk,5)
            if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
            if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
            if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
            if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
            if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
            if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
            if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

            dunb01 = (sgn*radb+uc-ac             )*qm01 &
                    +(beta*vp2*sx*rh             )*qm02 &
                    +(beta*vp2*sy*rh             )*qm03 &
                    +(beta*vp2*sz*rh             )*qm04 &
                    -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
            dunb02 = (0.0d0                      )*qm01 &
                    +(sgn*radb+uc-ac             )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sx*rrh                     )*qm05
            dunb03 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(sgn*radb+uc-ac             )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sy*rrh                     )*qm05
            dunb04 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(sgn*radb+uc-ac             )*qm04 &
                    +(sz*rrh                     )*qm05
            dunb05 = (0.0d0                      )*qm01 &
                    +(vp2*sx*rh                  )*qm02 &
                    +(vp2*sy*rh                  )*qm03 &
                    +(vp2*sz*rh                  )*qm04 &
                    +(sgn*radb+vp2/c2*(uc-ac)    )*qm05
            if (neqns>=6) dunb06 = (sgn*radb+uc-ac)*qm06
            if (neqns>=7) dunb07 = (sgn*radb+uc-ac)*qm07
            if (neqns>=8) dunb08 = (sgn*radb+uc-ac)*qm08
            if (neqns>=9) dunb09 = (sgn*radb+uc-ac)*qm09
            if (neqns>=10)dunb10 = (sgn*radb+uc-ac)*qm10
            if (neqns>=11)dunb11 = (sgn*radb+uc-ac)*qm11
            if (neqns>=12)dunb12 = (sgn*radb+uc-ac)*qm12

            dunb01 = 0.5D0*dunb01
            dunb02 = 0.5D0*dunb02
            dunb03 = 0.5D0*dunb03
            dunb04 = 0.5D0*dunb04
            dunb05 = 0.5D0*dunb05
            if (neqns>=6) dunb06 = 0.5D0*dunb06
            if (neqns>=7) dunb07 = 0.5D0*dunb07
            if (neqns>=8) dunb08 = 0.5D0*dunb08
            if (neqns>=9) dunb09 = 0.5D0*dunb09
            if (neqns>=10)dunb10 = 0.5D0*dunb10
            if (neqns>=11)dunb11 = 0.5D0*dunb11
            if (neqns>=12)dunb12 = 0.5D0*dunb12

            !JACOBI C+ (i,j,k-1)
            if (i2dimension>0) then
              dunc01 = 0.0D0
              dunc02 = 0.0D0
              dunc03 = 0.0D0
              dunc04 = 0.0D0
              dunc05 = 0.0D0
              if (neqns>=6) dunc06 = 0.0D0
              if (neqns>=7) dunc07 = 0.0D0
              if (neqns>=8) dunc08 = 0.0D0
              if (neqns>=9) dunc09 = 0.0D0
              if (neqns>=10)dunc10 = 0.0D0
              if (neqns>=11)dunc11 = 0.0D0
              if (neqns>=12)dunc12 = 0.0D0
            else
              ii = i
              jj = j
              kk = km

              sx = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
              sy = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
              sz = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)

              uc = sx*u   +sy*v   +sz*w
              ac = sx*urot+sy*vrot+sz*wrot

              sxyz2= sx*sx+sy*sy+sz*sz
              sxyz = dsqrt(sxyz2)

              ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

              eig1 = uc-ac
              eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
              eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

              eig1 = dabs(eig1)
              eig2 = dabs(eig2)
              eig3 = dabs(eig3)

              radc = 1.01D0*dmax1(eig1,eig2,eig3)

              qm01 = DQSTAR(ii,jj,kk,1)
              qm02 = DQSTAR(ii,jj,kk,2)
              qm03 = DQSTAR(ii,jj,kk,3)
              qm04 = DQSTAR(ii,jj,kk,4)
              qm05 = DQSTAR(ii,jj,kk,5)
              if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
              if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
              if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
              if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
              if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
              if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
              if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

              dunc01 = (sgn*radc+uc-ac             )*qm01 &
                      +(beta*vp2*sx*rh             )*qm02 &
                      +(beta*vp2*sy*rh             )*qm03 &
                      +(beta*vp2*sz*rh             )*qm04 &
                      -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
              dunc02 = (0.0d0                      )*qm01 &
                      +(sgn*radc+uc-ac             )*qm02 &
                      +(0.0d0                      )*qm03 &
                      +(0.0d0                      )*qm04 &
                      +(sx*rrh                     )*qm05
              dunc03 = (0.0d0                      )*qm01 &
                      +(0.0d0                      )*qm02 &
                      +(sgn*radc+uc-ac             )*qm03 &
                      +(0.0d0                      )*qm04 &
                      +(sy*rrh                     )*qm05
              dunc04 = (0.0d0                      )*qm01 &
                      +(0.0d0                      )*qm02 &
                      +(0.0d0                      )*qm03 &
                      +(sgn*radc+uc-ac             )*qm04 &
                      +(sz*rrh                     )*qm05
              dunc05 = (0.0d0                      )*qm01 &
                      +(vp2*sx*rh                  )*qm02 &
                      +(vp2*sy*rh                  )*qm03 &
                      +(vp2*sz*rh                  )*qm04 &
                      +(sgn*radc+vp2/c2*(uc-ac)    )*qm05
              if (neqns>=6) dunc06 = (sgn*radc+uc-ac)*qm06
              if (neqns>=7) dunc07 = (sgn*radc+uc-ac)*qm07
              if (neqns>=8) dunc08 = (sgn*radc+uc-ac)*qm08
              if (neqns>=9) dunc09 = (sgn*radc+uc-ac)*qm09
              if (neqns>=10)dunc10 = (sgn*radc+uc-ac)*qm10
              if (neqns>=11)dunc11 = (sgn*radc+uc-ac)*qm11
              if (neqns>=12)dunc12 = (sgn*radc+uc-ac)*qm12

              dunc01 = 0.5D0*dunc01
              dunc02 = 0.5D0*dunc02
              dunc03 = 0.5D0*dunc03
              dunc04 = 0.5D0*dunc04
              dunc05 = 0.5D0*dunc05
              if (neqns>=6) dunc06 = 0.5D0*dunc06
              if (neqns>=7) dunc07 = 0.5D0*dunc07
              if (neqns>=8) dunc08 = 0.5D0*dunc08
              if (neqns>=9) dunc09 = 0.5D0*dunc09
              if (neqns>=10)dunc10 = 0.5D0*dunc10
              if (neqns>=11)dunc11 = 0.5D0*dunc11
              if (neqns>=12)dunc12 = 0.5D0*dunc12
            endif

            dinv = 1.0D0/(1.0D0+h*(rada+radb+radc))

            DQSTAR(i,j,k,1) = dinv*(BLK%RHSS(i,j,k,1)+h*(duna01+dunb01+dunc01))
            DQSTAR(i,j,k,2) = dinv*(BLK%RHSS(i,j,k,2)+h*(duna02+dunb02+dunc02))
            DQSTAR(i,j,k,3) = dinv*(BLK%RHSS(i,j,k,3)+h*(duna03+dunb03+dunc03))
            DQSTAR(i,j,k,4) = dinv*(BLK%RHSS(i,j,k,4)+h*(duna04+dunb04+dunc04))
            DQSTAR(i,j,k,5) = dinv*(BLK%RHSS(i,j,k,5)+h*(duna05+dunb05+dunc05))
            if (neqns>=6) DQSTAR(i,j,k,6) = dinv*(BLK%RHSS(i,j,k, 6)+h*(duna06+dunb06+dunc06))
            if (neqns>=7) DQSTAR(i,j,k,7) = dinv*(BLK%RHSS(i,j,k, 7)+h*(duna07+dunb07+dunc07))
            if (neqns>=8) DQSTAR(i,j,k,8) = dinv*(BLK%RHSS(i,j,k, 8)+h*(duna08+dunb08+dunc08))
            if (neqns>=9) DQSTAR(i,j,k,9) = dinv*(BLK%RHSS(i,j,k, 9)+h*(duna09+dunb09+dunc09))
            if (neqns>=10)DQSTAR(i,j,k,10)= dinv*(BLK%RHSS(i,j,k,10)+h*(duna10+dunb10+dunc10))
            if (neqns>=11)DQSTAR(i,j,k,11)= dinv*(BLK%RHSS(i,j,k,11)+h*(duna11+dunb11+dunc11))
            if (neqns>=12)DQSTAR(i,j,k,12)= dinv*(BLK%RHSS(i,j,k,12)+h*(duna12+dunb12+dunc12))
          enddo
        enddo
      enddo


    !Backward sweep
      do k=kkend,kksta,-1
        do mm=mmin,mmax
          m = mmax-mm+mmin
          idiv1= m-jjend
          idiv2= m-jjsta
          ista = max0(idiv1,iista)
          iend = min0(idiv2,iiend)
          do i=ista,iend
            j = m-i
            ip= i+1
            jp= j+1
            kp= k+1

            if (MASK(i,j,k)/=1) cycle

            sgn = -1.0d0

            dt = BLK%DTLC(i,j,k)

            bb = 1.0d0+1.5d0*dt/physical_time_step
            bb = 1.0d0/bb

            h = bb*dt

            rh= BLK%DNST(i,j,k)
            u = BLK%ULCT(i,j,k)
            v = BLK%VLCT(i,j,k)
            w = BLK%WLCT(i,j,k)

            hst = BLK%SHCP(i,j,k)*CPREF
            hsp = BLK%DHPR(i,j,k)

            rrh = 1.0d0/rh

            c  = BLK%CSOS(i,j,k)
            vp = dmin1(BLK%VPRE(i,j,k),c)

            c2 = c*c
            vp2= vp*vp

            beta = (1.0d0-rh*hsp)/(rh*hst)

            rx = BLK%XCRD(i,j,k)-frame_origin(1)
            ry = BLK%YCRD(i,j,k)-frame_origin(2)
            rz = BLK%ZCRD(i,j,k)-frame_origin(3)

            urot = oy*rz-oz*ry
            vrot = oz*rx-ox*rz
            wrot = ox*ry-oy*rx

            !JACOBI A- (i+1,j,k)
            ii = ip
            jj = j
            kk = k

            sx = BLK%DXSX(i,j,k)/BLK%RJCB(i,j,k)
            sy = BLK%DXSY(i,j,k)/BLK%RJCB(i,j,k)
            sz = BLK%DXSZ(i,j,k)/BLK%RJCB(i,j,k)

            uc = sx*u   +sy*v   +sz*w
            ac = sx*urot+sy*vrot+sz*wrot

            sxyz2= sx*sx+sy*sy+sz*sz
            sxyz = dsqrt(sxyz2)

            ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

            eig1 = uc-ac
            eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
            eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

            eig1 = dabs(eig1)
            eig2 = dabs(eig2)
            eig3 = dabs(eig3)

            rada = 1.01D0*dmax1(eig1,eig2,eig3)

            qm01 = BLK%RHSS(ii,jj,kk,1)
            qm02 = BLK%RHSS(ii,jj,kk,2)
            qm03 = BLK%RHSS(ii,jj,kk,3)
            qm04 = BLK%RHSS(ii,jj,kk,4)
            qm05 = BLK%RHSS(ii,jj,kk,5)
            if (neqns>=6) qm06 = BLK%RHSS(ii,jj,kk,6)
            if (neqns>=7) qm07 = BLK%RHSS(ii,jj,kk,7)
            if (neqns>=8) qm08 = BLK%RHSS(ii,jj,kk,8)
            if (neqns>=9) qm09 = BLK%RHSS(ii,jj,kk,9)
            if (neqns>=10)qm10 = BLK%RHSS(ii,jj,kk,10)
            if (neqns>=11)qm11 = BLK%RHSS(ii,jj,kk,11)
            if (neqns>=12)qm12 = BLK%RHSS(ii,jj,kk,12)

            duna01 = (sgn*rada+uc-ac             )*qm01 &
                    +(beta*vp2*sx*rh             )*qm02 &
                    +(beta*vp2*sy*rh             )*qm03 &
                    +(beta*vp2*sz*rh             )*qm04 &
                    -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
            duna02 = (0.0d0                      )*qm01 &
                    +(sgn*rada+uc-ac             )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sx*rrh                     )*qm05
            duna03 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(sgn*rada+uc-ac             )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sy*rrh                     )*qm05
            duna04 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(sgn*rada+uc-ac             )*qm04 &
                    +(sz*rrh                     )*qm05
            duna05 = (0.0d0                      )*qm01 &
                    +(vp2*sx*rh                  )*qm02 &
                    +(vp2*sy*rh                  )*qm03 &
                    +(vp2*sz*rh                  )*qm04 &
                    +(sgn*rada+vp2/c2*(uc-ac)    )*qm05
            if (neqns>=6) duna06 = (sgn*rada+uc-ac)*qm06
            if (neqns>=7) duna07 = (sgn*rada+uc-ac)*qm07
            if (neqns>=8) duna08 = (sgn*rada+uc-ac)*qm08
            if (neqns>=9) duna09 = (sgn*rada+uc-ac)*qm09
            if (neqns>=10)duna10 = (sgn*rada+uc-ac)*qm10
            if (neqns>=11)duna11 = (sgn*rada+uc-ac)*qm11
            if (neqns>=12)duna12 = (sgn*rada+uc-ac)*qm12

            duna01 = 0.5D0*duna01
            duna02 = 0.5D0*duna02
            duna03 = 0.5D0*duna03
            duna04 = 0.5D0*duna04
            duna05 = 0.5D0*duna05
            if (neqns>=6) duna06 = 0.5D0*duna06
            if (neqns>=7) duna07 = 0.5D0*duna07
            if (neqns>=8) duna08 = 0.5D0*duna08
            if (neqns>=9) duna09 = 0.5D0*duna09
            if (neqns>=10)duna10 = 0.5D0*duna10
            if (neqns>=11)duna11 = 0.5D0*duna11
            if (neqns>=12)duna12 = 0.5D0*duna12

            !JACOBI B- (i,j+1,k)
            ii = i
            jj = jp
            kk = k

            sx = BLK%DETX(i,j,k)/BLK%RJCB(i,j,k)
            sy = BLK%DETY(i,j,k)/BLK%RJCB(i,j,k)
            sz = BLK%DETZ(i,j,k)/BLK%RJCB(i,j,k)

            uc = sx*u   +sy*v   +sz*w
            ac = sx*urot+sy*vrot+sz*wrot

            sxyz2= sx*sx+sy*sy+sz*sz
            sxyz = dsqrt(sxyz2)

            ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

            eig1 = uc-ac
            eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
            eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

            eig1 = dabs(eig1)
            eig2 = dabs(eig2)
            eig3 = dabs(eig3)

            radb = 1.01D0*dmax1(eig1,eig2,eig3)

            qm01 = BLK%RHSS(ii,jj,kk,1)
            qm02 = BLK%RHSS(ii,jj,kk,2)
            qm03 = BLK%RHSS(ii,jj,kk,3)
            qm04 = BLK%RHSS(ii,jj,kk,4)
            qm05 = BLK%RHSS(ii,jj,kk,5)
            if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
            if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
            if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
            if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
            if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
            if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
            if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

            dunb01 = (sgn*radb+uc-ac             )*qm01 &
                    +(beta*vp2*sx*rh             )*qm02 &
                    +(beta*vp2*sy*rh             )*qm03 &
                    +(beta*vp2*sz*rh             )*qm04 &
                    -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
            dunb02 = (0.0d0                      )*qm01 &
                    +(sgn*radb+uc-ac             )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sx*rrh                     )*qm05
            dunb03 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(sgn*radb+uc-ac             )*qm03 &
                    +(0.0d0                      )*qm04 &
                    +(sy*rrh                     )*qm05
            dunb04 = (0.0d0                      )*qm01 &
                    +(0.0d0                      )*qm02 &
                    +(0.0d0                      )*qm03 &
                    +(sgn*radb+uc-ac             )*qm04 &
                    +(sz*rrh                     )*qm05
            dunb05 = (0.0d0                      )*qm01 &
                    +(vp2*sx*rh                  )*qm02 &
                    +(vp2*sy*rh                  )*qm03 &
                    +(vp2*sz*rh                  )*qm04 &
                    +(sgn*radb+vp2/c2*(uc-ac)    )*qm05
            if (neqns>=6) dunb06 = (sgn*radb+uc-ac)*qm06
            if (neqns>=7) dunb07 = (sgn*radb+uc-ac)*qm07
            if (neqns>=8) dunb08 = (sgn*radb+uc-ac)*qm08
            if (neqns>=9) dunb09 = (sgn*radb+uc-ac)*qm09
            if (neqns>=10)dunb10 = (sgn*radb+uc-ac)*qm10
            if (neqns>=11)dunb11 = (sgn*radb+uc-ac)*qm11
            if (neqns>=12)dunb12 = (sgn*radb+uc-ac)*qm12

            dunb01 = 0.5D0*dunb01
            dunb02 = 0.5D0*dunb02
            dunb03 = 0.5D0*dunb03
            dunb04 = 0.5D0*dunb04
            dunb05 = 0.5D0*dunb05
            if (neqns>=6) dunb06 = 0.5D0*dunb06
            if (neqns>=7) dunb07 = 0.5D0*dunb07
            if (neqns>=8) dunb08 = 0.5D0*dunb08
            if (neqns>=9) dunb09 = 0.5D0*dunb09
            if (neqns>=10)dunb10 = 0.5D0*dunb10
            if (neqns>=11)dunb11 = 0.5D0*dunb11
            if (neqns>=12)dunb12 = 0.5D0*dunb12

            !JACOBI C- (i,j,k+1)
            if (i2dimension>0) then
              dunc01 = 0.0D0
              dunc02 = 0.0D0
              dunc03 = 0.0D0
              dunc04 = 0.0D0
              dunc05 = 0.0D0
              if (neqns>=6) dunc06 = 0.0D0
              if (neqns>=7) dunc07 = 0.0D0
              if (neqns>=8) dunc08 = 0.0D0
              if (neqns>=9) dunc09 = 0.0D0
              if (neqns>=10)dunc10 = 0.0D0
              if (neqns>=11)dunc11 = 0.0D0
              if (neqns>=12)dunc12 = 0.0D0
            else
              ii = i
              jj = j
              kk = kp

              sx = BLK%DZTX(i,j,k)/BLK%RJCB(i,j,k)
              sy = BLK%DZTY(i,j,k)/BLK%RJCB(i,j,k)
              sz = BLK%DZTZ(i,j,k)/BLK%RJCB(i,j,k)

              uc = sx*u   +sy*v   +sz*w
              ac = sx*urot+sy*vrot+sz*wrot

              sxyz2= sx*sx+sy*sy+sz*sz
              sxyz = dsqrt(sxyz2)

              ggd = dsqrt((uc-ac)*(uc-ac)*(1.0d0-vp2/c2)*(1.0d0-vp2/c2)+4.0d0*vp2*sxyz2)

              eig1 = uc-ac
              eig2 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)+ggd)
              eig3 = 0.5d0*((uc-ac)*(1.0d0+vp2/c2)-ggd)

              eig1 = dabs(eig1)
              eig2 = dabs(eig2)
              eig3 = dabs(eig3)

              radc = 1.01D0*dmax1(eig1,eig2,eig3)

              qm01 = BLK%RHSS(ii,jj,kk,1)
              qm02 = BLK%RHSS(ii,jj,kk,2)
              qm03 = BLK%RHSS(ii,jj,kk,3)
              qm04 = BLK%RHSS(ii,jj,kk,4)
              qm05 = BLK%RHSS(ii,jj,kk,5)
              if (neqns>=6) qm06 = DQSTAR(ii,jj,kk,6)
              if (neqns>=7) qm07 = DQSTAR(ii,jj,kk,7)
              if (neqns>=8) qm08 = DQSTAR(ii,jj,kk,8)
              if (neqns>=9) qm09 = DQSTAR(ii,jj,kk,9)
              if (neqns>=10)qm10 = DQSTAR(ii,jj,kk,10)
              if (neqns>=11)qm11 = DQSTAR(ii,jj,kk,11)
              if (neqns>=12)qm12 = DQSTAR(ii,jj,kk,12)

              dunc01 = (sgn*radc+uc-ac             )*qm01 &
                      +(beta*vp2*sx*rh             )*qm02 &
                      +(beta*vp2*sy*rh             )*qm03 &
                      +(beta*vp2*sz*rh             )*qm04 &
                      -(beta*(1.0d0-vp2/c2)*(uc-ac))*qm05
              dunc02 = (0.0d0                      )*qm01 &
                      +(sgn*radc+uc-ac             )*qm02 &
                      +(0.0d0                      )*qm03 &
                      +(0.0d0                      )*qm04 &
                      +(sx*rrh                     )*qm05
              dunc03 = (0.0d0                      )*qm01 &
                      +(0.0d0                      )*qm02 &
                      +(sgn*radc+uc-ac             )*qm03 &
                      +(0.0d0                      )*qm04 &
                      +(sy*rrh                     )*qm05
              dunc04 = (0.0d0                      )*qm01 &
                      +(0.0d0                      )*qm02 &
                      +(0.0d0                      )*qm03 &
                      +(sgn*radc+uc-ac             )*qm04 &
                      +(sz*rrh                     )*qm05
              dunc05 = (0.0d0                      )*qm01 &
                      +(vp2*sx*rh                  )*qm02 &
                      +(vp2*sy*rh                  )*qm03 &
                      +(vp2*sz*rh                  )*qm04 &
                      +(sgn*radc+vp2/c2*(uc-ac)    )*qm05
              if (neqns>=6) dunc06 = (sgn*radc+uc-ac)*qm06
              if (neqns>=7) dunc07 = (sgn*radc+uc-ac)*qm07
              if (neqns>=8) dunc08 = (sgn*radc+uc-ac)*qm08
              if (neqns>=9) dunc09 = (sgn*radc+uc-ac)*qm09
              if (neqns>=10)dunc10 = (sgn*radc+uc-ac)*qm10
              if (neqns>=11)dunc11 = (sgn*radc+uc-ac)*qm11
              if (neqns>=12)dunc12 = (sgn*radc+uc-ac)*qm12

              dunc01 = 0.5D0*dunc01
              dunc02 = 0.5D0*dunc02
              dunc03 = 0.5D0*dunc03
              dunc04 = 0.5D0*dunc04
              dunc05 = 0.5D0*dunc05
              if (neqns>=6) dunc06 = 0.5D0*dunc06
              if (neqns>=7) dunc07 = 0.5D0*dunc07
              if (neqns>=8) dunc08 = 0.5D0*dunc08
              if (neqns>=9) dunc09 = 0.5D0*dunc09
              if (neqns>=10)dunc10 = 0.5D0*dunc10
              if (neqns>=11)dunc11 = 0.5D0*dunc11
              if (neqns>=12)dunc12 = 0.5D0*dunc12
            endif

            dinv = 1.0D0/(1.0D0+h*(rada+radb+radc))

            BLK%RHSS(i,j,k,1) = DQSTAR(i,j,k,1)-dinv*h*(duna01+dunb01+dunc01)
            BLK%RHSS(i,j,k,2) = DQSTAR(i,j,k,2)-dinv*h*(duna02+dunb02+dunc02)
            BLK%RHSS(i,j,k,3) = DQSTAR(i,j,k,3)-dinv*h*(duna03+dunb03+dunc03)
            BLK%RHSS(i,j,k,4) = DQSTAR(i,j,k,4)-dinv*h*(duna04+dunb04+dunc04)
            BLK%RHSS(i,j,k,5) = DQSTAR(i,j,k,5)-dinv*h*(duna05+dunb05+dunc05)
            if (neqns>=6) BLK%RHSS(i,j,k,6) = DQSTAR(i,j,k, 6)-dinv*h*(duna06+dunb06+dunc06)
            if (neqns>=7) BLK%RHSS(i,j,k,7) = DQSTAR(i,j,k, 7)-dinv*h*(duna07+dunb07+dunc07)
            if (neqns>=8) BLK%RHSS(i,j,k,8) = DQSTAR(i,j,k, 8)-dinv*h*(duna08+dunb08+dunc08)
            if (neqns>=9) BLK%RHSS(i,j,k,9) = DQSTAR(i,j,k, 9)-dinv*h*(duna09+dunb09+dunc09)
            if (neqns>=10)BLK%RHSS(i,j,k,10)= DQSTAR(i,j,k,10)-dinv*h*(duna10+dunb10+dunc10)
            if (neqns>=11)BLK%RHSS(i,j,k,11)= DQSTAR(i,j,k,11)-dinv*h*(duna11+dunb11+dunc11)
            if (neqns>=12)BLK%RHSS(i,j,k,12)= DQSTAR(i,j,k,12)-dinv*h*(duna12+dunb12+dunc12)
          enddo
        enddo
      enddo

      if (i2dimension>0) then
        do k=BLK%ksta_wh,BLK%kend_wh
          do j=BLK%jsta_wh,BLK%jend_wh
            do i=BLK%ista_wh,BLK%iend_wh
              BLK%RHSS(i,j,k,4) = 0.0d0
            enddo
          enddo
        enddo
      endif

      deallocate(DQSTAR)

    end subroutine lu_sgs_3d
  !-----------------------------------------------------------------------------
  end subroutine lu_sgs

end module mdl_implicit_sgs
