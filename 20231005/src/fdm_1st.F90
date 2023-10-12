module mdl_fdm1
  use mdl_band
  implicit none

  integer, private, parameter :: align_mem_access = 1 !<0>:simple <1>:longer vector length

contains
!------------------------------------------------------------------------------!
!
!F : Input function
!
!DF: Output function
!
!idir
!   <1>:xsi <2>:eta <3>:zta
!
!iperi
!   <1>:periodic <else>:bounded
!   If 'iperi==1', scheme designated by isc(3) is applied for all points.
!
!icompact
!   <1>:Compact scheme <else>:explicit scheme
!
!Choose Scheme
!----------- icompact==1 -------------------------------------------------------
! isc(1)[Boundary Point 1]
!      <1>:2nd-order-explicit <2>:2nd-order-compact
!      <3>:4th-order-compact  <4>:6th-order-compact
!      <else>:1st-order-explicit
! isc(2)[Boundary Point 2]
!      <1>:4th-order-compact  <2>:6th-order-compact
!      <else>:2nd-order-explicit
! isc(3)[Interior Points]
!      <1>:6th-order-compact <else>:4th-order-compact
!      <2>:4th-order-compact-optimized
! isc(4)[Boundary Point N-1]
!      <1>:4th-order-compact  <2>:6th-order-compact
!      <else>:2nd-order-explicit
! isc(5)[Boundary Point N]
!      <1>:2nd-order-explicit <2>:2nd-order-compact
!      <3>:4th-order-compact  <4>:6th-order-compact
!      <else>:1st-order-explicit
!----------- icompact/=1 -------------------------------------------------------
! isc(1)[Boundary Point 1]
!      <1>:2nd-order-explicit <2>:4th-order-explicit
!      <else>:1st-order-explicit
! isc(2)[Boundary Point 2]
!      <1>:4th-order-explicit
!      <else>:2nd-order-explicit
! isc(3)[Interior Points]
!      <1>:4th-order-explicit <else>:2nd-order-explicit
! isc(4)[Boundary Point N-1]
!      <1>:4th-order-explicit
!      <else>:2nd-order-explicit
! isc(5)[Boundary Point N]
!      <1>:2nd-order-explicit <2>:4th-order-explicit
!      <else>:1st-order-explicit
!-------------------------------------------------------------------------------!
  subroutine calc_fdm_1st(F,DF,idir,iperi,icompact,isc,imax,jmax,kmax)
    real(8), dimension(:,:,:), intent(in) :: F
    real(8), dimension(:,:,:), intent(out):: DF
    integer, intent(in) :: idir,iperi,icompact
    integer, dimension(:), intent(in) :: isc
    integer, intent(in) :: imax,jmax,kmax

    integer :: i,j,k
    integer :: ii,jj,iimax,jjmax

    real(8), dimension(:,:), allocatable ::  F_ALN
    real(8), dimension(:,:), allocatable :: DF_ALN

    if (idir==1) then
      iimax = jmax*kmax
      jjmax = imax
    else if (idir==2) then
      iimax = imax*kmax
      jjmax = jmax
    else if (idir==3) then
      iimax = imax*jmax
      jjmax = kmax
    endif

    allocate( F_ALN(iimax,jjmax))
    allocate(DF_ALN(iimax,jjmax))

    if (idir==1) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = j+jmax*(k-1)
              jj = i
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = jj
            j = mod(ii-1,jmax)+1
            k = (ii-1)/jmax+1
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    else if (idir==2) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(k-1)
              jj = j
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = jj
            k = (ii-1)/imax+1
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    else if (idir==3) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(j-1)
              jj = k
              F_ALN(ii,jj) = F(i,j,k)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = (ii-1)/imax+1
            k = jj
            F_ALN(ii,jj) = F(i,j,k)
          enddo
        enddo
      endif

    endif


    call calc_fdm_1st_align(F_ALN,DF_ALN,iperi,icompact,isc,iimax,jjmax)


    if (idir==1) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = j+jmax*(k-1)
              jj = i
              DF(i,j,k) = DF_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = jj
            j = mod(ii-1,jmax)+1
            k = (ii-1)/jmax+1
            DF(i,j,k) = DF_ALN(ii,jj)
          enddo
        enddo
      endif

    else if (idir==2) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(k-1)
              jj = j
              DF(i,j,k) = DF_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = jj
            k = (ii-1)/imax+1
            DF(i,j,k) = DF_ALN(ii,jj)
          enddo
        enddo
      endif

    else if (idir==3) then

      if (align_mem_access==0) then
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              ii = i+imax*(j-1)
              jj = k
              DF(i,j,k) = DF_ALN(ii,jj)
            enddo
          enddo
        enddo
      else if (align_mem_access==1) then
        do jj=1,jjmax
          do ii=1,iimax
            i = mod(ii-1,imax)+1
            j = (ii-1)/imax+1
            k = jj
            DF(i,j,k) = DF_ALN(ii,jj)
          enddo
        enddo
      endif

    endif


    deallocate( F_ALN)
    deallocate(DF_ALN)

  end subroutine calc_fdm_1st
!-------------------------------------------------------------------------------
  subroutine calc_fdm_1st_align(F,DF,iperi,icompact,isc,iimax,jjmax)
    real(8), dimension(:,:), intent(in) :: F
    real(8), dimension(:,:), intent(out):: DF
    integer, intent(in) :: iperi,icompact
    integer, dimension(:), intent(in) :: isc
    integer, intent(in) :: iimax,jjmax

    integer :: ii,jj
    real(8), dimension(:,:), allocatable :: aa,bb,cc,xx
    real(8) :: alpha_ic,a_ic,b_ic
    real(8) :: alpha_1,a1,b1,c1,d1,e1,f1
    real(8) :: alpha_n,an,bn,cn,dn,en,fn
    real(8) :: alpha_21,alpha_22,a2,b2,c2,d2,e2,f2
    real(8) :: alpha_m1,alpha_m2,am,bm,cm,dm,em,fm

    if (jjmax<5) return

    call set_coefficients

    allocate(aa(1    ,jjmax))
    allocate(bb(1    ,jjmax))
    allocate(cc(1    ,jjmax))
    allocate(xx(iimax,jjmax))

    !Set L.H.S
    if (iperi==1) then
      do jj=1,jjmax
        aa(1,jj) = alpha_ic
        bb(1,jj) = 1.0d0
        cc(1,jj) = alpha_ic
      enddo
    else
      aa(1,1) = 0.0d0
      bb(1,1) = 1.0d0
      cc(1,1) = alpha_1
      aa(1,2) = alpha_21
      bb(1,2) = 1.0d0
      cc(1,2) = alpha_22
      do jj=3,jjmax-2
        aa(1,jj) = alpha_ic
        bb(1,jj) = 1.0d0
        cc(1,jj) = alpha_ic
      enddo
      aa(1,jjmax-1) = alpha_m1
      bb(1,jjmax-1) = 1.0d0
      cc(1,jjmax-1) = alpha_m2
      aa(1,jjmax  ) = alpha_n
      bb(1,jjmax  ) = 1.0d0
      cc(1,jjmax  ) = 0.0d0
    endif

    call set_rhs

    if (icompact==1) then
      DF = 0.0d0
      if (iperi==1) then
        call tridiag1_cyc(aa,bb,cc,xx,DF,iimax,jjmax)
      else
        call tridiag1(aa,bb,cc,xx,DF,iimax,jjmax)
      endif
    else
      DF = xx
    endif

    deallocate(aa)
    deallocate(bb)
    deallocate(cc)
    deallocate(xx)

  contains
!-------------------------------------------------------------------------------
    subroutine set_rhs

      !Interior
      do jj=3,jjmax-2
        do ii=1,iimax
          xx(ii,jj) = a_ic*(F(ii,jj+1)-F(ii,jj-1)) &
                     +b_ic*(F(ii,jj+2)-F(ii,jj-2))
        enddo
      enddo

      !Near boundary
      if (iperi==1) then

        do ii=1,iimax
          xx(ii,1) = a_ic*(F(ii,2)-F(ii,jjmax  )) &
                    +b_ic*(F(ii,3)-F(ii,jjmax-1))
          xx(ii,2) = a_ic*(F(ii,3)-F(ii,1      )) &
                    +b_ic*(F(ii,4)-F(ii,jjmax  ))

          xx(ii,jjmax  ) = a_ic*(F(ii,1    )-F(ii,jjmax-1)) &
                          +b_ic*(F(ii,2    )-F(ii,jjmax-2))
          xx(ii,jjmax-1) = a_ic*(F(ii,jjmax)-F(ii,jjmax-2)) &
                          +b_ic*(F(ii,1    )-F(ii,jjmax-3))
        enddo

      else

        do ii=1,iimax
          xx(ii,1) = a1*F(ii,1) + b1*F(ii,2) &
                    +c1*F(ii,3) + d1*F(ii,4) &
                    +e1*F(ii,5) + f1*F(ii,6)
          xx(ii,2) = a2*F(ii,1) + b2*F(ii,2) &
                    +c2*F(ii,3) + d2*F(ii,4) &
                    +e2*F(ii,5) + f2*F(ii,6)

          xx(ii,jjmax-1) =-am*F(ii,jjmax  ) - bm*F(ii,jjmax-1) &
                          -cm*F(ii,jjmax-2) - dm*F(ii,jjmax-3) &
                          -em*F(ii,jjmax-4) - fm*F(ii,jjmax-5)
          xx(ii,jjmax  ) =-an*F(ii,jjmax  ) - bn*F(ii,jjmax-1) &
                          -cn*F(ii,jjmax-2) - dn*F(ii,jjmax-3) &
                          -en*F(ii,jjmax-4) - fn*F(ii,jjmax-5)
        enddo

      endif

    end subroutine set_rhs
!-------------------------------------------------------------------------------
    subroutine set_coefficients

      !Boundary 1
      if (icompact==1) then !---------------------------------------------------
        if (isc(1)==1) then
          !2nd-order-explicit
          alpha_1 = 0.0d0
          a1 =-3.0d0/2.0d0
          b1 = 2.0d0
          c1 =-1.0d0/2.0d0
          d1 = 0.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        else if (isc(1)==2) then
          !2nd-order-compact
          alpha_1 = 1.0d0
          a1 =-2.0d0
          b1 = 2.0d0
          c1 = 0.0d0
          d1 = 0.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        else if (isc(1)==3) then
          !4th-order-compact
          alpha_1 = 3.0d0
          a1 =-17.0d0/6.0d0
          b1 = 3.0d0/2.0d0
          c1 = 3.0d0/2.0d0
          d1 =-1.0d0/6.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        else if (isc(1)==4) then
          !6th-order-compact
          alpha_1 = 5.0d0
          a1 =-197.0d0/60.0d0
          b1 =-5.0d0/12.0d0
          c1 = 5.0d0
          d1 =-5.0d0/3.0d0
          e1 = 5.0d0/12.0d0
          f1 =-1.0d0/20.0d0
        else
          !1st-order-explicit
          alpha_1 = 0.0d0
          a1 =-1.0d0
          b1 = 1.0d0
          c1 = 0.0d0
          d1 = 0.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        endif
      else !--------------------------------------------------------------------
        if (isc(1)==1) then
          !2nd-order-explicit
          alpha_1 = 0.0d0
          a1 =-3.0d0/2.0d0
          b1 = 2.0d0
          c1 =-1.0d0/2.0d0
          d1 = 0.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        else if (isc(1)==2) then
          !4th-order-explicit
          alpha_1 = 0.0d0
          a1 =-25.0d0/12.0d0
          b1 = 4.0d0
          c1 =-3.0d0
          d1 = 4.0d0/3.0d0
          e1 =-1.0d0/4.0d0
          f1 = 0.0d0
        else
          !1st-order-explicit
          alpha_1 = 0.0d0
          a1 =-1.0d0
          b1 = 1.0d0
          c1 = 0.0d0
          d1 = 0.0d0
          e1 = 0.0d0
          f1 = 0.0d0
        endif
      endif !-------------------------------------------------------------------

      !Boundary 2
      if (icompact==1) then !---------------------------------------------------
        if (isc(2)==1) then
          !4th-order-compact
          alpha_21 = 1.0d0/4.0d0
          alpha_22 = 1.0d0/4.0d0
          a2 =-3.0d0/4.0d0
          b2 = 0.0d0
          c2 = 3.0d0/4.0d0
          d2 = 0.0d0
          e2 = 0.0d0
          f2 = 0.0d0
        else if (isc(2)==2) then
          !6th-order-compact
          alpha_21 = 2.0d0/11.0d0
          alpha_22 = 2.0d0/11.0d0
          a2 =-20.0d0/33.0d0
          b2 =-35.0d0/132.0d0
          c2 = 34.0d0/33.0d0
          d2 =-7.0d0/33.0d0
          e2 = 2.0d0/33.0d0
          f2 =-1.0d0/132.0d0
        else
          !2nd-order-explicit
          alpha_21 = 0.0d0
          alpha_22 = 0.0d0
          a2 =-1.0d0/2.0d0
          b2 = 0.0d0
          c2 = 1.0d0/2.0d0
          d2 = 0.0d0
          e2 = 0.0d0
          f2 = 0.0d0
        endif
      else !--------------------------------------------------------------------
        if (isc(2)==1) then
          !4th-order-explicit
          alpha_21 = 0.0d0
          alpha_22 = 0.0d0
          a2 =-1.0d0/4.0d0
          b2 =-5.0d0/6.0d0
          c2 = 3.0d0/2.0d0
          d2 =-1.0d0/2.0d0
          e2 = 1.0d0/12.0d0
          f2 = 0.0d0
        else
          !2nd-order-explicit
          alpha_21 = 0.0d0
          alpha_22 = 0.0d0
          a2 =-1.0d0/2.0d0
          b2 = 0.0d0
          c2 = 1.0d0/2.0d0
          d2 = 0.0d0
          e2 = 0.0d0
          f2 = 0.0d0
        endif
      endif !-------------------------------------------------------------------

      !Interior
      if (icompact==1) then !---------------------------------------------------
        if (isc(3)==1) then
          !6th-order-compact
          alpha_ic = 1.0d0/3.0d0
          a_ic = 7.0d0/9.0d0
          b_ic = 1.0d0/36.0d0
        else if (isc(3)==2) then
          !4th-order-compact-optimized
          alpha_ic = 0.381365d0
          a_ic = 1.5875767d0/2.0d0
          b_ic = 0.1751533d0/4.0d0
        else
          !4th-order-compact
          alpha_ic = 1.0d0/4.0d0
          a_ic = 3.0d0/4.0d0
          b_ic = 0.0d0
        endif
      else !--------------------------------------------------------------------
        if (isc(3)==1) then
          !4th-order-explicit
          alpha_ic = 0.0d0
          a_ic = 2.0d0/3.0d0
          b_ic =-1.0d0/12.0d0
        else
          !2nd-order-explicit
          alpha_ic = 0.0d0
          a_ic = 1.0d0/2.0d0
          b_ic = 0.0d0
        endif
      endif !-------------------------------------------------------------------

      !Boundary N-1
      if (icompact==1) then !---------------------------------------------------
        if (isc(4)==1) then
          !4th-order-compact
          alpha_m1 = 1.0d0/4.0d0
          alpha_m2 = 1.0d0/4.0d0
          am =-3.0d0/4.0d0
          bm = 0.0d0
          cm = 3.0d0/4.0d0
          dm = 0.0d0
          em = 0.0d0
          fm = 0.0d0
        else if (isc(4)==2) then
          !6th-order-compact
          alpha_m1 = 2.0d0/11.0d0
          alpha_m2 = 2.0d0/11.0d0
          am =-20.0d0/33.0d0
          bm =-35.0d0/132.0d0
          cm = 34.0d0/33.0d0
          dm =-7.0d0/33.0d0
          em = 2.0d0/33.0d0
          fm =-1.0d0/132.0d0
        else
          !2nd-order-explicit
          alpha_m1 = 0.0d0
          alpha_m2 = 0.0d0
          am =-1.0d0/2.0d0
          bm = 0.0d0
          cm = 1.0d0/2.0d0
          dm = 0.0d0
          em = 0.0d0
          fm = 0.0d0
        endif
      else !--------------------------------------------------------------------
        if (isc(4)==1) then
          !4th-order-explicit
          alpha_m1 = 0.0d0
          alpha_m2 = 0.0d0
          am =-1.0d0/4.0d0
          bm =-5.0d0/6.0d0
          cm = 3.0d0/2.0d0
          dm =-1.0d0/2.0d0
          em = 1.0d0/12.0d0
          fm = 0.0d0
        else
          !2nd-order-explicit
          alpha_m1 = 0.0d0
          alpha_m2 = 0.0d0
          am =-1.0d0/2.0d0
          bm = 0.0d0
          cm = 1.0d0/2.0d0
          dm = 0.0d0
          em = 0.0d0
          fm = 0.0d0
        endif
      endif !-------------------------------------------------------------------

      !Boundary N
      if (icompact==1) then !---------------------------------------------------
        if (isc(5)==1) then
          !2nd-order-explicit
          alpha_n = 0.0d0
          an =-3.0d0/2.0d0
          bn = 2.0d0
          cn =-1.0d0/2.0d0
          dn = 0.0d0
          en = 0.0d0
          fn = 0.0d0
        else if (isc(5)==2) then
          !2nd-order-compact
          alpha_n = 1.0d0
          an =-2.0d0
          bn = 2.0d0
          cn = 0.0d0
          dn = 0.0d0
          en = 0.0d0
          fn = 0.0d0
        else if (isc(5)==3) then
          !4th-order-compact
          alpha_n = 3.0d0
          an =-17.0d0/6.0d0
          bn = 3.0d0/2.0d0
          cn = 3.0d0/2.0d0
          dn =-1.0d0/6.0d0
          en = 0.0d0
          fn = 0.0d0
        else if (isc(5)==4) then
          !6th-order-compact
          alpha_n = 5.0d0
          an =-197.0d0/60.0d0
          bn =-5.0d0/12.0d0
          cn = 5.0d0
          dn =-5.0d0/3.0d0
          en = 5.0d0/12.0d0
          fn =-1.0d0/20.0d0
        else
          !1st-order-explicit
          alpha_n = 0.0d0
          an =-1.0d0
          bn = 1.0d0
          cn = 0.0d0
          dn = 0.0d0
          en = 0.0d0
          fn = 0.0d0
        endif
      else !--------------------------------------------------------------------
        if (isc(5)==1) then
          !2nd-order-explicit
          alpha_n = 0.0d0
          an =-3.0d0/2.0d0
          bn = 2.0d0
          cn =-1.0d0/2.0d0
          dn = 0.0d0
          en = 0.0d0
          fn = 0.0d0
        else if (isc(5)==2) then
          !4th-order-explicit
          alpha_n = 0.0d0
          an =-25.0d0/12.0d0
          bn = 4.0d0
          cn =-3.0d0
          dn = 4.0d0/3.0d0
          en =-1.0d0/4.0d0
          fn = 0.0d0
        else
          !1st-order-explicit
          alpha_n = 0.0d0
          an =-1.0d0
          bn = 1.0d0
          cn = 0.0d0
          dn = 0.0d0
          en = 0.0d0
          fn = 0.0d0
        endif
      endif !-------------------------------------------------------------------

    end subroutine set_coefficients

  end subroutine calc_fdm_1st_align

end module mdl_fdm1
