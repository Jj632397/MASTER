module mdl_bound_wall
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_initmp
  implicit none

  integer, private, parameter :: nbound = 7
  integer, private, parameter, dimension(nbound) :: id_dms_wall = (/1,1,2,3,3,4,4/) !DomainID
  integer, private, parameter, dimension(nbound) :: id_bds_wall = (/3,4,4,3,4,1,2/) !FaceID

contains

  subroutine bound_wall(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: is_blk,js_blk,ks_blk,isplt_blk,jsplt_blk,ksplt_blk
    integer :: i,j,k,l,n

    real(8) :: rh,u,v,w,ei,p,t,ak,cfr,ys
    real(8) :: a_ne,b_ne,c_ne,d_ne,r_ne

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain, &
                                                   is_blk,js_blk,ks_blk)

    !4th-order
    a_ne = 48.0d0
    b_ne =-36.0d0
    c_ne = 16.0d0
    d_ne = -3.0d0

    r_ne = 25.0d0

    !1st-order
    ! a_ne = 1.0d0
    ! b_ne = 0.0d0
    ! c_ne = 0.0d0
    ! d_ne = 0.0d0

    ! r_ne = 1.0d0

    do n=1,nbound

      if (id_domain/=id_dms_wall(n)) cycle

      if (id_bds_wall(n)==1.and.is_blk==1) then

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            p  =( a_ne*BLK%PRSS(BLK%ista+1,j,k) &
                 +b_ne*BLK%PRSS(BLK%ista+2,j,k) &
                 +c_ne*BLK%PRSS(BLK%ista+3,j,k) &
                 +d_ne*BLK%PRSS(BLK%ista+4,j,k) )/r_ne

            t  =( a_ne*BLK%TMPR(BLK%ista+1,j,k) &
                 +b_ne*BLK%TMPR(BLK%ista+2,j,k) &
                 +c_ne*BLK%TMPR(BLK%ista+3,j,k) &
                 +d_ne*BLK%TMPR(BLK%ista+4,j,k) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(BLK%ista,j,k,1) = rh
!            BLK%QCS1(BLK%ista,j,k,2) = rh*u
!            BLK%QCS1(BLK%ista,j,k,3) = rh*v
!            BLK%QCS1(BLK%ista,j,k,4) = rh*w
!            BLK%QCS1(BLK%ista,j,k,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(BLK%ista,j,k,1) = t
            BLK%QCS1(BLK%ista,j,k,2) = u
            BLK%QCS1(BLK%ista,j,k,3) = v
            BLK%QCS1(BLK%ista,j,k,4) = w
            BLK%QCS1(BLK%ista,j,k,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              ys =( a_ne*BLK%YSPC(BLK%ista+1,j,k,l) &
                   +b_ne*BLK%YSPC(BLK%ista+2,j,k,l) &
                   +c_ne*BLK%YSPC(BLK%ista+3,j,k,l) &
                   +d_ne*BLK%YSPC(BLK%ista+4,j,k,l) )/r_ne
!              BLK%QCS1(BLK%ista,j,k,5+l) = rh*ys
              BLK%QCS1(BLK%ista,j,k,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              BLK%QCS1(BLK%ista,j,k,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      else if (id_bds_wall(n)==2.and.is_blk==isplt_blk) then

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            p  =( a_ne*BLK%PRSS(BLK%iend-1,j,k) &
                 +b_ne*BLK%PRSS(BLK%iend-2,j,k) &
                 +c_ne*BLK%PRSS(BLK%iend-3,j,k) &
                 +d_ne*BLK%PRSS(BLK%iend-4,j,k) )/r_ne

            t  =( a_ne*BLK%TMPR(BLK%iend-1,j,k) &
                 +b_ne*BLK%TMPR(BLK%iend-2,j,k) &
                 +c_ne*BLK%TMPR(BLK%iend-3,j,k) &
                 +d_ne*BLK%TMPR(BLK%iend-4,j,k) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(BLK%iend,j,k,1) = rh
!            BLK%QCS1(BLK%iend,j,k,2) = rh*u
!            BLK%QCS1(BLK%iend,j,k,3) = rh*v
!            BLK%QCS1(BLK%iend,j,k,4) = rh*w
!            BLK%QCS1(BLK%iend,j,k,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(BLK%iend,j,k,1) = t
            BLK%QCS1(BLK%iend,j,k,2) = u
            BLK%QCS1(BLK%iend,j,k,3) = v
            BLK%QCS1(BLK%iend,j,k,4) = w
            BLK%QCS1(BLK%iend,j,k,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              ys =( a_ne*BLK%YSPC(BLK%iend-1,j,k,l) &
                   +b_ne*BLK%YSPC(BLK%iend-2,j,k,l) &
                   +c_ne*BLK%YSPC(BLK%iend-3,j,k,l) &
                   +d_ne*BLK%YSPC(BLK%iend-4,j,k,l) )/r_ne
!              BLK%QCS1(BLK%iend,j,k,5+l) = rh*ys
              BLK%QCS1(BLK%iend,j,k,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              BLK%QCS1(BLK%iend,j,k,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      else if (id_bds_wall(n)==3.and.js_blk==1) then

        do k=BLK%ksta,BLK%kend
          do i=BLK%ista,BLK%iend
            p  =( a_ne*BLK%PRSS(i,BLK%jsta+1,k) &
                 +b_ne*BLK%PRSS(i,BLK%jsta+2,k) &
                 +c_ne*BLK%PRSS(i,BLK%jsta+3,k) &
                 +d_ne*BLK%PRSS(i,BLK%jsta+4,k) )/r_ne

            t  =( a_ne*BLK%TMPR(i,BLK%jsta+1,k) &
                 +b_ne*BLK%TMPR(i,BLK%jsta+2,k) &
                 +c_ne*BLK%TMPR(i,BLK%jsta+3,k) &
                 +d_ne*BLK%TMPR(i,BLK%jsta+4,k) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(i,BLK%jsta,k,1) = rh
!            BLK%QCS1(i,BLK%jsta,k,2) = rh*u
!            BLK%QCS1(i,BLK%jsta,k,3) = rh*v
!            BLK%QCS1(i,BLK%jsta,k,4) = rh*w
!            BLK%QCS1(i,BLK%jsta,k,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(i,BLK%jsta,k,1) = t
            BLK%QCS1(i,BLK%jsta,k,2) = u
            BLK%QCS1(i,BLK%jsta,k,3) = v
            BLK%QCS1(i,BLK%jsta,k,4) = w
            BLK%QCS1(i,BLK%jsta,k,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              ys =( a_ne*BLK%YSPC(i,BLK%jsta+1,k,l) &
                   +b_ne*BLK%YSPC(i,BLK%jsta+2,k,l) &
                   +c_ne*BLK%YSPC(i,BLK%jsta+3,k,l) &
                   +d_ne*BLK%YSPC(i,BLK%jsta+4,k,l) )/r_ne
!              BLK%QCS1(i,BLK%jsta,k,5+l) = rh*ys
              BLK%QCS1(i,BLK%jsta,k,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,BLK%jsta,k,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      else if (id_bds_wall(n)==4.and.js_blk==jsplt_blk) then

        do k=BLK%ksta,BLK%kend
          do i=BLK%ista,BLK%iend
            p  =( a_ne*BLK%PRSS(i,BLK%jend-1,k) &
                 +b_ne*BLK%PRSS(i,BLK%jend-2,k) &
                 +c_ne*BLK%PRSS(i,BLK%jend-3,k) &
                 +d_ne*BLK%PRSS(i,BLK%jend-4,k) )/r_ne

            t  =( a_ne*BLK%TMPR(i,BLK%jend-1,k) &
                 +b_ne*BLK%TMPR(i,BLK%jend-2,k) &
                 +c_ne*BLK%TMPR(i,BLK%jend-3,k) &
                 +d_ne*BLK%TMPR(i,BLK%jend-4,k) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(i,BLK%jend,k,1) = rh
!            BLK%QCS1(i,BLK%jend,k,2) = rh*u
!            BLK%QCS1(i,BLK%jend,k,3) = rh*v
!            BLK%QCS1(i,BLK%jend,k,4) = rh*w
!            BLK%QCS1(i,BLK%jend,k,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(i,BLK%jend,k,1) = t
            BLK%QCS1(i,BLK%jend,k,2) = u
            BLK%QCS1(i,BLK%jend,k,3) = v
            BLK%QCS1(i,BLK%jend,k,4) = w
            BLK%QCS1(i,BLK%jend,k,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              ys =( a_ne*BLK%YSPC(i,BLK%jend-1,k,l) &
                   +b_ne*BLK%YSPC(i,BLK%jend-2,k,l) &
                   +c_ne*BLK%YSPC(i,BLK%jend-3,k,l) &
                   +d_ne*BLK%YSPC(i,BLK%jend-4,k,l) )/r_ne
!              BLK%QCS1(i,BLK%jend,k,5+l) = rh*ys
              BLK%QCS1(i,BLK%jend,k,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,BLK%jend,k,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      else if (id_bds_wall(n)==5.and.ks_blk==1) then

        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            p  =( a_ne*BLK%PRSS(i,j,BLK%ksta+1) &
                 +b_ne*BLK%PRSS(i,j,BLK%ksta+2) &
                 +c_ne*BLK%PRSS(i,j,BLK%ksta+3) &
                 +d_ne*BLK%PRSS(i,j,BLK%ksta+4) )/r_ne

            t  =( a_ne*BLK%TMPR(i,j,BLK%ksta+1) &
                 +b_ne*BLK%TMPR(i,j,BLK%ksta+2) &
                 +c_ne*BLK%TMPR(i,j,BLK%ksta+3) &
                 +d_ne*BLK%TMPR(i,j,BLK%ksta+4) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(i,j,BLK%ksta,1) = rh
!            BLK%QCS1(i,j,BLK%ksta,2) = rh*u
!            BLK%QCS1(i,j,BLK%ksta,3) = rh*v
!            BLK%QCS1(i,j,BLK%ksta,4) = rh*w
!            BLK%QCS1(i,j,BLK%ksta,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(i,j,BLK%ksta,1) = t
            BLK%QCS1(i,j,BLK%ksta,2) = u
            BLK%QCS1(i,j,BLK%ksta,3) = v
            BLK%QCS1(i,j,BLK%ksta,4) = w
            BLK%QCS1(i,j,BLK%ksta,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              ys =( a_ne*BLK%YSPC(i,j,BLK%ksta+1,l) &
                   +b_ne*BLK%YSPC(i,j,BLK%ksta+2,l) &
                   +c_ne*BLK%YSPC(i,j,BLK%ksta+3,l) &
                   +d_ne*BLK%YSPC(i,j,BLK%ksta+4,l) )/r_ne
!              BLK%QCS1(i,j,BLK%ksta,5+l) = rh*ys
              BLK%QCS1(i,j,BLK%ksta,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,j,BLK%ksta,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      else if (id_bds_wall(n)==6.and.ks_blk==ksplt_blk) then

        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            p  =( a_ne*BLK%PRSS(i,j,BLK%kend-1) &
                 +b_ne*BLK%PRSS(i,j,BLK%kend-2) &
                 +c_ne*BLK%PRSS(i,j,BLK%kend-3) &
                 +d_ne*BLK%PRSS(i,j,BLK%kend-4) )/r_ne

            t  =( a_ne*BLK%TMPR(i,j,BLK%kend-1) &
                 +b_ne*BLK%TMPR(i,j,BLK%kend-2) &
                 +c_ne*BLK%TMPR(i,j,BLK%kend-3) &
                 +d_ne*BLK%TMPR(i,j,BLK%kend-4) )/r_ne

            ak = cp_i/(cp_i-cr_i)

            rh = p/(cr_i*t)
            u  = 0.0d0
            v  = 0.0d0
            w  = 0.0d0
            ei = p/(rh*(ak-1.0d0))

!            BLK%QCS1(i,j,BLK%kend,1) = rh
!            BLK%QCS1(i,j,BLK%kend,2) = rh*u
!            BLK%QCS1(i,j,BLK%kend,3) = rh*v
!            BLK%QCS1(i,j,BLK%kend,4) = rh*w
!            BLK%QCS1(i,j,BLK%kend,5) = rh*(ei+0.5d0*(u*u+v*v+w*w))
            BLK%QCS1(i,j,BLK%kend,1) = t
            BLK%QCS1(i,j,BLK%kend,2) = u
            BLK%QCS1(i,j,BLK%kend,3) = v
            BLK%QCS1(i,j,BLK%kend,4) = w
            BLK%QCS1(i,j,BLK%kend,5) = p
          enddo
        enddo

        do l=1,neqns-5
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              ys =( a_ne*BLK%YSPC(i,j,BLK%kend-1,l) &
                   +b_ne*BLK%YSPC(i,j,BLK%kend-2,l) &
                   +c_ne*BLK%YSPC(i,j,BLK%kend-3,l) &
                   +d_ne*BLK%YSPC(i,j,BLK%kend-4,l) )/r_ne
!              BLK%QCS1(i,j,BLK%kend,5+l) = rh*ys
              BLK%QCS1(i,j,BLK%kend,5+l) = ys
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do j=BLK%jsta,BLK%jend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,j,BLK%kend,neqn_tbl_sa) = 0.0d0
            enddo
          enddo
        endif

      endif

    enddo

  end subroutine bound_wall

end module mdl_bound_wall
