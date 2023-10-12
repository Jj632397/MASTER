module mdl_bound_inflow_outflow_tpipe
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_initmp
  use mdl_look_up_table
  use mdl_table
  implicit none

contains

  subroutine bound_inflow_outflow_tpipe(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: is_blk,js_blk,ks_blk,isplt_blk,jsplt_blk,ksplt_blk
    integer :: i,j,k,l

    integer :: ivalid

    real(8) :: u,v,w,p,t
    real(8) :: a_ne,b_ne,c_ne,d_ne

    real(8) :: tp_pos, tp_ppf

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

    !inflow boundary of main passage
    if (id_domain==1.and.is_blk==1) then

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            BLK%PRSS(BLK%ista,j,k) =( a_ne*BLK%PRSS(BLK%ista+1,j,k) &
                                     +b_ne*BLK%PRSS(BLK%ista+2,j,k) &
                                     +c_ne*BLK%PRSS(BLK%ista+3,j,k) &
                                     +d_ne*BLK%PRSS(BLK%ista+4,j,k) )/25.0d0
          enddo
        enddo

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            p = BLK%PRSS(BLK%ista,j,k)
            t = TTUNI_MP/TTREF

          tp_pos = BLK%YCRD(BLK%ista,j,k) - (-0.5d0)
          tp_ppf = 4.0d0*tp_pos*(1.0d0-tp_pos)     ! Poiseuille flow
          
          u  = UUUNI_MP/UUREF*tp_ppf
!            u  = UUUNI_MP/UUREF
            v  = 0.0d0
            w  = 0.0d0

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
              BLK%QCS1(BLK%ista,j,k,5+l) = 0.0d0
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              BLK%QCS1(BLK%ista,j,k,neqn_tbl_sa) = tkin_tbl_sa*VSUNI_MP/RHUNI_MP
            enddo
          enddo
        endif

        if (iact_passive_scalar>0) then
          do k=BLK%ksta,BLK%kend
            do j=BLK%jsta,BLK%jend
              BLK%QCS1(BLK%ista,j,k,neqn_passive_scalar) = 0.0D0
              BLK%QCS1(BLK%ista,j,k,neqn_passive_scalar+1) = 0.0D0
            enddo
          enddo
        endif

    endif


    !inflow boundary of sub passage
    if (id_domain==4.and.js_blk==1) then

        do k=BLK%ksta,BLK%kend
          do i=BLK%ista,BLK%iend
            BLK%PRSS(i,BLK%jsta,k) =( a_ne*BLK%PRSS(i,BLK%jsta+1,k) &
                                     +b_ne*BLK%PRSS(i,BLK%jsta+2,k) &
                                     +c_ne*BLK%PRSS(i,BLK%jsta+3,k) &
                                     +d_ne*BLK%PRSS(i,BLK%jsta+4,k) )/25.0d0
          enddo
        enddo

        do k=BLK%ksta,BLK%kend
          do i=BLK%ista,BLK%iend
            p = BLK%PRSS(i,BLK%jsta,k)
            t = TTUNI_SP/TTREF

          tp_pos = BLK%XCRD(i,BLK%jsta,k) - (-0.5d0)
          tp_ppf = 4.0d0*tp_pos*(1.0d0-tp_pos)     ! Poiseuille flow

            u  = 0.0d0
            v  = UUUNI_SP/UUREF*tp_ppf
!            v  = UUUNI_SP/UUREF
            w  = 0.0d0

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
              BLK%QCS1(i,BLK%jsta,k,5+l) = 0.0d0
            enddo
          enddo
        enddo

        if (iact_tbl_sa>0) then
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,BLK%jsta,k,neqn_tbl_sa) = tkin_tbl_sa*VSUNI_SP/RHUNI_SP
            enddo
          enddo
        endif

        if (iact_passive_scalar>0) then
          do k=BLK%ksta,BLK%kend
            do i=BLK%ista,BLK%iend
              BLK%QCS1(i,BLK%jsta,k,neqn_passive_scalar) = SCALR_SP
              BLK%QCS1(i,BLK%jsta,k,neqn_passive_scalar+1) = 0.0D0
            enddo
          enddo
        endif

    endif


    !outflow boundary of main passage
    if (id_domain==3.and.is_blk==isplt_blk) then

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            BLK%PRSS(BLK%iend,j,k) = PPBCK/PPREF
            BLK%TMPR(BLK%iend,j,k) =( a_ne*BLK%TMPR(BLK%iend-1,j,k) &
                                     +b_ne*BLK%TMPR(BLK%iend-2,j,k) &
                                     +c_ne*BLK%TMPR(BLK%iend-3,j,k) &
                                     +d_ne*BLK%TMPR(BLK%iend-4,j,k) )/25.0d0
            BLK%ULCT(BLK%iend,j,k) =( a_ne*BLK%ULCT(BLK%iend-1,j,k) &
                                     +b_ne*BLK%ULCT(BLK%iend-2,j,k) &
                                     +c_ne*BLK%ULCT(BLK%iend-3,j,k) &
                                     +d_ne*BLK%ULCT(BLK%iend-4,j,k) )/25.0d0
            BLK%VLCT(BLK%iend,j,k) =( a_ne*BLK%VLCT(BLK%iend-1,j,k) &
                                     +b_ne*BLK%VLCT(BLK%iend-2,j,k) &
                                     +c_ne*BLK%VLCT(BLK%iend-3,j,k) &
                                     +d_ne*BLK%VLCT(BLK%iend-4,j,k) )/25.0d0
            BLK%WLCT(BLK%iend,j,k) =( a_ne*BLK%WLCT(BLK%iend-1,j,k) &
                                     +b_ne*BLK%WLCT(BLK%iend-2,j,k) &
                                     +c_ne*BLK%WLCT(BLK%iend-3,j,k) &
                                     +d_ne*BLK%WLCT(BLK%iend-4,j,k) )/25.0d0
          enddo
        enddo

        do k=BLK%ksta,BLK%kend
          do j=BLK%jsta,BLK%jend
            p = BLK%PRSS(BLK%iend,j,k)
            t = BLK%TMPR(BLK%iend,j,k)

            u = BLK%ULCT(BLK%iend,j,k)
            v = BLK%VLCT(BLK%iend,j,k)
            w = BLK%WLCT(BLK%iend,j,k)

            if (u<0.0d0) u = -1.0d0*u

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
              BLK%QCS1(BLK%iend,j,k,5+l) =( a_ne*BLK%QCS1(BLK%iend-1,j,k,5+l) &
                                           +b_ne*BLK%QCS1(BLK%iend-2,j,k,5+l) &
                                           +c_ne*BLK%QCS1(BLK%iend-3,j,k,5+l) &
                                           +d_ne*BLK%QCS1(BLK%iend-4,j,k,5+l) )/25.0d0
            enddo
          enddo
        enddo


    endif

  end subroutine bound_inflow_outflow_tpipe

end module mdl_bound_inflow_outflow_tpipe
