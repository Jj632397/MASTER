module mdl_result
    use mdl_decompo
    use mdl_mpisub_sbsp
    use mdl_param
    use mdl_refs
    use mdl_block
    use mdl_table
    use mdl_initmp
    implicit none

contains

  subroutine result(BLK)
    type(block), intent(inout) :: BLK

    integer :: i,j,k,l
    real(8) :: rh,rht,rhp,hst,hsp,c2

    integer :: ivalid

    if (i2dimension>0) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            BLK%QCS1(i,j,k,4) = 0.0d0
          enddo
        enddo
      enddo
    endif

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          BLK%TMPR(i,j,k) = BLK%QCS1(i,j,k,1)
          BLK%ULCT(i,j,k) = BLK%QCS1(i,j,k,2)
          BLK%VLCT(i,j,k) = BLK%QCS1(i,j,k,3)
          BLK%WLCT(i,j,k) = BLK%QCS1(i,j,k,4)
          BLK%PRSS(i,j,k) = BLK%QCS1(i,j,k,5)
          BLK%VSBV(i,j,k) = 0.0d0
        enddo
      enddo
    enddo

    do l=1,neqns-5
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            BLK%YSPC(i,j,k,l) = BLK%QCS1(i,j,k,5+l)
          enddo
        enddo
      enddo
    enddo

    if (iact_tbl_sa>0) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            BLK%YSPC(i,j,k,neqn_tbl_sa-5) = dmax1(BLK%YSPC(i,j,k,neqn_tbl_sa-5),0.0d0)
          enddo
        enddo
      enddo
    endif

    call look_up_table_calc_func(TBL_TP_to_DNST,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%DNST,RHREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_HTPS,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%HTPS,TTREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_DRTM,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%DRTM,RHREF/TTREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_DRPR,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%DRPR,RHREF/PPREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_SHCP,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%SHCP,CPREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_DHPR,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%DHPR,TTREF/PPREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_VSSV,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%VSSV,VSREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    call look_up_table_calc_func(TBL_TP_to_CKPV,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    BLK%CKPV,CKREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)


    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          rh  = BLK%DNST(i,j,k)
          rht = BLK%DRTM(i,j,k)
          rhp = BLK%DRPR(i,j,k)
          hst = BLK%SHCP(i,j,k)*CPREF
          hsp = BLK%DHPR(i,j,k)

          c2 = rh*hst/(rh*rhp*hst+rht*(1.0d0-rh*hsp))

          BLK%CSOS(i,j,k) = dsqrt(c2)
        enddo
      enddo
    enddo

  end subroutine result

end module mdl_result
