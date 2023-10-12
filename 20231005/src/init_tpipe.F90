module mdl_init_tpipe
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_block
  use mdl_param
  use mdl_refs
  use mdl_initmp
  implicit none

contains

  subroutine init_tpipe(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: i,j,k,l

    real(8) :: u,v,w,p,t,tk

    real(8) :: tp_pos, tp_ppf
    
    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    if (id_domain==1.or.id_domain==2.or.id_domain==3) then

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend

            p = PPBCK/PPREF
            t = TTUNI_MP/TTREF
          tp_pos = BLK%YCRD(i,j,k) - (-0.5d0)
          tp_ppf = 4.0d0*tp_pos*(1.0d0-tp_pos)     ! Poiseuille flow
          
            u = UUUNI_MP/UUREF*tp_ppf
!            u = UUUNI_MP/UUREF
            v = 0.0d0
            w = 0.0d0

            if (iact_tbl_sa>0) then
              tk = tkin_tbl_sa*VSUNI_MP/RHUNI_MP
            endif

            BLK%QCS1(i,j,k,1) = t
            BLK%QCS1(i,j,k,2) = u
            BLK%QCS1(i,j,k,3) = v
            BLK%QCS1(i,j,k,4) = w
            BLK%QCS1(i,j,k,5) = p
            if (iact_tbl_sa>0) then
              BLK%QCS1(i,j,k,neqn_tbl_sa) = tk
            endif
            
            if (iact_passive_scalar>0) then
              BLK%QCS1(i,j,k,neqn_passive_scalar) = 0.0D0
              BLK%QCS1(i,j,k,neqn_passive_scalar+1) = 0.0D0
            endif

          enddo
        enddo
      enddo

    endif


    if (id_domain==4) then

      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend

            p = PPBCK/PPREF
            t = TTUNI_SP/TTREF
            
          tp_pos = BLK%XCRD(i,j,k) - (-0.5d0)
          tp_ppf = 4.0d0*tp_pos*(1.0d0-tp_pos)     ! Poiseuille flow
            
            u = 0.0d0
            v = UUUNI_SP/UUREF *tp_ppf
!            v = UUUNI_SP/UUREF
            w = 0.0d0

            if (iact_tbl_sa>0) then
              tk = 5.0d0*VSUNI_SP/RHUNI_SP
            endif

            BLK%QCS1(i,j,k,1) = t
            BLK%QCS1(i,j,k,2) = u
            BLK%QCS1(i,j,k,3) = v
            BLK%QCS1(i,j,k,4) = w
            BLK%QCS1(i,j,k,5) = p
            if (iact_tbl_sa>0) then
              BLK%QCS1(i,j,k,neqn_tbl_sa) = tk
            endif

            if (iact_passive_scalar>0) then
              BLK%QCS1(i,j,k,neqn_passive_scalar) = SCALR_SP
              BLK%QCS1(i,j,k,neqn_passive_scalar+1) = 0.0D0
            endif

          enddo
        enddo
      enddo

    endif

    do l=1,neqns
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            BLK%QCS2(i,j,k,l) = BLK%QCS1(i,j,k,l)
            BLK%QCS3(i,j,k,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

  end subroutine init_tpipe

end module mdl_init_tpipe
