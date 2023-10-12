module mdl_schs_source
  use mdl_input_manager
  use mdl_mpisub_sbsp
  use mdl_decompo
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_fdm1
  use mdl_halo
  use mdl_steps
  use mdl_table
  implicit none

  real(8), public, parameter :: SCHS1 = 1.0d3

  real(8), save, public, allocatable, dimension(:) :: SCHS_PKT
  real(8), save, public, allocatable, dimension(:) :: SCHS_PKK
  
contains
!-------------------------------------------------------------------------------
  subroutine schs_source(BLK,RHS)
    type(block), intent(inout) :: BLK
    real(8), intent(inout), dimension(:,:,:,:) :: RHS
    
    integer :: i,j,k
    integer :: mode
    integer :: ivalid
    real(8) :: cfdim, srct, cklno, dckin
    integer :: n, iFLG_SCHS
    real(8) :: logkw,logkh
    real(8),parameter :: logrf=86.0d0,gasco=8.3144626d0,psi=-450000d0
    real(8) :: k_ref,k_int,ohmns,rho,tmp,t_ref
    real(8) :: a1_ne,b1_ne,c1_ne,d1_ne,e1_ne,f1_ne,g1_ne
    real(8) :: a2_ne,b2_ne,c2_ne,d2_ne,e2_ne,f2_ne,g2_ne

    
    ! call look_up_table_calc_func(TBL_TP_to_DNST,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
    ! BLK%DNST,RHREF,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

    !mode=0 is directly getting knet from LES code
    !mode=1 is getting knet fron Ni.txt
    mode=0
    if(mode==0) then 

      call look_up_table_calc_func(TBL_TP_to_DIELEC,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
      BLK%DIEL,1.0d0,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)

      a1_ne = -4.098d0
      b1_ne = -3245.2d0
      c1_ne = 223630.0d0
      d1_ne = -39984000.0d0
      e1_ne = 13.957d0
      f1_ne = -1262.3d0
      g1_ne = 856410.0d0

      a2_ne = -4.8d0
      b2_ne = 22.4d0
      c2_ne = -14.8d0

      k_ref=86.0d0
      t_ref=380d0+273.15d0

  !! FURUKAWA
  !! MEMO--------------
  !! This roop is not vectorized due to the inner loop. Need the implementation for vectorizing.
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend

            rho=BLK%DNST(i,j,k)*RHREF*1.0D-3
            tmp=BLK%TMPR(i,j,k)*TTREF
            logkh=a2_ne*rho**2+b2_ne*rho+c2_ne
            logkw=(a1_ne+b1_ne/tmp+c1_ne/(tmp**2)+d1_ne/(tmp**3)) &
                +(e1_ne+f1_ne/tmp+g1_ne/(tmp**2))*DLOG10(rho)
            ohmns=10**logkw/dsqrt(10**logkh)
            k_int=dexp((BLK%DIEL(i,j,k)-1)/(2*BLK%DIEL(i,j,k)+1)*psi/(gasco*t_ref)+k_ref)
            BLK%KNET(i,j,k)=k_int*(ohmns**2)
            !k_net=8.22d0

            srct = BLK%KNET(i,j,k)*BLK%DNST(i,j,k)*RHREF*BLK%YSPC(i,j,k,1)
            cfdim= XYREF/(RHREF*UUREF)
            
            RHS(i,j,k,neqn_passive_scalar) = RHS(i,j,k,neqn_passive_scalar) - BLK%RJCB(i,j,k)*cfdim*srct
            RHS(i,j,k,neqn_passive_scalar+1) = RHS(i,j,k,neqn_passive_scalar+1) + BLK%RJCB(i,j,k)*cfdim*srct
          enddo
        enddo
      enddo
    else if(mode==1) then
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend

            iFLG_SCHS = 0
            if(SCHS_PKT(7)<=BLK%TMPR(i,j,k)*TTREF) then
              BLK%KNET(i,j,k) = exp(SCHS_PKK(7))
            else if(SCHS_PKT(1)>=BLK%TMPR(i,j,k)*TTREF) then
              BLK%KNET(i,j,k) = exp(SCHS_PKK(1))
            else 
              do n=2,7
                if(SCHS_PKT(n)>BLK%TMPR(i,j,k)*TTREF.and.iFLG_SCHS==0) then
                  dckin = (SCHS_PKK(n)-SCHS_PKK(n-1))/(SCHS_PKT(n)-SCHS_PKT(n-1))
                  cklno = SCHS_PKK(n-1) + dckin*(BLK%TMPR(i,j,k)*TTREF-SCHS_PKT(n-1))
                  BLK%KNET(i,j,k) = exp(cklno)
                  iFLG_SCHS = 1
                end if
              end do
            end if

            srct = BLK%KNET(i,j,k)*BLK%DNST(i,j,k)*RHREF*BLK%YSPC(i,j,k,1)
            cfdim= XYREF/(RHREF*UUREF)
            
            RHS(i,j,k,neqn_passive_scalar) = RHS(i,j,k,neqn_passive_scalar) - BLK%RJCB(i,j,k)*cfdim*srct
            RHS(i,j,k,neqn_passive_scalar+1) = RHS(i,j,k,neqn_passive_scalar+1) + BLK%RJCB(i,j,k)*cfdim*srct
          enddo
        enddo
      enddo
    end if

  end subroutine schs_source

  subroutine schs_init
    integer :: i
    integer :: myrank,nprocs,iproc
    integer :: ierror

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
#else
    myrank = 0
    nprocs = 1
#endif
    allocate(SCHS_PKT(7))
    allocate(SCHS_PKK(7))
    
    do iproc=0,nprocs-1
      if (myrank==iproc) then
        
        open(11,file='./data/reaction_rate/Ni.txt',form='formatted')
         read(11,*) (SCHS_PKT(i),i=1,7)
         read(11,*) (SCHS_PKK(i),i=1,7)
      endif
    end do
    
  end subroutine schs_init

  subroutine schs_final
    deallocate(SCHS_PKT)
    deallocate(SCHS_PKK)
  end subroutine schs_final

end module mdl_schs_source
