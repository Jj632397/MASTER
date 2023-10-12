program main
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_steps
  use mdl_readgrid
  use mdl_halo
  use mdl_table
  use mdl_metric
  use mdl_initmp
  use mdl_result
  use mdl_time_delta
  use mdl_dfw
  use mdl_output
  use mdl_restart
  use mdl_explicit
  use mdl_implicit
  use mdl_exchange_domain_halo
  use mdl_init_tpipe
  use mdl_mean_and_rms
  use mdl_output_tpipe
  use mdl_schs_source  !! FURUSAWA
  implicit none

  integer :: ierror
  integer :: myrank,id_block
  integer :: isend,nstep
  integer :: imax,jmax,kmax
  integer :: mode_time_integral
  type(block) :: BLK

#ifdef _USE_MPI
  call MPI_Init(ierror)
#endif

  call input_manager('data/input.txt')

  call decompo_init

  call mpisub_sbsp_init

  call mpisub_sbsp_get_myrank_world(myrank)

  call Get_BlockID_From_Rank(myrank,id_block)

  call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

  call block_allocate(BLK,imax,jmax,kmax)

  call exchange_domain_halo_init

  call halo_init

  call table_init

  call schs_init !! FURUSAWA

  call initmp

  call readgrid(BLK)

  call metric(BLK)

  call halo_exchange_metric(BLK)

  call distance_from_wall(BLK)

  call init_tpipe(BLK)

  call restart(0,BLK)

  call result(BLK)

  call bound(BLK)

  call result(BLK)

  call preconditioning_velocity_parameter(BLK)

  call time_delta(BLK)

  if (myrank==0) then
    open(32,file='data/cond.txt',form='formatted')
    write(32,2000) XYREF
    write(32,2001) UUREF
    write(32,2002) RHREF
    write(32,2003) TTREF
    write(32,2004) PPREF
    write(32,2005) VSREF
    write(32,2006) CPREF
    write(32,2007) DTREF
    write(32,2008) RYNLD
    close(32)
  endif

  call steps_init

  call output_grid(BLK)
  call output_grid_planes(BLK)

  call output_grid_block_with_halo(BLK)

  call output_xsix(BLK)

  call output_xsix_block_with_halo(BLK)

  call output_phys(BLK)

  call output_phys_block_with_halo(BLK)

  call input_time_integral(mode_time_integral)

  !Main loop
  do while (1.eq.1)

    if (mode_time_integral==0) then
      call explicit(BLK)
    else if (mode_time_integral==1) then
      call implicit(BLK)
    endif

    !****************************************************
    !Increment step count (Physical time iteration)     *
    call steps_increment !                              *
    !Output subroutine can be called hereafter.         *
    !****************************************************

    call steps_get_step_now(nstep)

    if (mode_time_integral==0) then
      if (myrank==0) then
        open(33,file='data/progress.txt',form='formatted')
        write(33,*) nstep
        close(33)
      endif
    endif

    call mean_and_rms(BLK)

    call output_tpipe_inflow_outflow(BLK)

    if (nstep.eq.10) then
      call output_phys(BLK)
    endif

    if (nstep.eq.100) then
      call output_phys(BLK)
    endif

!    if (nstep.eq.12600) then
!      call output_phys(BLK)
!    endif

    if ((mod(nstep,1000).eq.0).and.(nstep.gt.0)) then
      call output_phys(BLK)
      call output_phys_planes(BLK)
    endif

    call output_l2norm(BLK)

    if ((mod(nstep,20000).eq.0).and.(nstep.gt.0)) then
      call restart(1,BLK)
    endif

    call steps_isend(isend)
    if (isend.eq.1) exit
  enddo

  call steps_final

  call table_final

  call halo_final

  call schs_final !! FURUSAWA

  call block_deallocate(BLK)

  call exchange_domain_halo_final
  call mpisub_sbsp_final
  call decompo_final

#ifdef _USE_MPI
  call MPI_Finalize(ierror)
#endif

2000  format('XYREF      = ',e13.6,'[m       ]')
2001  format('UUREF      = ',e13.6,'[m/s     ]')
2002  format('RHREF      = ',e13.6,'[kg/m^3  ]')
2003  format('TTREF      = ',e13.6,'[K       ]')
2004  format('PPREF      = ',e13.6,'[Pa      ]')
2005  format('VSREF      = ',e13.6,'[Pa*s    ]')
2006  format('CPREF      = ',e13.6,'[J/(kg*K)]')
2007  format('DTREF      = ',e13.6,'[s       ]')
2008  format('Re         = ',e13.6)
2009  format('Mach       = ',e13.6)
end program main
