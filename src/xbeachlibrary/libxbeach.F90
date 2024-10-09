module libxbeach_module
   use iso_c_binding,        only: c_int, c_char
   use params,               only: parameters
   use spaceparamsdef,       only: spacepars
   use timestep_module,      only: timepars
   use ship_module,          only: ship
   use vegetation_module,    only: veggie


   implicit none
   private
   public executestep, outputext, final, init, getversion
   public par, tpar, s, sglobal
#ifdef USEMPI
   public slocal
#endif
   save

   type(parameters)                     :: par
   type(timepars)                       :: tpar
   type(spacepars), pointer             :: s
   type(spacepars), target              :: sglobal
   type(ship), dimension(:), pointer    :: sh

   integer                              :: n,it,error
   real*8                               :: tbegin

#ifdef USEMPI
   type(spacepars), target              :: slocal
   real*8                               :: t0,t01
   logical                              :: toall = .true.
   logical                              :: end_program
   integer                              :: nxbak, nybak
#endif


   !startinit

   !-----------------------------------------------------------------------------!
   ! Initialize program                                                          !
   !-----------------------------------------------------------------------------!


contains
   integer(c_int) function init()
      use params,               only: params_inio, all_input
      use spaceparams,          only: space_alloc_scalars, ranges_init
      use xmpi_module
      use initialize_module,    only: drifter_init, wave_init, sed_init, flow_init, discharge_init, hotstart_init_1, hotstart_init_2
      use initialize_module,    only: setbathy_init, grid_bathy
      use readtide_module,      only: readtide
      use readwind_module,      only: readwind
      use timestep_module,      only: timestep_init
      use logging_module,       only: writelog
      use groundwaterflow,      only: gw_init
      use logging_module,       only: start_logfiles, writelog_startup
      use means_module,         only: means_init
      use output_module,        only: output_init, output
      use ship_module,          only: ship_init
      use nonh_module,          only: nonh_init
      use vegetation_module,    only: veggie_init
      use paramsconst
      use rainfall_module
#ifdef USEMPI
      use xmpi_module
      use spaceparams
      use logging_module,       only: writelog_mpi
      use params,               only: distribute_par

      integer, dimension(12)               :: info
      character(256)                       :: line
      integer                              :: rank,i
#endif

      error   = 0
      n = 0

      ! setup of MPI
#ifdef USEMPI
      s=>slocal
      call xmpi_initialize
      call xmpi_barrier(toall)
      t0 = MPI_Wtime()
#endif

      ! create log files
      call start_logfiles(error)

      ! set starting time and date
      call cpu_time(tbegin)

      ! show statup message
      call writelog_startup()

      !-----------------------------------------------------------------------------!
      ! Initialize simulation                                                       !
      !-----------------------------------------------------------------------------!

      ! initialize time counter
      it      = 0

      ! read input from params.txt
      params_inio = .false.
      call all_input(par)

      ! allocate space scalars
      call space_alloc_scalars(sglobal)
      s => sglobal

      ! read grid and bathymetry
      call grid_bathy(s,par)

      ! distribute grid over processors
#ifdef USEMPI
      call xmpi_determine_processor_grid(s%nx,s%ny,par%mpiboundary,par%mmpi,par%nmpi,par%cyclic,error)
#if 0
      ! print information about the neighbours of the processes
      info                      = 0
      info(1)                   = xmpi_orank
      info(2)                   = xmpi_rank
      info(3)                   = xmpi_prow
      info(4)                   = xmpi_pcol
      info(5)                   = xmpi_left
      info(6)                   = xmpi_right
      info(7)                   = xmpi_top
      info(8)                   = xmpi_bot
      if(xmpi_isleft)  info(9)  = 1
      if(xmpi_isright) info(10) = 1
      if(xmpi_istop)   info(11) = 1
      if(xmpi_isbot)   info(12) = 1

      do i=5,8
         if(info(i) .eq. MPI_PROC_NULL) then
            info(i) = -99
         endif
      enddo

      call writelog("ls"," "," ranks and neigbours (-99 means: no neighbour):")
      call writelog("ls"," ",' ')
      call writelog("ls"," ","  orank rank pcol prow left right  top  bot isleft isright istop isbot")

      do rank = 0,xmpi_osize-1
         if (rank .ne. xmpi_omaster) then
            call xmpi_send(rank,xmpi_imaster,info)
            if (xmaster) then
               write(line,'(i5,i5,i5,i5,i5,i6,i5,i5,i7,i8,i6,i6)') info
            endif
            call writelog("ls"," ",trim(line))
         endif
      enddo
      call writelog("ls"," ",' ')
#endif
      call writelog_mpi(par%mpiboundary,error)
#endif

      ! initialize timestep
      call timestep_init(par, tpar)

      if (xmaster) then

         call writelog('ls','','Initializing .....')
      endif
      if (xmaster) then
         if(par%hotstart==1) then
            call hotstart_init_1(s,par)
         endif
      endif
      call setbathy_init      (s,par)
      ! initialize physics
      call readtide           (s,par)
      call readwind           (s,par)
      call flow_init          (s,par)
      call discharge_init     (s,par)
      call drifter_init       (s,par)
      call wave_init          (s,par)
      call gw_init            (s,par)
      call rainfall_init      (s,par)
      ! TODO, fix ordening of arguments....
      call sed_init           (s,par)
      call ship_init          (s,par,sh)   ! always need to call initialise in order
      ! to reserve memory on MPI subprocesses.
      ! Note: if par%ships==0 then don't allocate
      ! and read stuff for sh structures
      call veggie_init         (s,par)
      !
      if (xmaster) then
         if(par%hotstart==1) then
            call hotstart_init_2(s,par)
         endif
      endif
#ifdef USEMPI
      call distribute_par(par)
      s => slocal
      !
      ! here an hack to ensure that sglobal is populated, also on
      ! the not-(o)master processes, just to get valid addresses.
      !
      if (.not. xmaster .and. .not. xomaster) then
         !nxbak = sglobal%nx
         !nybak = sglobal%ny
         !sglobal%nx=0
         !sglobal%ny=0
         call space_alloc_arrays_dummies(sglobal,par)
         !sglobal%nx = nxbak
         !sglobal%ny = nybak
      endif
      call space_distribute_space (sglobal,s,par)
#endif

      call ranges_init(s)

      ! nonh_init does not always need to be called
      if (par%wavemodel==WAVEMODEL_NONH) call nonh_init(s,par)

      ! initialize output
      call means_init             (sglobal,s,par     )

      call output_init            (sglobal,s,par,tpar)


      ! store first timestep
      ! from this point on, xomaster will hang in subroutine output
      ! until a broadcast .true. is received
      call output(sglobal,s,par,tpar)
      init = 0
   end function init

   integer(c_int) function outputext()
      use output_module,        only: output, output_error
      ! store first timestep
      call output(sglobal,s,par,tpar,.false.)
      if(error==0) then
         outputext = 0
      elseif(error==1) then
         call output_error(s,sglobal,par,tpar)
         outputext = 1
      endif

   end function outputext
   !-----------------------------------------------------------------------------!
   ! Start simulation                                                            !
   !-----------------------------------------------------------------------------!

   !_____________________________________________________________________________

   integer(c_int) function executestep(dt)
      use loopcounters_module,  only: execute_counter
      use xmpi_module,          only: xcompute
      use drifter_module,       only: drifter
      use flow_timestep_module, only: flow
      use boundaryconditions,   only: wave_bc, flow_bc
      use morphevolution,       only: bed_update, setbathy_update, transus
      use wave_timestep_module, only: wave
      use timestep_module,      only: timestep, outputtimes_update
      use groundwaterflow,      only: gw_bc, gwflow
      use ship_module,          only: shipwave
      use vegetation_module,    only: vegatt
      use wetcells_module,      only: compute_wetcells
      use output_module,        only: log_progress
      use paramsconst
      use rainfall_module
      use logging_module
#ifdef USEMPI
      use xmpi_module
#endif

      real*8, optional :: dt

#ifdef USEMPI
      if (execute_counter .eq. 1) then
         ! exclude first pass from time measurement
         call xmpi_barrier
         t01 = MPI_Wtime()
      endif
#endif
      execute_counter = execute_counter + 1

      executestep = -1

      ! determine timestep
      if(xcompute) then

         ! determine this time step's wet points
         call compute_wetcells(s,par)
         !
         ! determine time step
         call timestep(s,par,tpar,it,dt=dt,ierr=error)
         call outputtimes_update(par, tpar)

         ! update log
         call log_progress(par)

         if (error==0) then
            !
            ! Boundary conditions
            call wave_bc        (sglobal,s,par)
            if (par%gwflow==1)                                      call gw_bc          (s,par)
            if ((par%flow==1).or.(par%wavemodel==WAVEMODEL_NONH))   call flow_bc        (s,par)
            !
            ! Compute timestep
            if (par%ships==1)                                       call shipwave       (s,par,sh)
            if (par%swave==1)                                       call wave           (s,par)
            if (par%vegetation==1)                                  call vegatt         (s,par)
            if (par%gwflow==1)                                      call gwflow         (s,par)
            if ((par%flow==1).or.(par%wavemodel==WAVEMODEL_NONH))   call flow           (s,par)
            if (par%ndrifter>0)                                     call drifter        (s,par)
            if (par%sedtrans==1)                                    call transus        (s,par)
            !
            ! Bed level update
            if ((par%morphology==1).and.(.not. par%setbathy==1)) call bed_update(s,par)
            if (par%setbathy==1)     call setbathy_update(s, par)
         endif
      endif

      n = n + 1
      executestep = 0
   end function executestep
   !_____________________________________________________________________________


   integer(c_int) function final()
      use logging_module,       only: writelog_finalize
      use xmpi_module


      !-----------------------------------------------------------------------------!
      ! Finalize simulation                                                         !
      !-----------------------------------------------------------------------------!

#ifdef USEMPI
      end_program = .true.
      call xmpi_send_sleep(xmpi_imaster,xmpi_omaster) ! wake up omaster
      call xmpi_bcast(end_program,toall)
      call xmpi_barrier(toall)
      call writelog_finalize(tbegin,n,par%t,par%nx,par%ny,t0,t01)
      call xmpi_finalize
#else
      call writelog_finalize(tbegin,n,par%t,par%nx,par%ny)
#endif
      final = 0
   end function final

   subroutine getversion(version)
      character(kind=c_char,len=*),intent(inout) :: version

      include 'version.def'
      include 'version.dat'

      version = trim(Build_Revision)
   end subroutine

end module libxbeach_module
