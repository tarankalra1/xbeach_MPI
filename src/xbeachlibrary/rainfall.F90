module rainfall_module
   implicit none
   save

   logical,save,private                     :: constantRainfall

contains

   subroutine rainfall_init(s,par)
      use params
      use spaceparams
      use filefunctions
      use logging_module
      use interp

      type(spacepars)                     :: s
      type(parameters)                    :: par

      integer                             :: i,j,it
      integer                             :: ier,fid,dummy

      if(.not. xmaster) return

      ! allocate global rainfall rate variable
      allocate(s%rainfallrate(s%nx+1,s%ny+1))

      if (par%rainfall==1) then
         ! read from file?
         if (par%rainfallratefile==' ') then
            constantRainfall = .true.
            allocate(s%rainfallinput(0,0,0))
            allocate(s%trainfallinput(0))
            s%rainfallrate = par%rainfallrate/1000.d0/3600.d0  ! convert from [mm/hr] to [m/s]
         else
            constantRainfall = .false.
            allocate(s%rainfallinput(s%nx+1,s%ny+1,par%nrainfallrate))
            allocate(s%trainfallinput(par%nrainfallrate))
            ! start file read
            fid = create_new_fid()
            ! call check_file_exist(par%rainfallratefile) ! already done at all_input subroutine level
            open (fid,file=par%rainfallratefile)
            do it=1,par%nrainfallrate
               read(fid,*,iostat=ier)s%trainfallinput(it)
               if (ier .ne. 0) then
                  call report_file_read_error(par%rainfallratefile)
               endif
               do j=1,s%ny+1
                  read(fid,*,iostat=ier)(s%rainfallinput(i,j,it),i=1,s%nx+1)
                  if (ier .ne. 0) then
                     call report_file_read_error(par%rainfallratefile)
                  endif
               enddo
            enddo
            close(fid)
            s%rainfallinput = s%rainfallinput/1000.d0/3600.d0  ! convert source file from [mm/hr] to [m/s]
            ! Interpolate initial rainfall
            do j=1,s%ny+1
               do i=1,s%nx+1
                  call LINEAR_INTERP(s%trainfallinput,s%rainfallinput(i,j,:),par%nrainfallrate, &
                  0.d0,s%rainfallrate(i,j),dummy)
               enddo
            enddo
         endif
      else
         ! give MPI bcast a memory address
         allocate(s%rainfallinput(0,0,0))
         allocate(s%trainfallinput(0))
         constantRainfall = .true.
         s%rainfallrate = 0.d0
      endif
   end subroutine rainfall_init

   subroutine rainfall_update(s, par)

      use params
      use spaceparams
      use interp

      implicit none

      type(spacepars)                     :: s
      type(parameters)                    :: par

      integer                             :: i,j,dummy

#ifdef USEMPI
      if (par%t<=par%dt) then
         call xmpi_bcast(constantRainfall)
      endif
#endif

      if (.not. constantRainfall) then
         ! interpolate from time series
         do j=1,s%ny+1
            do i=1,s%nx+1
               call LINEAR_INTERP(s%trainfallinput,s%rainfallinput(i,j,:),par%nrainfallrate, &
               par%t,s%rainfallrate(i,j),dummy)
            enddo
         enddo
      endif

   end subroutine rainfall_update

end module rainfall_module
