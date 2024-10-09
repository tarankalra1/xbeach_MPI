module wave_stationary_directions_module
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
   ! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
   ! Jaap van Thiel de Vries, Robert McCall                                  !
   !                                                                         !
   ! d.roelvink@unesco-ihe.org                                               !
   ! UNESCO-IHE Institute for Water Education                                !
   ! P.O. Box 3015                                                           !
   ! 2601 DA Delft                                                           !
   ! The Netherlands                                                         !
   !                                                                         !
   ! This library is free software; you can redistribute it and/or           !
   ! modify it under the terms of the GNU Lesser General Public              !
   ! License as published by the Free Software Foundation; either            !
   ! version 2.1 of the License, or (at your option) any later version.      !
   !                                                                         !
   ! This library is distributed in the hope that it will be useful,         !
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
   ! Lesser General Public License for more details.                         !
   !                                                                         !
   ! You should have received a copy of the GNU Lesser General Public        !
   ! License along with this library; if not, write to the Free Software     !
   ! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
   ! USA                                                                     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none
   save

contains

   subroutine wave_stationary_directions(s,par,callType)
      use params
      use spaceparams
      use roelvink_module
      use wave_functions_module
      use xmpi_module
      use logging_module
      use paramsconst

      IMPLICIT NONE

      type(spacepars), target     :: s
      type(parameters)            :: par
      integer                     :: callType

      integer,parameter           :: callTypeStationary = 0
      integer,parameter           :: callTypeDirections = 1

      integer                     :: i,imax,i1
      integer                     :: j
      integer                     :: itheta,iter,itermax
      integer                     :: ntheta_local
      integer,save                :: scheme_local,scheme_local_yadvec
      real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy
      real*8 , dimension(:,:),allocatable,save    :: uorb
      real*8 , dimension(:,:)  ,allocatable,save  :: sinh2kh
      real*8 , dimension(:,:)  ,allocatable,save  :: hhwlocal
      real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr,dder
      real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
      real*8 , dimension(:,:,:),allocatable,save  :: ee_local
      real*8 , dimension(:),allocatable,save      :: Hprev,thetaprev,gammafac_i,RH
      real*8 , dimension(:,:),allocatable,save    :: eeprev, Herr_field
      integer , dimension(:,:),allocatable,save   :: wetmask
      real*8                                      :: Herr,thetaerr,thetaerr_cyc,dtw
      real*8                                      :: dtheta_local
      logical                                     :: stopiterate
      !
      ! set dimension sizes
      select case (callType)
       case (callTypeStationary)
         ntheta_local = s%ntheta
       case (callTypeDirections)
         ntheta_local = s%ntheta_s
      end select
      !
      ! Allocate
      if (.not. allocated(xadvec)) then
         allocate(xadvec    (s%nx+1,s%ny+1,ntheta_local))
         allocate(yadvec    (s%nx+1,s%ny+1,ntheta_local))
         allocate(thetaadvec(s%nx+1,s%ny+1,ntheta_local))
         allocate(dd        (s%nx+1,s%ny+1,ntheta_local))
         allocate(dder      (s%nx+1,s%ny+1,ntheta_local))
         allocate(ee_local  (s%nx+1,s%ny+1,ntheta_local))
         allocate(dhdx(s%nx+1,s%ny+1))
         allocate(dhdy(s%nx+1,s%ny+1))
         allocate(dudx(s%nx+1,s%ny+1))
         allocate(dudy(s%nx+1,s%ny+1))
         allocate(dvdx(s%nx+1,s%ny+1))
         allocate(dvdy(s%nx+1,s%ny+1))
         allocate(uorb(s%nx+1,s%ny+1))
         allocate(sinh2kh(s%nx+1,s%ny+1))
         allocate(hhwlocal(s%nx+1,s%ny+1))
         allocate(gammafac_i(s%ny+1))
         allocate(RH(s%ny+1))
         allocate(Hprev(s%ny+1))
         allocate(eeprev(s%ny+1,ntheta_local))
         allocate(Herr_field(s%ny+1,ntheta_local))
         allocate(wetmask(s%ny+1, ntheta_local))
         allocate(thetaprev(s%ny+1))
         !
         if(par%scheme==SCHEME_WARMBEAM) then
            scheme_local = SCHEME_UPWIND_2 ! note, this is already stated as warning in log file during read of params.txt
         else
            scheme_local = par%scheme
         endif
         if(scheme_local .ne. SCHEME_UPWIND_1) then
            scheme_local_yadvec = SCHEME_UPWIND_1
            call writelog('lws','(a)', 'Warning: y-advection in stationary wave solver set to upwind_1 for convergence')
         else
            scheme_local_yadvec = scheme_local
         endif
      endif
      !
      ! Allocate additional variables for stationary wave computation
      if ((.not. allocated(xradvec)) .and. (callType==callTypeStationary)) then
         allocate(xradvec(s%nx+1,s%ny+1,ntheta_local))
         allocate(yradvec(s%nx+1,s%ny+1,ntheta_local))
         allocate(thetaradvec(s%nx+1,s%ny+1,ntheta_local))
         allocate(drr (s%nx+1,s%ny+1,ntheta_local))
      endif
      !
      ! Initialize
      xadvec      = 0.0d0
      yadvec      = 0.0d0
      thetaadvec  = 0.0d0
      dd          = 0.0d0
      dder        = 0.0d0
      dhdx        = 0.0d0
      dhdy        = 0.0d0
      dudx        = 0.0d0
      dudy        = 0.0d0
      dvdx        = 0.0d0
      dvdy        = 0.0d0
      sinh2kh     = 0.0d0
      uorb        = 0.0d0
      !
      if (callType==callTypeStationary) then
         xradvec     = 0.0d0
         yradvec     = 0.0d0
         thetaradvec = 0.0d0
         drr         = 0.0d0
      endif
      !
      ! Set switch for pointer variables (not strictly defined as Fortran pointers)
      select case (callType)
       case (callTypeStationary)
         hhwlocal = s%hhw
         ee_local = s%ee
         dtheta_local = s%dtheta
       case (callTypeDirections)
         hhwlocal = s%hhws
         ee_local = s%ee_s
         dtheta_local = s%dtheta_s
      end select
      !
      ! Slopes of water depth
      call slope2D(hhwlocal,s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy,s%wete)
      !
      ! Dano: limit slopes used in refraction to avoid unrealistic refraction speeds
      dhdx=sign(1.d0,dhdx)*min(abs(dhdx),0.1d0)
      dhdy=sign(1.d0,dhdy)*min(abs(dhdy),0.1d0)
      !
      ! slopes of flow (in case of WCI)
      if (par%wci==1) then
         select case (callType)
          case(callTypeStationary)
            call slope2D(s%u,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy,s%wete)
            call slope2D(s%v,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy,s%wete)
          case (callTypeDirections)
            call slope2D(s%uws,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy,s%wete)
            call slope2D(s%vws,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy,s%wete)
         end select
      else
         dudx = 0.d0
         dudy = 0.d0
         dvdx = 0.d0
         dvdy = 0.d0
      endif
      ! wwvv these slope routines are in wave_timestep, and are
      !   MPI-aware
      !
      ! Calculate once sinh(2kh)
      where(s%wete==1 .and. 2*hhwlocal*s%k<=3000.d0)
         sinh2kh=sinh(min(2*s%k*hhwlocal,10.0d0))
      elsewhere
         sinh2kh = 3000.d0
      endwhere
      !
      ! all dry cells have zero energy
      do itheta=1,ntheta_local
         where(s%wete==0)
            ee_local(:,:,itheta) = 0.d0
         endwhere
      enddo
      where(s%wete==0)
         s%E=0.d0
         s%H=0.d0
      endwhere
      !
      ! wave directions
      select case (callType)
       case(callTypeStationary)
         if (par%snells==0) then
            forall (i=1:s%nx+1,j=1:s%ny+1,s%wete(i,j)==1)
               s%thetamean(i,j)=(sum(ee_local(i,j,:)*s%thet(i,j,:))/ntheta_local) / &
               (max(sum(ee_local(i,j,:)),0.00001d0)/ntheta_local)
            endforall
         else
            s%thetamean=asin(sin(s%theta0-s%alfaz(1,1))*s%c/s%c(1,1))+s%alfaz(1,1)
            s%costh(:,:,1)=cos(s%thetamean-s%alfaz)
            s%sinth(:,:,1)=sin(s%thetamean-s%alfaz)
         endif
       case(callTypeDirections)
         forall (i=1:s%nx+1,j=1:s%ny+1,s%wete(i,j)==1)
            s%thetamean(i,j)=(sum(ee_local(i,j,:)*s%thet_s(i,j,:))/ntheta_local) / &
            (max(sum(ee_local(i,j,:)),0.00001d0)/ntheta_local)
         endforall
      end select
      !
      ! split wave velocities in wave grid directions theta
      select case (callType)
       case(callTypeStationary)
         call compute_wave_direction_velocities(s,par,0,dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh)
       case(callTypeDirections)
         call compute_wave_direction_velocities(s,par,2,dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh)
      end select
      !
      ! Initialize at the boundary
      s%E(1,:)=max(sum(ee_local(1,:,:),2)*dtheta_local,0.d0)
      s%H(1,:)=sqrt(s%E(1,:)/par%rhog8)
      if (callType==callTypeStationary) s%R(1,:)=max(sum(s%rr(1,:,:),2)*dtheta_local,0.0d0)
      !
      !
      ! Initialize loop over grid
      !
      ! write to screen that waves are updated
      if(xmaster) then
         select case (callType)
          case(callTypeStationary)
            call writelog('ls','(a,f0.2,a)','Computing wave transformation at t = ',par%t,' s')
          case(callTypeDirections)
            call writelog('ls','(a,f0.2,a)','Computing wave directions at t = ',par%t,' s')
         end select
      endif
      call progress_indicator(.true.,0.d0,5.d0,2.d0)
      imax=s%nx
      itermax = 0
      !
      !
      ! Start of loop of the grid
      do i=2,imax
         call progress_indicator(.false.,dble(i)/imax*100,5.d0,2.d0)
         !
         ! determine internal time step
         select case (callType)
          case(callTypeStationary)
            dtw=        .5d0*minval(s%dsu(i-1:i+1,jmin_ee:jmax_ee))/max(1.0d-6,maxval(abs(s%cgx(i-1:i+1,jmin_ee:jmax_ee,:))))
            dtw=min(dtw,.5d0*minval(s%dnv(i,jmin_ee:jmax_ee))      /max(1.0d-6,maxval(abs(s%cgy(i,jmin_ee:jmax_ee,:)      ))))
            dtw=min(dtw,.5d0*dtheta_local                          /max(1.0d-6,maxval(abs(s%ctheta(i,jmin_ee:jmax_ee,:)   ))))
          case(callTypeDirections)
            dtw=        .5d0*minval(s%dsu(i-1:i+1,jmin_ee:jmax_ee))/max(1.0d-6,maxval(abs(s%cgx_s(i-1:i+1,jmin_ee:jmax_ee,:))))
            dtw=min(dtw,.5d0*minval(s%dnv(i,jmin_ee:jmax_ee))      /max(1.0d-6,maxval(abs(s%cgy_s(i,jmin_ee:jmax_ee,:)      ))))
            dtw=min(dtw,.5d0*dtheta_local                          /max(1.0d-6,maxval(abs(s%ctheta_s(i,jmin_ee:jmax_ee,:)   ))))
         end select

         !
         !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
         call xmpi_allreduce(dtw,MPI_MIN)
#endif
         !
         ! Reduced internal time step in case of WCI
         if(par%wci==1) then
            dtw = min(dtw,par%dt)
         endif
         !
         ! Initialize iteration loop per gridline
         Herr=huge(1.d0)
         thetaerr = 2*par%px
         iter=0
         stopiterate=.false.
         !
         ! Start of iteration loop per gridline
         do while (stopiterate .eqv. .false.)
            iter=iter+1
            Hprev=s%H(i,:)
            eeprev=ee_local(i,:,:)
            Herr_field = 0
            wetmask = 0
            thetaprev = s%thetamean(i,:)
            !
            ! transform energy to wave action
            i1=max(i-2,1)
            do itheta=1,ntheta_local
               where(s%wete(i1:i+1,:)==1)
                  ee_local(i1:i+1,:,itheta) = ee_local(i1:i+1,:,itheta)/s%sigm(i1:i+1,:)
               endwhere
            enddo
            !
            ! Upwind Euler timestep propagation
            !
            ! M-Direction
            select case (callType)
             case(callTypeStationary)
               if (i==2) then
                  call advecxho(ee_local(i-1:i+1,:,:),s%cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
                  2,s%ny,ntheta_local,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1, &
                  s%wete(i-1:i+1,:),dtw,s%dsz)
               else
                  ! note: at the moment Warm-Beam is currently not a valid combination (automatic change to upwind_2 in params.F90)
                  call advecxho(ee_local(i-2:i+1,:,:),s%cgx(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
                  3,s%ny,ntheta_local,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),scheme_local, &
                  s%wete(i-2:i+1,:),dtw,s%dsz)
               endif
             case(callTypeDirections)
               if (i==2) then
                  call advecxho(ee_local(i-1:i+1,:,:),s%cgx_s(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
                  2,s%ny,ntheta_local,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1, &
                  s%wete(i-1:i+1,:),dtw,s%dsz)
               else
                  ! note: at the moment Warm-Beam is currently not a valid combination (automatic change to upwind_2 in params.F90)
                  call advecxho(ee_local(i-2:i+1,:,:),s%cgx_s(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
                  3,s%ny,ntheta_local,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),scheme_local, &
                  s%wete(i-2:i+1,:),dtw,s%dsz)
               endif
            end select
            !
            ! N-Direction
            if (s%ny>0) then
               select case (callType)
                case(callTypeStationary)
                  call advecyho(ee_local(i,:,:),s%cgy(i,:,:),yadvec(i,:,:),  &
                  0,s%ny,ntheta_local,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),scheme_local_yadvec,s%wete(i,:),dtw,s%dnz)
                case(callTypeDirections)
                  call advecyho(ee_local(i,:,:),s%cgy_s(i,:,:),yadvec(i,:,:),  &
                  0,s%ny,ntheta_local,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),scheme_local_yadvec,s%wete(i,:),dtw,s%dnz)
               end select
            endif
            !
            ! Theta-Direction
            select case (callType)
             case(callTypeStationary)
               call advecthetaho(ee_local(i,:,:),s%ctheta(i,:,:),thetaadvec(i,:,:),0,s%ny,ntheta_local, &
               dtheta_local,scheme_local,s%wete(i,:))
             case(callTypeDirections)
               call advecthetaho(ee_local(i,:,:),s%ctheta_s(i,:,:),thetaadvec(i,:,:),0,s%ny,ntheta_local, &
               dtheta_local,scheme_local,s%wete(i,:))
            end select
            !
            ! update wave action
            ee_local(i,:,:)=ee_local(i,:,:)-dtw*(xadvec(i,:,:) + yadvec(i,:,:) + thetaadvec(i,:,:))
            !
            ! transform back to wave energy
            do itheta=1,ntheta_local
               where(s%wete(i1:i+1,:)==1)
                  ee_local(i1:i+1,:,itheta) = ee_local(i1:i+1,:,itheta)*s%sigm(i1:i+1,:)
                  ee_local(i1:i+1,:,itheta)=max(ee_local(i1:i+1,:,itheta),0.0d0)
               endwhere
            enddo
#ifdef USEMPI
            call xmpi_shift(ee_local(i-1:i+1,:,:),SHIFT_Y_R,1,2)
            call xmpi_shift(ee_local(i-1:i+1,:,:),SHIFT_Y_L,3,4)
#endif
            !
            ! Energy integrated over wave directions,Hrms
            where(s%wete(i,:)==1)
               s%E(i,:)=sum(ee_local(i,:,:),2)*dtheta_local
               s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            endwhere
            !
            ! adjust wave energy and height for gammax
            gammafac_i = (par%gammax*hhwlocal(i,:))/max(s%H(i,:),par%eps/1000) ! inverse of gammafac = s%H(i,:)/(par%gammax*hhwlocal(i,:))
            do j=1,s%ny+1
               if (s%wete(i,j)==1 .and. gammafac_i(j)<1.d0) then
                  do itheta=1,ntheta_local
                     ee_local(i,j,itheta)=ee_local(i,j,itheta)*gammafac_i(j)**2
                  enddo
                  s%H(i,j)=s%H(i,j)*gammafac_i(j)
                  s%E(i,j)=s%E(i,j)*gammafac_i(j)**2
               endif
            enddo
            !
            ! redo wave directions (Dano: not for Snellius; Robert: also always for callTypeDirections)
            select case (callType)
             case(callTypeStationary)
               if (par%snells==0) then
                  where(s%wete(i,:)==1)
                     s%thetamean(i,:) = (sum(ee_local(i,:,:)*s%thet(i,:,:),2)/ntheta_local)/ &
                     (max(sum(ee_local(i,:,:),2),0.000010d0)/ntheta_local)
                  endwhere
               endif
             case(callTypeDirections)
               where(s%wete(i,:)==1)
                  s%thetamean(i,:) = (sum(ee_local(i,:,:)*s%thet_s(i,:,:),2)/ntheta_local)/ &
                  (max(sum(ee_local(i,:,:),2),0.000010d0)/ntheta_local)
               endwhere
            end select
            !
            !
            ! Compute wave dissipation
            !
            ! Dissipation by breaking
            select case(par%break)
             case(BREAK_ROELVINK1,BREAK_ROELVINK2,BREAK_ROELVINK_DALY)
               call roelvink       (par,s,i)
             case(BREAK_BALDOCK)
               call baldock        (par,s,i)
             case(BREAK_JANSSEN)
               call janssen_battjes(par,s,i)
            end select
            !
            ! Dissipation by bed friction
            where(s%fw(i,:)>0.d0 .and. s%wete(i,:)==1 .and. hhwlocal(i,:)<=par%fwcutoff)
               uorb(i,:)=par%px*s%H(i,:)/par%Trep/sinh(min(max(s%k(i,:),0.01d0)*hhwlocal(i,:),10.0d0))
               s%Df(i,:)=0.28d0*par%rho*s%fw(i,:)*uorb(i,:)**3 ! Robert: note in wave-instationary the coefficient is 2/(3pi)??
            elsewhere
               s%Df(i,:) = 0.d0
            end where
            !
            ! Dissipation by vegetation computed in vegetation module
            !
            !
            ! Distribution of dissipation over directions and frequencies
            do itheta=1,ntheta_local
               where(s%wete(i,:)==1)
                  ! just dissipation due to breaking, used for roller energy balance (in callTypeStationary only)
                  dder(i,:,itheta)=ee_local(i,:,itheta)*s%D(i,:)/max(s%E(i,:),0.00001d0)
                  ! Then all short wave energy dissipation, including bed friction and vegetation, used to reduce wave height
                  ! Robert: note Dveg needs update for callTypeDirections, because based on only instantaneous H
                  dd(i,:,itheta)=dder(i,:,itheta) + ee_local(i,:,itheta)*(s%Df(i,:)+s%Dveg(i,:))/max(s%E(i,:),0.00001d0)
               elsewhere
                  dder(i,:,itheta)=0.d0
                  dd(i,:,itheta)=0.d0
               endwhere
            end do
            !
            !
            ! Calculate roller energy balance (only needed for callTypeStationary)
            if (callType==callTypeStationary .and. par%roller==1) then
               !
               ! M-direction
               if (i==2) then
                  call advecxho(s%rr(i-1:i+1,:,:),s%cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),    &
                  2,s%ny,ntheta_local,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1, &
                  s%wete(i-1:i+1,:),dtw,s%dsz)
               else
                  ! note: at the moment Warm-Beam is currently not a valid combination (automatic change to upwind_2 in params.F90)
                  call advecxho(s%rr(i-2:i+1,:,:),s%cx(i-2:i+1,:,:),xradvec(i-2:i+1,:,:),    &
                  3,s%ny,ntheta_local,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),scheme_local, &
                  s%wete(i-2:i+1,:),dtw,s%dsz)
               endif
               !
               ! N-Direction
               if (s%ny>0) then
                  call advecyho(s%rr(i,:,:),s%cy(i,:,:),yradvec(i,:,:),   &
                  0,s%ny,ntheta_local,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),scheme_local_yadvec,s%wete(i,:),dtw,s%dnz)
               endif
               !
               ! Theta-Direction
               call advecthetaho(s%rr(i,:,:),s%ctheta(i,:,:),thetaradvec(i,:,:), &
               0,s%ny,ntheta_local,dtheta_local,scheme_local,s%wete(i,:))

               s%rr(i,:,:)=s%rr(i,:,:)-dtw*(xradvec(i,:,:)+yradvec(i,:,:)+thetaradvec(i,:,:))
               s%rr(i,:,:)=max(s%rr(i,:,:),0.0d0)
#ifdef USEMPI
               call xmpi_shift(s%rr(i-1:i,:,:),SHIFT_Y_R,1,2)
               call xmpi_shift(s%rr(i-1:i,:,:),SHIFT_Y_L,3,4)
#endif
            endif
            !
            !
            ! Compute energy source and sink function
            !
            ! update internal timestep
            do j=jmin_ee,jmax_ee
               if (s%wete(i,j)==1) then
                  do itheta=1,ntheta_local
                     if(ee_local(i,j,itheta)>0.d0 .and. dd(i,j,itheta)>0.d0) then
                        dtw=min(dtw,.5*ee_local(i,j,itheta)/dd(i,j,itheta))
                     endif
                  enddo
               endif
            enddo
            ! Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
            call xmpi_allreduce(dtw,MPI_MIN)
#endif
            !
            ! Timestep update of source and sink
            do j=jmin_ee,jmax_ee
               if(s%wete(i,j)==1) then
                  do itheta=1,ntheta_local
                     ee_local(i,j,itheta)=ee_local(i,j,itheta)-dtw*dd(i,j,itheta)
                     ee_local(i,j,itheta)=max(ee_local(i,j,itheta),0.0d0)
                     if(callType==callTypeStationary) then
                        if (par%roller==1) then  !Christophe
                           drr(i,j,itheta) = 2*par%g*par%beta*max(s%rr(i,j,itheta),0.0d0)/&
                           sqrt(s%cx(i,j,itheta)**2 +s%cy(i,j,itheta)**2)
                           s%rr(i,j,itheta) = s%rr(i,j,itheta)+dtw*dder(i,j,itheta)-dtw*drr(i,j,itheta)
                           s%rr(i,j,itheta) = max(s%rr(i,j,itheta),0.0d0)
                        else
                           drr(i,j,itheta)= 0.0d0
                           s%rr(i,j,itheta)= 0.0d0
                        end if
                     endif
                  enddo
               else
                  ee_local(i,j,:)=0.0d0
                  if(callType==callTypeStationary) then
                     s%rr(i,j,:)=0.0d0
                     drr(i,j,:)=0.0d0
                  endif
               endif
            enddo
            !
            ! Lateral boundary conditions
            if (xmpi_isleft .and. s%ny>0) then
               do itheta=1,ntheta_local
                  ee_local(i,1,itheta)=ee_local(i,2,itheta)
                  if(callType==callTypeStationary) s%rr(i,1,itheta)=s%rr(i,2,itheta)
               enddo
               s%k(:,1)=s%k(:,2)
               s%sigm(:,1)=s%sigm(:,2)
            endif
            if (xmpi_isright .and. s%ny>0) then
               do itheta=1,ntheta_local
                  ee_local(i,s%ny+1,itheta)=ee_local(i,s%ny,itheta)
                  if(callType==callTypeStationary) s%rr(i,s%ny+1,itheta)=s%rr(i,s%ny,itheta)
               end do
               s%k(:,s%ny+1)=s%k(:,s%ny)
               s%sigm(:,s%ny+1)=s%sigm(:,s%ny)
            endif
            !
            !
            ! Redo wave directions (Dano: not for Snellius; Robert: also always for callTypeDirections)
            select case (callType)
             case(callTypeStationary)
               if (par%snells==0) then
                  where(s%wete(i,:)==1)
                     s%thetamean(i,:) = (sum(ee_local(i,:,:)*s%thet(i,:,:),2)/ntheta_local)/ &
                     (max(sum(ee_local(i,:,:),2),0.000010d0)/ntheta_local)
                  endwhere
               endif
             case(callTypeDirections)
               where(s%wete(i,:)==1)
                  s%thetamean(i,:) = (sum(ee_local(i,:,:)*s%thet_s(i,:,:),2)/ntheta_local)/ &
                  (max(sum(ee_local(i,:,:),2),0.000010d0)/ntheta_local)
               endwhere
            end select
            !
            ! forward-copy wave directions on dry cells in case they flood before the next wave update
            where(s%wete(i,:)==0)
               s%thetamean(i,:) = s%thetamean(i-1,:)
            endwhere
            !
            !
            ! Energy integrated over wave directions,Hrms
            where(s%wete(i,:)==1)
               s%E(i,:)=sum(ee_local(i,:,:),2)*dtheta_local
               s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            elsewhere
               s%E(i,:)=0.d0
               s%H(i,:)=0.d0
            endwhere
            !
            ! Roller energy
            if (callType==callTypeStationary) then
               where(s%wete(i,:)==1)
                  s%R(i,:)=sum(s%rr(i,:,:),2)*dtheta_local
                  s%DR(i,:)=sum(drr(i,:,:),2)*dtheta_local
               elsewhere
                  s%R(i,:)=0.d0
                  s%DR(i,:)=0.d0
               endwhere
            endif
            !
            ! Correct roller for gammax
            ! adjust wave energy and height for gammax
            if(callType==callTypeStationary .and. par%rollergammax==1) then
               RH = sqrt(s%R(i,:)/par%rhog8)
               gammafac_i = (par%gammax*hhwlocal(i,:))/max(RH,par%eps/1000) ! inverse of gammafac = s%H(i,:)/(par%gammax*hhwlocal(i,:))
               do j=1,s%ny+1
                  if (s%wete(i,j)==1 .and. gammafac_i(j)<1.d0) then
                     do itheta=1,ntheta_local
                        s%rr(i,j,itheta)=s%rr(i,j,itheta)*gammafac_i(j)**2
                     enddo
                     s%R(i,j)=RH(j)*gammafac_i(j)
                  endif
               enddo
            endif
            !
            !
            ! Wave height and direction difference with last iteration
            select case (par%truncationType)
             case (TRUNCATIONTYPE_ABSOLUTE)
               Herr=maxval(abs(Hprev(jmin_ee:jmax_ee)-s%H(i,jmin_ee:jmax_ee)),mask=s%wete(i,jmin_ee:jmax_ee)==1)
             case (TRUNCATIONTYPE_RELATIVE)
               Herr=maxval(abs(Hprev(jmin_ee:jmax_ee)-s%H(i,jmin_ee:jmax_ee))/s%H(i,jmin_ee:jmax_ee),mask=s%wete(i,jmin_ee:jmax_ee)==1)
             case (TRUNCATIONTYPE_ABSOLUTE_EE)
               do itheta=1,ntheta_local
                  wetmask(jmin_ee:jmax_ee,itheta) = s%wete(i,jmin_ee:jmax_ee)
               enddo
               Herr_field(jmin_ee:jmax_ee,:) = sqrt(abs(eeprev(jmin_ee:jmax_ee,:)-ee_local(i,jmin_ee:jmax_ee,:))/par%rhog8)
               Herr = 0.d0
               do itheta=1,ntheta_local
                  Herr = max(Herr,maxval(Herr_field(jmin_ee:jmax_ee,itheta)))
               enddo
             case (TRUNCATIONTYPE_RELATIVE_EE)
               do itheta=1,ntheta_local
                  wetmask(jmin_ee:jmax_ee,itheta) = s%wete(i,jmin_ee:jmax_ee)
               enddo
               Herr_field(jmin_ee:jmax_ee,:) = abs(eeprev(jmin_ee:jmax_ee,:)-ee_local(i,jmin_ee:jmax_ee,:))/max(ee_local(i,jmin_ee:jmax_ee,:),0.00001d0)
               Herr = 0.d0
               do itheta=1,ntheta_local
                  Herr = max(Herr,maxval(Herr_field(jmin_ee:jmax_ee,itheta)))
               enddo
            end select
            thetaerr = 0.d0
            do j=jmin_ee,jmax_ee
               thetaerr_cyc = minval( (/ abs(thetaprev(j)-(s%thetamean(i,j)-2*par%px)) , &
               abs(thetaprev(j)-(s%thetamean(i,j)         )) , &
               abs(thetaprev(j)-(s%thetamean(i,j)+2*par%px)) /) )
               thetaerr = max(thetaerr,thetaerr_cyc)
            enddo
            !
            ! communicate error across MPI domains
#ifdef USEMPI
            call xmpi_allreduce(Herr,MPI_MAX)
            call xmpi_allreduce(thetaerr,MPI_MAX)
#endif
            !
            !
            ! Stopping criteria
            if (iter<par%maxiter) then
               if (Herr<par%maxerror .and. thetaerr<par%maxerror_angle) then
                  stopiterate=.true.
               endif
            else
               stopiterate=.true.
               if(xmaster) then
                  if (Herr>=par%maxerror) then
                     call writelog('lsw','(a,i4,a,i4,a,f5.4)','Wave propagation row ',i,', iteration ',iter,', H error: ',min(Herr,0.99))
                  endif
                  if (thetaerr>=par%maxerror_angle) then
                     call writelog('lsw','(a,i4,a,i4,a,f5.4)','Wave propagation row ',i,', iteration ', &
                     iter,', dir error: ',min(thetaerr/par%px*180,0.99))
                  endif
               endif
            endif
            !
         enddo ! end interation loop per gridline
         !
         ! Store maximum iterations update
         itermax = max(itermax,iter)
      enddo ! end loop over grid 2,nx
      !
      !
      ! write summary to log
      call writelog('ls','(a,i4)','Maximum number of iterations: ',itermax)
      !
      !
      ! store local energy "pointer" variable to s structure
      select case (callType)
       case(callTypeStationary)
         s%ee= ee_local
       case(callTypeDirections)
         s%ee_s = ee_local
      end select
      !
      !
      ! Boundary conditions at nx+1
      s%E(s%nx+1,:)    = s%E(s%nx,:)
      s%H(s%nx+1,:)    = s%H(s%nx,:)
      s%k(s%nx+1,:)    = s%k(s%nx,:)
      s%sigm(s%nx+1,:) = s%sigm(s%nx,:)
      s%cg(s%nx+1,:)   = s%cg(s%nx,:)
      s%c(s%nx+1,:)    = s%c(s%nx,:)
      select case (callType)
       case(callTypeStationary)
         s%ee(s%nx+1,:,:) = s%ee(s%nx,:,:)
         s%rr(s%nx+1,:,:) = s%rr(s%nx,:,:)
         s%R(s%nx+1,:)    = s%R(s%nx,:)
         s%DR(s%nx+1,:)   = s%DR(s%nx,:)
         s%thet(s%nx+1,:,:) = s%thet(s%nx,:,:)
       case(callTypeDirections)
         s%ee_s(s%nx+1,:,:) = s%ee_s(s%nx,:,:)
         s%thet_s(s%nx+1,:,:) = s%thet_s(s%nx,:,:)
      end select
      !
      !
      ! Compute radiation stress and forcing terms (callType=callTypeStationary only)
      if (callType==callTypeStationary) then
         !
         ! wave forces
         call compute_wave_forces(s)
         !
         ! orbital velocity
         uorb=par%px*s%H/par%Trep/sinh(min(max(s%k,0.01d0)*s%hhw,10.0d0))!%ML this caused dzb changes that are strange
         s%urms=uorb/sqrt(2.d0)
         !
         ! Stokes drift and orbital velocities
         call compute_stokes_drift(s,par)
      endif
      !
   end subroutine wave_stationary_directions
end module wave_stationary_directions_module
