#ifndef MP_3PTS
      subroutine exchange_3d_tile (Istr,Iend,Jstr,Jend, A)
#else
      subroutine exchange_3d_3pts_tile (Istr,Iend,Jstr,Jend, A)
#endif
!
! Set periodic boundary conditions (if any) for a three-dimensional
! field A of RHO-, U-, V- or PSI-type. This file is designed to
! generate five different subroutines, by redefining (via CPP) the
! name of the subroutine exchange_2d_tile above and the starting
! indices ISTART = [Istr for U-,PSI-type; IstrR for V-,RHO-type]
! and JSTART = [Jstr for V-,PSI-type; JstrR for U-,RHO-type] below,
! as well as macro KSTART for the vertical RHO- and W-types. See
! also mounting file exchange.F
!
      implicit none
#include "param.h"
#include "scalars.h"
      integer Npts,ipts,jpts
#ifndef MP_3PTS
      parameter (Npts=2)
#else
      parameter (Npts=3)
#endif
      real A(GLOBAL_2D_ARRAY,KSTART:N)
      integer Istr,Iend,Jstr,Jend, i,j,k
#if defined MPI && defined OPENMP
      real hybrid_A3d(GLOBAL_2D_ARRAY,0:N)
      integer hI0,hI1,hJ0,hJ1
      common /hybrid_exchange_3d/ hybrid_A3d
#endif
!
#include "compute_auxiliary_bounds.h"
!
!$acc kernels if(compute_on_device) default(present)  
#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE JSTART,JendR
# endif
# ifdef MPI
      if (NP_XI.eq.1) then
# endif
        if (WESTERN_EDGE) then
          do k=KSTART,N
            do j=J_RANGE
              do ipts=1,Npts
                A(Lm+ipts,j,k)=A(ipts,j,k)
	      enddo
            enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do k=KSTART,N
            do j=J_RANGE
	      do ipts=1,Npts
                A(ipts-Npts,j,k)=A(Lm+ipts-Npts,j,k)
              enddo
            enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef J_RANGE
#endif            /* EW_PERIODIC */

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE ISTART,IendR
# endif
# ifdef MPI
      if (NP_ETA.eq.1) then
# endif
        if (SOUTHERN_EDGE) then
          do k=KSTART,N
            do i=I_RANGE
	      do jpts=1,Npts
                A(i,Mm+jpts,k)=A(i,jpts,k)
	      enddo
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do k=KSTART,N
            do i=I_RANGE
	      do jpts=1,Npts
                A(i,jpts-Npts,k)=A(i,Mm+jpts-Npts,k)
	      enddo
            enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef I_RANGE
#endif               /* NS_PERIODIC */

#if defined EW_PERIODIC && defined NS_PERIODIC
# ifdef MPI
      if (NP_XI.eq.1 .and. NP_ETA.eq.1) then
# endif
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
          do k=KSTART,N
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(Lm+ipts,Mm+jpts,k)=A(ipts,jpts,k)
	      enddo
	    enddo
          enddo
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          do k=KSTART,N
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(ipts-Npts,Mm+jpts,k)=A(Lm+ipts-Npts,jpts,k)
	      enddo
	    enddo
          enddo
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
          do k=KSTART,N
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(Lm+ipts,jpts-Npts,k)=A(ipts,Mm+jpts-Npts,k)
	      enddo
	    enddo
          enddo
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
          do k=KSTART,N
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(ipts-Npts,jpts-Npts,k)=A(Lm+ipts-Npts,Mm+jpts-Npts,k)
	      enddo
	    enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
#endif
!$acc end kernels		   
#ifdef MPI
# ifdef OPENMP
      if (hybrid_master_exchange) then
        k=N-KSTART+1
#  ifndef MP_3PTS
#   ifndef MP_M3FAST_SEDLAYERS
        call MessPass3D_tile (Istr,Iend,Jstr,Jend, A,k)
#   else
        call MessPass3D_sl_tile (Istr,Iend,Jstr,Jend, A,k)
#   endif
#  else
#   ifndef MP_M3FAST_SEDLAYERS
        call MessPass3D_3pts_tile (Istr,Iend,Jstr,Jend, A,k)
#   else
        call MessPass3D_3pts_sl_tile (Istr,Iend,Jstr,Jend, A,k)
#   endif
#  endif
      else
      hI0=Istr
      hI1=Iend
      hJ0=Jstr
      hJ1=Jend
      if (WESTERN_EDGE) hI0=Istr-Npts
      if (EASTERN_EDGE) hI1=Iend+Npts
      if (SOUTHERN_EDGE) hJ0=Jstr-Npts
      if (NORTHERN_EDGE) hJ1=Jend+Npts

      if (Istr.eq.1) then
        do k=KSTART,N
          do j=hJ0,hJ1
            do i=1,min(Npts,Iend)
              hybrid_A3d(i,j,k)=A(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Iend.eq.Lmmpi) then
        do k=KSTART,N
          do j=hJ0,hJ1
            do i=max(Istr,Lmmpi-Npts+1),Lmmpi
              hybrid_A3d(i,j,k)=A(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Jstr.eq.1) then
        do k=KSTART,N
          do j=1,min(Npts,Jend)
            do i=hI0,hI1
              hybrid_A3d(i,j,k)=A(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Jend.eq.Mmmpi) then
        do k=KSTART,N
          do j=max(Jstr,Mmmpi-Npts+1),Mmmpi
            do i=hI0,hI1
              hybrid_A3d(i,j,k)=A(i,j,k)
            enddo
          enddo
        enddo
      endif
C$OMP BARRIER
C$OMP MASTER
# endif
      k=N-KSTART+1
# ifndef MP_3PTS
#   ifndef MP_M3FAST_SEDLAYERS
#    ifdef OPENMP
      call MessPass3D_tile (1,Lmmpi,1,Mmmpi,
     &                     hybrid_A3d(START_2D_ARRAY,KSTART),k)
#    else
      call MessPass3D_tile (Istr,Iend,Jstr,Jend,  A,k)
#    endif
#   else
#    ifdef OPENMP
      call MessPass3D_sl_tile (1,Lmmpi,1,Mmmpi,
     &                     hybrid_A3d(START_2D_ARRAY,KSTART),k)
#    else
      call MessPass3D_sl_tile (Istr,Iend,Jstr,Jend,  A,k)
#    endif
#   endif      
#  else
!!  MP_3PTS        
#   ifndef MP_M3FAST_SEDLAYERS       
#    ifdef OPENMP
      call MessPass3D_3pts_tile (1,Lmmpi,1,Mmmpi,
     &                     hybrid_A3d(START_2D_ARRAY,KSTART),k)
#    else
      call MessPass3D_3pts_tile (Istr,Iend,Jstr,Jend,  A,k)
#    endif
#   else
#    ifdef OPENMP
      call MessPass3D_3pts_sl_tile (1,Lmmpi,1,Mmmpi,
     &                     hybrid_A3d(START_2D_ARRAY,KSTART),k)
#    else
      call MessPass3D_3pts_sl_tile (Istr,Iend,Jstr,Jend,  A,k)      
#    endif
#   endif
#endif   
# ifdef OPENMP
C$OMP END MASTER
C$OMP BARRIER
      hI0=max(1-Npts,Istr-Npts)
      hI1=min(Lmmpi+Npts,Iend+Npts)
      hJ0=max(1-Npts,Jstr-Npts)
      hJ1=min(Mmmpi+Npts,Jend+Npts)
      if (Istr.eq.1 .and. WEST_INTER2) then
        do k=KSTART,N
          do j=hJ0,hJ1
            do i=1-Npts,0
              A(i,j,k)=hybrid_A3d(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Iend.eq.Lmmpi .and. EAST_INTER2) then
        do k=KSTART,N
          do j=hJ0,hJ1
            do i=Lmmpi+1,Lmmpi+Npts
              A(i,j,k)=hybrid_A3d(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Jstr.eq.1 .and. SOUTH_INTER2) then
        do k=KSTART,N
          do j=1-Npts,0
            do i=hI0,hI1
              A(i,j,k)=hybrid_A3d(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (Jend.eq.Mmmpi .and. NORTH_INTER2) then
        do k=KSTART,N
          do j=Mmmpi+1,Mmmpi+Npts
            do i=hI0,hI1
              A(i,j,k)=hybrid_A3d(i,j,k)
            enddo
          enddo
        enddo
      endif
C$OMP BARRIER
      endif
# endif
#   ifdef  BAND_DEBUG          
      chkbandname='none'
#   endif
#endif
#if defined OPENMP && defined OPENACC
      if (.not.SOUTHERN_EDGE) then
!$acc update host(A(:,Jstr:Jstr+Npts-1,KSTART:N))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update host(A(:,Jend-Npts+1:Jend,KSTART:N))
      endif
C$OMP BARRIER
      if (.not.SOUTHERN_EDGE) then
!$acc update device(A(:,Jstr-Npts:Jstr-1,KSTART:N))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update device(A(:,Jend+1:Jend+Npts,KSTART:N))
      endif
#endif
      return
      end

#ifndef MP_3PTS
# define MP_3PTS
# include "exchange_3d_tile.h"
# undef MP_3PTS
#endif
