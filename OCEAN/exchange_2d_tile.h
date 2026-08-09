#if (!defined MP_3PTS) && (!defined MP_1PTS)
      subroutine exchange_2d_tile (Istr,Iend,Jstr,Jend, A)
#elif defined MP_1PTS
      subroutine exchange_2d_1pts_tile (Istr,Iend,Jstr,Jend, A)
#else
      subroutine exchange_2d_3pts_tile (Istr,Iend,Jstr,Jend, A)
#endif
!
! Set periodic boundary conditions (if any) for a two-dimensional
! field A of ZETA-, U-, V- or PSI-type. This file is designed to
! generate four different subroutines, by redefining (via CPP) the
! name of the subroutine exchange_2d_tile above and the starting
! indices ISTART = [Istr for U-,PSI-type; IstrR for V-,ZETA-type] 
! and JSTART = [Jstr for V-,PSI-type; JstrR for U-,ZETA-type]
! below. See also mounting file exchange.F  
!
      implicit none
#include "param.h"
#include "scalars.h"
      integer Npts,ipts,jpts
#if (!defined MP_3PTS) && (!defined MP_1PTS)
      parameter (Npts=2)
#elif defined MP_1PTS
      parameter (Npts=1)
# else
      parameter (Npts=3)
# endif
      real A(GLOBAL_2D_ARRAY)
      integer Istr,Iend,Jstr,Jend, i,j
#if defined MPI && defined OPENMP
      real hybrid_A2d(GLOBAL_2D_ARRAY)
      integer hI0,hI1,hJ0,hJ1
      common /hybrid_exchange_2d/ hybrid_A2d
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
          do j=J_RANGE
            do ipts=1,Npts
              A(Lm+ipts,j)=A(ipts,j)
	    enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_RANGE
            do ipts=1,Npts
              A(ipts-Npts,j)=A(Lm+ipts-Npts,j)
            enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef J_RANGE
#endif

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
          do i=I_RANGE
	    do jpts=1,Npts
              A(i,Mm+jpts)=A(i,jpts)
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_RANGE
	    do jpts=1,Npts
              A(i,jpts-Npts)=A(i,Mm+jpts-Npts)
	    enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
# ifdef MPI
      if (NP_XI.eq.1 .and. NP_ETA.eq.1) then
# endif
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
	  do jpts=1,Npts
	    do ipts=1,Npts
	      A(Lm+ipts,Mm+jpts)=A(ipts,jpts)
	    enddo
	  enddo
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          do jpts=1,Npts
	    do ipts=1,Npts
	      A(ipts-Npts,Mm+jpts)=A(Lm+ipts-Npts,jpts)
	    enddo
	  enddo
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
	  do jpts=1,Npts
	    do ipts=1,Npts
	      A(Lm+ipts,jpts-Npts)=A(ipts,Mm+jpts-Npts)
	    enddo
	  enddo
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
	  do jpts=1,Npts
	    do ipts=1,Npts
	      A(ipts-Npts,jpts-Npts)=A(Lm+ipts-Npts,Mm+jpts-Npts)
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
#  if (!defined MP_3PTS) && (!defined MP_1PTS)
        call MessPass2D_tile (Istr,Iend,Jstr,Jend, A)
#  elif defined MP_1PTS
        call MessPass2D_1pts_tile (Istr,Iend,Jstr,Jend, A)
#  else
        call MessPass2D_3pts_tile (Istr,Iend,Jstr,Jend, A)
#  endif
      else
!
! Assemble the rank boundary in shared storage.  A can be either a
! shared model field or a thread-private scratch field, so passing the
! master thread's A directly to MessPass is not correct.
!
      hI0=Istr
      hI1=Iend
      hJ0=Jstr
      hJ1=Jend
      if (WESTERN_EDGE) hI0=Istr-Npts
      if (EASTERN_EDGE) hI1=Iend+Npts
      if (SOUTHERN_EDGE) hJ0=Jstr-Npts
      if (NORTHERN_EDGE) hJ1=Jend+Npts

      if (Istr.eq.1) then
        do j=hJ0,hJ1
          do i=1,min(Npts,Iend)
            hybrid_A2d(i,j)=A(i,j)
          enddo
        enddo
      endif
      if (Iend.eq.Lmmpi) then
        do j=hJ0,hJ1
          do i=max(Istr,Lmmpi-Npts+1),Lmmpi
            hybrid_A2d(i,j)=A(i,j)
          enddo
        enddo
      endif
      if (Jstr.eq.1) then
        do j=1,min(Npts,Jend)
          do i=hI0,hI1
            hybrid_A2d(i,j)=A(i,j)
          enddo
        enddo
      endif
      if (Jend.eq.Mmmpi) then
        do j=max(Jstr,Mmmpi-Npts+1),Mmmpi
          do i=hI0,hI1
            hybrid_A2d(i,j)=A(i,j)
          enddo
        enddo
      endif
C$OMP BARRIER
C$OMP MASTER
# endif
#if (!defined MP_3PTS) && (!defined MP_1PTS)
# ifdef OPENMP
      call MessPass2D_tile (1,Lmmpi,1,Mmmpi, hybrid_A2d)
# else
      call MessPass2D_tile (Istr,Iend,Jstr,Jend,  A)
# endif
#elif defined MP_1PTS
# ifdef OPENMP
      call MessPass2D_1pts_tile (1,Lmmpi,1,Mmmpi, hybrid_A2d)
# else
      call MessPass2D_1pts_tile (Istr,Iend,Jstr,Jend,  A)
# endif
# else
# ifdef OPENMP
      call MessPass2D_3pts_tile (1,Lmmpi,1,Mmmpi, hybrid_A2d)
# else
      call MessPass2D_3pts_tile (Istr,Iend,Jstr,Jend,  A)
# endif
# endif
# ifdef OPENMP
C$OMP END MASTER
C$OMP BARRIER
!
! Copy received rank halos back to every participating field.  The
! Adjacent tiles can copy the same corner values; all sources are the
! already-completed shared staging field.
!
      hI0=max(1-Npts,Istr-Npts)
      hI1=min(Lmmpi+Npts,Iend+Npts)
      hJ0=max(1-Npts,Jstr-Npts)
      hJ1=min(Mmmpi+Npts,Jend+Npts)
      if (Istr.eq.1 .and. WEST_INTER2) then
        do j=hJ0,hJ1
          do i=1-Npts,0
            A(i,j)=hybrid_A2d(i,j)
          enddo
        enddo
      endif
      if (Iend.eq.Lmmpi .and. EAST_INTER2) then
        do j=hJ0,hJ1
          do i=Lmmpi+1,Lmmpi+Npts
            A(i,j)=hybrid_A2d(i,j)
          enddo
        enddo
      endif
      if (Jstr.eq.1 .and. SOUTH_INTER2) then
        do j=1-Npts,0
          do i=hI0,hI1
            A(i,j)=hybrid_A2d(i,j)
          enddo
        enddo
      endif
      if (Jend.eq.Mmmpi .and. NORTH_INTER2) then
        do j=Mmmpi+1,Mmmpi+Npts
          do i=hI0,hI1
            A(i,j)=hybrid_A2d(i,j)
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
!$acc update host(A(:,Jstr:Jstr+Npts-1))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update host(A(:,Jend-Npts+1:Jend))
      endif
C$OMP BARRIER
      if (.not.SOUTHERN_EDGE) then
!$acc update device(A(:,Jstr-Npts:Jstr-1))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update device(A(:,Jend+1:Jend+Npts))
      endif
#endif
           
      return
      end

#if (!defined MP_3PTS) && (!defined MP_1PTS)
#  define MP_3PTS
#  include "exchange_2d_tile.h"
#  undef MP_3PTS
#ifndef MP_1PTS
#define MP_1PTS
#include "exchange_2d_tile.h"
#undef MP_1PTS
#endif
# endif
