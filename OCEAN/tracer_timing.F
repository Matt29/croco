!======================================================================
! Optional wall-clock timing for the 3-D tracer kernels.
!
! The storage is indexed by OpenMP thread.  A thread therefore only
! updates its own cache line and no critical section is needed in the
! time-stepping path.  Enable with CPP key TRACER_TIMING.
!======================================================================
#include "cppdefs.h"

#ifdef TRACER_TIMING
      subroutine tracer_timing_init
      implicit none
# include "param.h"
# include "tracer_timing.h"
      integer slot,trd
      real*8 tracer_seconds(8,tracer_timing_nslots,0:NPP-1)
      integer tracer_calls(16,tracer_timing_nslots,0:NPP-1)
      common /tracer_timing_data/ tracer_seconds,tracer_calls
      save /tracer_timing_data/

      do trd=0,NPP-1
        do slot=1,tracer_timing_nslots
          tracer_seconds(1,slot,trd)=0.D0
          tracer_calls(1,slot,trd)=0
        enddo
      enddo
      return
      end


      subroutine tracer_timing_start(clock_count,clock_rate,clock_max)
      implicit none
      integer*8 clock_count,clock_rate,clock_max

      call system_clock(clock_count,clock_rate,clock_max)
      return
      end


      subroutine tracer_timing_add(slot,trd,clock_start,clock_rate,
     &                             clock_max)
      implicit none
# include "param.h"
# include "tracer_timing.h"
      integer slot,trd
      integer*8 clock_start,clock_rate,clock_max,clock_stop
      real*8 elapsed
      real*8 tracer_seconds(8,tracer_timing_nslots,0:NPP-1)
      integer tracer_calls(16,tracer_timing_nslots,0:NPP-1)
      common /tracer_timing_data/ tracer_seconds,tracer_calls
      save /tracer_timing_data/

      call system_clock(clock_stop,clock_rate,clock_max)
      if (clock_rate.le.0) return
      if (clock_stop.ge.clock_start) then
        elapsed=dble(clock_stop-clock_start)/dble(clock_rate)
      else
        elapsed=(dble(clock_max)-dble(clock_start)+dble(clock_stop)
     &                                        +1.D0)/dble(clock_rate)
      endif
!
!     Each OpenMP worker owns one value of trd.  The padded first
!     dimensions put data from different threads on separate cache lines.
!
      tracer_seconds(1,slot,trd)=tracer_seconds(1,slot,trd)+elapsed
      tracer_calls(1,slot,trd)=tracer_calls(1,slot,trd)+1
      return
      end


      subroutine tracer_timing_report
      implicit none
# include "param.h"
# include "tracer_timing.h"
# include "scalars.h"
# ifdef MPI
#  include "ncscrum.h"
# endif
      integer slot,trd,calls,active_threads
      real*8 tracer_seconds(8,tracer_timing_nslots,0:NPP-1)
      real*8 total_seconds,max_seconds,average_seconds
      real*8 seconds_per_call,max_over_average
      integer tracer_calls(16,tracer_timing_nslots,0:NPP-1)
      common /tracer_timing_data/ tracer_seconds,tracer_calls
      save /tracer_timing_data/

      do slot=1,tracer_timing_nslots
        total_seconds=0.D0
        max_seconds=0.D0
        calls=0
        active_threads=0
        do trd=0,NPP-1
          total_seconds=total_seconds+tracer_seconds(1,slot,trd)
          calls=calls+tracer_calls(1,slot,trd)
          if (tracer_calls(1,slot,trd).gt.0) then
            active_threads=active_threads+1
            max_seconds=max(max_seconds,tracer_seconds(1,slot,trd))
          endif
        enddo
        average_seconds=0.D0
        seconds_per_call=0.D0
        max_over_average=0.D0
        if (active_threads.gt.0) then
          average_seconds=total_seconds/dble(active_threads)
          if (average_seconds.gt.0.D0) then
            max_over_average=max_seconds/average_seconds
          endif
        endif
        if (calls.gt.0) then
          seconds_per_call=total_seconds/dble(calls)
        endif
# ifdef MPI
        write(stdout,1000)
     &    'TRACER_TIMING rank=',mynode,'kernel=',slot,
     &    'threads=',active_threads,'calls=',calls,
     &    'thread_seconds=',total_seconds,'max_seconds=',max_seconds,
     &    'avg_seconds=',average_seconds,
     &    'seconds_per_call=',seconds_per_call,
     &    'max_over_avg=',max_over_average
# else
        write(stdout,1001)
     &    'TRACER_TIMING kernel=',slot,
     &    'threads=',active_threads,'calls=',calls,
     &    'thread_seconds=',total_seconds,'max_seconds=',max_seconds,
     &    'avg_seconds=',average_seconds,
     &    'seconds_per_call=',seconds_per_call,
     &    'max_over_avg=',max_over_average
# endif
      enddo
      write(stdout,'(A)')
     &  'TRACER_TIMING kernels: 1=pre_step3d 2=step3d_t 3=t3dmix'
 1000 format(A,I5,1X,A,I2,1X,A,I3,1X,A,I10,1X,A,F12.5,
     &       1X,A,F12.5,1X,A,F12.5,1X,A,ES12.4,1X,A,F8.4)
 1001 format(A,I2,1X,A,I3,1X,A,I10,1X,A,F12.5,
     &       1X,A,F12.5,1X,A,F12.5,1X,A,ES12.4,1X,A,F8.4)
      return
      end
#endif
