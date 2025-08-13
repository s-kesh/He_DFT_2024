!--------------------------------------------------------------------
!--                      Subroutine TIMER                         ---
!--------------------------------------------------------------------
!
! This subroutine gives the time in seconds (REAL*4) measured from
! the beginning of the run.
!
! Version for GNU compilers.
!
subroutine timer(secs)
    real (kind=4) :: secs
        real (kind=8) :: t

            call cpu_time(t)
                secs = real(t, kind=4)
                    
                        return
                        end subroutine timer
