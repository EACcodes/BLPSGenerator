subroutine c_wrtlog(msg)

  implicit none

 include 'mpif.h'

!#if defined HAVE_CONFIG_H
!#include "config.h"
!#endif

  character(len=*), intent(in) :: msg;
  INTEGER :: me, ierr
  INTEGER,save :: master = 0
  
  me = 0

  !Determine who I am
  call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)

      if(me/=master) RETURN

      ! only master node write to invKS_log file
      ! open file  ...
      open(unit=2222,access='append',action='write',status='unknown', file='invKS_log',form='formatted');
      
      IF (msg(1:13) .eq. 'WELCOMEMSG_BF') THEN
        WRITE(2222,*)'---------------------------------------------------------'
        write(2222,*)'      Revised version based on ABINIT v6.61              '
        write(2222,*)'                                                         '
        write(2222,*)'-    This program inverts Kohn-Sham equation for        -'
        write(2222,*)'-    a given electron density.                          -'
        write(2222,*)'-    A BFGS and parallel version                        -'
        write(2222,*)'-    Author: Chen Huang 2011                            -'
        write(2222,*)'---------------------------------------------------------'
      ELSE
        write (2222, '(A)') trim(msg);
      ENDIF
      

      CLOSE (2222)

end subroutine
