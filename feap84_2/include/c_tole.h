
                   ! CONTACT TOLERANCES COMMON for contact pair

      real*8          zero,tlipen,tlopen,tlouts
      common /c_tole/
     &  zero,      ! Tol. to replace zero
     &  tlipen,    ! Tol. Initial PENetration
     &  tlopen,    ! Tol. OPENing of gap
     &  tlouts     ! Tol. OUT of Segment for slavenode
