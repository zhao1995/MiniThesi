
!-----[--+---------+---------+---------+---------+---------+---------+-]
!     Modification log                                  Date(dd/mm/year)
!     1. Increase dimension to 1024 to force loops on        03/10/2011
!        long arrays

      real*8          hr
      integer                  mr
      common /comblk/ hr(1024),mr(1024)
