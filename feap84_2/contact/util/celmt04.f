c$Id:$
      subroutine celmt04 (ndm,ndf,x,u,
     &    csw,npair,cs01,cs02,cm01,cm02,cp0,ix1,ix2,cm1,cm2,ch1,ch2,ch3,
     &    w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact DRIVER # 1

c      Purpose: Management of specific contact formulation

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         csw     - Contact switch
c         npair   - # of current pair
c         cs01(*) - Contact surfaces control data for surface 1
c         cs02(*) - Contact surfaces control data for surface 2
c         cm01(*) - Contact material control data for surface 1
c         cm02(*) - Contact material control data for surface 2
c         cp0(*)  - Contactpair control data
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         cm1(*)  - Contact materials data storage for surface 1
c         cm2(*)  - Contact materials data storage for surface 2

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c                 - Data exchange with main program subroutine calls
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'c_tole.h'
      include  'iofile.h'
      include  'print.h'

      character w1(*)*(*),w3(*)*(*)
      integer   csw,npair,ix1(dnope1,*),ix2(dnope2,*),ndm,ndf
      real*8    cs01(nr0,n0c1:*),cs02(nr0,n0c1:*),cm01(nr0,n0c2:*)
      real*8    cm02(nr0,n0c2:*), cp0(nr0,n0c3:*),cm1(*),cm2(*)
      real*8    ch1(lh1,*),ch2(lh1,*),ch3(lh3,*),x(ndm,*),u(ndf,*)

      save

      call cdebug0 ('    celmt04',csw)

      if ((csw.ne.0).and.(csw.ne.200).and.(csw.ne.400)) write (*,3001)

3001  format (' WARNING - dummy contact element CELMT04 called')

      end
