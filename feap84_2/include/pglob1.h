
      real*8          units
      integer                  g2type,gdtype,gtdof
      common /pglob1/ units(3),g2type,gdtype,gtdof

      logical         gtypfl,gdeffl,gomgfl,gtdofl
      common /pglob1/ gtypfl,gdeffl,gomgfl,gtdofl

      real*8          gomega   ,gomex   ,gomev
      common /pglob2/ gomega(3),gomex(3),gomev(3)

      real*8          gquadn,gaugm
      common /pglob2/ gquadn,gaugm

      integer         geqnum,gneq,gpart
      common /pglob3/ geqnum,gneq,gpart
