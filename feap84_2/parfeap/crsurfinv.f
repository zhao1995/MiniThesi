c$Id:$
      subroutine crsurfinv(pairn,cs0,cp0, icsn)

      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'cdata.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    pairn,fac
      real*8     cs0(nr0,n0c1:nc01,*), cp0(nr0,n0c3:nc03,*)

      integer    i, surf, ofs,neps,dnope,nope,icsn(*)

c     Loop over pairs


      do i = 1,2
        surf = nint(cp0(i+1,0,pairn))

c       Find information

        ofs  = nint(cs0(2,-1,surf))-1
        neps = nint(cs0(3,-1,surf))
        dnope= nint(cs0(4,-1,surf))
        nope = nint(cs0(2,0,surf))

c       Mark slave side with negative values
        if(i.eq.1) then
          fac = -1
        else
          fac = 1
        endif

        call crsetinv(pairn, mr(np(133)+ofs),dnope,nope,neps, icsn,fac)
      end do ! i



      end

      subroutine crsetinv(np, ics,dnope,nope,neps, icsinv,fac)

      implicit   none

      integer    np, dnope,nope,neps,fac
      integer    ics(dnope,neps), icsinv(*)

      integer    ns,is, nn


      do ns = 1,neps
        do is = 1,nope
          nn = ics(is,ns)
          if(nn.gt.0) then
            icsinv(nn) = np*fac
          endif
        end do ! is
      end do ! ns

      end
