c$Id:$
      subroutine mdo (n, ia,ja, max, v,l, head,last,next, mark, flag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  mdo finds a minimum degree ordering of rows and
c               columns of a symmetric matrix m stored in (ia,ja,a)
c               format.

c     Additional parameters

c        max  - declared dimension of one-dimensional arrays v and l;
c               max must be at least  n+2k,  where k is number of
c               nonzeroes in strict upper triangle of m

c        v    - integer one-dimensional work array;  dimension = max

c        l    - integer one-dimensional work array;  dimension = max

c        head - integer one-dimensional work array;  dimension = n

c        last - integer one-dimensional array used to return permutation
c               of rows and columns of m corresponding to minimum
c               degree ordering;  dimension = n

c        next - integer one-dimensional array used to return inverse of
c               permutation returned in last;  dimension = n

c        mark - integer one-dimensional work array (may be same as v);
c               dimension = n

c        flag - integer error flag;  values and their meanings are -
c                 0      no errors detected
c                 11n+1  insufficient storage in mdo

c     Definitions of internal parameters

c     ---------+-----------------------------------------------------
c     v(s)     | value field of list entry
c     ---------+-----------------------------------------------------
c     l(s)     | link field of list entry  (0 => end of list)
c     ---------+-----------------------------------------------------
c     l(vi)    | pointer to element list of uneliminated vertex vi
c     ---------+-----------------------------------------------------
c     l(ej)    | pointer to boundary list of active element ej
c     ---------+-----------------------------------------------------
c     head(d)  | vj => vj head of d-list d
c              |  0 => no vertex in d-list d

c              |                  vi uneliminated vertex
c              |          vi in ek           |       vi not in ek
c     ---------+-----------------------------+-----------------------
c     next(vi) | undefined but nonnegative   | vj => vj next in d-list
c              |                             |  0 => vi tail of d-list
c     ---------+-----------------------------+-----------------------
c     last(vi) | (not set until mdp)         | -d => vi head of d-list d
c              |-vk => compute degree        | vj => vj last in d-list
c              | ej => vi prototype of ej    |  0 => vi not in  d-list
c              |  0 => do not compute degree |
c     ---------+-----------------------------+-----------------------
c     mark(vi) | mark(vk)                    | nonneg tag < mark(vk)

c              |                   vi eliminated vertex
c              |      ei active element      |           otherwise
c     ---------+-----------------------------+-----------------------
c     next(vi) | -j => vi was j-th vertex    | -j => vi was j-th vertex
c              |       to be eliminated      |       to be eliminated
c     ---------+-----------------------------+-----------------------
c     last(vi) |  m => size of ei = m        | undefined
c     ---------+-----------------------------+-----------------------
c     mark(vi) | -m => overlap count of ei   | undefined
c              |       with ek = m           |
c              | otherwise nonnegative tag   |
c              |       < mark(vk)            |

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,max,flag,tag,dmin,vk,ek,tail, k
      integer   ia(*),ja(*),v(*),l(*),head(*),last(*),next(*),mark(*)

      equivalence  (vk,ek)

      save

c     Initialization

      tag = 0
      call  mdi (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return

      k    = 0
      dmin = 1

c     While  k < n  do

      do while (k.lt.n)

c       Search for vertex of minimum degree

        do while(head(dmin).le.0)
          dmin = dmin + 1
        end do ! while

c       Remove vertex vk of minimum degree from degree list

        vk         = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin

c       Number vertex vk, adjust tag, and tag vk

        k        =  k + 1
        next(vk) = -k
        last(ek) =  dmin - 1
        tag      =  tag + last(ek)
        mark(vk) =  tag

c       Form element ek from uneliminated neighbors of vk

        call  mdm (vk,tail, v,l, last,next, mark)

c       Purge inactive elements and do mass elimination

        call  mdp (k,ek,tail, v,l, head,last,next, mark)

c       Update degrees of uneliminated vertices in ek

        call  mdu (ek,dmin, v,l, head,last,next, mark)

      end do ! while

c     Generate inverse permutation from permutation

      do k = 1,n
        next(k)       = -next(k)
        last(next(k)) =  k
      end do ! k

      end

      subroutine snf (n,p,ip,ia,ja,a,d,iju,ju,iu,u,il,jl, flag)

c-----[--.----+----.----+----.-----------------------------------------]
c     Definitions of internal parameters (during k-th row elimination)

c       (d(i),i=k,n) contains k-th row of u (expanded)

c       il(i) points to first nonzero element in columns k,...,n of
c         row i of u

c       jl contains lists of rows to be added to uneliminated rows --
c         i ge k => jl(i) is first row to be added to row i
c         i lt k => jl(i) is row following row i in some list of rows
c         in either case, jl(i) = 0 indicates end of a list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'print.h'

      integer   n, flag, vj, i,j,k, jmin,jmax, mu, nexti
      integer   p(*),ip(*), ia(*),ja(*), iju(*),ju(*),iu(*),il(*), jl(*)
      real*4    etime, tt, tary(2)
      real*8    dk, ukidi, a(*), d(*), u(*)

      save

c     Initialization

      do k = 1,n
        d(k)  = 0.d0
        jl(k) = 0
      end do ! k

c     For each row k

      do k = 1,n

        if(ior.lt.0 .and. mod(k,2000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3001) k,tary
        endif

c       Initialize k-th row with elements nonzero in row p(k) of m

        do j = ia(p(k)),ia(p(k)+1) - 1
          vj = ip(ja(j))
          if (k.le.vj)  d(vj) = a(j)
        end do ! j

c       Modify k-th row by adding in rows i with u(i,k) ne 0

        dk = d(k)
        i  = jl(k)
        do while (i.gt.0)
          nexti = jl(i)

c         Compute multiplier and update diagonal element

          ukidi    = -u(il(i)) * d(i)
          dk       =  dk + ukidi * u(il(i))
          u(il(i)) =  ukidi

c         Add multiple of row i to k-th row ...

          jmin = il(i)   + 1
          jmax = iu(i+1) - 1
          if (jmin.le.jmax) then
            mu = iju(i) - iu(i)
            do j = jmin,jmax
              d(ju(mu+j)) = d(ju(mu+j)) + ukidi * u(j)
            end do ! j

c           And add i to row list for next nonzero entry

            il(i) = jmin
            j     = ju(mu+jmin)
            jl(i) = jl(j)
            jl(j) = i
          end if
          i = nexti
        end do ! i

c       Check for zero pivot and save diagonal element

        if (dk.eq.0.d0)  go to 100
        d(k) = 1.d0 / dk

c       Save nonzero entries in k-th row of u ...

        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin.le.jmax) then
          mu = iju(k) - jmin
          do j = jmin,jmax
            u(j)        = d(ju(mu+j))
            d(ju(mu+j)) = 0.d0
          end do ! j

c         And add k to row list for first nonzero entry in k-th row

          il(k)           = jmin
          jl(k)           = jl(ju(mu+jmin))
          jl(ju(mu+jmin)) = k
        endif
      end do ! k

      flag = 0
      return

c     Error -- insufficient storage for u

      flag = 7*n + 1
      return

c     Error -- zero pivot

 100  flag = 8*n + k

c     Format

3001  format( ' SNF: Equation = ',i10,'  Time: CPU = ',f12.2,
     &        ' System = ',f12.2)

      end

      subroutine sdrv (n,p,ip,ia,ja,a,b,z,isp,rsp,esp,path, flag)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  SDRV solves sparse symmetric positive definite systems
c               of linear equations.  The solution process is divided
c               into three stages --

c              ssf - coefficient matrix m is factored symbolically to
c                    determine where fillin will occur during numeric
c                    factorization.

c              snf - m is factored numerically into product ut-d-u,
c                    where d is diagonal and u is unit upper triangular.

c              sns - linear system  mx = b  is solved using ut-d-u
c                    factorization from snf.

c     Description:

c       For several systems with same coefficient matrix, ssf and snf
c       need be done only once (for first system);  then sns is done
c       once for each additional right-hand side.  For several systems
c       whose coefficient matrices have same nonzero structure, ssf
c       need be done only once (for first system);  then snf and sns
c       are done once for each additional system.

c     Storage of sparse matrices

c       The nonzero entries of matrix m are stored row-by-row in the
c       array a.  To identify individual nonzero entries in each row,
c       we need to know in which column each entry lies.  These column
c       indices are stored in array ja;  i.e., if  a(k) = m(i,j),  then
c       ja(k) = j.  To identify individual rows, we need to know where
c       each row starts.  These row pointers are stored in array ia;
c       i.e., if m(i,j) is first nonzero entry (stored) in i-th row
c       and  a(k) = m(i,j), then ia(i) = k. Moreover, ia(n+1) points to
c       first location following last element in last row.
c       Thus, number of entries in i-th row is  ia(i+1) - ia(i),
c       nonzero entries in i-th row are stored consecutively in

c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),

c       and corresponding column indices are stored consecutively in

c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).

c       Since coefficient matrix is symmetric, only nonzero entries
c       in upper triangle need be stored, for example, matrix

c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )

c       could be stored as

c            | 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia | 1  4  5  8 12 14
c         ja | 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a | 1  2  3  4  2  5  6  3  6  7  8  8  9

c       or (symmetrically) as

c            | 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia | 1  4  5  7  9 10
c         ja | 1  3  4  2  3  4  4  5  5
c          a | 1  2  3  4  5  6  7  8  9          .


c     Reordering rows and columns of m

c       A symmetric permutation of rows and columns of coefficient
c       matrix m (e.g., which reduces fillin or enhances numerical
c       stability) must be specified.  The solution z is returned in the
c       original order.

c       To specify trivial ordering (i.e., identity permutation),
c       set  p(i) = ip(i) = i, i=1,...,n.  in this case, p & ip can be
c       same array.

c       If a nontrivial ordering (i.e., not identity permutation) is
c       specified & m stored symmetrically (i.e., not both m(i,j) and
c       m(j,i) are stored for i ne j), then odrv should be called (with
c       path = 3 or 5) to symmetrically reorder (ia,ja,a) before calling
c       sdrv.  This is to ensure that if m(i,j) will be in upper
c       triangle of m with respect to new ordering, then m(i,j) is
c       stored in row i (thus m(j,i) is not stored);  whereas if m(i,j)
c       will be in strict lower triangle of m, then m(j,i) is stored in
c       row j (and thus m(i,j) is not stored).


c     Parameters

c       n    - number of variables/equations

c       p    - integer one-dimensional array specifying a permutation of
c              rows and columns of m;  dimension = n

c       ip   - integer one-dimensional array containing inverse of the
c              permutation specified in p; i.e., ip(p(i)) = i, i=1,..,n;
c              dimension = n

c       ia   - integer one-dimensional array containing pointers to
c              delimit rows in ja and a;  dimension = n+1

c       ja   - integer one-dimensional array containing column indices
c              corresponding to elements of a;  dimension = number of
c              nonzero entries in m stored

c       a    - real one-dimensional array containing nonzero entries in
c              coefficient matrix m, stored by rows;  dimension =
c              number of nonzero entries in m stored

c       b    - real one-dimensional array containing right-hand side b;
c              b and z can be same array;  dimension = n

c       z    - real one-dimensional array containing solution x;  z and
c              b can be same array;  dimension = n

c       isp  - integer one-dimensional array used for working storage;
c              isp and rsp should be equivalenced;

c       rsp  - real one-dimensional array used for working storage;  rsp
c              and isp should be equivalenced;

c       esp  - integer variable;  if sufficient storage was available to
c              perform symbolic factorization (ssf), then esp is set to
c              amount of excess storage provided (negative if
c              insufficient storage was available to perform numeric
c              factorization (snf))

c       path - integer path specification; values & their meanings are:
c                1  perform ssf, snf, and sns
c                2  perform snf and sns (isp/rsp is assumed to have been
c                     set up in an earlier call to sdrv (for ssf))
c                3  perform sns only (isp/rsp is assumed to have been
c                     set up in earlier call to sdrv (for ssf and snf))
c                4  perform ssf
c                5  perform ssf and snf
c                6  perform snf only (isp/rsp is assumed to have been
c                     set up in an earlier call to sdrv (for ssf))

c       flag - integer error flag;  values and their meanings are -
c                  0     no errors detected
c                 2n+k   duplicate entry in a  --  row = k
c                 6n+k   insufficient storage in ssf  --  row = k
c                 7n+1   insufficient storage in snf
c                 8n+k   zero pivot  --  row = k
c                10n+1   insufficient storage in sdrv
c                11n+1   illegal path specification
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   n, esp,  path,  flag, il,iu,iju, jl,ju
      integer   d, u, tmp
      integer   p(*), ip(*),  ia(*), ja(*),  isp(*)
      real*8    a(*),  b(*),  z(*),  rsp(*)

      save

c     Allocate storage and factor m symbolically to determine fill-in

      iju = 1
      iu  = iju + n
      jl  = iu  + n + 1
      ju  = jl  + n

      if ((path-1) * (path-4) * (path-5) .eq. 0) then
        fp(1) = np(111)
        fp(2) = n + fp(1)

        call ssf (n, p,ip, ia,ja, isp(iju),isp(ju),isp(iu),
     &            mr(fp(1)),mr(fp(2)),isp(jl), flag)

        if (flag.ne.0)  go to 100
      endif

c     Allocate storage for real arrays

      il   = ju     + isp(iju+(n-1))
      tmp  = ((il-1)+(ipr-1)) / ipr  +  1
      d    = tmp    + n
      u    = d      + n

c     Set required storage for problem

      esp  = 3*n + 1 + isp(iju+(n-1)) + ipr*(2*n + isp(iu+n))

c     Factor numerically

      if ((path-1) * (path-2) * (path-5) * (path-6) .eq. 0) then

        call snf (n, p,ip, ia,ja,a,rsp(d), isp(iju),isp(ju),isp(iu),
     &            rsp(u), isp(il),isp(jl), flag)
        if (flag.ne.0)  go to 100
      endif

c     Solve system of linear equations  mx = b

      if ((path-1) * (path-2) * (path-3) .eq. 0) then

        call sns (n, p, rsp(d),isp(iju),isp(ju),isp(iu),rsp(u), z,b,
     &            rsp(tmp))

      endif

c     Error detected in ssf, snf, or sns

 100  return

      end

      subroutine mdm (vk,tail, v,l, last,next, mark)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Minimum degree computations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   tag,tail, s,ls,vk,vs,es, b,lb,vb, blp,blpmax
      integer   v(*),l(*), last(*),next(*), mark(*)

      equivalence  (vs, es)

      save

c     Initialize tag and list of uneliminated neighbors

      tag  = mark(vk)
      tail = vk

c     For each vertex/element vs/es in element list of vk

      ls = l(vk)
      do while (ls.gt.0)
        s  = ls
        ls = l(s)
        vs = v(s)

c       If vs is uneliminated vertex, then tag and append to list of
c          uneliminated neighbors

        if (next(vs).ge.0) then
          mark(vs) = tag
          l(tail)  = s
          tail     = s

c       If es is active element, then ...
c           for each vertex vb in boundary list of element es

        else

          lb     = l(es)
          blpmax = last(es)
          do blp = 1,blpmax
            b  = lb
            lb = l(b)
            vb = v(b)

c           If vb is untagged vertex, then tag and append to list of
c              uneliminated neighbors

            if (mark(vb).lt.tag) then
              mark(vb) = tag
              l(tail)  = b
              tail     = b
            end if
          end do ! blp

c         Mark es inactive

          mark(es) = tag
        end if
      end do ! while

c     Terminate list of uneliminated neighbors

      l(tail) = 0

      end

      subroutine mdi (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Minimum degree

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,max,tag, flag, sfs, vi,dvi, vj,j
      integer   ia(*),ja(*), v(*),l(*), head(*),last(*),next(*),mark(*)

      save

c     Initialize degrees, element lists, and degree lists

      do vi = 1,n
        mark(vi) = 1
        l(vi)    = 0
        head(vi) = 0
      end do ! vi
      sfs = n + 1

c     Create nonzero structure

c     For each nonzero entry a(vi,vj) in strict upper triangle

      do vi=1,n
        do j = ia(vi),ia(vi+1)-1
          vj = ja(j)
          if (vi.lt.vj) then
            if (sfs.ge.max)  go to 101

c           Enter vj in element list for vi

            mark(vi) = mark(vi) + 1
            v(sfs)   = vj
            l(sfs)   = l(vi)
            l(vi)    = sfs
            sfs      = sfs + 1

c           Enter vi in element list for vj

            mark(vj) = mark(vj) + 1
            v(sfs)   = vi
            l(sfs)   = l(vj)
            l(vj)    = sfs
            sfs      = sfs + 1
          endif
        end do ! j
      end do ! vi

c     Create degree lists and initialize mark vector

      do vi=1,n
        dvi = mark(vi)
        next(vi)  =  head(dvi)
        head(dvi) =  vi
        last(vi)  = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        mark(vi)  =  tag
      end do ! vi

      return

c     Error -- insufficient storage

 101  flag = 9*n + vi

      end

      subroutine mdp (k,ek,tail, v,l, head,last,next, mark)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Minimum degree numerical

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   k,ek,tail,tag,free, i,li,vi,lvi,evi,s,ls,es,ilp,ilpmax
      integer   v(*), l(*),  head(*), last(*), next(*), mark(*)

      save

c     Initialize tag

      tag = mark(ek)

c     For each vertex vi in ek

      li = ek
      ilpmax = last(ek)
      do ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)

c       Remove vi from degree list

        if (last(vi).ne.0) then
          if (last(vi).lt.0)  then
            head(-last(vi)) = next(vi)
          elseif(last(vi).gt.0) then
            next(last(vi)) = next(vi)
          endif
          if (next(vi).gt.0)  last(next(vi)) = last(vi)
        endif

c       Remove inactive items from element list of vi

        ls = vi
   4    s  = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).ge.tag) then
            free = ls
            l(s) = l(ls)
            ls   = s
          endif
        go to 4

c       If vi is interior vertex, then remove from list and eliminate

   6    lvi = l(vi)
        if (lvi.eq.0) then
          l(i)     =  l(li)
          li       =  i
          k        =  k+1
          next(vi) = -k
          last(ek) =  last(ek) - 1

c       Else ...
c         classify vertex vi

        else
          if (l(lvi).ne.0)  go to 9
          evi = v(lvi)
          if (next(evi).ge.0)  go to 9

c         If vi is prototype vertex, then mark as such, initialize
c           overlap count for corresponding element, and move vi to end
c           of boundary list

          if (mark(evi).ge.0) then
            last(vi)  =  evi
            mark(evi) = -1
            l(tail)   =  li
            tail      =  li
            l(i)      =  l(li)
            li        =  i

c         Else if vi is duplicate vertex, then mark as such and adjust
c           overlap count for corresponding element

          else
            last(vi)  =  0
            mark(evi) =  mark(evi) - 1
          endif
          go to 10

c         Else mark vi to compute degree

   9      last(vi)    = -ek

c         Insert ek in element list of vi

  10      v(free)     = ek
          l(free)     = l(vi)
          l(vi)       = free
        endif
      end do ! ilp

c     Terminate boundary list

      l(tail) = 0

      end

      subroutine sns (n, p, d, iju,ju,iu,u, z, b, tmp)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sparse solution

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'endata.h'
      include  'print.h'

      integer   n, j,k, mu, p(*),  iju(*), ju(*), iu(*)
      real*4    etime, tt, tary(2)
      real*8    tmpk, d(*), u(*),  z(*), b(*),  tmp(*)

      save

c     Set tmp to permuted b

      do k = 1,n
        tmp(k) = b(p(k))
      end do ! k

c     Solve  U^t * D * y = b  by forward substitution & compute energy

      aengy = 0.d0
      do k = 1,n
        tmpk = tmp(k)
        mu   = iju(k) - iu(k)
        do j = iu(k),iu(k+1)-1
          tmp(ju(mu+j)) = tmp(ju(mu+j)) + u(j) * tmpk
        end do ! j
        tmp(k) = tmpk * d(k)
        aengy  = aengy + tmpk*tmp(k)
        if(ior.lt.0 .and. mod(k,5000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3001) 'Fwd.',k,tary
        endif
      end do ! k

c     Solve  U x = y  by back substitution

      do k = n,1,-1
        mu = iju(k) - iu(k)
        do j = iu(k),iu(k+1) - 1
          tmp(k) = tmp(k) + u(j) * tmp(ju(mu+j))
        end do ! j
        z(p(k)) = tmp(k)
        if(ior.lt.0 .and. mod(k,5000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3001) 'Back',k,tary
        endif
      end do ! k

c     Format

3001  format( ' SNS: ',a4,'  Eq. = ',i10,'  Time: CPU = ',f12.2,
     &        ' System = ',f12.2)
      end

      subroutine mdu (ek,dmin, v,l, head,last,next, mark)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sparse solution factorizations

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ek,dmin
      integer   tag, vi,evi,dvi, s,vs,es, b,vb, i,ilp,ilpmax,blp,blpmax
      integer   v(*), l(*),  head(*), last(*), next(*), mark(*)

      equivalence  (vs, es)

      save

c     Initialize tag

      tag = mark(ek) - last(ek)

c     For each vertex vi in ek

      i      = ek
      ilpmax = last(ek)
      do 10 ilp = 1,ilpmax
        i  = l(i)
        vi = v(i)

c       If vi neither prototype nor duplicate vertex, then merge
c         elements to compute degree

        if(last(vi).ne.0) then
          if(last(vi).lt.0) then
            tag = tag + 1
            dvi = last(ek)

c           For each vertex/element vs/es in element list of vi

            s = l(vi)
   2        s = l(s)
            if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3

c           If vs is uneliminated vertex, then tag and adjust degree

            mark(vs) = tag
            dvi      = dvi + 1
            go to 5

c           If es is active element, then expand
c             check for outmatched vertex

   3        if (mark(es).lt.0)  go to 6

c           For each vertex vb in es

            b      = es
            blpmax = last(es)
            do blp = 1,blpmax
              b  = l(b)
              vb = v(b)

c             If vb is untagged, then tag and adjust degree

              if (mark(vb).lt.tag) then
                mark(vb) = tag
                dvi      = dvi + 1
              endif
            end do ! blp

   5        go to 2

c           Else if vi outmatched vertex, then adjust overlaps but do
c             not compute degree

   6        last(vi) = 0
            mark(es) = mark(es) - 1
   7        s = l(s)
            if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7

c         Else if vi is prototype vertex, then calculate degree by
c           inclusion/exclusion and reset overlap count

          elseif(last(vi).gt.0) then
            evi = last(vi)
            dvi = last(ek) + last(evi) + mark(evi)
            mark(evi) = 0
          endif

c         Insert vi in appropriate degree list

   9      next(vi)  =  head(dvi)
          head(dvi) =  vi
          last(vi)  = -dvi
          if (next(vi).gt.0)  last(next(vi)) = vi
          dmin      =  min(dmin,dvi)
        endif

  10  continue

      end

      subroutine  sro (n, ip, ia,ja,a, q, r, dflag)

c-----[--.----+----.----+----.-----------------------------------------]
c  Additional parameters

c    q     - integer one-dimensional work array;  dimension = n

c    r     - integer one-dimensional work array;  dimension = number of
c            nonzero entries in upper triangle of m

c    dflag - logical variable;  if dflag = .true., then store nonzero
c            diagonal elements at beginning of row
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   dflag
      integer   n, i,j,k, jak,jmin,jmax, ilast
      integer   ip(*),  ia(*), ja(*),  q(*), r(*)
      real*8    a(*),  ak

      save

c     Phase 1: Find row in which to store each nonzero

      do i = 1,n
        q(i) = 0
      end do ! i

c     For each nonzero element a(j)

      do i = 1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        do j = jmin,jmax

c         Find row (=r(j)) & column (=ja(j)) in which to store a(j)

          k = ja(j)
          if (ip(k).lt.ip(i)) then
            ja(j) = i
          else
            k     = i
          endif
          r(j) = k

c         Increment count of nonzeroes (=q(r(j)) in that row

          q(k) = q(k) + 1
        end do ! j
      end do ! i

c     Phase 2: Find new ia and permutation to apply to (ja,a) determine
c              pointers to delimit rows in permuted (ja,a)

      do i = 1,n
        ia(i+1) = ia(i) + q(i)
        q(i)    = ia(i+1)
      end do ! i

c     Determine where each (ja(j),a(j)) is stored in permuted (ja,a)
c       for each nonzero element (in reverse order)

      ilast = 0
      jmin  = ia(1)
      jmax  = ia(n+1) - 1
      do j = jmax,jmin,-1
        i = r(j)

c       Put (off-diagonal) nonzero in last unused location in row

        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast) then
          q(i) = q(i) - 1
          r(j) = q(i)

c       Put diagonal nonzero at beginning of row

        else
          r(j)  = ia(i)
          ilast = i
        endif

      end do ! j

c     Phase 3: Permute (ja,a) to upper triangular form
c             (wrt new ordering)

      do j = jmin,jmax
        do while(r(j).ne.j)
          k     = r(j)
          r(j)  = r(k)
          r(k)  = k

          jak   = ja(k)
          ja(k) = ja(j)
          ja(j) = jak

          ak    = a(k)
          a(k)  = a(j)
          a(j)  = ak
        end do ! while
      end do ! j

      end

      subroutine ssf (n,p,ip,ia,ja,iju,ju,iu,q,mark,jl, flag)

c-----[--.----+----.----+----.-----------------------------------------]
c     Definitions of internal parameters (during k-th row elimination)

c       q contains an ordered linked list representation of nonzero
c         structure of k-th row of u --
c           q(k) is first column with a nonzero entry
c           q(i) is next column with a nonzero entry after column i
c         in either case, q(i) = n+1 indicates end of list

c       jl contains lists of rows to be merged into uneliminated rows --
c           i ge k => jl(i) is first row to be merged into row i
c           i lt k => jl(i) is row following row i in some list of
c           rows in either case, jl(i) = 0 indicates end of a list

c       mark(i) is last row stored in ju for which u(mark(i),i) ne 0

c       jumin and juptr are indices in ju of first and last
c         elements in last row saved in ju

c       luk is number of nonzero entries in k-th row
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'print.h'

      logical   clique
      integer   p(*),ip(*), ia(*),ja(*), iju(*),ju(*),iu(*),q(*)
      integer   mark(*), jl(*), n,flag, i,j,k,m
      integer   tag, vj, qm, jmin,jmax,jumin,juptr, lui,luk,lmax
      real*4    etime, tt, tary(2)

      save

c     Initialization

      jumin = 1
      juptr = 0
      iu(1) = 1
      do k = 1,n
        mark(k) = 0
        jl(k)   = 0
      end do ! k

c     For each row k

      do k = 1,n

        if(ior.lt.0 .and. mod(k,5000).eq.0 .and. prnt) then
          tt = etime(tary)
          write(*,3001) k,tary
        endif
        luk  = 0
        q(k) = n + 1

        tag = mark(k)
        if (jl(k).ne.0) then
          clique = jl(jl(k)).eq.0
        else
          clique = .false.
        endif

c       Initialize nonzero structure of k-th row to row p(k) of m

        jmin = ia(p(k))
        jmax = ia(p(k)+1) - 1
        do j = jmin,jmax
          vj = ip(ja(j))
          if (vj.gt.k) then

            qm = k
   2        m  = qm
            qm = q(m)
            if (qm.lt.vj)  go to 2
            if (qm.eq.vj)  go to 102
            luk   = luk + 1
            q(m)  = vj
            q(vj) = qm
            if (mark(vj).ne.tag)  clique = .false.
          endif

        end do ! j

c       If exactly one row is to be merged into k-th row and there is
c       a nonzero entry in every column in that row in which there is
c       nonzero entry in row p(k) of m, then do not compute fill-in,
c       just use column indices for row which was to have been merged

        if (clique) then
          iju(k) = iju(jl(k)) + 1
          luk    = iu(jl(k)+1) - (iu(jl(k))+1)
          go to 17
        endif

c       Modify nonzero structure of k-th row by computing fill-in
c       for each row i to be merged in

        lmax   = 0
        iju(k) = juptr

        i = k
        do while(jl(i).ne.0)
          i = jl(i)

c         Merge row i into k-th row

          lui  = iu(i+1) - (iu(i)+1)
          jmin = iju(i) +  1
          jmax = iju(i) + lui
          qm   = k

          do j = jmin,jmax
            vj = ju(j)
   7        m  = qm
            qm = q(m)
            if (qm.lt.vj) go to 7
            if (qm.ne.vj) then
              luk   = luk + 1
              q(m)  = vj
              q(vj) = qm
              qm    = vj
            endif
          end do ! j

c         Remember length and position in ju of longest row merged

          if (lui.gt.lmax) then
            lmax   = lui
            iju(k) = jmin
          endif

        end do ! while

c       If k-th row is same length as longest row merged,
c       then use column indices for that row

        if (luk.ne.lmax) then

c         If tail of last row saved in ju is same as head
c         of k-th row, then overlap two sets of column indices --
c         search last row saved for first nonzero entry in k-th row

          i = q(k)
          if (jumin.gt.juptr)  go to 12
          do jmin=jumin,juptr
            if     (ju(jmin)-i .gt. 0) then
              go to 12
            elseif (ju(jmin)-i .eq. 0) then
              go to 13
            endif
          end do ! jmin

  12      go to 15

c         ... and then test whether tail matches head of k-th row

  13      iju(k) = jmin
          do j = jmin,juptr
            if (ju(j).ne.i)  go to 15
            i = q(i)
            if (i.gt.n)  go to 17
          end do ! j

          juptr = jmin - 1

c         save nonzero structure of k-th row in ju

  15      i     = k
          jumin = juptr +  1
          juptr = juptr + luk
          do j = jumin,juptr
            i       = q(i)
            ju(j)   = i
            mark(i) = k
          end do ! j
          iju(k) = jumin
        endif

c       Add k to row list for first nonzero element in k-th row

  17    if (luk.gt.1) then
          i     = ju(iju(k))
          jl(k) = jl(i)
          jl(i) = k
        endif

        iu(k+1) = iu(k) + luk
      end do ! k

      flag = 0
      return

c     Error -- duplicate entry in a

 102  flag = 2*n + p(k)
      return

c     Format

3001  format( ' SSF: Equation = ',i10,'  Time: CPU = ',f12.2,
     &        ' System = ',f12.2)

      end

      subroutine sderr(iflag,neq)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sparse solution error driver

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   iflag,neq,error

      save

c     Error messages for smpak-solution driver

      error = iflag/neq

      if     (error.eq.2) then
        write(iow,1000) mod(iflag,neq)
      elseif (error.eq.6)  then
        write(iow,1010)
      elseif (error.eq.7)  then
        write(iow,1020)
      elseif (error.eq.8)  then
        write(iow,1030) mod(iflag,neq)
      elseif (error.eq.10) then
        write(iow,1040)
      elseif (error.eq.11) then
        write(iow,1050)
      endif
      write (*,*)
      write (*,*) '**** ERROR in SMPAK (SDRV) ****'
      call plstop()

c     Formats

1000  format('**** SMPAK-ERROR (SDRV):',
     &           ' Duplicate entry in A, line',i5,' ****')

1010  format('**** SMPAK-ERROR (SDRV):',
     &           ' Insufficient memory in SSF ****')

1020  format('**** SMPAK-ERROR (SDRV):',
     &           ' Insufficient memory in SNF ****')

1030  format('**** SMPAK-ERROR (SDRV):',
     &           ' Zero pivot in line',i5,' ****')

1040  format('**** SMPAK-ERROR (SDRV):',
     &           ' Insufficient memory in SDRV ****')

1050  format('**** SMPAK-ERROR (SDRV):',
     &           ' Illegal path specification ****')

      end

      subroutine odrv (n, ia,ja,a, p,ip, esp,isp, path, flag)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  ODRV finds a minimum degree ordering of rows and
c               columns of a symmetric matrix m stored in (ia,ja,a)
c               format (see below).  For reordered matrix, work
c               and storage required to perform Gaussian elimination is
c               (usually) significantly reduced.

c     Description:
c       If only nonzero entries in upper triangle of m are being
c       stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c       with diagonal entries placed first in each row.  This is to
c       ensure that if m(i,j) will be in upper triangle of m with
c       respect to new ordering, then m(i,j) is stored in row i (and
c       thus m(j,i) is not stored);  whereas if m(i,j) will be in the
c       strict lower triangle of m, then m(j,i) is stored in row j (and
c       thus m(i,j) is not stored).

c     Storage of sparse matrices

c       The nonzero entries of matrix m are stored row-by-row in the
c       array a.  To identify individual nonzero entries in each row,
c       we need to know in which column each entry lies.  These column
c       indices are stored in array ja;  i.e., if  a(k) = m(i,j),  then
c       ja(k) = j.  To identify individual rows, we need to know where
c       each row starts.  These row pointers are stored in array ia;
c       i.e., if m(i,j) is first nonzero entry (stored) in i-th row
c       and  a(k) = m(i,j), then  ia(i) = k.  Moreover, ia(n+1) points
c       to first location following last element in last row.
c       Thus, number of entries in i-th row is  ia(i+1) - ia(i),
c       nonzero entries in i-th row are stored consecutively in

c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),

c       and corresponding column indices are stored consecutively in

c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).

c       Since coefficient matrix is symmetric, only nonzero entries
c       in upper triangle need be stored.  For example, matrix

c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )

c       could be stored as

c            | 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia | 1  4  5  8 12 14
c         ja | 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a | 1  2  3  4  2  5  6  3  6  7  8  8  9

c       or (symmetrically) as

c            | 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia | 1  4  5  7  9 10
c         ja | 1  3  4  2  3  4  4  5  5
c          a | 1  2  3  4  5  6  7  8  9          .


c     Parameters

c        n    - order of matrix

c        ia   - integer one-dimensional array contains pointer to
c               delimit rows in ja and a;  dimension = n+1

c        ja   - integer one-dimensional array containing column indices
c               corresponding to elements of a;  dimension = number of
c               nonzero entries in (the upper triangle of) m

c        a    - real one-dimensional array containing nonzero entries in
c               (the upper triangle of) m, stored by rows;  dimension =
c               number of nonzero entries in (the upper triangle of) m

c        p    - integer one-dimensional array returns permutation
c               of rows & columns of m corresponding to minimum
c               degree ordering;  dimension = n

c        ip   - integer one-dimensional array used to return inverse of
c               permutation returned in p;  dimension = n

c        esp  - declared dimension of one-dimensional array isp;  esp
c               must be >= 3n+4k,  where k is number of nonzeroes
c               in strict upper triangle of m

c        isp  - integer one-dimensional array used for working storage;
c               dimension = esp

c        path - integer path specification; values & their meanings are:
c                 1  find minimum degree ordering only
c                 2  find minimum degree ordering & reorder symmetrical
c                      stored matrix (used when only nonzero entries in
c                      upper triangle of m are being stored)
c                 3  reorder symmetrically stored matrix as specified by
c                      input permutation (used when ordering has already
c                      been determined and only nonzero entries in
c                      upper triangle of m are being stored)
c                 4  as 2 but put diagonal entries at start of each row
c                 5  as 3 but put diagonal entries at start of each row

c        flag - integer error flag;  values and their meanings are -
c                   0    no errors detected
c                  9n+k  insufficient storage in md
c                 10n+1  insufficient storage in odrv
c                 11n+1  illegal path specification
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   dflag
      integer   n,esp,path,flag, v,l,head,max,next, tmp,q
      integer   ia(*), ja(*),  p(*), ip(*),  isp(*)
      real*8    a(*)

      save

c     Initialize error flag and validate path specification

      flag = 0
      if (path.ge.1 .or. 5.ge.path) then

c       Allocate storage and find minimum degree ordering

        if ((path-1) * (path-2) * (path-4) .eq. 0) then
          max = (esp-n)/2
          v    = 1
          l    = v     +  max
          head = l     +  max
          next = head  +  n
          if (max.lt.n)  go to 110

          call  mdo (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip,
     &               isp(v), flag)
          if (flag.ne.0)  go to 100
        endif

c       Allocate storage and symmetrically reorder matrix

        if ((path-2) * (path-3) * (path-4) * (path-5) .eq. 0) then
          tmp = (esp+1) -      n
          q   = tmp     - (ia(n+1)-1)
          if (q.lt.1)  go to 110

          dflag = path.eq.4 .or. path.eq.5
          call sro (n, ip, ia,ja,a, isp(tmp), isp(q), dflag)

        endif

c     Error -- illegal path specified

      else
        flag = 11*n + 1
      endif
      return

c     Error -- error detected in mdo

 100  return

c     Error -- insufficient storage

 110  flag = 10*n + 1

      end

      subroutine oderr(iflag)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sparse solution error messages

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   iflag

      save

c     Error mesages for smpak ordering driver

      if     (iflag.eq. 9)  then
        write(iow,1000)
      elseif (iflag.eq.10) then
        write(iow,1010)
      elseif (iflag.eq.11) then
        write(iow,1020)
      endif
      write(*,*)
      write(*,*) '**** ERROR in SMPAK (odrv) ****'
      call plstop()

1000  format('**** SMPAK-ERROR (ODRV):',
     &           ' Insufficient memory in MD ****')

1010  format('**** SMPAK-ERROR (ODRV):',
     &           ' Insufficient memory in ODRV ****')

1020  format('**** SMPAK-ERROR (ODRV):',
     &           ' Illegal path specification ****')

      end
