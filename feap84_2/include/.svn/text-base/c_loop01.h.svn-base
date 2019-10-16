c     Loop over all surface 1 elements

      do ke = 1,neps1
        nel1 = ke
        do kn = 1,nope1
          nod1 = kn

c         Skip if node just checked

          if ((nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c           Search for closest master node (global search

            kset=nel1
            if (nod1.eq.2) then
              kset = neps1 + 1
            endif

