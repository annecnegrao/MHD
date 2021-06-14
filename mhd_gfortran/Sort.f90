!===================================================================
! From Leonard J. Moss of SLAC:
! Here's a hybrid quicksort I wrote a number of years ago. It's
! based on suggestions in Knuth, volume 3, and performs much better
! than a pure quicksort on short or partially ordered input arrays.  
!===================================================================
!
!     sort -- sort, real input, iseq output
!
!
!     input:  n     integer
!             values  real
!
!     output: iseq integer (dimension n)
!
! This routine performs an in-memory sort of the first n elements of
! array values, returning into array iseq the indices of elements of
! values arranged in ascending order. Thus,
!
!    values(iseq(1)) will be the smallest number in array values;
!    values(iseq(n)) will be the largest number in values.
!
! The original data is not physically rearranged. The original order
! of equal input values is not necessarily preserved.
!
!===================================================================
!
! sortrx uses a hybrid quicksort algorithm, based on several
! suggestions in Knuth, volume 3, section 5.2.2.  in particular, the
! "pivot key" [my term] for dividing each subsequence is chosen to be
! the median of the first, last, and middle values of the subsequence;
! and the quicksort is cut off when a subsequence has 9 or fewer
! elements, and a straight insertion sort of the entire array is done
! at the end. The result is comparable to a pure insertion sort for
! very short arrays, and very fast for very large arrays (of order 12
! micro-sec/element on the 3081k for arrays of 10k elements). It is
! also not subject to the poor performance of the pure quicksort on
! partially ordered data.
!
! created:  15 jul 1986  Len Moss
!
!===================================================================
      subroutine sort(n,values,iseq) 
      integer   n,iseq(n)
      real      values(n)
 
      integer   lstk(31),rstk(31),istk
      integer   l,r,i,j,p,iseqp,iseqt
      real      valuesp
 
!     quicksort cutoff
!
!     Quit quicksort-ing when a subsequence contains m or fewer
!     elements and finish off at end with straight insertion sort.
!     According to Knuth, v.3, the optimum value of m is around 9.
 
      integer   m
      parameter (m=9)
 
!===================================================================
!
!     make initial guess for iseq
 
      do i=1,n
         iseq(i)=i
      enddo
 
!     if array is short, skip quicksort and go directly to
!     the straight insertion sort.
 
      if (n.le.m) goto 900
 
!===================================================================
!
!     quicksort
!
!     the "qn:"s correspond roughly to steps in algorithm q,
!     knuth, v.3, pp.116-117, modified to select the median
!     of the first, last, and middle elements as the "pivot
!     key" (in knuth's notation, "k").  also modified to leave
!     data in place and produce an iseq array.  to simplify
!     comments, let values[i]=values(iseq(i)).
 
! q1: initialize
      istk=0
      l=1
      r=n
 
  200 continue
 
! q2: sort the subsequence values[l]..values[r].
!
!     at this point, values[l] <= values[m] <= values[r] for all l < l,
!     r > r, and l <= m <= r.  (first time through, there is no
!     values for l < l or r > r.)
 
      i=l
      j=r
 
! q2.5: select pivot key
!
!     let the pivot, p, be the midpoint of this subsequence,
!     p=(l+r)/2; then rearrange iseq(l), iseq(p), and iseq(r)
!     so the corresponding values values are in increasing order.
!     the pivot key, valuesp, is then values[p].
 
      p=(l+r)/2
      iseqp=iseq(p)
      valuesp=values(iseqp)
 
      if (values(iseq(l)) .gt. valuesp) then
         iseq(p)=iseq(l)
         iseq(l)=iseqp
         iseqp=iseq(p)
         valuesp=values(iseqp)
      endif
 
      if (valuesp .gt. values(iseq(r))) then
         if (values(iseq(l)) .gt. values(iseq(r))) then
            iseq(p)=iseq(l)
            iseq(l)=iseq(r)
         else
            iseq(p)=iseq(r)
         endif
         iseq(r)=iseqp
         iseqp=iseq(p)
         valuesp=values(iseqp)
      endif
 
!     now we swap values between the right and left sides and/or
!     move valuesp until all smaller values are on the left and all
!     larger values are on the right.  neither the left or right
!     side will be internally ordered yet; however, valuesp will be
!     in its final position.
 
  300 continue
 
! q3: search for datum on left >= valuesp
!
!     at this point, values[l] <= valuesp.  we can therefore start scanning
!     up from l, looking for a value >= valuesp (this scan is guaranteed
!     to terminate since we initially placed valuesp near the middle of
!     the subsequence).
 
         i=i+1
         if (values(iseq(i)).lt.valuesp) goto 300
 
  400 continue
 
! q4: search for datum on right <= valuesp
!
!     at this point, values[r] >= valuesp.  we can therefore start scanning
!     down from r, looking for a value <= valuesp (this scan is guaranteed
!     to terminate since we initially placed valuesp near the middle of
!     the subsequence).
 
         j=j-1
         if (values(iseq(j)).gt.valuesp) goto 400
 
! q5: have the two scans collided?
 
      if (i.lt.j) then
 
! q6: no, interchange values[i] <--> values[j] and continue
 
         iseqt=iseq(i)
         iseq(i)=iseq(j)
         iseq(j)=iseqt
         goto 300
      else
 
! q7: yes, select next subsequence to sort
!
!     at this point, i >= j and values[l] <= values[i] == valuesp <= values[r],
!     for all l <= l < i and j < r <= r.  if both subsequences are
!     more than m elements long, push the longer one on the stack and
!     go back to quicksort the shorter; if only one is more than m
!     elements long, go back and quicksort it; otherwise, pop a
!     subsequence off the stack and quicksort it.
 
         if (r-j .ge. i-l .and. i-l .gt. m) then
            istk=istk+1
            lstk(istk)=j+1
            rstk(istk)=r
            r=i-1
         else if (i-l .gt. r-j .and. r-j .gt. m) then
            istk=istk+1
            lstk(istk)=l
            rstk(istk)=i-1
            l=j+1
         else if (r-j .gt. m) then
            l=j+1
         else if (i-l .gt. m) then
            r=i-1
         else
! q8: pop the stack, or terminate quicksort if empty
            if (istk.lt.1) goto 900
            l=lstk(istk)
            r=rstk(istk)
            istk=istk-1
         endif
         goto 200
      endif
 
  900 continue
 
!===================================================================
!
! q9: straight insertion sort
 
      do 950 i=2,n
         if (values(iseq(i-1)) .gt. values(iseq(i))) then
            iseqp=iseq(i)
            valuesp=values(iseqp)
            p=i-1
  920       continue
               iseq(p+1) = iseq(p)
               p=p-1
               if (p.gt.0) then
                  if (values(iseq(p)).gt.valuesp) goto 920
               endif
            iseq(p+1) = iseqp
         endif
  950    continue
 
!===================================================================
!
!     all done
      return
      end