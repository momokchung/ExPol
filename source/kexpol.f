c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kexpol  --  exch polarization term assignment   ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kexpol" assigns the constant prefactor and damping alpha
c     for exchange polarization interactions and processes any new or
c     changed values for these parameters
c
c
      subroutine kexpol
      use atomid
      use atoms
      use inform
      use iounit
      use kexpl
      use keys
      use mpole
      use polpot
      use expol
      use sizes
      implicit none
      integer i,k,ii
      integer ia,ic,next
      real*8 ppr,apr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing exchange polarization parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'EXCHIND ') then
            k = 0
            ppr = 0.0d0
            apr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  ppr,apr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Polarization',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Size',11x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (k .le. maxclass) then
                  peppre(k) = ppr
                  pepdmp(k) = apr
                  if (.not. silent) then
                     write (iout,30)  k,ppr,apr
   30                format (6x,i6,7x,2f15.4,f15.3)
                  end if
               else
                  write (iout,40)
   40             format (/,' KEXPOL  --  Too many Exch Polarization',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(prepep))  deallocate (prepep)
      if (allocated(dmppep))  deallocate (dmppep)
      allocate (prepep(n))
      allocate (dmppep(n))
c
c     assign the repulsion size, alpha and valence parameters 
c
      do i = 1, n
         prepep(i) = 0.0d0
         dmppep(i) = 0.0d0
         ic = class(i)
         if (ic .ne. 0) then
            prepep(i) = peppre(ic)
            dmppep(i) = pepdmp(ic)
         end if
      end do
c
c     process keywords containing atom specific exchange polarization
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'EXCHIND ') then
            ia = 0
            ppr = 0.0d0
            apr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,ppr,apr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Polarization Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Size',12x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,ppr,apr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
               prepep(ia) = ppr
               dmppep(ia) = apr
            end if
   70       continue
         end if
      end do
c
c     condense exch polarization sites to the list of multipole sites
c
      nexpol = 0
      do ii = 1, npole
         i = ipole(ii)
         if (prepep(i) .ne. 0)  nexpol = nexpol + 1
         prepep(ii) = prepep(i)
         dmppep(ii) = dmppep(i)
      end do
c
c     turn on the exchange polarization potential if used
c
      if (nexpol .ne. 0)  use_expol = .true.
      return
      end