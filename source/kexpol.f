c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
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
      real*8 kpr,ppr,apr
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
         if (keyword(1:8) .eq. 'EXCHPOL ') then
            k = 0
            kpr = 0.0d0
            ppr = 0.0d0
            apr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  kpr,ppr,apr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Polarization',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Spring',11x,'Size',
     &                       8x,'Damp'/)
               end if
               if (k .le. maxclass) then
                  pepk(k) = kpr
                  peppre(k) = ppr
                  pepdmp(k) = apr
                  if (.not. silent) then
                     write (iout,30)  k,kpr,ppr,apr
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
      if (allocated(kpep))  deallocate (kpep)
      if (allocated(prepep))  deallocate (prepep)
      if (allocated(dmppep))  deallocate (dmppep)
      allocate (kpep(n))
      allocate (prepep(n))
      allocate (dmppep(n))
c
c     assign the spring constant, prefactor and alpha parameters
c
      do i = 1, n
         kpep(i) = 0.0d0
         prepep(i) = 0.0d0
         dmppep(i) = 0.0d0
         ic = class(i)
         if (ic .ne. 0) then
            kpep(i) = pepk(ic)
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
         if (keyword(1:8) .eq. 'EXCHPOL ') then
            ia = 0
            kpr = 0.0d0
            ppr = 0.0d0
            apr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,kpr,ppr,apr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Polarization Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Spring',12x,'Size',
     &                       8x,'Damp'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,kpr,ppr,apr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
               kpep(ia) = kpr
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
         if (kpep(i) .ne. 0)  nexpol = nexpol + 1
         kpep(ii) = kpep(i)
         prepep(ii) = prepep(i)
         dmppep(ii) = dmppep(i)
      end do
c
c     turn on the exchange polarization potential if used
c
      if (nexpol .ne. 0)  use_expol = .true.
      return
      end
