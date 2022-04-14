c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine expolar  --  Atom overlap for exch polarization ##
c     ##                                                             ##
c     #################################################################
c
c
c     "expolar" calculates the total sum of overlap squared between 
c     atoms that is used in exchange polarization calculations
c
c     literature reference:
c
c
      subroutine expolar (overlap2)
      use atoms
      use bound
      use cell
      use couple
      use expol
      use mpole
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 s2
      real*8, allocatable :: pscale(:)
      real*8 overlap2(*)
      character*6 mode
c
c
c     zero out the value of the overlap squared at each site
c
      do ii = 1, npole
         overlap2(ii) = 0.0d0
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     find the overlap2
c
      do ii = 1, npole-1
         i = ipole(ii)
         sizi = prepep(ii)
         alphai = dmppep(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               sizk = prepep(kk)
               alphak = dmppep(kk)
               sizik = sizi*sizk
               call dampep (r,sizik,alphai,alphak,s2)
               overlap2(ii) = overlap2(ii) + s2*s2*pscale(k)
               overlap2(kk) = overlap2(kk) + s2*s2*pscale(k)
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
