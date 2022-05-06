c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine expolar3  --  ExchPol variable polarizability  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "expolar3" calculates the variable polarizability due to exchange
c     polarization
c
c     literature reference:
c
c
      subroutine expolar3 (polscale,invpolscale)
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
      real*8 det
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 springi,springk
      real*8 s2,ds2
      real*8 p33i, p33k
      real*8 kS2i(3,3)
      real*8 kS2k(3,3)
      real*8 ps(3,3)
      real*8, allocatable :: pscale(:)
      real*8 polscale(3,3,npole)
      real*8 invpolscale(3,3,npole)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     set polscale and invpolscale to the identity matrix
c
      do i = 1, npole
         polscale(1,1,i) = 1.0d0
         polscale(2,1,i) = 0.0d0
         polscale(3,1,i) = 0.0d0
         polscale(1,2,i) = 0.0d0
         polscale(2,2,i) = 1.0d0
         polscale(3,2,i) = 0.0d0
         polscale(1,3,i) = 0.0d0
         polscale(2,3,i) = 0.0d0
         polscale(3,3,i) = 1.0d0
         invpolscale(1,1,i) = 1.0d0
         invpolscale(2,1,i) = 0.0d0
         invpolscale(3,1,i) = 0.0d0
         invpolscale(1,2,i) = 0.0d0
         invpolscale(2,2,i) = 1.0d0
         invpolscale(3,2,i) = 0.0d0
         invpolscale(1,3,i) = 0.0d0
         invpolscale(2,3,i) = 0.0d0
         invpolscale(3,3,i) = 1.0d0
      end do
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     find the variable polarizability
c
      do ii = 1, npole-1
         i = ipole(ii)
         springi = kpep(ii)
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
               springk = kpep(kk)
               sizk = prepep(kk)
               alphak = dmppep(kk)
               sizik = sizi*sizk
               call dampep (r,sizik,alphai,alphak,s2,ds2)
               p33i = springi*s2*pscale(k)
               p33k = springk*s2*pscale(k)
               call exrotate(i,k,p33i,p33k,kS2i,kS2k)
               do j = 1, 3
                  do m = 1, 3
                     polscale(j,m,ii) = polscale(j,m,ii) + kS2i(j,m)
                     polscale(j,m,kk) = polscale(j,m,kk) + kS2k(j,m)
                  end do
               end do
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
c     invert polscale matrix
c
      do ii = 1, npole
         do j = 1, 3
            do m = 1, 3
               ps(j,m) = polscale(j,m,ii)
            end do
         end do
         det = ps(1,1)*(ps(2,2)*ps(3,3) - ps(3,2)*ps(2,3))
     &       - ps(1,2)*(ps(2,1)*ps(3,3) - ps(2,3)*ps(3,1))
     &       + ps(1,3)*(ps(2,1)*ps(3,2) - ps(2,2)*ps(3,1))
         invpolscale(1,1,ii) = (ps(2,2)*ps(3,3)-ps(3,2)*ps(2,3))/det
         invpolscale(1,2,ii) = (ps(1,3)*ps(3,2)-ps(1,2)*ps(3,3))/det
         invpolscale(1,3,ii) = (ps(1,2)*ps(2,3)-ps(1,3)*ps(2,2))/det
         invpolscale(2,1,ii) = (ps(2,3)*ps(3,1)-ps(2,1)*ps(3,3))/det
         invpolscale(2,2,ii) = (ps(1,1)*ps(3,3)-ps(1,3)*ps(3,1))/det
         invpolscale(2,3,ii) = (ps(2,1)*ps(1,3)-ps(1,1)*ps(2,3))/det
         invpolscale(3,1,ii) = (ps(2,1)*ps(3,2)-ps(3,1)*ps(2,2))/det
         invpolscale(3,2,ii) = (ps(3,1)*ps(1,2)-ps(1,1)*ps(3,2))/det
         invpolscale(3,3,ii) = (ps(1,1)*ps(2,2)-ps(2,1)*ps(1,2))/det
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine exrotate  --  rotate polarizability matrix ##
c     ##                                                        ##
c     ############################################################
c
c
c     "exrotate" rotates and inverts the variable polarizability tensor
c     due to exchange polarization
c
c
      subroutine exrotate (i,k,p33i,p33k,kS2i,kS2k)
      use atoms
      use math
      use mpole
      use polpot
      implicit none
      integer i,ii,k,kk
      integer j,l,m,o
      real*8 r,dot
      real*8 eps,angle
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 dx,dy,dz
      real*8 p33i,p33k
      real*8 ai(3,3)
      real*8 ak(3,3)
      real*8 kS2i(3,3)
      real*8 kS2k(3,3)
c
c
c     use the identity matrix as the default rotation matrix
c
      ai(1,1) = 1.0d0
      ai(2,1) = 0.0d0
      ai(3,1) = 0.0d0
      ai(1,2) = 0.0d0
      ai(2,2) = 1.0d0
      ai(3,2) = 0.0d0
      ai(1,3) = 0.0d0
      ai(2,3) = 0.0d0
      ai(3,3) = 1.0d0
      if (use_expol) then
c
c     get coordinates and frame definition for the multipole site
c
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xk = x(k)
         yk = y(k)
         zk = z(k)
c
c     get the rotation matrix elements
c
         dx = xk - xi
         dy = yk - yi
         dz = zk - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         ai(1,3) = dx / r
         ai(2,3) = dy / r
         ai(3,3) = dz / r
         dx = 1.0d0
         dy = 0.0d0
         dz = 0.0d0
         dot = ai(1,3)
         eps = 0.707d0
         if (abs(dot) .gt. eps) then
            dx = 0.0d0
            dy = 1.0d0
            dot = ai(2,3)
         end if
         dx = dx - dot*ai(1,3)
         dy = dy - dot*ai(2,3)
         dz = dz - dot*ai(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         ai(1,1) = dx / r
         ai(2,1) = dy / r
         ai(3,1) = dz / r
         ai(1,2) = ai(3,1)*ai(2,3) - ai(2,1)*ai(3,3)
         ai(2,2) = ai(1,1)*ai(3,3) - ai(3,1)*ai(1,3)
         ai(3,2) = ai(2,1)*ai(1,3) - ai(1,1)*ai(2,3)
         ak = ai
         ak(1,2) = -ak(1,2)
         ak(2,2) = -ak(2,2)
         ak(3,2) = -ak(3,2)
         ak(1,3) = -ak(1,3)
         ak(2,3) = -ak(2,3)
         ak(3,3) = -ak(3,3)
c
c     apply R^T from left and R from right to rotate kS2 matrix
c
         do j = 1, 3
            do l = 1, 3
               kS2i(j,l) = p33i*ai(j,3)*ai(l,3)
               kS2k(j,l) = p33k*ak(j,3)*ak(l,3)
            end do
         end do
      else
         do j = 1, 3
            do l = 1, 3
               kS2i(j,l) = 0.0d0
               kS2k(j,l) = 0.0d0
            end do
         end do
      end if
      return
      end