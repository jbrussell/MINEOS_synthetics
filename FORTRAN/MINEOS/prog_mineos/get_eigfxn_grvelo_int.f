      program draw_eigfxn_grvelo_int
c
c     mineos_q performs the modal q integration on a mineos output file
c     it leaves the eigenfunctions untouched, outputting an ascii file
c     which contains the mode header, including the calculated Q
c
c     will integrate both spheroidal and toroidal modes
c
c     mineos_q expects a raw mineos output file containing:
c       two record header
c       then 1 record for every mode, depending on type
c        spheroidal modes: nn, ll, w, q, gv, u(1:knots), up(1:knots),
c                          v(1:knots), vp(1:knots), 
c                          phi(1:knots), phip(1:knots)
c        radial modes:     nn, ll, w, q, gv, u(1:knots), up(1:knots)
c        toroidal modes:   nn, ll, w, q, gv, w(1:knots), wp(1:knots)
c    
c     03/24/87 spheroidal, radial, and toroidal calculation checked 
c              against q values from Sailor and Dziewonski, 1978       
c     11/04/88 rechecked calculations - error ~1.e-3 for certain modes
c     07/02/90 test added for Q calculation in mineos
c     01/09/91 check for overlapping eigenfrequency added (thanks, Jim)
c     04/22/91 minor bugs cleaned up (concerning anisotropic prem, 
c              overlapping eigenfrequency test,gravity scalar on ocean 
c              floor for oceanic model) pp
c     07/20/91 updated to handle anisotropic models.  Calculated q is 
c              still equivelent isotropic Q,as in PREM
c              Tested by comparing Q results from an isotropic model 
c              run here as both isotropic and anisotropic.
c              Occasionally Q values differ by 2 part in 10**8 
c              (i.e. 272.43210 vs 272.43212), which seems OK
c        01/92 w,h1,gv from mineos are now real*8 and nvec 
c              is increased by 3
c        04/92 option to overwrite modal Q's if already calculated
c              in mineos
c
c      JBR 05/2022: This version uses the normalization integral to estimate I_1
c                   rather than the vertical eigenfunction at the surface. This
c                   mainly so that we don't have to worry about the water layer.
c                   Equations 3 and 4 in Lin et al. (2012)
c
c
      include 'parameter.h'
c
      real*8 intg(nknot),x(nknot)
c  JBR addition 5/2022
      real*8 intg2(nknot)
      real*8 dnj, dr
      real*8 tsum
      real*8 u_norm, v_norm
      integer*4 knot_surf
c  JBR end addition
      real*8 u, up, v , vp, wl, wp, f
      real*8 ff, mm, kk, scale, pi, rn, rj, third
      real*8 fl1(0:maxll), fl2(0:maxll),fl3(0:maxll)
      real*8 mu(nknot),muq(nknot)
      real*8 kappa(nknot),kapq(nknot)
      
      real*8 wd,h1d,gvd
c
      real*4 w,h1,gv
      real*4 bulkq(nknot),shearq(nknot),rq(nknot)
      real*4 mbq(nknot),bbq(nknot),msq(nknot),bsq(nknot)
      real*4 rad(nknot)
      real*4 interple, wlast
      real*4 abuf(maxbyte5+3),radius(nknot)
      real*4 dn(nknot),alpha(nknot),beta(nknot),qa(nknot),qb(nknot)
      real*4 alphah,betah,eta
      real*4 temp(nknot,5)
      real*4 old(maxn, maxll)
c
      integer*4 nvec
      integer*4 nn, ll, knot, ifirst, obs
c
      character*1 comp(3),ans,dtest,set
      character*256 q_file,m_file,o_file
c
      logical lq
c
      common/ablk/nn,ll,wd,h1d,gvd,buf(nknot5)
      equivalence(nn,abuf)
c
      data comp /'S','T','S'/
      data pi/3.1415926535d0/
      data rn/6371000.d0/
c JBR - begin modify (5/3/2022)
C      data old /maxold*0.0,maxold*0.0/
      data old /maxold*0./
c JBR - end modify
      data lq /.false./
      data ifirst /0/
      data qflag /999999.0/
c
cad added 1/10/18 and also 7/28/2021
      character*250 outfile
C  JBR    parameter (maxper=10)
      parameter (maxper=200)
      real*8 period(maxper)
      real*8 diffsave(maxper),persave(maxper)
      real*8 uz0(maxper),grvelo(maxper)
      real*8 phvelo(maxper),vz0(maxper)
      real*8 i1(maxper)
cad end of additions
      third = 1.d0/3.d0
      fthird = 4.d0/3.d0
c
c     open mode file and read header records - 
c     beginning of loop over mode files
c
      print*,' Enter name of input mode file'
      read(*,'(a)')m_file 
      if (m_file .eq. ' ') then
        stop
      endif
      open(unit=2,file=m_file,form='unformatted',access='sequential')
csun     +    recl=36000)
      read(2) jcom, wmin, wmax, lmin, lmax, wgrav
      read(2) knot, nic, noc, ifanis, tref, (rad(i), i=1,knot),
     &    (dn(i),i=1,knot), (alpha(i),i=1,knot), (beta(i),i=1,knot),
     &    ((temp(i,j),i=1, knot), j = 1, 5)

      ifirst = ifirst + 1

cpp <<
c
      print*, ' mode type = ',jcom,' lmin = ',lmin,' lmax = ', lmax
      print*, ' wmin = ', wmin,' wmax = ', wmax
      print*, ' knots = ',knot

cad added 7/28/2021
      write(6,*)'Number of periods?'
      read(5,*)numper
      if(numper.gt.maxper) stop 'numper too large'
      
      do i=1,numper
        write(6,*)'Enter desired period (s) for ',i
        read(5,*)period(i)
        diffsave(i)=9999.
      enddo

cad added 1/10/18
c      write(6,*)'Enter desired n and l'
c      read(5,*)nwant,lwant
      
      write(6,*)'Name of output file?'
C  JBR    read(5,*)outfile
      read(*,'(a)')outfile 
      open(1,file=outfile)
      
      write(6,*)'Surface knot? (used to normalize eigenfunctions)'
      read(5,*)knot_surf
cad end of additions

      nocor = knot - noc
c
c     zero knot+1 radius BEFORE resetting knots for toroidal case
c        jbg 6/92
c
      rad(knot+1) = 0.0
      rad(1)=1.0
c
      if (jcom .eq. 2) then
        knot = nocor
      endif
ccc
c
c     set up integration parameters
c
      do 20 j = lmin, lmax
        f = dble(j)
        fl1(j) = f * (f + 1.d0)
        fl2(j) = dsqrt(fl1(j))
        fl3(j) = (f - 1.d0) * (f + 2.d0)
 20   continue
c
c     scale is the normalization for the eigenfunctions - 
c     checked and rechecked 11/5/88
c
      scale = 1.d0/(rn*dsqrt(rn*pi*6.6723d-11)*5515.d0)
c
      if (jcom .eq. 3) then
        nvec=5*knot+8
      else
        nvec=2*knot+8
      endif
      k4 = 4*knot
      k3 = 3*knot
      k2 = 2*knot
      k1 = knot-1
c
c     begin reading mode files
c
 30   continue
      read(2,end=999)(abuf(i),i=1,nvec)
c      print*,'ll= ',ll,' nn= ',nn,' w= ',wd,' h1= ',h1d,' gv= ',gvd
c
c     look for spurious mode
c
c.... put w,q,gv in real*4 array
      w  = sngl(wd)
      h1 = sngl(h1d)
      gv = sngl(gvd)
cad added 7/29/2021
      cvel=(w*rad(knot))/(ll+0.5)/1000.
cad end of addition
cad added 1/10/18
c      if(nn.eq.nwant.and.ll.eq.lwant)then
c        write(6,*)'Found a match. T = ',2.*pi/w,' sec'
c      endif
cad end of addition      
c      print*,'ll= ',ll,' nn= ',nn,' w= ',w,' h1= ',h1,' gv= ',gv
      if (nn .lt. 0) then
        print *, ' apparent error in mode calculation'
        print *, nn,' ',comp(jcom),' ',ll
        if (nn.lt.-1) then
          print*,'n=-10: skipping next mode as well'
          read(2,end=999)
        end if
        go to 30 
      else if (nn .gt. 0) then
        wlast = old(nn,ll+1)
        if (abs(w - wlast) .lt. 1.e-7*w) then
          print*,' skipping overlapping eigenfrequency'
          print *, nn,' ',comp(jcom),' ',ll
          go to 30 
        endif 
      endif
c
c     check to see if this mode has been read already
c
      if (old(nn+1,ll+1) .ne. w) then
        old(nn+1,ll+1) = w
ccc      
c
c         equations for Q calculation taken from p. 266-267 
c         of Woodhouse, 1980
c
c         calculate compressional and shear energy densities
c         NORMALIZATIONS:
c         up, vp, and wp are normalized by 1/rn through (rj/rn) 
c                 in the equations
c         v, vp, w, and wp all carry  an implicit sqrt(ll(ll+1)) 
c                 normalization
c
cad added 7/28/2021
          T=2.*pi/w
cad end of addition	  
          if (jcom .eq. 3.and.nn.eq.0) then
            do j = knot_surf, knot_surf
              rj = dble(rad(j))
              u = dble(buf(j))*scale
              up = dble(buf(j + knot))*scale
              v = dble(buf(j + k2))*scale
              vp = dble(buf(j + k3))*scale
              
C   JBR added 5/2022
              do l = 1, knot
                rj = dble(rad(l))
                dnj = dble(dn(l))
                u_norm = dble(buf(l))*scale / u
C                v_norm = dble(buf(l + k2))*scale / v
                v_norm = dble(buf(l + k2))*scale / u
                intg2(l) = dnj*(u_norm**2 + v_norm**2)*rj**2
C          print *, dnj,' ',rj,' ',u,' ',v
C          print *, intg2(j)
              enddo
c         form the integrand - two-point trapezoidal integration
            tsum = 0.d0
            do k = 1, knot-1
              dr = rad(k+1) - rad(k)
              tsum = tsum + 0.5d0*(intg2(k)+intg2(k+1))*dr
C              print *, tsum, 0.5d0*(intg2(k)+intg2(k+1))*dr
            enddo
c   JBR end of addition       
              
          do i=1,numper
            d=abs(T-period(i))
        if(d.lt.diffsave(i))then
          diffsave(i)=d
          persave(i)=T
          uz0(i)=u 
          vz0(i)=v
          grvelo(i)=gvd
          phvelo(i)=cvel
          i1(i) = tsum
        endif
          enddo 
          
            enddo        
          elseif (jcom .eq. 1.and.nn.eq.0) then
            do j= knot, knot
              rj = dble(rad(j))
              u = dble(buf(j))*scale
              up = dble(buf(j + knot))*scale
              v = 0.d0
              vp = 0.d0
          do i=1,numper
            d=abs(T-period(i))
        if(d.lt.diffsave(i))then
          diffsave(i)=d
          persave(i)=T
          uz0(i)=u 
          vz0(i)=v
          grvelo(i)=gvd
          phvelo(i)=cvel
        endif
          enddo 
            enddo
          elseif (jcom .eq. 2.and.nn.eq.0) then
        write(6,*)'not sure if works for T b/c of rad indexing. Stop.'
        stop
            do j = 1, knot
              rj = dble(rad(j))
              wl = dble(buf(j))*scale
              wp = dble(buf(j + knot))*scale
              if(nn.eq.nwant.and.ll.eq.lwant)then
            write(1,"(3e15.5)")rj,wl,wp*(rj/rn)		
          endif
              mm = ((fl3(ll) * wl**2) + ((rj/rn)*wp - wl)**2)
c              if (nn.eq.34 .and. ll.eq.1) then
c                 print*,nn,ll,mm,mu(j),muq(j)
c              end if
              intg(j) = (mu(j)/muq(j))*mm
            enddo
          endif
 35       continue
c

      else
        print*, 
     &   ' mode ',nn,comp(jcom),ll,' overlaps with previously read mode'
      endif
      goto 30
c
c     exit this program
c
999   continue
      close(2)
c
c     formats and other such garbage
c
cad added 1/10/18 and 7/28/2021
      do i=1,numper
        if(diffsave(i).le.5.)then
          write(1,"(7e15.5)")period(i),persave(i),
     +     uz0(i),grvelo(i),phvelo(i),vz0(i),i1(i)
        else
           write(6,*)'diffsave too large for ',period(i)
        endif
      enddo
      close(1)
cad end of additions

110   format(f8.0,3f9.0)
120   format(2i6,5f15.5)
129   format(i4)
130   format(3f11.0)
      end
c
c
c
      subroutine interpol(n1, n2, x, y, m, b)
c
c     computes the coefficients for linear interpolation
c     y = mx + b
c
c     inputs:
c       n1:      lower bound for interpolation
c       n2:      upper bound for interpolation
c       x(n):    points at which the function is evaluated
c       y(n):    function to be interpolated
c     outputs:
c       m(n):    slopes of lines
c       b(n):    intercepts
c
      save
      parameter (n=1000)
      real x(n), y(n)
      real b(n), m(n)
c
      if ((n2-n1) .gt. n) then
        print*,' array limits exceeded in interpl'
        stop
      endif
      do i = n1, n2-1
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        if (dx .eq. 0.) then
          m(i) = 999.0
        else
          m(i) = dy/dx
        endif
        b(i) = y(i) - m(i)*x(i)
      end do
      return
      end
c
c
c
      real function interple(n1, n2, x, dx, xlast, y, m, b)
c
c     given the coefficients for linear interpolation
c     this routine calculates y for an input x
c
c     inputs:
c       n1:      lower bound
c       n2:      upper bound
c       x(n):    array of x-values
c       dx:      point a which the function is to be evaluated
c       y(n):    function to be interpolated
c       m(n-1):  slopes
c       b(n-1):  intercepts
c     returned
c       y:       interpolated value
c
      parameter (n=1000)
      real x(n), dx, y(n)
      real b(n), m(n), xlast
c
      if ((n2-n1) .gt. n) then
        print*,' array limits exceeded in interpl'
        stop
      endif
c
      do i = n1, n2
        if (dx .eq. x(i)) then
          if (dx .eq. x(i+1)) then
            if (xlast .eq. 0.) then
              interple = y(i+1)
              return
            elseif (xlast .lt. x(i)) then
              interple = y(i)
              return
            else
              interple = y(i+1)
              return
            endif
          else
            interple = y(i)
            return
          endif
        elseif ((dx .gt. x(i)) .and. (dx .lt. x(i+1))) then
          if (m(i) .ge. 999.0) then
            if (xlast .lt. dx) then
              interple = y(i)
            else
              interple = y(i+1)
            endif
          else
            interple = m(i)*dx + b(i)
          endif
          return
        endif
      end do
20    continue
c
c     outside array bounds - extrapolate
c
      if (dx .lt. x(n1)) then
        interple = m(n1)*dx + b(n1)
      elseif (dx .gt. x(n2)) then
        interple = m(n2)*dx + b(n2)
      else
        print*,' error in interpolation'
      endif
      return
      end

