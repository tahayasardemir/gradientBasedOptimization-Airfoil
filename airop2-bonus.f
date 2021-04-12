c******************************************************************************
c   Gradient based airfoil shape optimization coupled with a panel code      **
c   AE443 - Computational Aerodynamics                                       **
c   Prof.Dr. Ismail H. Tuncer                                                **
c******************************************************************************
      program AIRFOIL_OPTIMIZATION
      parameter(mx=301)
      real xbase(mx),ybase(mx),cptarget(mx),cp(mx),dOdV(44),dstepnm_1
      real ybase2(mx)
      data designobj,designobj_min/1.,0.1E-4/ dpert/4.13e-05/
      data nspot/0/

c..Basic variables
c  np         : # of panels defining the airfoil surface (even number)
c  npp        : # of points defining the airfoil surface, np+1
c  nopmax     : Max # of optimization steps
c  nspot      : Counter for spot calls
c  xbase()    : Baseline airfoil coordinates (panel coordinates)
c  ybase()    :  xbase is fixed, ybase is updated along the optimization steps
c  cptarget() : target pressure distribution (at the mid panel locations)
c  cp()       : Pressure distribution of the current/designed airfoil
c  dOdV(22)   : Unit gradient vector components

c..Read the target Cp distr. and the baseline airfoil profile
      call INPUT(mx,np,nopmax,xbase,ybase,cptarget)
c..Start the optimization loop
      nop = 0
      dstep = 0.002

      open(55,file='objective.dat')
      open(10,file='grad.dat')
      open(65,file='line.dat')
      DO WHILE( designobj .gt. designobj_min  .and.  nop .lt. nopmax )
         nop = nop+1
c..Evaluate the objective function for the baseline/current airfoil
         call SPOT(mx,np,nspot,xbase,ybase,cp)
         designobj = OBJECTIVE(mx,np,xbase,cptarget,cp)
         print*, nop,designobj
         write(55,*) nop, designobj
         write(65,*) nop, dstep

c..Output intermediate airfoil profile and Cp distributions
         if( mod(nop,10) .eq. 0 .or. nop .eq. 1 )
     >       call CPOUT(mx,np,nop,xbase,ybase,cp)

c..Evaluate the Gradient vector
        call GRAD(mx,np,nspot,designobj,dpert,xbase,ybase,cptarget,dOdV)
        write(10,'(i3,e15.7)') (n,dodv(n),n=1,22)

c..Do a line-search along the gradient vector for the max stepsize
       do i=1,np+1
       ybase2(i)=ybase(i)
       enddo

       call NEWFOIL(mx,np,dstep,dOdV,xbase,ybase)
       call SPOT(mx,np,nspot,xbase,ybase,cp)
       dobj_n = OBJECTIVE(mx,np,xbase,cptarget,cp)

c       print *, designobj, dobj_n

       if (dobj_n.gt.designobj) then
          do while(dobj_n.gt.designobj)

          if (dstep.le.dpert) then
          do i=1,np+1
          ybase(i) = ybase2(i)
          enddo

          dstep = dpert
          call NEWFOIL(mx,np,dstep,dOdV,xbase,ybase)
          call SPOT(mx,np,nspot,xbase,ybase,cp)
          dobj_n = OBJECTIVE(mx,np,xbase,cptarget,cp)
          if(dobj_n.gt.designobj) goto 40
          exit

          else
          dstep = dstep*0.9
          endif

          do i=1,np+1
          ybase(i) = ybase2(i)
          enddo

          call NEWFOIL(mx,np,dstep,dOdV,xbase,ybase)
          call SPOT(mx,np,nspot,xbase,ybase,cp)
          dobj_n = OBJECTIVE(mx,np,xbase,cptarget,cp)

          print *, "loop", designobj, dobj_n

          enddo
       endif

      ENDDO


c..Output the final airfoil profile and Cp distributions
40    call CPOUT(mx,np,nop,xbase,ybase,cp)

      print *, nspot
      write(55,*) nspot
      stop
      end

c------------------------------------------------------------------------
      subroutine INPUT(mx,np,nopmax,xbase,ybase,cptarget)
      real xbase(mx),ybase(mx),cptarget(mx)
      character*32 filename, line
      logical ok

      print*,' '
      filename='cptarget.dat'
      print*,' Enter filename for TARGET Cp  [cptarget.dat]: '
      read (*, '(a)') line
      if( line .ne. '') filename = line
      inquire(FILE=filename,EXIST=ok)
      if( .not. ok ) then
        print*, filename,'does NOT EXIST!'
        stop
      endif
      open (1, file=filename,status='old')
      do i=1,mx
         read (1,*,end=10) xmid, cptarget(i)
      enddo
   10 close(1)
      npcp = i-1

      print*,' '
      filename='baseline.dat'
      print*,' Enter filename for BASELINE airfoil [baseline.dat]: '
      read (*, '(a)') line
      if( line .ne. '') filename = line
      inquire(FILE=filename,EXIST=ok)
      if( .not. ok ) then
        print*,' ERROR:', filename,'does NOT EXIST!'
        stop
      endif
      open (1, file=filename,status='old')
      do i=1,mx
         read(1,*,end=20) xbase(i),ybase(i)
      enddo
   20 close(1)
      np = i-2

      if( npcp .ne. np) then
        print*,' ERROR: # of panels are NOT the same:',npcp,"<->",np
        stop
      else
        print*,' # of panels:', np
      endif

      print*,' Enter the # of optimization steps : '
      read (*,*) nopmax
      print*,' '

      return
      end

c------------------------------------------------------------------------
      subroutine GRAD(mx,np,nspot,designobj,dpert,xbase,ybase,cptarget,
     >                dOdV)
c..Compute the gradient vector and its magnitude
      real xbase(mx),ybase(mx),cptarget(mx),cp(mx),ypert(mx),dOdV(42)

      nle = np/2 + 1
      npp = np+1
      dOdVmag = 0.
      do ndv = 1,21 !..perturb the lower surface with bump functions
        nfun = ndv
        call PERTURB(mx,npp,nfun,2,nle,dpert,xbase,ybase, ypert)
        call SPOT(mx,np,nspot,xbase,ypert,cp)
        dOdV(ndv)= (OBJECTIVE(mx,np,xbase,cptarget,cp)-designobj)/dpert
        dOdVmag  = dOdVmag + dOdV(ndv)**2
      enddo

      do ndv = 22,42   !..perturb the upper surface with bump functions
        nfun = ndv-21
        call PERTURB(mx,npp,nfun,nle,np,dpert,xbase,ybase, ypert)
        call SPOT(mx,np,nspot,xbase,ypert,cp)
        dOdV(ndv) =(OBJECTIVE(mx,np,xbase,cptarget,cp)-designobj)/dpert
        dOdVmag   = dOdVmag + dOdV(ndv)**2
      enddo
      dOdVmag  = sqrt(dOdVmag)

c..Evaluate the unit gradient vector
      do n = 1,42
        dOdV(n) = - dOdV(n)/dOdVmag   !..minus sign is used for minimization
      enddo

      return
      end

c------------------------------------------------------------------------
      subroutine PERTURB(mx,npp,nfun,ns,ne,delta,xbase,ybase,ypert)
      real xbase(mx),ybase(mx),ypert(mx)
      do i=1,npp
        ypert(i) = ybase(i)
      enddo
      do n=ns,ne
        ypert(n) = ypert(n) + delta*BUMP_HH(nfun,xbase(n))
c       ypert(n) = ypert(n) + delta*BUMP_Hicks(nfun,xbase(n))
      enddo
      return
      end

c-------------------BUMP_Poly-------------------------------------------
      function  BUMP_Poly(n,x)
      real pow(22)
      data pow/0.2,0.3,0.4,0.5,0.7,1.,2.,3.,4.,6.,12.,
     >12.,6.,4.,3.,2.,1.,0.7,0.5,0.4,0.3,0.2/
c..Define the bump functions
      if ( n .lt. 45 ) then
         BUMP_Poly = -4*(x**pow(n)-0.5)**2 + 1.
      else
         print*,' ERROR: No such a BUMP_Poly function:',n
         stop
      endif
      return
      end



c-----------------Generalized BUMP_Hicks-Henne -------------------------
c..Generalized Hicks-Henne functions,  Ref: D. mavripilis, AIAA J. Vol. 48, No.12 Dec 2010
c..  Bump_HH = sin(pi*x**(ln(0.5)/ln(xbump))**wbump
c..            wbump = width of the bump = 4

c..the locations of bump functions are defined for each 0.05 of chord line
c..In other words, we have 20 design variables (from 0.05 to 0.95) for
c..both upper and lower chord and 2 more additional variables (at 0.01
c..and 0.99)

      function BUMP_HH(n,x)
      real xbump(21)
      data pi/3.1415927/ wbump/4.0/
      
c..Define the bump functions here
      do k=1,21
c..apply the bump functions near the leading edge
         if(k.eq.1) then
           xbump(k)=0.01
         else if (k.lt.21) then
           xbump(k)=0.05*(k-1)
         else
           xbump(k)=0.97
         end if
      end do
      
      if (n.le.21) then
         BUMP_HH=(sin(pi*(x**(alog(0.5)/alog(xbump(n))))))**wbump
      else
         print*,' ERROR: Generalized Bump Hicks-Henne is not defined',n
      end if
      return
      end


c-------------------BUMP_Hicks-Henne-------------------------------------
      function  BUMP_Hicks(n,x)
      real pow(10)
      data pow/3*0, 0.5757166, 0.7564708, 1., 1.356915,
     +              1.943358,  3.106283,  6.578813 /
      data pi/3.1415927/ dv12/1.0/

c..Define the bump functions
      if ( n .eq. 1)  then
         BUMP_Hicks = 10*(1-x)*x**dv12/exp(40*x)
      elseif ( n .eq. 2)  then
         BUMP_Hicks = 10*(1-x)*x**dv12/exp(20*x)
      elseif ( n .eq. 3)  then
         BUMP_Hicks = sqrt(x)*(1.-x)/exp(3*x)
      elseif ( n .lt. 11 ) then
         BUMP_Hicks = sin(pi*x**pow(n))**5
c..   elseif ( n .eq. 11 ) then     !..moves the trailing edge, not used
c..      BUMP_Hicks = x**10
      else
         print*,' ERROR: No such a BUMP_Hicks-Henne function:',n
         stop
      endif
      return
      end

c------------------------------------------------------------------------
      subroutine CPOUT(mx,np,io,xbase,ybase,cp)
      real xbase(mx),ybase(mx),cp(mx)
      character filename*16,filename2*16,string*8,ext*3

      write(string,'(f5.3)') float(io)/1000
      read(string,'(2x,a3)') ext
      filename = 'cp-'//ext//'.dat'
      open(2,file=filename)
      do n = 1,np
         xm = 0.5*(xbase(n)+xbase(n+1))
         write(2,'(3E13.5)') xm, cp(n)
      enddo
      close(2)

      filename2 = 'af-'//ext//'.dat'
      open(3,file=filename2)
      do n=1,np+1
         write(3,'(2E13.5)') xbase(n), ybase(n)
      enddo
      close(3)
      return
      end

c------------------------------------------------------------------------
      function OBJECTIVE(mx,np,xbase,cptarget,cp)
      real xbase(mx),cptarget(mx),cp(mx)
      
      OBJECTIVE = 0.0
      
      do i=1,np
      dx = xbase(i+1)-xbase(i)
      OBJECTIVE = OBJECTIVE + abs((cptarget(i)-cp(i))*dx)
      enddo


      return
      end

c..New Airfoil with Bumped
      subroutine NEWFOIL(mx,np,dstep,dOdV,xbase,ybase)
      real xbase(mx),ybase(mx),dOdV(42),dstep
      integer nle, np
      nle = np/2 + 1
      
      do i=2,nle
      do ndv = 1,21
         ybase(i)=ybase(i)+dstep*dOdV(ndv)*BUMP_HH(ndv,xbase(i))
      enddo
      enddo

      do i=nle+1,np+1
      do ndv = 22,42
         ybase(i)=ybase(i)+dstep*dOdV(ndv)*BUMP_HH(ndv-21,xbase(i))
      enddo
      enddo

      return
      end


c..Panel Method Subroutine
      subroutine SPOT(mx,np,nspot,x,y,cp)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx),cp(mx)
      real Sinf(mx,mx),Ginf(mx,mx),A(mx+1,mx+2)
      real Sigma(mx),Gamma
      integer np, nspot

      npanel = np
c..Calculate Panel Angles and Midpoints
      call CALCULATE(mx,npanel,x,y,xm,ym,costhe,sinthe)
c..Setup the influence coefficients
      call INFCOEF(mx,npanel,x,y,xm,ym,costhe,sinthe,Sinf,Ginf,A)
c..Solve the system of equations
      call GAUSS(mx,npanel,A)

c..Retrieve the solution
      neqn = npanel+1
      do i = 1, npanel
         Sigma(i) = A(i,neqn+1)
      enddo
      Gamma       = A(neqn,neqn+1)

c..Compute Cp
      call LOADS(mx,npanel,x,y,xm,ym,costhe,sinthe,
     >           Gamma,Sigma,Ginf,Sinf,cp)

      nspot = nspot + 1

      return
      end

C----------------------------------------------------------------------
      subroutine CALCULATE(mx,npanel,x,y,xm,ym,costhe,sinthe)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)

c---- Compute the panel properties
      do i = 1, npanel
c-------Control points at the mid-panel locations.
        xm(i) = 0.5*(x(i+1) + x(i))
        ym(i) = 0.5*(y(i+1) + y(i))
c----Panel size and slopes/angles
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        plen = sqrt(dx**2 + dy**2)
        sinthe(i) = dy/plen
        costhe(i) = dx/plen
      enddo

      return
      end

C-----------------------------------------------------------------------
      subroutine INFCOEF(mx,npanel,x,y,xm,ym,costhe,sinthe,Sinf,Ginf,A)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      real Sinf(mx,mx),Ginf(mx,mx),A(mx+1,mx+2)

       neqn   = npanel + 1
       pi     = acos(-1.0)
       pi2inv = 1./(2*pi)

c..Initialize coefs.
       do i = 1, neqn
       do j = 1, neqn
         A(i,j) = 0.0
       enddo
       enddo

       DO i = 1,npanel
c..Find contribution of jth panel on the Vn of ith.
       DO j = 1,npanel
          if ( j .eq. i ) then
             flog = 0.0
             ftan = pi
          else
             dxj  = xm(i) - x(j)
             dxjp = xm(i) - x(j+1)
             dyj  = ym(i) - y(j)
             dyjp = ym(i) - y(j+1)
             flog = 0.5*alog((dxjp**2 + dyjp**2)/(dxj**2  + dyj**2))
             ftan = atan2( dyjp*dxj - dxjp*dyj , dxjp*dxj + dyjp*dyj )
          endif
          ctimtj = costhe(i) * costhe(j) + sinthe(i) * sinthe(j)
          stimtj = sinthe(i) * costhe(j) - costhe(i) * sinthe(j)
c..Store Ginf and Sinf values  to be used in Cp computations
          Sinf(i,j) = pi2inv * ( ftan * ctimtj + flog * stimtj )
          Ginf(i,j) = pi2inv * ( flog * ctimtj - ftan * stimtj )

          A(i,j)    = Sinf(i,j)                !..Source contribution
          A(i,neqn) = A(i,neqn) + Ginf(i,j)    !..Vortex contribution

c..Find the contribution of jth panel on the Vt of 1st and the last panels
          if (( i .eq. 1 ) .or. ( i .eq. npanel )) then
             A(neqn,j)    = A(neqn,j)    - Ginf(i,j)
             A(neqn,neqn) = A(neqn,neqn) + Sinf(i,j)
          endif
       ENDDO
c..Fill in the RHS
       A(i,neqn+1) = sinthe(i)
       ENDDO

c..Fill in the RHS  with Vt_inf of 1st and the last panels
       A(neqn,neqn+1) = -(costhe(1) + costhe(npanel))

       return
       end

c..Gauss Solver for Solution Matrix
      subroutine GAUSS(mx,npanel,A)
c..Performs Gaussian elimination on matrix A
      real A(mx+1, mx+2)

      neqn = npanel+1
      krhs = neqn + 1

c..Do full matrix elimination sequence.
      do i = 2, neqn
         im = i - 1
c..Eliminate the (i-1)th unknown from i-th through neqn-th equations.
         do j = i, neqn
            r = A(j,im) / A(im,im)
            do k = i, krhs
               A(j,k) = A(j,k) - r * A(im,k)
            enddo
         enddo
      enddo

c..Back subtitution.
      k = krhs
      A(neqn,k) = A(neqn,k) / A(neqn,neqn)
      do l = 2, neqn
         i = neqn + 1 - l
         do j = i+1, neqn
            A(i,k) = A(i,k) - A(i,j) * A(j,k)
         enddo
         A(i,k) = A(i,k) / A(i,i)
      enddo

      return
      end

c..Calculate Cp
      subroutine LOADS(mx,npanel,x,y,xm,ym,costhe,sinthe,
     +                 Gamma,Sigma,Ginf,Sinf,cp)
      real x(mx),y(mx),xm(mx),ym(mx),costhe(mx),sinthe(mx)
      real Sigma(mx),Sinf(mx,mx),Ginf(mx,mx),cp(mx)

c..Find V_tangential and Cp at the mid-point of i-th panel.
c..use influence coefficients already stored in Sinf and Ginf arrays
      do i = 1, npanel
c..Free Stream contribution
         Vtan = costhe(i)
c..Add the contribution of j-th panel.
         do j = 1, npanel
           Vtan = Vtan - Sigma(j)*Ginf(i,j) + Gamma*Sinf(i,j)
         enddo
         cp(i) = 1.0 - Vtan**2
      enddo

      return
      end
