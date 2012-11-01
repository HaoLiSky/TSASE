!c    EAM version f:   June 1996  Modified 4/98 by Graeme Henkelman.
!c    Potential 6
!c      This routine evaluates the energy and force for a LEPS
!c      potential acting between three atoms in 1-D.  Only the
!c      middle atom can move.  The leftmost atom is assumed to
!c      be located at 0  and the rightmost one at rl (an input parameter).                                                                     
!c      The three atoms are treated as being of one component,
!c      i.e. there is only one call to gagafe for the whole system.
!c      The first atom in the input configuration is assumed to be
!c      the movable atom.  So, rAB=x1=RA(1) and
!c                             rBC=rl-RA(1).
!c      The second atom has the harmonic oscillator coordinate as
!c      x coordinate.
!c      Use with forces routine in forcesAllch1call.f
!c      Need to set   indTers = .true.                          
!c                                        
!C     SUBROUTINE TO COMPUTE FORCE and ENERGY for three atoms
!c     which can form only one bond and are coupled linearly to
!c     a harmonic oscillator.                     
!c                                                                                         
!      SUBROUTINE LEPSPHO(RA1,FA1,UTOT)
      SUBROUTINE FORCE(RA1,FA1,UTOT)

      parameter (MAXPOTPAR = 20)   !max # of parameters for each int.
      parameter (MAXATOMS = 2)
      parameter (MAXCOO = 3*MAXATOMS)

      implicit real*8 (a-h,o-z)

      common /unscale/RALOCAL(MAXCOO)
      common /cputim/cput0,cputas,cputil,cputig,cpudel,cputr2,cputbl,
     &               cpudef,cputpb,cpupbs,cpucpr,cpu372,cpu370,cpu373
     &              ,cpu380,cpu381,cpu382,cpu383,cpu375,cpu377,cpu376
     &              ,cpu3761

      DIMENSION RA1(MAXCOO),FA1(MAXCOO),Potpar(MAXPOTPAR),UTOT(1)
      dimension rd(3),ap1(3),r0(3),alp(3),dd(3),exch(3),exchp(3)

      Potpar(1) = 1.05
      Potpar(2) = 1.80
      Potpar(3) = 1.05                                  
      Potpar(4) = 0.742       
      Potpar(5) = 0.742
      Potpar(6) = 0.742
      Potpar(7) = 1.942
      Potpar(8) = 1.942
      Potpar(9) = 1.942
      Potpar(10) = 4.746        
      Potpar(11) = 4.746      
      Potpar(12) = 3.445
      Potpar(13) = 3.742
      Potpar(14) = 0.090
      Potpar(15) = 1.50 
      Potpar(16) = 1.0
      Potpar(17) = 1.0

!c   ------------------------------------------------------------------------                                                                                    
!c   Interpret potential parameters:
      do i=1,3
        ap1(i) = potpar(i)
        r0(i) = potpar (3+i)
        alp(i) = potpar(6+i)
        dd(i) = potpar(9+i)
      enddo
      rl=potpar(13)
      sprc=potpar(14)
      couplingc=potpar(15)

      xHO=RA1(4)

!c  scaling factors for debugging:
      fQ=potpar(16)
      fJ=potpar(17)

!c          ---------------------------------------------------------
!c  Find distances between the atom pairs:  rAB,  rBC  and rAC
!c  store in rd(1), rd(2), and rd(3), resp.
        rd(1)=RA1(1)
        rd(2)=rl-RA1(1)
        rd(3)=rd(1)+rd(2)

!c  Evaluate potential:            
!c   first sum over the Coulomb interactions:
      sumcoul=0.0
      do i=1,3
        coul=cQ(rd(i),r0(i),alp(i),dd(i),ap1(i))
        sumcoul=sumcoul+coul               
      enddo
      sumcoul=sumcoul*fQ

!c   then sum over the exchange terms:
!c    first diagonal terms: 
      sumexch=0
      do i=1,3
         exch(i)=xJ(rd(i),r0(i),alp(i),dd(i),ap1(i))
         sumexch=sumexch+(exch(i))**2
      enddo

!c    now the cross terms:
!c                      AB+BC+AC
      sumexch=sumexch-(exch(1)*exch(2)+exch(2)*exch(3)+exch(1)*exch(3))
      sumexch=sqrt(sumexch)
      sumexch=sumexch*fJ
      UTOT(1)=sumcoul-sumexch

!c ----------------------------------------------------------------------
!c    Now add the HO term and the interaction term:

      HOterm=2.*sprc*couplingc**2*(RA1(1)-(0.5*rl-1.3*xHO/couplingc))**2
      UTOT(1)=UTOT(1)+HOterm

!c ----------------------------------------------------------------------
!c   Calculate the force:            

      do i=1,3
        exchp(i)=xJp(rd(i),r0(i),alp(i),dd(i),ap1(i))
      enddo

       FA1(1)=-fQ*(cQp(rd(1),r0(1),alp(1),dd(1),ap1(1))-
     &             cQp(rd(2),r0(2),alp(2),dd(2),ap1(2)))
       cQp1=cQp(rd(1),r0(1),alp(1),dd(1),ap1(1))
       cQp2=cQp(rd(2),r0(2),alp(2),dd(2),ap1(2))

       FA1(1)=FA1(1)-fJ*(
     &      -0.5/sumexch*(2.*exchp(1)*exch(1)-2.*exchp(2)*exch(2)-
     &      exchp(1)*exch(2)+exch(1)*exchp(2)+exchp(2)*exch(3)-
     &      exchp(1)*exch(3)) )

       FA1(1)=FA1(1)-4.*sprc*couplingc**2*
     &                 (RA1(1)-(0.5*rl-1.3*xHO/couplingc))

       FA1(4)=-5.2*sprc*couplingc*(RA1(1)-(0.5*rl-1.3*xHO/couplingc))

       FA1(2)=0.0
       FA1(3)=0.0
       FA1(5)=0.0
       FA1(6)=0.0

      RETURN
      end

!c -----------------------------------------------------------------------
!c   The Coulomb function:

      function  cQ(r,r0,alpha,dd,ap1)
      implicit real*8 (a-h,o-z)

      cQ=0.5*dd*(1.5*exp(-2.*alpha*(r-r0))-exp(-alpha*(r-r0)))/ap1

      return
      end

!c -----------------------------------------------------------------------
!c   The derivative of the Coulomb function:

      function  cQp(r,r0,alpha,dd,ap1)
      implicit real*8 (a-h,o-z)

      cQp=0.5*dd*(-3.*alpha*exp(-2.*alpha*(r-r0))+
     &    alpha*exp(-alpha*(r-r0)))/ap1

      return
      end

!c -----------------------------------------------------------------------
!c   The Exchange function:

      function  xJ(r,r0,alpha,dd,ap1)
      implicit real*8 (a-h,o-z)

      xJ=0.25*dd*(exp(-2.*alpha*(r-r0))-6.*exp(-alpha*(r-r0)))/ap1

      return
      end

!c -----------------------------------------------------------------------
!c   The derivative of the Exchange function:

      function  xJp(r,r0,alpha,dd,ap1)
      implicit real*8 (a-h,o-z)

      xJp=0.25*dd*(-2.0*alpha*exp(-2.*alpha*(r-r0)) +
     &             6.0*alpha*exp(-alpha*(r-r0))) /ap1

      return
      end
