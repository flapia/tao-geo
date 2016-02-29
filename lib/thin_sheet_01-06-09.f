CC **********************************************************************
CC  Set of subroutines to calculate:
CC   -The velocity field of a thin sheet (subroutine VELOCITY_FIELD)
CC      from the boundary conditions (subroutine READ_BC)
CC   -Thickening of the sheet (subroutine THICKEN)
CC **********************************************************************


CCC======================================================================
CC  SUBROUTINE  velvis
CC  Calculation of velocity field of a thin sheet, from a given 
CC      lithospheric structure and viscosity. 
CC	viscosity is velocity depending -> iteration
CC
CC   ix=0,...,n.  iy=0,...,m.
CC   x positiu cap a la dreta, y positiu cap dalt.
CC   vel [m/s]	Array of the velocity field: u1, v1, u2, v2, ... unn, vnn.
CC			From west to east and south to north.
CC   vel(nincogn)=(u,v)(nincogn)
CC   velold : velocity field at the first iteration
CC   vel = alfa * vel + (1-alfa) * velold    (alfa=.5 is used)
CC   ep_limit [s-1]  Imposed limit for calculated vert.str.rate 'epuntzz'.
CC
CC  CALCULO EL NOU CAMP DE VELOCITATS SENSE TENIR EN COMPTE L'ANTERIOR
CC  A L'ENTRAR A LA SUBRUTINA INICIALITZO EL VECTOR DE VELOCITATS, vel.
CC  Si la deformacio es molt gran -> Disminueixo el Dt
CC  Com m‚s petits siguin el Dx i Dy -> m‚s petit ha de ser el Dt
CC     (strain maxim) = Dt*(strain rate maxim) < quotadef

      SUBROUTINE velvis (TEMPSMa, Dt, AX, BY, n, m, nn, 
     +                  nitermax, tallmax, alfa, g, 
     +                  visTer, vis, vissup, visinf, 
     +			elevation, sediment, s, GL, u, v, epuntzz,
     +			RHOH2O, rosed, roc, RHOAST,
     +			roalfa, TISOTER, Tmoho, ZASTH, BCfile)

CC  TEMPSMa[Ma]	If <> 0 then uses the previous velocity field.
CC  Dt [s]   	Time increment, if the amount of deformation is too 
CC			large then Dt will be reduced to fit 'quotadef'
CC  AX [m]	domain length in x
CC  BY [m]	domain length in y
CC  n, m	Number of nodes (in x,y directions) minus 1: column index 
CC			runs from 0 (west) to n (east)
CC  nn		Number of nodes (n+1)*(m+1)
CC  nitermax	max. number of iterations for convergence between 
CC			viscosity and velocity (if exceeded => 
CC			divergence).
CC			If 0 => viscosity doesn't depend on strain rate.
CC  tallmax	Criterium of convergence. If SUM(Dvel/vel) < tallmax =>
CC			convergence. (p.e., =0.01).
CC  alfa	0 => explicit ;  alfa=1 => implicit. 0.5 recommended
CC  g [m/s2] 	Gravity acceleration at the surface of the study planet
CC  visTer [Pa]	Array of thermal term of viscosity: 
CC			1/2 * stregth[Pa*m] / layer.thickness[m]
CC  vis [Pa*s]	Array of effective viscosity of the layer:
CC			visTer / ref.strain.rate
CC			Returns the calculated viscosity (a function of 
CC			the new strain rate).
CC  vissup [Pa*s]  Imposed upper limit for calculated viscosity 'vis'. (~10^25 Pa*s)
CC  visinf [Pa*s]  Imposed lower limit for calculated viscosity 'vis'. (~10^22 Pa*s)
CC  elevation [m]  Array of topography (positive above sea level).
CC  sediment [m]   Array of sediment thickness.	
CC  s [m]	Array of crustal thickness (needed for the lateral 
CC			pressure variations). Does not include sediments.
CC  GL [m]	Array of lithospheric thickness (needed for the lateral 
CC			pressure variations). Includes sedims. and crust.
CC  u, v [m/s]	Arrays of velocity in x, y directions in the previous time step.
CC			Sorted from west to east and south to north.
CC			Returns the new calculated velocity.
CC  epuntzz [1/s]  Array returning the vertical strain rate.
CC  RHOH2O [kg/m3] Density of sea water.
CC  rosed [kg/m3]  Density of sediments.	
CC  roc	[kg/m3]	Density of the crust.
CC  RHOAST [kg/m3] Density of asthenosphere.	
CC			of pressure.
CC  roalfa [1/K]   Volumetric thermal expansion coeff. of lith. mantle (3.5e-5 K-1).
CC  TISOTER [K]	Temperature at which the lithosphere has density RHOAST (1573 K = 1300 C).
CC			rom = RHOAST * (1+roalfa/2*(TISOTER-Tmoho))
CC  Tmoho [K]	Array of temperature distribution at Moho.
CC  ZASTH [m]	Depth of compensation to calculate lateral variations (e.g., 300e3 m)
CC  BCfile [char*80]  Boundary conditions file name 

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       CHARACTER*84 TITOL_BC, BCfile
       PARAMETER (pi=3.1415926535897932D0, FACTEMP=3.1536D7, 
     +		quotadef=8.D-2, DP_sup=10000.D6, DP_inf=-10000.D6)
c       DOUBLE PRECISION, allocatable :: b(:), PRESSIO(:)
c       DOUBLE PRECISION, allocatable :: a(:,:) 
       DIMENSION u(nn),v(nn),average_Pressure(nn),vis(nn),visTer(nn),
     +           sediment(nn),elevation(nn),s(nn),GL(nn),epuntzz(nn),
     +           Tmoho(nn)
CC  Parametres per a les Condicions de Contorn
       PARAMETER (vr=1.D-10,IBC=2,IFAULT=1)   

       PRINT*,'====== SUBROUTINE velvis  ======'

c       allocate( b(nincogn), PRESSIO(nn) ,STAT=ist1)
c       allocate( a(n1,n2) ,STAT=ist2)
       GLmin=500.D3
       GLmax=0.D0
       elevmin=500.D3
       elevmax=-500.D3
       vismax=0.D0
       vismin=1.D35
       sedmin=500.D3
       sedmax=0.D0
       Pmin=500.D30
       Pmax=0.D0
       rommin=500.D30
       rommax=0.D0

       Dx=AX/n
       Dy=BY/m
       DO 132 kixy=1,nn
	  sed=sediment(kixy)
	  hm=GL(kixy)-s(kixy)-sed
	  ha=ZASTH-GL(kixy)
	  water=ABS(elevation(kixy))
	  rom=RHOAST*(1+((roalfa/2.D0)*(TISOTER-Tmoho(kixy))))
	  if(elevation(kixy).gt.0.D0) water=0.D0
	  Fw=(g*RHOH2O*water*water)/2.D0
	  Fsed=(g*RHOH2O*water*sed)+((g*rosed*sed*sed)/2.D0)
     	  Fc=(g*(RHOH2O*water+rosed*sed)*s(kixy))+
     +			((g*roc*s(kixy)*s(kixy))/2.D0)
	  Fm=(g*(RHOH2O*water+rosed*sed+roc*s(kixy))*hm)+
     +			(g*rom*hm*hm/2.D0)
	  Fa=(g*(RHOH2O*water+rosed*sed+roc*s(kixy)+rom*hm)*ha)+
     +			(g*RHOAST*ha*ha/2.D0)
     	  average_Pressure(kixy)=(Fw+Fsed+Fc+Fm+Fa)/ZASTH
              GLmin=MIN(GLmin,GL(kixy))
              GLmax=MAX(GLmax,GL(kixy))
	      elevmin=MIN(elevmin,elevation(kixy))
	      elevmax=MAX(elevmax,elevation(kixy))
	      vismin=MIN(vismin,vis(kixy))
	      vismax=MAX(vismax,vis(kixy))
	      Pmin=MIN(Pmin,average_Pressure(kixy))
	      Pmax=MAX(Pmax,average_Pressure(kixy))
	      rommin=MIN(rommin,rom)
	      rommax=MAX(rommax,rom)
	      sedmin=MIN(sedmin,sed)
	      sedmax=MAX(sedmax,sed)
 132   CONTINUE
       IF(TEMPSMa.EQ.0) WRITE(6,364) DP_inf/1.D6, DP_sup/1.D6
 364   FORMAT(/2X,'Delta Presion limits: ',F11.2,' MPa < DPresion <',
     +		F11.2,' MPa')
       
       WRITE(6,365) GLmin/1.D3,GLmax/1.D3,elevmin,elevmax,sedmin,sedmax,
     +		     rommin,rommax,ZASTH/1.D3,Pmin,Pmax,vismin,vismax
 365   FORMAT(/2X,'LITHOSPHERIC THICKNESS: minimum =',F9.4,' km,',
     +		'   maximum =',F9.4,' km'/
     +	2X,'ELEVATION: minima =',F9.2,' m,   maxima =',F9.2,' m'/
     +	2X,'SEDIMENTS: minim =',F9.2,' m,    maxim =',F9.2,' m'/
     +	2X,'DENSITY OF LITHOSPHERIC MANTLE: minima =',F9.2,
     +	' km/m3,    maxima =',F9.2,' km/m3'/2X, 
     + 'PRESION, FORCE OF THE LITHOSPHERIC COLUMN OF DEPTH: ',F8.2,' km'
     +  /5X,'minima =',1P,E12.4,' N/m2,   maxima =',1P,E12.4,' N/m2'/
     +	2X,'VISCOSITY: minima =',1P,E11.3,' Pa.s,   maxima =',1P,E11.3,
     +	' Pa.s'/)

	nbanda=4*(n+1)+7 
	nincogn=2*nn

       CALL VELOCITY_FIELD (TEMPSMa, Dt, AX, BY, n, m, nn,
     +			vissup, visinf, visTer, vis, nincogn,
     +			tallmax, alfa, nitermax, average_Pressure, u, v,
     +			BCfile, nbanda)
      CALL vertical_strain_rate (AX, BY, n, m, nn, u, v, epuntzz)

       RETURN
       END
              
CC ********************************************************************
CC ********************************************************************

       SUBROUTINE VELOCITY_FIELD (TEMPSMa, Dt, AX, BY, n, m, nn,
     +			vissup, visinf, visTer, vis, nincogn,
     +			tallmax, alfa, nitermax, average_Pressure, u, v,
     +			BCfile, nbanda)
C     +                  nincogn)
C     nbanda)
C     nitermax)
Ctallmax, alfa)
C     +                  visTer, vis, vissup, visinf,
C     +			average_Pressure, u, v, BCfile)

CC  TEMPSMa[Ma]	If != 0 then uses the previous velocity field.
CC  Dt [s]   	Time increment, if the amount of deformation is too 
CC			large then Dt will be reduced to fit 'quotadef'
CC  AX [m]	domain length in x
CC  BY [m]	domain length in y
CC  n, m	Number of nodes (in x,y directions) minus 1: column index 
CC			runs from 0 (west) to n (east)
CC  nn		Number of nodes (n+1)*(m+1)
CC  nincogn	2*nn
CC  nbanda	4(n+1)+7 
CC  nitermax	max. number of iterations for convergence between 
CC			viscosity and velocity (if exceeded => 
CC			divergence).
CC			If 0 => viscosity doesn't depend on strain rate.
CC  tallmax	Criterium of convergence. If SUM(Dvel/vel) < tallmax =>
CC			convergence. (p.e., =0.01).
CC  alfa	0 => explicit ;  alfa=1 => implicit. 0.5 recommended
CC  visTer [Pa]	Array of thermal term of viscosity: 
CC			1/2 * stregth[Pa*m] / layer.thickness[m]
CC  vis [Pa*s]	Array of effective viscosity of the layer:
CC			visTer / ref.strain.rate
CC			Returns the calculated viscosity (a function of 
CC			the new strain rate).
CC  vissup [Pa*s]  Imposed upper limit for calculated viscosity 'vis'. (~10^25 Pa*s)
CC  visinf [Pa*s]  Imposed lower limit for calculated viscosity 'vis'. (~10^22 Pa*s)
CC average_Pressure  [Pa=N/m2] Integral of (g*rho(z)*z*dz)/L between 0 and ZASTH, 
CC			where L is the thickness of the thin sheet, and
CC			ZASTH is the level of compensation.
CC  u, v [m/s]	Arrays of velocity in x, y directions in the previous time step.
CC			Sorted from west to east and south to north.
CC			Returns the new calculated velocity.
CC  BCfile [char*84]  Boundary conditions file name 

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       CHARACTER*84 TITOL_BC, BCfile
       PARAMETER (pi=3.1415926535897932D0, FACTEMP=3.1536D7, 
     +		quotadef=8.D-2, DP_sup=10000.D6, DP_inf=-10000.D6)
C       DOUBLE PRECISION, allocatable :: b(:), PRESSIO(:)
c       DOUBLE PRECISION, allocatable :: a(:,:) 
       DIMENSION a(nincogn,nbanda), u(nn), v(nn), b(nincogn),
     +		average_Pressure(nn), visTer(nn), vis(nn), visold(nn),
     +		vel(nincogn), velold(nincogn)

CC  Parametres per a les Condicions de Contorn
       PARAMETER (vr=1.D-10,IBC=2,IFAULT=1)   

       PRINT*,'====== SUBROUTINE VELOCITY_FIELD  ======'
c       PRINT*,'TEMPSMa, Dt, AX, BY, n, m: ', TEMPSMa, Dt, AX, BY, n, m
c       PRINT*,'nn, nincogn, nbanda:',nn, nincogn,nbanda
c       PRINT*,'vissup,visinf,tallmax,alfa:',vissup,visinf,tallmax,alfa
        
c       allocate( b(nincogn), PRESSIO(nn) ,STAT=ist1)
c       allocate( a(n1,n2) ,STAT=ist2)

       tall=0.D0
       vismax=0.D0
       vismin=1.D35
       Dx=AX/n
       Dy=BY/m
       kvel=1
       DO 132 kixy=1,nn
          visold(kixy)=vis(kixy)
          vel(kvel)=0.D0
          vel(kvel+1)=0.D0
          velold(kvel)=u(kixy)
          velold(kvel+1)=v(kixy)
          kvel=kvel+2
132   CONTINUE
C	  PRINT*,'1 u:',u(1),'1 v:',v(1)
C	  PRINT*,'nn u:',u(nn),'nn v:',v(nn)
C	  PRINT*,'1 visTer:',visTer(1),'nn visTer:',visTer(nn)
C	  PRINT*,'1 vis:',vis(1),'  aver_Pres :',average_Pressure(1)
C	  PRINT*,'nn vis:',vis(nn),'  aver_Pres :',average_Pressure(nn)

CC *************** CONSTRUCCIO DE LA MATRIU *********************
       Li=2*(n+1)+3
       Ls=Li
       Ld=Li+1
       kp=Ld+2*(n+1) 
       kn=Ld-2*(n+1)
       
       niter=0
 111   CONTINUE 
       niter=niter+1
         do 5 i=1,nincogn
            b(i)=0.D0
            do 5 j=1,nbanda
               a(i,j)=0.D0
5        continue

CC  ----------  Condicions de Contorn --------------------
CC  POL DE ROTACIO: 1:NUVEL-1A.  2:Argus. 3:RM2.  4:P071.  5:Mckenzie
		IF(IBC.EQ.1) THEN
			dlonpol=-20.6D0
			dlatpol=21.0D0
			omegada=0.12D0
		ENDIF
		IF(IBC.EQ.2) THEN
			dlonpol=-20.3D0
			dlatpol=18.8D0
			omegada=0.104D0
		ENDIF
		IF(IBC.EQ.3) THEN
			dlonpol=-21.19D0
			dlatpol=25.23D0
			omegada=0.1D0
		ENDIF
		IF(IBC.EQ.4) THEN
			dlonpol=-23.5D0
			dlatpol=29.2D0
			omegada=0.14D0
		ENDIF
		IF(IBC.EQ.5) THEN
			dlonpol=-28.2D0
			dlatpol=22.7D0
			omegada=0.27D0
		ENDIF
C       CALL fixCCSud (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada)
C       CALL fixCCWest (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada,vis,nn)
C       CALL fixCCEst (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
C     +			dlonpol,dlatpol,omegada,IFAULT)
C       CALL fixCCNort (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
       
      CALL READ_BC (TEMPSMa,niter,n,m,nn,nincogn,nbanda,Dx,Dy,vis,
     +			a,b,TITOL_BC, BCfile)

CC  ---------- NO POL ----------------------------------------
C       CALL fixCCSud (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
C       CALL fixCCWest (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
C       CALL fixCCEst (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
C       CALL fixCCNort (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
C  ------------------------------------------------------
C   PUNTS INTERIORS
	DP_max=-1.D100
	DP_min=1.D100
	Dvis_max=-1.D100
	Dvis_min=1.D100
       DO 9 iy=1,m-1
           DO 10 ix=1,n-1
               k=ix+1+iy*(n+1)
                GLv=vis(k)
                G1=vis(k+1)
                G2=vis(k-1)
                G3=vis(k+n+1)
                G4=vis(k-n-1)
C  ---------------------------------------------------------------------
               GLvx=(G1-G2)/(2.D0*Dx)
               GLvy=(G3-G4)/(2.D0*Dy)
	       Dvis_max=MAX(Dvis_max,GLvx,GLvy)
	       Dvis_min=MIN(Dvis_min,GLvx,GLvy)
	       DP_x=(average_Pressure(k+1)-average_Pressure(k-1))/
     +	       		(2.D0*Dx)
	       DP_y=(average_Pressure(k+n+1)-average_Pressure(k-n-1))/
     +	       		(2.D0*Dy)
CC	LIMIT DEL GRADIENT DE LA average_Pressure
		DP_x=MIN(DP_x,DP_sup)
	       	DP_x=MAX(DP_x,DP_inf)
		DP_y=MIN(DP_y,DP_sup)
	       	DP_y=MAX(DP_y,DP_inf)
               DP_max=MAX(DP_max,DP_x,DP_y)
               DP_min=MIN(DP_min,DP_x,DP_y)
	       Tx=-1.D0*DP_x
	       Ty=-1.D0*DP_y
               l=2*ix+1+2*iy*(n+1)
               DIAGONAL1=-(8.D0*GLv)/(Dx*Dx)-(2.D0*GLv)/(Dy*Dy)
               a(l,Ld+2)=((4.D0*GLv)/(Dx*Dx)+(2.D0*GLvx)/Dx)/DIAGONAL1
               a(l,kp)=(GLv/(Dy*Dy)+GLvy/(2.D0*Dy))/DIAGONAL1
               a(l,Ld)=1.D0
               a(l,kn)=(GLv/(Dy*Dy)-GLvy/(2.D0*Dy))/DIAGONAL1
               a(l,Ld-2)=((4.D0*GLv)/(Dx*Dx)-(2.D0*GLvx)/Dx)/DIAGONAL1
               a(l,kp+3)=((3.D0*GLv)/(4.D0*Dx*Dy))/DIAGONAL1
               a(l,kn+3)=-a(l,kp+3)
               a(l,kp-1)=-a(l,kp+3)
               a(l,kn-1)=a(l,kp+3)
               a(l,kp+1)=(GLvx/Dy)/DIAGONAL1
               a(l,kn+1)=-a(l,kp+1)
               a(l,Ld+3)=(GLvy/(2.D0*Dx))/DIAGONAL1
               a(l,Ld-1)=-a(l,Ld+3)
               b(l)=Tx/DIAGONAL1
               l=2*ix+2+2*iy*(n+1)
               DIAGONAL2=-2.D0*GLv*(1.D0/(Dx*Dx)+4.D0/(Dy*Dy))
               a(l,Ld+2)=(GLv/(Dx*Dx)+GLvx/(2.D0*Dx))/DIAGONAL2
               a(l,kp)=((4.D0*GLv)/(Dy*Dy)+(2.D0*GLvy)/Dy)/DIAGONAL2
               a(l,Ld)=1.D0
               a(l,kn)=((4.D0*GLv)/(Dy*Dy)-(2.D0*GLvy)/Dy)/DIAGONAL2
               a(l,Ld-2)=(GLv/(Dx*Dx)-GLvx/(2.D0*Dx))/DIAGONAL2
               a(l,kp+1)=((3.D0*GLv)/(4.D0*Dx*Dy))/DIAGONAL2
               a(l,kn+1)=-a(l,kp+1)
               a(l,kp-3)=-a(l,kp+1)
               a(l,kn-3)=a(l,kp+1)
               a(l,Ld+1)=(GLvy/Dx)/DIAGONAL2
               a(l,Ld-3)=-a(l,Ld+1)
               a(l,kp-1)=(GLvx/(2.D0*Dy))/DIAGONAL2
               a(l,kn-1)=-a(l,kp-1)
               b(l)=Ty/DIAGONAL2
10        CONTINUE 
 9     CONTINUE 

       CALL sistbanda(a,b,nincogn,nbanda,Li,vel)
CC  Boundary conditions 
C	    DO 557 ix=0,n
C	    	kvel=2*ix+1
C 557		WRITE(6,565) ix,0,1,vel(kvel),vel(kvel+1)		
C	    DO 558 iy=1,m-1
C	    	kvel=2*n+1+2*iy*(n+1)
C 558		WRITE(6,565) n,iy,4,0.0,0.0		
CC 558		WRITE(6,*) n,iy,4,vel(kvel),vel(kvel+1)		
C	    DO 559 ix=n,0,-1
C	    	kvel=2*ix+1+2*m*(n+1)
C 559		WRITE(6,565) ix,m,1,0.0,0.0		
CC 559		WRITE(6,565) ix,m,1,vel(kvel),vel(kvel+1)		
C	    DO 560 iy=m-1,1,-1
C	    	kvel=1+2*iy*(n+1)
C 560		WRITE(6,565) 0,iy,4,0.0,0.0				
C 560		WRITE(6,*) 0,iy,vel(kvel),vel(kvel+1)				
C 565	FORMAT(3I4,2G15.6)
      
           IF(TEMPSMa.EQ.0.D0.AND.niter.EQ.1) GOTO 19 
           IF(nitermax.LE.1) GOTO 19 
           DO 20 kvel=1,nincogn 
 20              vel(kvel)=alfa*vel(kvel)+(1.D0-alfa)*velold(kvel)
 19    CONTINUE
	epmax=0.D0

       IF(nitermax.EQ.0) THEN
	  PRINT*,'DELTA PRESSIO: minim =',DP_min/1.D6,'MPa (10**6 N/m2)'
	  PRINT*,'               maxim =',DP_max/1.D6,'MPa (10**6 N/m2)'  
	      PRINT*,' Dvis_min =',Dvis_min, ', Dvis_max =',Dvis_max 
              PRINT*,' NO ITERO. LA VISCOSITAT ES LA MATEIXA'
	      epmig=0
	      cepmig=0
              DO 42 iy=1,m-1
                  DO 42 ix=1,n-1
                      kxy=ix+1+iy*(n+1)
                      kd=2*ix+1+2*iy*(n+1)
                      kdp=kd+2*(n+1)
                      kdn=kd-2*(n+1)
                      ux=(vel(kd+2)-vel(kd-2))/(2.D0*Dx)
                      vx=(vel(kd+3)-vel(kd-1))/(2.D0*Dx)
                      vy=(vel(kdp+1)-vel(kdn+1))/(2.D0*Dy)
                      uy=(vel(kdp)-vel(kdn))/(2.D0*Dy)
                      epuntzz=-(ux+vy)
	    	 Esecond=((ux*ux+vy*vy+epuntzz*epuntzz)/2.D0)+
     +                            ((1.D0/4.D0)*(uy+vx)*(uy+vx))	
     	    	      epeffec=DSQRT(Esecond)
	    	      IF(ABS(epeffec).GT.epmax) epmax=ABS(epeffec)
	    	      epmig=epmig+ABS(epeffec)
		      cepmig=cepmig+1
 42          CONTINUE
	     epmig=epmig/cepmig
             GOTO 999
       ENDIF
       IF(niter.EQ.1) WRITE(6,60) visinf,vissup
 60    FORMAT(/' ITERATION',32X,'______ limited points for ______'/
     +		2X'NUMBER',6X,'Dumig',8X,'Dvismig',8X,
     +		'inf. viscosity',4X,'sup. viscosity'/
     +		42X,1P,G9.2,' Pa.s',4X,1P,G9.2,' Pa.s' ) 

       NVISSUP=0
       NVISINF=0
C ****************  TROBO LA NOVA VISCOSITAT *********************
       epmig=0
       Kepmig=0
       DO 25 iy=1,m-1
         DO 25 ix=1,n-1
	    kxy=ix+1+iy*(n+1)
            kd=2*ix+1+2*iy*(n+1)
            kdp=kd+2*(n+1)
            kdn=kd-2*(n+1)
            ux=(vel(kd+2)-vel(kd-2))/(2.D0*Dx)
            vx=(vel(kd+3)-vel(kd-1))/(2.D0*Dx)
            vy=(vel(kdp+1)-vel(kdn+1))/(2.D0*Dy)
            uy=(vel(kdp)-vel(kdn))/(2.D0*Dy)
            epuntzz=-(ux+vy)
C            Epunt=DSQRT(2.D0*(ux*ux+vy*vy+ux*vy+
C     +                            (1.D0/4.D0)*(uy+vx)*(uy+vx)))

CC  !!!  Controlar quin strain rate utilitzar per calcular la viscositat !!!
	    Esecond=((ux*ux+vy*vy+((ux+vy)*(ux+vy)))/2.D0)+
     +                            ((1.D0/4.D0)*(uy+vx)*(uy+vx))	
     	    epeffec=DSQRT(Esecond)
	    epmig=epmig+ABS(epeffec)
	    cepmig=cepmig+1
	    IF(ABS(epeffec).GT.epmax) epmax=ABS(epeffec)
C	    exy=(uy+vx)/2.D0
	    visnova=ABS(visTer(kxy)/epeffec)
            vis(kxy)=alfa*visnova+(1.D0-alfa)*visold(kxy)
            IF(vis(kxy).GT.vissup) THEN
            	vis(kxy)=vissup
           	NVISSUP=NVISSUP+1
            ENDIF
            IF(vis(kxy).LT.visinf) THEN
            	vis(kxy)=visinf
           	NVISINF=NVISINF+1
            ENDIF
 25    CONTINUE
	epmig=epmig/cepmig

      Dumig=0.D0
      Dumax=0.D0
      velmax=0.D0
      vismax=0.D0
      vismin=1.D35
      Dvismig=0.D0
      Dvismax=0.D0
      DO 29 iy=1,m-1
          DO 29 ix=1,n-1
              kd=ix+1+iy*(n+1)
              k=2*ix+1+2*iy*(n+1)
              velmod=DSQRT(vel(k)*vel(k)+vel(k+1)*vel(k+1))
              veloldm=DSQRT(velold(k)*velold(k)+velold(k+1)*velold(k+1))
              Du=(ABS(velmod-veloldm))/velmod  
              Dumig=Dumig+Du  
              Dumax=MAX(Du,Dumax)
              Dvis=(ABS(visold(kd)-vis(kd)))/vis(kd)
              Dvismig=Dvismig+Dvis
              Dvismax=MAX(Dvis,Dvismax) 
              vismax=MAX(vis(kd),vismax)
              vismin=MIN(vismin,vis(kd))
              velmax=MAX(velmod,velmax) 
C            WRITE(6,63) ix,iy,Du/velmod
C 63         FORMAT(10X,'(ix,iy)',2I4,' Du/velmod=',1F10.5)
29    CONTINUE
      Dvismig=Dvismig/((m-1)*(n-1)) 
      Dumig=Dumig/((m-1)*(n-1))
      IF(velmax.EQ.0.D0.OR.vismax.EQ.0.D0) THEN
              PRINT*,'             EL CAMP DE VELOCITATS ES NUL'
              PRINT*,' Terme maxim de la velocitat (velmax) =',velmax,
     +           '     Terme maxim de la viscositat (vismax) =',vismax
              GOTO 999 
      ENDIF
C       tall=Dumax/velmax
C       tallvis=Dvismax/vismax       
       tall=Dumig
       tallvis=Dvismig
       WRITE(6,65) niter,tall,tallvis,NVISINF,NVISSUP 
 65    FORMAT(2X,I3,5X,F10.5,4X,F10.5,12X,I4,14X,I4)  
       if(niter.ge.nitermax) THEN    
             PRINT*,'  MAXIMUM NUMBER OF ITERATIONS'
             GOTO 999
       ENDIF
       IF(tall.GT.tallmax) then
          k=1
          DO 27 kxy=1,nn
               visold(kxy)=vis(kxy)
               velold(k)=vel(k)
               velold(k+1)=vel(k+1)
               k=k+2
27        CONTINUE
          GOTO 111
       ENDIF               
       PRINT*,'     '
       PRINT*,'      VELOCITY CONVERGED'
999    CONTINUE 


CC    FIXO LA VISCOSITAT I VERTICAL STRAIN RATE A LES VORES 
             DO 35 ix=1,n-1
                 kxys=ix+1
                 kxyn=ix+1+m*(n+1)
                 vis(kxys)=vis(kxys+n+1)
                 vis(kxyn)=vis(kxyn-n-1)
 35          CONTINUE
             DO 37 iy=0,m
                 kxyw=1+iy*(n+1)
                 kxye=n+1+iy*(n+1)
                 vis(kxyw)=vis(kxyw+1)
                 vis(kxye)=vis(kxye-1)
 37          CONTINUE

C ----------------------------------------------------------------------
      strainmx=epmax*Dt
      WRITE(6,61) TITOL_BC,vismin,vismax,epmax,epmig,strainmx,epmig*Dt
 61   FORMAT(/4X,'Boundary Conditions: ',A80,/
     +    6X,'CALCULATED THE VELOCITY FIELD:'/
     +    6X,'VISCOSITY: minimum =',1P,G12.4,' Pa.s,   maximum =',
     +    	1P,G12.4,' Pa.s'/
     +	6X,'EFFECTIVE STRAIN RATE   maximum:',1P,G10.2,' s-1,  medium:',
     +		1P,G10.2,' s-1'/6X,
     +	'EFFECTIVE STRAIN   maximum:',1P,G10.2,'    medium:',1P,G10.2/)
     	 
          IF(strainmx.GT.quotadef) THEN
               Dt=(quotadef/epmax)-5.D4
               Dtany=Dt/FACTEMP
               PRINT*,'LA DEFORMACIO ES MASSA GRAN -> DISMINUEIXO '
               PRINT*,'L`INCREMENT DE TEMPS  Dt =',Dtany,' anys'
               strainmx=epmax*Dt
               PRINT*,'         NOU STRAIN MAXIM :',strainmx
          ENDIF

CC	Divideixo el vector vel (velocitat) en dos vectors:
CC	u (velocitat en x) i v (velocitat en y).
	kvel=1
	DO 95 iy=0,m
	    DO 95 ix=0,n
		kxy=ix+1+iy*(n+1)
		u(kxy)=vel(kvel)
		v(kxy)=vel(kvel+1)
		kvel=kvel+2
 95	CONTINUE

C	deallocate(b)
      
       RETURN
       END

       
CC ********************************************************************
CC ********************************************************************
CC		SUBROUTINE READ_BC
CC   Read the Boundary Conditions from the file: BC.in
CC   First line: Title
CC	ix,iy,ITBC,t1,t2 
CC	ITBC: Boundary Condition type.
CC	t1,t2: Boundary Condition 1 and 2.

CC  ITBC=1	velocity fixed ->  vel_x=t1  and  vel_y=t2 
CC  ITBC=12	stress xx and xy fixed  -> tau_xx=t1  and  tau_xy=t2 
CC  ITBC=13	stress xx and yy fixed  -> tau_xx=t1  and  tau_yy=t2 
CC  ITBC=23	stress xy and yy fixed  -> tau_xy=t1  and  tau_yy=t2 
CC  ITBC=4	free slip (vel normal=0, tau_xy=0) 

CC  Stress free -> stress normal i xy nuls.

      SUBROUTINE READ_BC (TEMPSMa,niter,n,m,nn,nincogn,nbanda,Dx,Dy,
     +				vis,a,b,TITOL_BC,BCfile)   
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a(nincogn,nbanda),b(nincogn),vis(nn)
      CHARACTER*84 TITOL_BC,BCfile
      PARAMETER (NBCmax=600,FACVEL=3.1536D10)

      OPEN(4,file=BCfile)
      NBC=2*(n+1)+2*(m-1)
      Li=2*(n+1)+3
      Ls=Li
      Ld=Li+1
      kp=Ld+2*(n+1)
      kn=Ld-2*(n+1)

       READ(4,1001) TITOL_BC
 1001  FORMAT(A80)
       IF(TEMPSMa.EQ.0.AND.niter.EQ.1) WRITE (6,3) TITOL_BC
 3     FORMAT (/2X,'READ BOUNDARY CONDITIONS: ',/4X,A80/
     +    4X,'(ix, iy, BC type:   u,v / tau_1,tau_2 )')
 33     FORMAT (1X,3I4,1P,2G14.6)
 
       NBC_read=0
       DO 7 inBC=1,NBCmax
           READ(4,*,IOSTAT=IOS) ix,iy,ITBC,t1,t2
           IF (IOS < 0) THEN
               GO TO 120
           ELSE IF (IOS > 0) THEN
               WRITE (*, 90) NDATA
 90            FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.')
               STOP
           ENDIF
C       		IF(TEMPSMa.EQ.0.AND.niter.EQ.1) 
C     +			WRITE (6,33) ix,iy,ITBC,t1,t2
           NDATA = NDATA + 1
	   NBC_read=NBC_read+1
	   leq1=2*ix+1+2*iy*(n+1)
	   leq2=2*ix+2+2*iy*(n+1)
	   kxy=ix+1+iy*(n+1)
	   kp=Ld+2*(n+1)
           kn=Ld-2*(n+1)
C	   PRINT*,'ix,iy',ix,iy,'kxy, leq1:',kxy,leq1
C	   PRINT*,'            ITBC,t1,t2:',ITBC,t1,t2
CCC ====================================================
CCCC	ITBC=1  fixo (u,v)
	    IF(ITBC.EQ.1) THEN
C	    	PRINT*,'A: ix,iy:',ix,iy,'  ITBC:',ITBC
	   	vel_x=t1
		vel_y=t2
             	a(leq1,Ld)=1.D0
             	b(leq1)=vel_x
             	a(leq2,Ld)=1.D0
             	b(leq2)=vel_y
		GOTO 99
	    ENDIF	   
CCC ====================================================
CCCC	ITBC=12  fixo:  tau_xx (eq.1) i  tau_xy (eq.2)
CCCC	ITBC=13  fixo:  tau_xx (eq.1) i  tau_yy (eq.2)
CCCC	ITBC=23  fixo:  tau_xy (eq.1) i  tau_yy (eq.2)
	    IF(ITBC.GT.11.AND.ITBC.LT.20) THEN
C	    	PRINT*,'B: ix,iy:',ix,iy,'  ITBC:',ITBC
	   	tau_xx=t1
		GOTO 411
	    ENDIF
 412	    CONTINUE
	    IF(ITBC.EQ.12) THEN
C	    	PRINT*,'C: ix,iy:',ix,iy,'  ITBC:',ITBC
		tau_xy=t2
		ieq=1
		GOTO 522
	    ENDIF
	    IF(ITBC.EQ.13) THEN
C	    	PRINT*,'D: ix,iy:',ix,iy,'  ITBC:',ITBC
		tau_yy=t2
		GOTO 813
	    ENDIF
	    IF(ITBC.EQ.23) THEN
C	    	PRINT*,'E: ix,iy:',ix,iy,'  ITBC:',ITBC
		tau_xy=t1
		ieq=0
		GOTO 522
	    ENDIF
 525	    CONTINUE
	    IF(ITBC.EQ.23) THEN
C	    	PRINT*,'F: ix,iy:',ix,iy,'  ITBC:',ITBC
		tau_yy=t2
		GOTO 813
	    ENDIF
CCC ====================================================
CCCC	ITBC=4  free slip (v norm=0, dv(tang)/dx(tang) =0)
	    IF(ITBC.EQ.4) THEN
		IF(iy.EQ.0) THEN
C	    	    PRINT*,'G: ix,iy:',ix,iy,'  ITBC:',ITBC
             	    a(leq1,kp)=1.D0
             	    a(leq1,Ld)=-1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(iy.EQ.m) THEN
C	    	    PRINT*,'H: ix,iy:',ix,iy,'  ITBC:',ITBC
             	    a(leq1,Ld)=1.D0
             	    a(leq1,kn)=-1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(ix.EQ.0) THEN
C	    	    PRINT*,'I: ix,iy:',ix,iy,'  ITBC:',ITBC
             	    a(leq1,Ld)=1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld+2)=1.D0
             	    a(leq2,Ld)=-1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
		IF(ix.EQ.n) THEN
C	    	    PRINT*,'J: ix,iy:',ix,iy,'  ITBC:',ITBC
             	    a(leq1,Ld)=1.D0
             	    b(leq1)=0.D0
             	    a(leq2,Ld)=1.D0
             	    a(leq2,Ld-2)=-1.D0
		    b(leq2)=0.D0
		    GOTO 99
	    	ENDIF		
	    ENDIF		

CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition   tau_xx:  Equation 1
 411		CONTINUE
		IF(ix.NE.0.AND.ix.NE.n) THEN	
		    a(leq1,Ld+2)=1.D0
		    a(leq1,Ld-2)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/vis(kxy)
		ENDIF
		IF(ix.EQ.0) THEN	
		    a(leq1,Ld+2)=1.D0
		    a(leq1,Ld)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/(2.D0*vis(kxy))
		ENDIF
		IF(ix.EQ.n) THEN	
		    a(leq1,Ld)=1.D0
		    a(leq1,Ld-2)=-1.D0
             	    b(leq1)=(Dx*tau_xx)/(2.D0*vis(kxy))
		ENDIF
		GOTO 412
CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition	tau_xy:  Equation 1 (ieq=0),	  Equation 2 (ieq=1)
 522		CONTINUE
		leq=leq1+ieq
		IF(iy.NE.0.AND.iy.NE.m) THEN
		    a(leq,kp-ieq)=1.D0/(2.D0*Dy)
		    a(leq,kn-ieq)=-1.D0/(2.D0*Dy)
		ENDIF
		IF(iy.EQ.0) THEN
		    a(leq,kp-ieq)=1.D0/Dy
		    a(leq,Ld-ieq)=-1.D0/Dy
		ENDIF
		IF(iy.EQ.m) THEN
		    a(leq,Ld-ieq)=1.D0/Dy
		    a(leq,kn-ieq)=-1.D0/Dy
		ENDIF
		IF(ix.NE.0.AND.ix.NE.n) THEN
		    a(leq,Ld+3-ieq)=1.D0/(2.D0*Dx)
		    a(leq,Ld-1-ieq)=-1.D0/(2.D0*Dx)
		ENDIF		    
		IF(ix.EQ.0) THEN
		    a(leq,Ld+3-ieq)=1.D0/Dx
		    a(leq,Ld+1-ieq)=-1.D0/Dx
		ENDIF
		IF(ix.EQ.n) THEN
		    a(leq,Ld+1-ieq)=1.D0/Dx
		    a(leq,Ld-1-ieq)=-1.D0/Dx
		ENDIF
             	b(leq)=tau_xy/vis(kxy)
	    	IF(ieq.EQ.1) GOTO 99
	    	IF(ITBC.EQ.23) GOTO 525
CCCC---------------------------------------------------------
CCCC---------------------------------------------------------
CC  Condition   tau_yy:  Equation 2
 813		CONTINUE
		IF(iy.NE.0.AND.iy.NE.m) THEN	
		    a(leq2,kp)=1.D0
		    a(leq2,kn)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/vis(kxy)
		ENDIF
		IF(iy.EQ.0) THEN	
		    a(leq2,kp)=1.D0
		    a(leq2,Ld)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/(2.D0*vis(kxy))
		ENDIF
		IF(iy.EQ.m) THEN	
		    a(leq2,Ld)=1.D0
		    a(leq2,kn)=-1.D0
             	    b(leq2)=(Dy*tau_yy)/(2.D0*vis(kxy))
		ENDIF
	        GOTO 99
CCCC---------------------------------------------------------
		
 99	CONTINUE
 7      CONTINUE

 100  CONTINUE
      IF(NBC_read.EQ.0) PRINT*,' FILE NO READ. number of data:',NBC_read
      WRITE(6,110) NBCmax
 110  FORMAT(/' Dataset contained more than ',
     &          I5,' data; Check the dimensions.')
 120  CONTINUE
C      WRITE(6,130) NBC_read
 130  FORMAT(/5X,'Reading of data completed:',I7,' data points.'/)
 	
      IF(NBC_read.NE.NBC) THEN
            WRITE(6,135) NBC_read,NBC
 135        FORMAT(' El numero de punts llegits:',I7,' no coincideix',
     &          'amb 2(n+1)+2(m-1)=',I7 //'PROGRAMA ATURAT')
            STOP
      ENDIF
      CLOSE(4)

      RETURN
      END








CC ********************************************************************
CC ********************************************************************
CC	Calculation of the vertical strain rate and control that it should be
CC		lower than ep_limit.

      SUBROUTINE vertical_strain_rate (AX, BY, n, m, nn, u, v, epuntzz) 

CC  epuntzz [1/s]  Array returning the vertical strain rate.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ep_limit=5.D-15)
      DIMENSION u(nn), v(nn), epuntzz(nn)

      WRITE(6,51) ep_limit
 51   FORMAT(/4X,'LIMIT OF THE VERTICAL STRAIN RATE: ',1P,G12.4,' s-1')
      Dx=AX/n
      Dy=BY/m
      epuntzz_min=1.D30
      epuntzz_max=-1.D30
      PRINT*,'AX,BY,n,m,nn,nn',AX,BY,n,m,nn,nn
     
       DO 42 iy=1,m-1
         DO 42 ix=1,n-1
	    	kxy=ix+1+iy*(n+1)
		ux=(u(kxy+1)-u(kxy-1))/(2.D0*Dx)
		vy=(v(kxy+n+1)-v(kxy-n-1))/(2.D0*Dy)
		epuntzz(kxy)=-(ux+vy)
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
 42   CONTINUE

CC    CONTINIUM VERTICAL STRAIN RATE TO THE BOUNDARIES
             DO 35 ix=1,n-1
                 kxys=ix+1
                 kxyn=ix+1+m*(n+1)
                 epuntzz(kxys)=epuntzz(kxys+n+1)
                 epuntzz(kxyn)=epuntzz(kxyn-n-1)
 35          CONTINUE
             DO 37 iy=0,m
                 kxyw=1+iy*(n+1)
                 kxye=n+1+iy*(n+1)
                 epuntzz(kxyw)=epuntzz(kxyw+1)
                 epuntzz(kxye)=epuntzz(kxye-1)
 37          CONTINUE
CC	LIMIT TO THE VERTICAL STRAIN RATE
          DO 125 kxy=1,nn
C  		limit the calculated vert.str.rate 'epuntzz'.
    		IF(epuntzz(kxy).GT.ep_limit) epuntzz(kxy)=ep_limit
            	IF(epuntzz(kxy).LT.(-1.D0*ep_limit)) 
     +			epuntzz(kxy)=-1.D0*ep_limit
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
 125	CONTINUE   
 
      WRITE(6,61) epuntzz_min,epuntzz_max
 61   FORMAT(4X,'VERTICAL STRAIN RATE:  minim=',1P,G12.4,
     +	  ' s-1,   maxim=',1P,G12.4,' s-1'/)

      RETURN
      END

CC **********************************************************************
CC **********************************************************************
CC  New thickness of a layer with velocity field (u,v) during a time 
CC	interval Dt(s)
CC	D(thickness)/Dt=
CC		thickness*epuntzz-(u*d(thickness)/dx+v*d(thickness)/dy)= 
CC  (u,v) (m/s) : array of horizontal velocities.

       SUBROUTINE THICKEN (Dt, n, m, nn, AX, BY, u, v, epuntzz, 
     +			       thickness)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION u(nn),v(nn), epuntzz(nn), 
     +		thickness(nn),thickness_old(nn)

      Dx=AX/n
      Dy=BY/m
      DO 53 kxy=1,nn
	 thickness_old(kxy)=thickness(kxy)
 53   CONTINUE 
	 
      DO 13 iy=1,m-1
         DO 13 ix=1,n-1
            kxy=ix+1+iy*(n+1)
            TERME1=thickness_old(kxy)*epuntzz(kxy)
CC  SISTEMA DE COORDENADES LAGRANGIA (seguint el material) -> TERME2 = 0
C           TERME2=0.D0
            sx=(thickness_old(kxy+1)-thickness_old(kxy-1))/(2.D0*Dx)
            sy=(thickness_old(kxy+n+1)-thickness_old(kxy-n-1))/(2.D0*Dy)
            TERME2=u(kxy)*sx+v(kxy)*sy
            thickness(kxy)=thickness_old(kxy)+Dt*(TERME1-TERME2)
	    IF(thickness(kxy).LT.0.0) thickness(kxy)=0.D0
13    CONTINUE 
CC - No lateral variations of the thickening on the boundaries 
CC     ix=0 / ix=n   d(thickness)/dx=0
      DO 5 iy=1,m-1
         kxy0=1+iy*(n+1)
         kxyn=n+1+iy*(n+1)
         thickness(kxy0)=thickness(kxy0+1)
         thickness(kxyn)=thickness(kxyn-1)
 5    CONTINUE
CC     iy=0 / iy=m   d(thickness)/dy=0  
      DO 7 ix=0,n
         kxy0=ix+1
         kxym=ix+1+m*(n+1)
         thickness(kxy0)=thickness(kxy0+n+1)
         thickness(kxym)=thickness(kxym-n-1)
 7    CONTINUE
CC --------------------------------------------------------------------
       	thickmin=200.D3
       	thickmax=0.D0
         DO 210 kxy=1,nn 
                thickmin=MIN(thickmin,thickness(kxy))
                thickmax=MAX(thickmax,thickness(kxy))
 210     CONTINUE      
      WRITE(6,65) thickmin,thickmax
 65   FORMAT(3X,'MINIMUM THICKNESS:  ',F10.2,' m',3X,
     +        'MAXIMUM THICKNESS: ',F10.2,' m'/)

      RETURN
      END
