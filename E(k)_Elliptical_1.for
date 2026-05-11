CDec$ FixedFormLineSize:132
CDec$ Real:8
	
	Program Ek_Elliptical
	
		Use DfIMSL
		Use DFLib
		Use MsFLib 
C=======================================================================
	Real(8):: Ec(10001),  Yn(10001), Hst, pi, eps, Xec(10001), DxF1(10001),DxF2(10001)   !DxF - flange tg(fi) inclanation
	Real(8):: X1(10001),X2(10001), Kw1(10001),Kw2(10001), KwAgg(10001),DxFAgg(10001), RadevRb(10001), a, X0, R0
	Real(8), Allocatable:: Break1(:), CSCoef1(:,:), Break2(:), CSCoef2(:,:), X1U(:), XeU(:)
	Real(8), Allocatable:: BreakRe(:),CSCoefRe(:,:), RR(:), XX(:)
	Character(4):: Sign, Mode*6, Aagg, Ragg*4
C=======================================================================
	a=0.0; R0=16.0 ! а - доля нормально эволюционирующих пятен (удлиняющихся поперёк) в общеим числе; 
	             !Х0 - Начальное соотношение радиусов кривизны Ra/Rb  при нулевой производной 
	Mode='Base' ! сдвиг в Х0 раз от соотношения 1:1 в сторону нормальных пятен (вытянутых поперёк) ShiftN Base
	Hst=0.0001; pi=DConst('pi'); eps=1.E-24
	Write(Aagg,'(f0.2)') a; Write(Ragg,'(f0.1)') R0
C=======================================================================
C=======================================================================
	SELECT CASE (Mode)

	CASE('Base')

	Ec=(/(Hst*(I-1), I=1,10001)/) !эксцентриситет e
	Xec=Ec**2 ! аргумент функций эллиптических интегралов =e**2 
	
	X1=1./Sqrt(1.-Xec) ! a>b 
	X2=Sqrt(1.-Xec) ! a<b
	X2(10001)=eps; Xec(10001)=eps 

	Kw1(1)=2.; Kw2(1)=2.
	
	Do I=2,10001
	Kw1(I)=0.5*pi*Xec(I)/(dElK(Xec(I))-dElE(Xec(I)))
	Kw2(I)=0.5*pi*Xec(I)/(DSqrt(1.-Xec(I))*(dElE(Xec(I))-(1.-Xec(I))*dElK(Xec(I))))
	End Do
C=======================================================================
	Open(205, File='Kw1&Kw2.dat')
	
	Do I=10001, 2, -1
	Write(205,*) real(X2(I)), real(Kw2(I)) 
	End Do
	
	Do I=1, 10001
	Write(205,*) real(X1(I)), real(Kw1(I))
	End Do
	
	Close(205)
C=================================================
C=================================================
	X1 = Log10(X1); X2 = Log10(X2) 
	Kw1=Log10(Kw1); Kw2=Log10(Kw2)
	
	Open(206, File='Kw1&Kw2_Log.dat')
	
	Do I=10000, 2, -1
	Write(206,*) real(X2(I)), real(Kw2(I)) 
	End Do
	
	Do I=1, 10000
	Write(206,*) real(X1(I)), real(Kw1(I))
	End Do
	!Write(206,*) 4.5, Log0.
	Close(206)
C=================================================
C=======================================================================
	!Open(201, File='E(k).dat')
	!Write(201,*) 'TITLE     = "E(k)"'; Write(201,*) 'VARIABLES = "V1"'; Write(201,*) '"V2"'; Write(201,*) 'ZONE T="ZONE 001"'
	!Write(201,*)  'I=', 10001,',',' J=1, K=1,F=BLOCK'; Write(201,*) 'DT=(SINGLE SINGLE )'
	!Do I=1, 10001
	!Write(201,*) Xn(I), Yn(I)
	!EndDo
	Close(201)
!	Open(202, File='E(k)_Ellips.dat')
	!Write(202,*) 'TITLE     = "E(k)"'; Write(202,*) 'VARIABLES = "V1"'; Write(202,*) '"V2"'; Write(202,*) 'ZONE T="ZONE 001"'
	!Write(202,*)  'I=', 10001,',',' J=1, K=1,F=BLOCK'; Write(202,*) 'DT=(SINGLE SINGLE )'
!	Do I=1, 10001
!	Write(202,*) Xn(I), Yn(I)-1.-(0.5*pi-1.)*DSqrt(1.-Xn(I)*Xn(I))-0.32*DLog(1.-Xn(I)+eps)*DSqrt(1.-Xn(I))*DExp(-1./(Xn(I)+eps))-
!     -                         (-0.0034+0.0034*DCos(6.2832*1.2*Xn(I))+0.004*Xn(I)**8)
!	EndDo
	Close(202)
!	Yn=Yn-(1.+0.5*(2.*(1.-Xn)-(1.-Xn)**2.)*(DLog(4.)-0.5*DLog(2.*(1.-Xn)-(1.-Xn)**2.+eps)-0.5))-
!     -                    0.064*(1+DCos(pi*Xn))-0.0096*(0.97-DCos(2.*pi*Xn))-
!     -                    ((Xn-0.667)*(1.-Dcos(6.2831852*Xn))*0.02+0.0052)*DExp(-1./DSqrt(1.-Xn+eps)-1./DSqrt(Xn+eps))*20.
	!Open(203, File='E(k)_Ln_Cos.dat')
	!Write(203,*) 'TITLE     = "E(k)"'; Write(203,*) 'VARIABLES = "V1"'; Write(203,*) '"V2"'; Write(203,*) 'ZONE T="ZONE 001"'
	!Write(203,*)  'I=', 10001,',',' J=1, K=1,F=BLOCK'; Write(203,*) 'DT=(SINGLE SINGLE )'
	!Do I=1, 10001
	!Write(203,*) Xn(I), Yn(I) 
	!EndDo
	Close(203)
C=======================================================================
C==========================1/cos^2: normal spots========================
	X1=0.;X2=0.; Xec=0.; DxF1=0.;DxF2=0.; Kw1=0.;Kw2=0.; RadevRb=1. ! начальное пятно считаем для а/b=1 во всех случаях
	
	X1=(/(1.+0.005*(I-1), I=1,10001)/); X2=1./X1 ! отношение пятен от 1(X0) до 51 и от 1/51 до 1
	
	Xec=1.-1./X1**2 ! Для всех случаев ! Строим сетку отношений пятен. По ней строим сетку производных для обоих случаев
					! Потом по новой сетке строим коэффициенты. Потом их складываем с помощью аппроксимации сплайном
	
	DO I=2, 10001
	DxF1(I)=-DSQRT((Xec(I)*dElK(Xec(I))/(dElK(Xec(I))-dElE(Xec(I)))-1.)*X1(I)**2-1.) !конформный контакт
	End DO !DxF1(1)=0.
C==========================cos: abnormal spots==========================	
	DO I=2, 10001
	DxF2(I)=-DSQRT(((Xec(I)*dElK(Xec(I))/(dElK(Xec(I))-dElE(Xec(I)))-1.)*X1(I)**2)**2-1.) ! произвольный (точечный) контакт
	END DO !DxF2(1)=0.
C=======================================================================
	Do I=2, 10001
	RadevRb(I)=(dElE(Xec(I))*X1(I)**2-dElK(Xec(I)))/(dElK(Xec(I))-dElE(Xec(I)))
	EndDo
C=======================================================================
	DFmax=Abs(DxF1(10001)); Hdxagg=DFmax/10000.
	ND2=MinLoc(DxF2,1,DxF2.GE.-DFmax)
	DxFAgg=-(/(Hdxagg*(I-1), I=1, 10001)/) ! универсальная сетка для суперпозиции коэффициентов

	Kw1(1)=2.; Kw2(1)=2.

	Do I=2,10001
	Kw1(I)=0.5*pi*Xec(I)/(dElK(Xec(I))-dElE(Xec(I)))
	Kw2(I)=0.5*pi*Xec(I)/(DSqrt(1.-Xec(I))*(dElE(Xec(I))-(1.-Xec(I))*dElK(Xec(I))))
	End Do
C==================Akima spline approximation===========================
	Allocate(Break1(10001),CSCoef1(4,10001), Break2(ND2), CSCoef2(4,ND2))

	CALL DCSAkm(10001, DxF1, Kw1, Break1, CSCoef1)
	CALL DCSAkm(ND2,DxF2(1:ND2),Kw2(1:ND2),Break2,CSCoef2)
C-----------------------------------------------------------------------
	DO I=1, 10001	
	KwAgg(I)=a*DCSVal(DxFAgg(I), 10000, Break1,CSCoef1)+(1.-a)*DCSVal(DxFAgg(I), ND2-1, Break2,CSCoef2)
	End Do

C=======================================================================
C=======================================================================
	CASE('ShiftN')

	X1=0.;X2=0.; Xec=0.; DxF1=0.;DxF2=0.; Kw1=2.;Kw2=2.
	Allocate(BreakRe(10001), CSCoefRe(4,10001), RR(10001), XX(10001))
	RR=1.; XX=1.
C=======================================================================
	XX=(/(1.+0.0025*(I-1), I=1,10001)/)
	Xec=1.-1./XX**2
	Do I=2, 10001
	RR(I)=(dElE(Xec(I))*XX(I)**2-dElK(Xec(I)))/(dElK(Xec(I))-dElE(Xec(I)))
	EndDo

	CALL DCSAkm(10001, RR, XX, BreakRe, CSCoefRe)  ! чтобы определить соотношение пятен Х0 по отношению радиусов К0
C=======================================================================
	X0=DCSVal(R0, 10000, BreakRe, CSCoefRe)
	NX0=Int(10000*(X0-1.)/25.)
	Allocate(X1U(NX0+1), XeU(10001)); X1U=0.; XeU=0.
	
	X1U=(/(X0-0.0025*(I-1), I=1, NX0+1)/) ! от X0 до 1
	X1= (/(X0+0.0025*(I-1), I=1, 10001)/) ! от Х0 и далее до 27
	X2= (/(1./(1.+0.0025*(I-1)), I=1, 10001)/) ! от 1 до ~1/25

	Xec=1.-1./X1**2                                                   ! Нижняя нормальная ветвь для Kw1
	XeU(1 : NX0+1) = (/(1.-1./X1U(I)**2, I=1, NX0+1)/)                ! Верхняя нормальная ветвь для Kw1
	XeU(NX0+2:10001)=(/(1.-X2(I-NX0-1)**2, I=NX0+2, 10001)/)          ! Верхняя абнормальная ветвь для Kw2 
C=======================================================================
	DO I=2, 10001
	DxF1(I)=-DSQRT((Xec(I)*dElK(Xec(I))/(dElK(Xec(I))-dElE(Xec(I)))-1.)*X1(I)**2/R0-1.) ! конформный контакт, нижняя ветвь
	END DO 

	DFmax=Abs(DxF1(10001)); Hdxagg=DFmax/10000.
	DxFAgg=-(/(Hdxagg*(I-1), I=1, 10001)/) ! универсальная сетка для суперпозиции коэффициентов
C-----------------------------------------------------------------------
	DO I=2, NX0+1    ! химера: отношение радиусов прямое, но меняется радиус Rb
	DxF2(I)=-DSQRT(R0**2/((XeU(I)*dElK(XeU(I))/(dElK(XeU(I))-dElE(XeU(I)))-1.)*X1U(I)**2)**2-1.) ! точечный контакт с регрессом отношения a:b от 2:1 к 1:1
	END DO

	DxF2(NX0+2)=-DSQRT((R0*X2(1)**2)**2-1.)
	DO I=NX0+3, 10001
	IF(DxF2(I-1).LT.-DFmax) EXIT
	DxF2(I)=-DSQRT(R0**2*((XeU(I)*dElK(XeU(I))/(dElK(XeU(I))-dElE(XeU(I)))-1.)/X2(I-NX0-1)**2)**2-1.) ! абнормальная часть верхней ветви, начинается с отношения 1:1, как обычно. Точечный контакт
	END DO
	ND2=I-1
C=======================================================================
	IF(R0.NE.1.) Kw1(1)=0.5*pi*Xec(1)/(dElK(Xec(1))-dElE(Xec(1)))
	Do I=2,10001 ! начинаем с I=2 на случай R0=1, когда эллиптические интегралы дают сбой на концах значений аргумента
	Kw1(I)=0.5*pi*Xec(I)/(dElK(Xec(I))-dElE(Xec(I)))
	End Do

	Do I=1, NX0+1
	Kw2(I)=0.5*pi*XeU(I)/(dElK(XeU(I))-dElE(XeU(I)))
	If(dIfNaN(Kw2(I))) Kw2(I)=2.0
	End Do
	Do I=NX0+2, ND2
	Kw2(I)=0.5*pi*XeU(I)/(DSqrt(1.-XeU(I))*(dElE(XeU(I))-(1.-XeU(I))*dElK(XeU(I))))
	If(dIfNaN(Kw2(I))) Kw2(I)=2.0
	End Do
C=======================================================================
	Allocate(Break1(10001),CSCoef1(4,10001), Break2(ND2), CSCoef2(4,ND2))

	CALL DCSDec(10001,DxF1, Kw1, 0, 0., 1, 0., Break1, CSCoef1)
	CALL DCSDec(ND2,DxF2(1:ND2),Kw2(1:ND2),0,0.,1,0.,Break2,CSCoef2)

	DO I=1, 10001	
	KwAgg(I)=a*DCSVal(DxFAgg(I), 10000, Break1,CSCoef1)+(1.-a)*DCSVal(DxFAgg(I), ND2-1, Break2,CSCoef2)
	End Do

	END SELECT
C=======================================================================
C=======================================================================
	Open (210, File='KwAgg_FlangeIncl_'//Aagg//'_'//Trim(Ragg)//'.dat')
	Do I=1, 10001
	WRITE(210, *) DxFAgg(I), KwAgg(I)
	End Do
	Close(210)

	Open (211, File='Kw1_FlangeIncl.dat')
	Do I=1, 10001
	WRITE(211, *) DxF1(I), Kw1(I)
	End Do
	Close(211)

	Open (212, File='Kw2_FlangeIncl.dat')
	Do I=1, ND2
	WRITE(212, *) DxF2(I), Kw2(I)
	End Do
	Close(212)

	Open (213, File='Ra%Rb_a%b.dat')
	!Do I=1, 10001
	write (213, *) '/'
	Do I=1, 2000
	WRITE(213, '(5x, A1,2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1)') '&', 
     &	  RadevRb((I-1)*5+1),',',RadevRb((I-1)*5+2),',',RadevRb((I-1)*5+3),',',RadevRb((I-1)*5+4),',',RadevRb(I*5),',' !(I)
	End do
	write (213, *) '/'
		WRITE(213, *) '     '
	write (213, *) '/'
	Do I=1, 2000
	WRITE(213, '(5x, A1,2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1, 2x, F0.16,A1)') '&',
     &      X1((I-1)*5+1),',', X1((I-1)*5+2),',', X1((I-1)*5+3),',', X1((I-1)*5+4),',', X1(I*5),','
	End do
	write (213, *) '/'
	!End Do
	Close(213)

	Open (214, File='LgRa%Rb_Lga%b.dat')
	Do I=1, 10001
	!WRITE(214, *) Log10(RadevRb(I)), Log10(X1(I))
	End Do
	Close(214)
C=======================================================================
	Print *, 'End of calculation'
	Pause 600
C=======================================================================
	End Program 