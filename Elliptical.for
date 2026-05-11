CDec$ FixedFormLineSize:132
	
	Program Elliptical_Asympt
	
		Use DfIMSL
C=======================================================================
	Real(8):: X, xN, Y, p
	Character(4):: Sign ! y=1+0.5*(2*(1-x)-(1-x)**2)*(log(4)-0.5*log(2*(1-x)-(1-x)**2)-0.5)
C=======================================================================
2	Print *, 'Enter first number'
	Accept *, X
1	Print *, 'Enter operation sign (+,-,"/",*,pwr,ln,log10,sin,cos,tg,ctg,arcsin,arccos,arcctg,arctg)'
	Accept *, Sign
	Print *, 'Enter second number'
	Accept *, Y
	Print *, X, Sign, Y
C=======================================
	Select case(Trim(Sign))
	Case('+')
	X=X+Y
	Case('-')
	X=X-Y
	Case('/',':')
		If(Y.ne.0.) Then
	X=X/Y
		Else
	Print *, 'Error in input data!'
		EndIf
	Case('*')
	X=X*Y
	Case('pwr')
		If(X.lt.0.or.Y.Lt.0) Then
	Print *, 'Error in data!'
		Else
	xN=X; p=1._8/Y
	X=dble(X)**dble(Y)

	X=X-(X**p-xN)/(p*X**(p-1._8))
		End If

	Case('ln','Ln','LN','Log','LOG','LOg','loG','lN','log')
		If(X.gt.0) Then
	X=Y*Log(X)
		Else
	Print *, 'Error in data!'
		End If
			 
	Case Default
	Print *, 'Error in data'
	End Select

	Print *, 'Result=', dble(X)
	Print *, 'Continue? Y/N/S'
	Accept *, Sign
	
	If(Trim(Sign).eq.'Y') Then
	Go to 1
	ElseIf(Trim(Sign).eq.'N') Then
	Print *, ''; Print *, 'New calculation...'
	Go To 2
	Else
	End If


	Print *, 'End of calculation'
C=======================================================================
	End Program 