!---------------------------------------------------------------------------------------------
! Subrotina para calcular o valor ideal de:
!1- dtcal (intervalo de tempo de calculo da propagação)
!2- ntr (numero de subtrechos em que sera dividida a propagação)
!---------------------------------------------------------------------------------------------

subroutine parcunge
use vars_main,only: nc,rugman,icodmusk,dtp,acel,adren,decl,qref,brio,srio,ntrmax,icellfp,dtfp, &
                  & cel,dt,ntrc,nfp
implicit none
integer ntr,ic
real coef,dtcal,dx      	!variaveis auxiliares	
real along,hflow,decliv  	!variaveis auxiliares

icodmusk=0 !por principio considera todas as propagacoes por muskingun cunge linear
dtcal=max(dtp/24.,600.) ! (1 hora dados diarios, 10' dados horarios)

do ic=1,nc !loop das células

	if(acel(ic)>0.1.and.adren(ic)/acel(ic)>1.5)then ! celulas com rio
		coef=1.67*decl(ic)**0.3/rugman**0.6
		cel(ic)=coef*(qref(ic)/brio(ic))**0.4
		along=srio(ic)*1000. !converte para metros 
		dx=0.5*cel(ic)*dtcal*(1+sqrt(1+3.068*qref(ic)/(2.*decl(ic)*brio(ic)*cel(ic)*cel(ic)*dtcal)))
		ntr=nint(along/dx) !calcula numero de trechos
		ntr=max(ntr,1)
				
		if(ntr > ntrmax)then 
			write(*,*) 'numero de sub-trechos maior que', ntrmax,' na rotina parcunge!!!'
			stop 'parou '
		endif
		
		dt(ic)=dtcal
		ntrc(ic)=ntr
	else !areas fonte (sinks), sem rio, nao interessam
		dt(ic)=0.0
		ntrc(ic)=0
	endif

enddo !fim do loop das células

if(nfp > 0) then ! existem celulas com planicie de alagamento	
    ! define o passo de tempo para o calculo na planicie (hunter et al. 
    dx=sqrt(acel(icellfp(1))) ! distancia entre celulas da planicie
    decliv=0.01/dx ! declividade minima entre celulas de planicie (admitindo 1 cm na diferenca de cota entre elemntos)
    hflow=2. ! maxima cota no elemento de planicie em metros
    dtfp=dx**2/4.*(2*rugman*decliv**0.5/hflow**1.667)
    dtfp=dtp
endif
	
return
end