!---------------------------------------------------------------------------------------------
! Subrotina que calcula os escoamentos superficial e subsuperficial
!---------------------------------------------------------------------------------------------
subroutine escoamentos(ic,iu,ib)
use vars_main,only: dtp,neta,csi,smax,ssub,lambdai,lambdam,mu,asat,ac,sssu,ssmax,pre,ec,dss, &
                  & dd,dsup,kss,lambda,lambdan,alpha,d1,d2,tanb,srad,srmax,tsub,dsub
implicit none
integer ic,iu,ib,k
real saux,sfc,dsaux,deltat,nmenos,aux,mumenos,dsmax,lambdax

! variables de entrada
! pre : precipitacao (mm)
! ec  : lamina interceptada
! sssu: armazenamento na camada superior do solo (mm)
! ssub: armazenamento total (mm)

! parametros
! kss : condutividade hidraulica da camada superior do solo (m/dia)
! tsub : transmissividade para o perfil do solo totalmente saturado (m2/dia)
! neta : parametro neta de brooks corey
! smax : armazenamento medio na camada superior do solo (mm)
! ac(i) e lambdai(i) histograma do índice topográfico da celula
! alpha: fator de anisotropia para condutividade hidráulica horizontal 
! csi  : percentual minimo de armazenamento subterraneo para início do escoamento base
! ac   : area de contribuicao media da celula
! d1   : espessura horizonte superior do solo
! d2   : espessura horizonte inferior do solo 

!variaveis de saida
!dsup : escoamento superficial
!dss  : escoamento subsuperficial
!dsub : escaomento subterraneo
!dd   : drenagem profunda
!ssub : armazenamento de todo o perfil do solo atualizado
!sssu : armazenamento na camada superior do solo (atualizado)

!variaveis internas
!nmenos : neta -1
!aux: auxiliar 
!asat: area saturada
!deltat	: intervalo de tempo em dias

deltat=dtp/86400.
nmenos = neta(iu) - 1
sfc=csi(ib,iu)*smax(ib,iu)
dsmax=smax(ib,iu)-sfc

! determina a area de contribuicao variavel usando o armazenamento subterraneo do inicio do intervalo
dsaux=ssub(ic,iu)-sfc
lambdax=lambdai(ic,1)	
if(dsaux > 0.) then
    aux=lambdam(ic,iu)*(dsmax/dsaux)**mu(ib,iu) 
	lambdax=min(aux,lambdax)
endif
if(lambdax >= lambdai(ic,1)) then
	asat(ic,iu)=0.
else
	asat(ic,iu)=1.
	do k=2,49
		if(lambdax >= lambdai(ic,k)) then
			asat(ic,iu)=ac(ic,k-1)+(lambdax-lambdai(ic,k-1))*(ac(ic,k)-ac(ic,k-1))/(lambdai(ic,k)-lambdai(ic,k-1))
			exit
		endif
	enddo
endif

! calcula fluxo nao saturado 
! estima o armazenamento ponderado na primeira camada de solo
! lamina de agua que infiltra na area nao saturada = (1-asat)*(pre-ec)
! armazenamento na area saturada = asat*ssmax 	
! se ec(ic,iu) < 0 - orvalho, pre-ec > ec, e soma a precipitacao
saux=sssu(ic,iu)+asat(ic,iu)*ssmax(ib,iu)+(1-asat(ic,iu))*(pre(ic)-ec(ic,iu))
if(saux <= 0.) then
	sssu(ic,iu) = 0.
	dss = 0.
	dd = 0.
    dsup=0.
else
    if(saux > ssmax(ib,iu)) then
        ! escoamento superficial por saturacao fora da area de contribuicao variavel
	    dsup=saux-ssmax(ib,iu)
	    saux=ssmax(ib,iu)
    else
        dsup=0.
    endif
	aux = (1+nmenos*(1000*kss(ib,iu)*deltat/ssmax(ib,iu))*(lambda(ic,iu)/lambdan(ic,iu)) &
        & *(saux/ssmax(ib,iu))**nmenos)**(1/nmenos)
	sssu(ic,iu)=saux/aux
	! calcula a drenagem total em mm
	aux=saux-sssu(ic,iu)
	! calcula escoamento subsuperficial
	dss = min(aux,alpha(ib,iu)*d1(ib)*tanb(ic)*aux/lambda(ic,iu))
	dd= aux-dss ! calcula a drenagem profunda no intervalo
	! substrai o efeito da area saturada para o calculo da transpiracao
    sssu(ic,iu)=sssu(ic,iu)-asat(ic,iu)*ssmax(ib,iu) 
	if(sssu(ic,iu) < 0.) then
	    dd = dd + sssu(ic,iu)
	    sssu(ic,iu)=0.
	endif
endif

! transfere agua para a camada 2
if(d2(ib) > 0.) then
	saux=srad(ic,iu)+dd
    if(saux > srmax(ib,iu)) then
        ! excesso volta para a camada superior
        sssu(ic,iu) = sssu(ic,iu) + saux - srmax(ib,iu) 
        if(sssu(ic,iu) > ssmax(ib,iu)) then
            dsup=dsup+sssu(ic,iu)-ssmax(ib,iu)
		    sssu(ic,iu)=ssmax(ib,iu)
        endif
        saux=srmax(ib,iu)
    endif
    if(saux <= 0.) then
        srad(ic,iu)=0.
        dd=0.
    else
        aux = (1+nmenos*(1000*kss(ib,iu)*deltat/srmax(ib,iu))*(saux/srmax(ib,iu))**nmenos)**(1/nmenos)
        srad(ic,iu)=saux/aux
        dd = saux-srad(ic,iu) ! drenagem profunda no intervalo em mm
    endif
endif

! calcula escoamento superficial do excesso da area nao saturada e da area de contribuicao variavel
dsup=dsup+(pre(ic)-ec(ic,iu))*asat(ic,iu)
dsaux=dsaux+dd ! adiciona a drenagem profunda ao armazenamento subterraneo	
! calcula o escoamento subterraneo
if(dsaux > 0.) then
	if(mu(ib,iu).eq.1.) then
		aux=exp(-1000*tsub(ib,iu)*tanb(ic)*deltat/(dsmax*lambda(ic,iu)))
		ssub(ic,iu)=dsaux*aux
	else
		mumenos=mu(ib,iu)-1.
		aux=1+mumenos*1000*tsub(ib,iu)*tanb(ic)*deltat*(dsaux/dsmax)**mumenos/(lambdam(ic,iu)*dsmax)
        if(aux > 0) then
		    ssub(ic,iu)=dsaux/aux**(1/mumenos)
		else
		    ssub(ic,iu)=0.
		endif
	endif
	dsub = dsaux-ssub(ic,iu)
	ssub(ic,iu) = ssub(ic,iu) + sfc
else
	ssub(ic,iu) = dsaux + sfc
	dsub=0.
endif
return
end