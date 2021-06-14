!---------------------------------------------------------------------------------------------
! Subrotina para calculo da transpiracao
! Usa modelo de Jarvis - A Simple Empirical Model of Root Water Uptake, 1989
!---------------------------------------------------------------------------------------------
subroutine transpiracao(ic,iu,ib)
use vars_main,only: prad,e0,sssu,srad,ssub,rc,d1,fr,ett,etp,ec,asat,e0ag,delta,gamma,ra,ssmax,secrbs,pi,&
&imes,iaf,secrjah,secrjal,srmax,d2,smax,omegac,cover
implicit none
integer ic,iu,ib
real rcx,etssu,etrad,etsub,fsss,frad,fsub,omega,a,fr1,fr2,fr3,scrit,secr,aux
real etbs

!variables locais

! etssu : evaporacao da camada superior do solo
! etsub	: evaporacao da camada inferior do solo
! rcx	: resistencia superficial
! a		: energia disponível para transpiracao 
! omega	: indice de estresse ponderado
! 
! calcula a quantidade de energia disponivel para evaporacao descontando a interceptacao
if(ec(ic,iu) < 0.) then
    ! se for orvalho, já foi considerado no calculo de escoamento
    ec(ic,iu)=0.
endif
ett(ic,iu) = ec(ic,iu) ! evaporacao total
if(asat(ic,iu) >= 1.) then ! toda a celula esta saturada
	etp(ic,iu)=0.0 ! nao ha transpiracao
    a = e0ag-ec(ic,iu) ! energia disponivel
	if(a > 0.) then
	    ett(ic,iu) = ett(ic,iu) + a 	! atualiza a evaporacao total
        ! atualiza armazenamento da terceira camada de solo, a evaporacao da area saturada sai dessa camada
        ssub(ic,iu) = ssub(ic,iu) - a
        if(ssub(ic,iu) < 0.) then
            ett(ic,iu)  = ett(ic,iu) + ssub(ic,iu) 
	        ssub(ic,iu) = 0.
        endif
    endif
else ! celula com fracao da area nao coberta pela agua
	!evaporacao da area saturada do bloco eh calculada como evaporacao da camada 3 pois ela controla a area saturada
	!terceira camada: evaporacao da area saturada + transpiracao
	etsub = asat(ic,iu)*max(e0ag-ec(ic,iu),0.)
	
	! transpiracao e extracao radicular
    a = e0-ec(ic,iu) ! energia disponivel
    if(a > 0.) then ! existe energia disponivel para transpiracao
	    ! evaporacao do solo nu
 	    etbs=a*(delta+gamma)/(delta+gamma*(1.+52./ra)) ! evaporacao potencial solo desnudo
 	    if(d1(ib) < 0.1) then ! camada menor que 10 cm
 	        aux=1.
 	    else
 	        aux=1.-0.3/d1(ib) ! acochambrado, visando reduzir a umidade a superficie quando a camada 1 eh espessa
 	    endif
 	    ! considero que o armazenamento da crosta de solo é proporcional a espessura
 	    scrit=aux*sssu(ic,iu)/(ssmax(ib,iu)*secrbs(iu)) 
	    if(scrit < 1.) then
            etbs=etbs*0.25*(1-cos(pi*scrit))**2. ! lee and pielke, 1991
	    endif
		
	    rcx=rc(iu)/iaf(iu,imes) !resistencia superficial da vegetacao
	    etp(ic,iu)=a*(delta+gamma)/(delta+gamma*(1.+rcx/ra)) ! transpiracao potencial
       
       ! define fatores de estresse do modelo de jarvis
        secr=min(1.,secrjah(iu)*etp(ic,iu)+secrjal(iu))
        secr=max(secr,0.)
   	    !fator de estresse devido a umidade de solo na camada superior do solo
        scrit=secr*ssmax(ib,iu)
        fsss=min(sssu(ic,iu)/scrit,1.)
        
	    if(prad(iu,imes) < d1(ib)) then
	        ! a profundidade radicular e menor que a espessura da camada superior: nao ha evaporacao na camadas inferiores
	        fr1 = 1.0
	        fr2 = 0.0
	        frad = 0.0
	        fr3 = 0.0
	        fsub = 0.0
        else
            !fator de estresse devido a umidade de solo da segunda camada
            scrit=max(0.001,secr*srmax(ib,iu)) ! armazenamento critico da segunda camada
    	    frad=min(srad(ic,iu)/scrit,1.)
	        fr1 = 1.0 - exp(-fr * d1(ib) / prad(iu,imes)) 
            if (prad(iu,imes) < d1(ib)+d2(ib)) then 
	            ! se a prof radicular for menor que d1+d2, o armazenamento disponivel para
	            ! evaporacao profunda eh dado pela segunda camada
	            fr2 = 1.0 - fr1
	            fr3 = 0. 
	            fsub= 0.
	        else
                !fator de estresse devido a umidade de solo na camada inferior do solo
                scrit=secr*smax(ib,iu) ! armazenamento critico da camada inferior
	            ! se prad >d1+d2, srad > 0 e a agua disponivel para transpiracao 
	            ! profunda sai tambem da terceira camada
    	        fsub=min(ssub(ic,iu)/scrit,1.)
    	        
		        fr2 = 1.0 - fr1 - exp(-fr * (d1(ib)+d2(ib)) / prad(iu,imes))
    	        fr3 = 1.0 - fr1 - fr2
	        endif
        endif
	    ! fator de compensacao ponderado
	    omega=fr1*fsss+fr2*frad+fr3*fsub
	    omega=max(omega,omegac)
		
	    !camada superior de solo: evaporacao solo desnudo + transpiracao na area nao saturada
	    etssu = (1.- asat(ic,iu))*(cover(iu,imes)*etp(ic,iu)*fr1*fsss/omega + (1.-cover(iu,imes))*etbs)
        
        !segunda camada: transpiracao na area nao saturada
		etrad = (1.- asat(ic,iu))*cover(iu,imes)*fr2*etp(ic,iu)*frad/omega

	    !evaporacao da area saturada do bloco eh calculada como evaporacao da camada 3 pois ela controla a area saturada
	    !terceira camada: evaporacao da area saturada + transpiracao
	    etsub = etsub + (1.- asat(ic,iu))*cover(iu,imes)*fr3*etp(ic,iu)*fsub/omega
    else ! na existe energia disponivel para transpiracao
        etssu=0.
        etrad=0.        
    endif
    		
    ! atualiza os armazenamentos
    ! camada superior do solo
    sssu(ic,iu) = sssu(ic,iu) - etssu
    if(sssu(ic,iu) < 0.) then
	    etssu = etssu + sssu(ic,iu)
	    sssu(ic,iu) = 0.
    endif
    
    ! segunda camada de solo
    srad(ic,iu) = srad(ic,iu) - etrad
    if(srad(ic,iu) < 0.) then
	    etrad = etrad + srad(ic,iu)
	    srad(ic,iu) = 0.
    endif

    ! terceira camada de solo
    ssub(ic,iu) = ssub(ic,iu) - etsub
    if(ssub(ic,iu) < 0.) then
        etsub  = etsub + ssub(ic,iu) 
	    ssub(ic,iu) = 0.
    endif
    
    ! evaporacao total
    ett(ic,iu) = ett(ic,iu) + etssu + etrad + etsub
		
endif

return
end