!---------------------------------------------------------------------------------------------
! Rotina para a propagação na rede de drenagem
!---------------------------------------------------------------------------------------------
subroutine rede	   
use vars_main, only: qm1,qm2,qj1,qj2,icalib,celjus,numsubst,isubst,qlido,it,qcontorm,qcontorj,ncdom,idren,ibac,qsup,qrfp,dtp,nc,dt,&
& qcel1,qcel2,icell,qsub,qssu,ntrc,pmsub,pmssu,pmsup,vrsub,vrssu,vrsup,icodmusk,pjsub,pjssu,pjsup
use vars_calib, only: qjcal,nqcal,iqcal 
implicit none
real qjx,qaux,vt
integer ntr,ntc !numero de subtrechos, numero de subintervalos
integer ntclast,lt,jus
integer k,is,ib,ic,ic2,ifp !contadores

ntclast=dtp/dt(nc) ! numero de intervalo de tempo para resolver mc da ultima celula e usado para propagar de uma celula a outra
!condições de contorno (cc) de jusante no primeiro trecho eh a vazao de jusante do intervalo de tempo anterior
qcontorm=0.0
qcontorj(1:nc,1)=qj2(1:nc)

!em cada intervalo de tempo o que era t+1 vira t
qm1=qm2
qj1=qj2
qm2=0.0
qj2=0.0
qcel1=qcel2

if(icalib == 1 .and. nqcal > 0) then ! esta calibrando e tem cc de montante
	! condicao de contorno de montante 
    do ib=1,nqcal
        jus=celjus(iqcal(ib))
        is=0
        if(numsubst > 0) then ! verifica se a vazão calculada na célula deve ser substituьda por hidrograma lido
            do k=1,numsubst 
                if(iqcal(ib) == isubst(k)) then ! localiza a celula cuja vazao deve ser substituida
		            is=k
                    if(qlido(is,it) < 0.0) is=0 ! nao existe vazao de substituicao
                    exit
                endif
            enddo
        endif
        if(is > 0) then ! a vazao deve ser substituida
            qm2(jus)=qlido(is,it) ! insere a condicao de contorno como vazao de montante
			do lt=1,ntclast+1 !substitui valores calculados por valores lidos
				qcontorm(jus,lt)=qlido(is,it-1)+(lt-1)*(qlido(is,it)-qlido(is,it-1))/ntclast
			enddo
        else
            qm2(jus)=qjcal(ib,ntclast+1,it) ! insere a condicao de contorno como vazao de montante
			do lt=1,ntclast+1
                !o contorno de montante da celula jus recebe a saida do hidrograma recem calculado					
				qcontorm(jus,lt)=qjcal(ib,lt,it) 
			enddo
        endif
    enddo  
endif


do ic2=1,ncdom ! loop sobre as celulas do dominio
    ic=idren(ic2)
	ib=ibac(ic)
	ifp=icell(ic)

	if(ifp > 0 ) then ! celula de planice, soma/substrai a vazao de troca rio - planicie
	    qsup(ic)=qsup(ic)+qrfp(ifp) ! a troca rio-planicie eh tratada como fluxo superficial na celula
		if(qsup(ic) < 0. ) then ! foi retirada agua a mais no rio
	        write(*,*) 'qsup < 0 em rede'
	        read(*,*)
	    endif
	endif
	qcel2(ic)=max(qsub(ic)+qssu(ic)+qsup(ic),0.000001)  !soma vazoes geradas na celula no tempo t+1

	if(ntrc(ic) > 0)then !celula com rio

		!as vazoes que sao geradas nas celulas com rio entram a montante da propria celula e sao propagadas
		qaux = qm2(ic)+qcel2(ic) ! vazao total de montante
		
		!atualiza proporções das vazoes de montante
		pmsub(ic)=(pmsub(ic)*qm2(ic)+qsub(ic))/qaux !vazao subterranea
		pmssu(ic)=(pmssu(ic)*qm2(ic)+qssu(ic))/qaux !vazao sub-superficial
		pmsup(ic)=(pmsup(ic)*qm2(ic)+qsup(ic))/qaux !vazão superficial
		
		qm2(ic) = qaux ! atualiza vazao total de montante
					
		!atualiza proporções no volume do rio
		vrsub(ic)=vrsub(ic)+pmsub(ic)*qm2(ic)*dtp
		vrssu(ic)=vrssu(ic)+pmssu(ic)*qm2(ic)*dtp
		vrsup(ic)=vrsup(ic)+pmsup(ic)*qm2(ic)*dtp	
		
		!condicao de contorno Muskingum-Cunge de montante: as vazoes geradas na celula no tempo t+1 (qcel2) 
		!sao interpoladas linearmente entre t e t+1 usando vazao no tempo t (qcel1)
		
		do lt=1,ntclast+1
			qcontorm(ic,lt)=qcontorm(ic,lt)+qcel1(ic)+(lt-1)*(qcel2(ic)-qcel1(ic))/ntclast
		enddo

		ntr=ntrc(ic) !numero de subtrechos mc na celula
		ntc=dtp/dt(ic)	!numero de subintervalos de tempo mc na celula
		if(icodmusk(ic)==0)then !muskingun cunge linear
			call musk(qjx,ntr,ntc,ntclast,ic)
		elseif(icodmusk(ic)==1)then !muskingun cunge nao linear
			call musk_nl(qjx,ntr,ntc,ntclast,ic) 
		endif
		qj2(ic)=qjx

		!atualiza proporções no volume do rio usando vazao propagada na celula
		vt=vrsub(ic)+vrssu(ic)+vrsup(ic)
		pjsub(ic)=vrsub(ic)/vt
		pjssu(ic)=vrssu(ic)/vt
		pjsup(ic)=vrsup(ic)/vt
		vrsub(ic)=max(0.0,vrsub(ic)-pjsub(ic)*qj2(ic)*dtp)	
		vrssu(ic)=max(0.0,vrssu(ic)-pjssu(ic)*qj2(ic)*dtp)	
		vrsup(ic)=max(0.0,vrsup(ic)-pjsup(ic)*qj2(ic)*dtp)	

		!as vazoes propagadas na celula entram como condicao de contorno de montante da celula que esta a jusante
		jus=celjus(ic)
		qaux=max(qm2(jus)+qj2(ic),0.000001) ! vazao total de montante da cel a jusante

        if(numsubst > 0) then ! verifica se a vazão calculada na célula deve ser substituьda por hidrograma lido
            is=0
		    do k=1,numsubst 
			    if(ic==isubst(k)) then ! localiza a celula cuja vazao deve ser substituida
			        is=k
			        if(qlido(is,it) < 0.0) is=0
				    exit
                endif
            enddo
            
            if(is <= 0) then ! a celula nao tem substituto, ou ha falha nos dados substitutos                 
				! nao faz nada, usa o valor calculado
				do lt=1,ntclast+1
                    !o contorno de montante da celula jus recebe a saida do hidrograma recem calculado					
					qcontorm(jus,lt)=qcontorm(jus,lt)+qcontorj(ic,lt) 
				enddo
            else ! a vazao deve ser substituida
				qaux=qlido(is,it)
				do lt=1,ntclast+1 !substitui valores calculados por valores lidos
					qcontorm(jus,lt)=qlido(is,it-1)+(lt-1)*(qlido(is,it)-qlido(is,it-1))/ntclast
				enddo
			endif
		else ! nao ha substituicao
			do lt=1,ntclast+1
			    !o contorno de montante da celula jusante recebe a saida do hidrograma recem calculado
				qcontorm(jus,lt)=qcontorm(jus,lt)+qcontorj(ic,lt) 
			enddo
		endif

		! atualiza proporções da celula a jusante
		pmsub(jus)=(pmsub(jus)*qm2(jus)+pjsub(ic)*qj2(ic))/qaux  
		pmssu(jus)=(pmssu(jus)*qm2(jus)+pjssu(ic)*qj2(ic))/qaux
		pmsup(jus)=(pmsup(jus)*qm2(jus)+pjsup(ic)*qj2(ic))/qaux
		qm2(jus)=qaux ! atualiza a vazao de montante da celula a jusante
		
	else  !celula sem rio
		
		qj2(ic)=qcel2(ic) ! vazao total a jusante da celula em t+1 eh a vazao da celula
		! proporcoes a jusante
		pjsub(ic)=qsub(ic)/qj2(ic)
		pjssu(ic)=qssu(ic)/qj2(ic)
		pjsup(ic)=qsup(ic)/qj2(ic)
		
		! as vazoes que sao geradas nas celulas fonte entram a montante do trecho de rio que estр a jusante
		jus=celjus(ic)
		qaux=max(qm2(jus)+qj2(ic),0.000001)

		!condiãoo de contorno muskingum-cunge distribuida em cada intervalo de tempo da cel de jusante
		do lt=1,ntclast+1
			qcontorm(jus,lt)=qcontorm(jus,lt)+(qcel1(ic)+(lt-1)*(qcel2(ic)-qcel1(ic))/ntclast)
		enddo

		!atualiza as proporções
		pmsub(jus)=(pmsub(jus)*qm2(jus)+qsub(ic))/qaux !vazão subterranea
		pmssu(jus)=(pmssu(jus)*qm2(jus)+qssu(ic))/qaux !vazão sub-superficial
		pmsup(jus)=(pmsup(jus)*qm2(jus)+qsup(ic))/qaux !vazão superficial
		!atualiza a vazao de montante da celula de jusante
		qm2(jus)=qaux
	endif
enddo

if(icalib == -1 .and. nqcal > 0) then ! flag para indicar que na primeira rodada da calibracao, 
    ! armazena as cc de montante da bacia a ser calibrada usando a cc de jusante da cel de montante
    do ib=1,nqcal
        qjcal(ib,1:ntclast+1,it)=qcontorj(iqcal(ib),1:ntclast+1)
    enddo
endif

return
end
