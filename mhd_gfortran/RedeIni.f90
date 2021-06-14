!---------------------------------------------------------------------------------------------
!esta rotina inicializa a rede de drenagem no primeiro passo de tempo
!quando nao existem condicoes iniciais
!---------------------------------------------------------------------------------------------
subroutine redeini
use vars_main,only: qj2,qm2,pmsub,pmssu,pmsup,qcel1,qcel2,qrioini,vrsub,vrssu,vrsup,qcontorm, &
                  & dtp,qcontorj,ncdom,idren,ntrc,qsup,qssu,qsub,celjus,pjsub,pjssu,pjsup, &
                  & srio,cel,dt,nc
implicit none
integer ntclast !numero de subintervalos
!proporcoes de origem
real qaux
integer lt,jus
integer ic,ic2 !contadores

! numero de intervalo de tempo para resolver mc em cada celula, fixado a partir 
! da ultima celula que tem o maior intervalo de tempo
ntclast=dtp/dt(nc)

!condicoes iniciais da vazao no rio
qj2=0.0
qm2=0.0
pmsub=0.0
pmssu=0.0
pmsup=0.0
qcel1=0.0
qcel2=0.0
qrioini=0.0
vrsub=0.0
vrssu=0.0
vrsup=0.0
qcontorm=0.
qcontorj=0.

do ic2=1,ncdom ! loop sobre as celulas do dominio
    ic=idren(ic2)
	if(ntrc(ic) > 0)then	!celula com rio
		!transforma proporcoes em vazoes		
		pmsub(ic)=pmsub(ic)*qm2(ic) !vazao subterranea
		pmssu(ic)=pmssu(ic)*qm2(ic) !vazao sub-superficial
		pmsup(ic)=pmsup(ic)*qm2(ic) !vazao superficial
			
		!acumula vazao de montante mais vazao na celula
		qm2(ic)=max(qm2(ic)+qsup(ic)+qssu(ic)+qsub(ic),0.000001)
			
		!atualiza proporcoes
		pmsub(ic)=(pmsub(ic)+qsub(ic))/qm2(ic)
		pmssu(ic)=(pmssu(ic)+qssu(ic))/qm2(ic)
		pmsup(ic)=(pmsup(ic)+qsup(ic))/qm2(ic)

		!condicao de contorno Muskingum-Cunge de montante
		do lt=1,ntclast+1
			qcontorm(ic,lt)=qcontorm(ic,lt)+qm2(ic)
		enddo

       ! propaga para jusante
		jus=celjus(ic)
		qj2(ic)=qm2(ic)	!vazao constante no trecho com condicao inicial

		!condicao de contorno Muskingum-Cunge de jusante
		do lt=1,ntclast+1
			qcontorj(ic,lt)=qcontorj(ic,lt)+qj2(ic)
		enddo
			
		!a proporcao que entra sai
		pjsub(ic)=pmsub(ic)
		pjssu(ic)=pmssu(ic)
		pjsup(ic)=pmsup(ic)
			
		! transforma proporcao em vazao
		pmsub(jus)=pmsub(jus)*qm2(jus) !vazao subterranea
		pmssu(jus)=pmssu(jus)*qm2(jus) !vazao sub-superficial
		pmsup(jus)=pmsup(jus)*qm2(jus) !vazao superficial

		!o que sai do rio de uma cel. vai pra outra
		qm2(jus)=max(qm2(jus)+qj2(ic),0.00001) 

		!atualiza proporcoes
		pmsub(jus)=(pmsub(jus)+pjsub(ic)*qj2(ic))/qm2(jus)  
		pmssu(jus)=(pmssu(jus)+pjssu(ic)*qj2(ic))/qm2(jus)
		pmsup(jus)=(pmsup(jus)+pjsup(ic)*qj2(ic))/qm2(jus)

		!tentativa de proporcoes
		qaux=(qm2(ic)+qj2(ic))/2
		vrsub(ic)=0.3333*(1000.*srio(ic)*qaux/cel(ic))
		vrssu(ic)=0.3333*(1000.*srio(ic)*qaux/cel(ic))
		vrsup(ic)=0.3333*(1000.*srio(ic)*qaux/cel(ic))

	else !celula fonte
		jus=celjus(ic)
		pmsub(jus)=pmsub(jus)*qm2(jus) !vazao subterranea
		pmssu(jus)=pmssu(jus)*qm2(jus) !vazao sub-superficial
		pmsup(jus)=pmsup(jus)*qm2(jus) !vazao superficial

		qm2(jus)=max(qm2(jus)+qsup(ic)+qssu(ic)+qsub(ic),0.000001)

		!condicao de contorno Muskingum-Cunge
		do lt=1,ntclast+1
			qcontorm(jus,lt)=qcontorm(jus,lt)+qm2(jus)
		enddo

		!atualiza proporcoes
		pmsub(jus)=(pmsub(jus)+qsub(ic))/qm2(jus)  
		pmssu(jus)=(pmssu(jus)+qssu(ic))/qm2(jus)
		pmsup(jus)=(pmsup(jus)+qsup(ic))/qm2(jus)

	endif
enddo
!fim das condicoes iniciais
return
end