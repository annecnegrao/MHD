!---------------------------------------------------------------------------------------------
!subrotina de previsao
!---------------------------------------------------------------------------------------------
subroutine previsao
use vars_main ,only: dir_dados,filccini,filprevfin,nc,ncdom,icdom,adren,idren,cdata,sssu,nu,nc, &
                    & ssub,srad,vsup,vssu,vsub,ta,qm2,qj2,ntrc,it,cdatag,ihidg,nhidg,vrsup, &
                    & nfp,afp,qrfp,tontem,iano,idia,imes,qrioini,qcel2,pmssu,pmsup, &
                    & vfp,nt,ntprev,ihoraini,dtp,ihora,jdia,idini,pmsub,vrsub,vrssu,ntrmax
implicit none
integer iu,ic,k,irec !contadores
real tempo ! tempo de simulacao em dias
real qprev(nhidg,ntprev) !armazena hidrogramas onde se deseja gravar nas celulas designadas
character cdatal*16
character*100 fmt !fromatacao


! dominio da simulacao
ncdom=nc
do ic=1,nc
    icdom(ic)=ic ! por enquanto, simula todo o dominio na rotina previsao
enddo

!vetor com a sequencia de celulas com areas drenada crescente usados na resolucao na rotina rede
call sort(nc,adren,idren)

! executa o modelo no modo simulacao ate o dia do inicio da previsao
call modelo

!guarda valores das variaveis antes de comecar o ciclo de previsao
open(filccini,file=dir_dados//'entrada/'//'ccini_'//cdata//'.bin', &
    & status='unknown',form='unformatted',access='direct',recl=4*nc*(ntrmax+1))
irec = 1
write(filccini,rec=irec) ((sssu(ic,iu),iu=1,nu),ic=1,nc)
irec = irec + 1
write(filccini,rec=irec) ((srad(ic,iu),iu=1,nu),ic=1,nc)
irec = irec + 1
write(filccini,rec=irec) ((ssub(ic,iu),iu=1,nu),ic=1,nc)
irec = irec + 1
!dos reservatorios subterraneo, subsuperficial e superficial
write(filccini,rec=irec) vsup ; irec = irec + 1
write(filccini,rec=irec) vssu ; irec = irec + 1
write(filccini,rec=irec) vsub ; irec = irec + 1
write(filccini,rec=irec) ta ; irec = irec + 1
write(filccini,rec=irec) qm2 ; irec = irec + 1
write(filccini,rec=irec) qj2 ; irec = irec + 1
write(filccini,rec=irec) ntrc ; irec = irec + 1
write(filccini,rec=irec) qrioini ; irec = irec + 1
write(filccini,rec=irec) qcel2 ; irec = irec + 1
write(filccini,rec=irec) pmsub ; irec = irec + 1
write(filccini,rec=irec) pmssu ; irec = irec + 1
write(filccini,rec=irec) pmsup ; irec = irec + 1
write(filccini,rec=irec) vrsub ; irec = irec + 1
write(filccini,rec=irec) vrssu ; irec = irec + 1
write(filccini,rec=irec) vrsup ; irec = irec + 1
if(nfp>0) then
    write(filccini,rec=irec) vfp ; irec = irec + 1
    write(filccini,rec=irec) afp ; irec = irec + 1
    write(filccini,rec=irec) qrfp ; irec = irec + 1
endif
close(filccini)

!chama a subrotina de assimilacao de dados para atualizar o modelo do ultimo passo de tempo
call atualiza

!ciclo de previsao 
!grava informacoes no arquivo diario
open(filprevfin,file=dir_dados//'saida/'//'prev_'//cdata//'.hig',status='unknown') ! arquivo com a previsao
write(filprevfin,'(a22,a10)') ' previsao iniciada em ',cdata	
write(filprevfin,*)

! grava condicao inicial
cdatal=cdatag(it)(1:2)//'/'//cdatag(it)(3:4)//'/'//cdatag(it)(5:8)//' '//cdatag(it)(9:10)//':00'
write(fmt,'(a,i2.2,a)')"(a16,",nhidg,"f10.3)"
write(filprevfin,fmt) cdatal,(qj2(ihidg(k)),k=1,nhidg)

do while (it < nt+ntprev)
    it=it+1
    tempo=ihoraini/24.+float(it-1)*dtp/86400
	ihora=nint((tempo-int(tempo))*24.) ! hora da simulacao
	jdia=idini+int(tempo) ! dia juliano do calendario 
			
	call caldat(jdia,imes,idia,iano)
	write(cdata(1:2),'(i2.2)') idia
    write(cdata(3:4),'(i2.2)') imes
    write(cdata(5:8),'(i4)') iano    
    write(cdata(9:10),'(i2.2)') ihora
	cdatag(it)=cdata
		
	!subrotinas de leitura e preparacao de dados meteorologicos e de chuva
	call lemetprev    ! le dados de chuva e meteorologicos para calcular evaporacao

	tontem=ta
	
	call celula

	call rede
	
	if(nfp>0) then
	    ! subrotina de alagamento
	    call floodplain
	endif
	
	!guarda dados para gravar hidrogramas em locais definidos no arquivo parfix.hig
	do k=1,nhidg 
		! ihidg(k) - numero da celula em que se deseja o hidrograma
		qprev(k,it-nt)=qj2(ihidg(k)) !qprev hidrogramas previtos nos locais desejados
	enddo
enddo !fim do loop de tempo

do it=nt+1,nt+ntprev
	!grava arquivo de hidrograma de previsao
	cdatal=cdatag(it)(1:2)//'/'//cdatag(it)(3:4)//'/'//cdatag(it)(5:8)//' '//cdatag(it)(9:10)//':00'
	write(filprevfin,fmt) cdatal,(qprev(k,it-nt),k=1,nhidg) 
enddo
close(filprevfin)

return
end
