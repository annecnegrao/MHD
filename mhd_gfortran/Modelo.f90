!---------------------------------------------------------------------------------------------
!esta subrotina comanda o loop do tempo do modelo hidrologico e chama as 
!rotinas de balanco e propagacao nas celulas e de propagacao na rede de drenagem
!---------------------------------------------------------------------------------------------
subroutine modelo
use vars_main,only: nmapas,ett,ec,etp,asat,itmap,ihoraini,dtp,jdia,idini,ihora,imes,idia, &
                  & iano,cdata,dir_dados,ncdom,icdom,ibac,sssu,ssmax,ssub,smax,srad,srmax, &
                  & vsub,qesp,acel,cb,nu,tkb,qsub,vssu,vsup,qsup,qssu,nfp,afp,vfp,zfp,zfpb, &
                  & qrfp,ta,filccini,nc,ntrmax,qm2,qj2,qcel2,pmsub,pmssu,pmsup,vrsub,vrssu, &
                  & vrsup,iniaju,qrioini,ntrc,it,nt,cdatag,cdata,tontem,pre,preall,taall, &
                  & tdew,tdewall,vv,vvall,patm,patmall,roc,rocall,icalib,iqobs,nb,iexut,qr, &
                  & nobs,nbvfo,nhidg,qrg,ihidg,filpmedia,filsol,prec,qbsub,pjsub,qbssu, &
                  & pjssu,qbsup,pjsup
implicit none 
integer,allocatable:: naux(:) ! numero de subtrechos no arquivo de condicoes iniciais 
real,allocatable:: qaux(:,:) ! vazoes qrioini no arquivo de condicoes iniciais 
integer imapas,ic,ic2,ib,iu,itr1,itr2 !contadores
integer jc,ju
real tempo ! tempo de simulacao em dias
real somad,somaar,pmeddia
real vbx,aux,dx
real qrb(nb,nt)	!armazena hidrogramas nos exutorios das sub-bacias
logical existe

! condicao inicial uso da terra
! itmap indica o numero do intervalo de tempo onde muda o mapa de uso da terra
! se a data inicial da previsao for posterior a ultima data do uso da terra, itmap <= 0
if(itmap(nmapas) <= 0) then
    ! forca a leitura do ultimo mapa de uso da terra
    call leuso(nmapas)
else
    ! procura o primeiro mapa de solo apos inicio da simulacao 
	do imapas=1,nmapas
		if(itmap(imapas) >= 0)then
			call leuso(imapas)
            exit
		endif
	enddo
endif

!zera evaporacao e area saturada
ett=0.0
ec=0.0
etp=0.0
asat=0.0

! condicao inicial de armazenamento do sistema: corresponde a data inicial -1
! cria data para a condicao inicial
tempo=ihoraini/24.-dtp/86400
if(tempo < 0) then
    tempo=tempo+1.
	jdia=idini+int(tempo)-1 ! dia juliano do calendario
else
	jdia=idini+int(tempo) ! dia juliano do calendario
endif

ihora=nint((tempo-int(tempo))*24.) ! hora da simulacao
call caldat(jdia,imes,idia,iano)	
write(cdata(1:2),'(i2.2)') idia
write(cdata(3:4),'(i2.2)') imes
write(cdata(5:8),'(i4)') iano    
write(cdata(9:10),'(i2.2)') ihora

!**********condicoes iniciais*********************************************	        
! verifica se existe arquivo com condicao inicial	
inquire(file=dir_dados//'entrada/ccini_'//cdata//'.bin',exist=existe)
if(.not. existe) then ! condicoes iniciais estimadas ja que nao existe arquivo inicial
	do ic2=1,ncdom
	    ic=icdom(ic2)
	    ib=ibac(ic)
	    ! do perfil de solo
	    !considera que solo esta na capacidade de campo no inicio da simulacao nas duas camadas
	    sssu(ic,1:nu)=0.4*ssmax(ib,1:nu)
	    ssub(ic,1:nu)=0.5*smax(ib,1:nu)
	    ! para definir a capacidade inicial da zona radicular usa o a espessura da terceira camada
	    srad(ic,1:nu)=0.3*srmax(ib,1:nu)
        ! dos reservatorios subterraneo, subsuperficial e superficial
	    vsub(ic)=qesp(ib)*acel(ic)*cb(ib) ! cb em segundos
	    vbx=vsub(ic)*exp(-dtp/tkb(ic))
	    qsub(ic)=(vsub(ic)-vbx)/dtp
	    vssu(ic)=0.0
	    vsup(ic)=0.0
		qsup(ic)=0.0
		qssu(ic)=0.0
	enddo
    ! da planicie de alagamento
	if(nfp > 0) then ! planicie de alagamento seca
	    afp=0.0
	    vfp=0.0
	    zfp=zfpb
	    qrfp=0.0
	endif
	!da rede de drenagem
	call redeini
	!de temperatura
	ta=28.0	
else ! condicoes iniciais dadas pelos valores armazenados do solo e na rede
    write(*,*) ' existe condicao inicial'
	open(filccini,file=dir_dados//'entrada/'//'ccini_'//cdata//'.bin', status='old',form='unformatted',access='direct')
	allocate(naux(nc),qaux(nc,ntrmax+1))
	read(filccini) ((sssu(ic,iu),iu=1,nu),ic=1,nc)
	read(filccini) ((srad(ic,iu),iu=1,nu),ic=1,nc)
	read(filccini) ((ssub(ic,iu),iu=1,nu),ic=1,nc)
	read(filccini) vsup
	read(filccini) vssu
	read(filccini) vsub
	read(filccini) ta
	read(filccini) qm2
	read(filccini) qj2
	read(filccini) naux
	read(filccini) qaux
	read(filccini) qcel2
	read(filccini) pmsub
	read(filccini) pmssu
	read(filccini) pmsup
    read(filccini) vrsub
    read(filccini) vrssu
    read(filccini) vrsup
    if(nfp>0) then
	    read(filccini) vfp
        read(filccini) afp
        read(filccini) qrfp
        zfp=zfpb
    endif
	close(filccini)
	iniaju=1 ! sem spin-up para inicio de comparacao
	! verifica se houve mudancas no  numero de subtrechos e interpola se necesario
	qrioini=qaux
    do ic=1,nc
        if(ntrc(ic) /= naux(ic)) then ! mudou o numero de subtrechos
            aux=(naux(ic)+1)/(ntrc(ic)+1)
            do itr1=2,ntrc(ic)+1
                dx=(itr1-1)*aux+1
                do itr2=2,naux(it)+1
                    if(dx <= itr2) then
                        qrioini(ic,itr1)=qaux(ic,itr2)+(dx-itr2)*(qaux(ic,itr2)-qaux(ic,itr2-1))
                        exit
                    endif
                enddo
            enddo
        endif
    enddo
    deallocate(naux,qaux)
endif
!****************************fim das condicoes iniciais*********************************
    
!	inicio do loop do tempo. 
it=0     	
do while (it < nt)
	it=it+1
    ! arruma as datas
    tempo=ihoraini/24.+float(it-1)*dtp/86400
	ihora=nint((tempo-int(tempo))*24.) ! hora da simulacao
	jdia=idini+int(tempo) ! dia juliano do calendario 
    ! arruma datas em format caracter para a simulacao
    call caldat(jdia,imes,idia,iano)
    write(cdata(1:2),'(i2.2)') idia
    write(cdata(3:4),'(i2.2)') imes
    write(cdata(5:8),'(i4)') iano
    write(cdata(9:10),'(i2.2)') ihora
    cdatag(it)=cdata

	do imapas=2,nmapas
		if(it == itmap(imapas))then
			call leuso(imapas)
			exit
		endif
	enddo
	! variaveis meteorologicas no intervalo
	tontem=ta
	
	do ic2=1,ncdom
	    ic=icdom(ic2)
	    pre(ic)=preall(ic,it) ! precipitacao no intervalo
        ta(ic)=taall(ic,it) ! temperatura do ar no intervalo
        tdew(ic)=tdewall(ic,it) ! temperatura ponto de orvalho no intervalo
        vv(ic)=vvall(ic,it) ! velocidade do vento no intervalo
        patm(ic)=patmall(ic,it) ! pressao atmosferica no intervalo
        roc(ic)=rocall(ic,it) ! radiacao global incidente no intervalo
    enddo
    
	!subrotina da celula
	call celula
	!subrotina da rede de drenagem
	call rede
	
	if(nfp>0) then
	    ! subrotina de alagamento
	    call floodplain ! verificar como funciona no dominio
	endif

	!armazena dados de vazao das celulas em que existe vazao observada 
	!qr(ib,it) vai ser comparado a qobs(ib,it) na rotina fobj
    !iqobs(ib) celulas que correspondem aos postos flu com dados
    if(icalib == 1) then ! esta calibrando
    	qr(nbvfo,it)=qj2(iqobs(nbvfo)) ! so armazena o posto que esta sendo calibrado
   	elseif(icalib == 0)then !esta simulando
   	    do ib=1,nb
		    !armazena vazoes das sub-bacias
		    qrb(ib,it)=qj2(iexut(ib)) ! iexut(ib) celula do exutorio da sub-bacia ib
        enddo   	
	    do ib=1,nobs
		    qr(ib,it)=qj2(iqobs(ib)) ! postos com vazoes observadas
	    enddo
		do ib=1,nhidg !guarda dados para gravar hidrogramas em locais definidos no arquivo parhig 
			!ihidg(k) eh o numero da celula em que se deseja o hidrograma
			qrg(ib,it)=qj2(ihidg(ib)) !qrg armazena os hidrogramas nos locais desejados
		enddo	    
	    	    
        somad=sum(pre(1:nc)*acel(1:nc))
        somaar=sum(acel(1:nc))
        pmeddia=somad/somaar
        write(filpmedia,'(a10,f10.2)') cdatag(it),pmeddia
	
	    !************** nosolo.hig ***************************************************************
	    !As linhas abaixo servem para gravar algumas variaveis detalhadas de um bloco e de
	    !uma celula. Os valores de jb e jc podem ser alterados, ju indica o bloco e jc a
	    !celula em que se desejam os dados. Os dados sao gravados num arquivo chamado nosolo.hig.
		ju=1
		jc=10
		write(filsol,'(i6,12f10.3)')it,pre(jc),ett(jc,ju),ec(jc,ju),etp(jc,ju),asat(jc,ju), &
                      & sssu(jc,ju),srad(jc,ju),ssub(jc,ju),qsub(jc),qssu(jc),qsup(jc)
	    !fim da saida de dados para o arquivo nosolo.hig
	    !******************************************************************************************
	    
	    do ic=1,nc
		    prec(ic)=prec(ic)+pre(ic) !acumula chuva na celula
	    enddo

		!armazena vazoes segundo a origem para uma celula
		! ihidg(nbvfo) celula em que se desejam os resultados de hidrograma separado por origem
		qbsub(it)=qj2(ihidg(nbvfo))*pjsub(ihidg(nbvfo))
		qbssu(it)=qbssu(it)+qj2(ihidg(nbvfo))*pjssu(ihidg(nbvfo))
		qbsup(it)=qbsup(it)+qj2(ihidg(nbvfo))*pjsup(ihidg(nbvfo))
	endif

enddo !fim do loop do tempo

return
end