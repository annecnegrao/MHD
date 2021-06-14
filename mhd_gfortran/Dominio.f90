!---------------------------------------------------------------------------------------------
!subrotina que eh chamada na opcao calibra e seleciona como dominio a bacia a ser calibrada, 
!para diminuir o tempo de calibracao
!---------------------------------------------------------------------------------------------
subroutine dominio
use vars_main, only: celjus,iexut,idren,icdom,nc,ncdom,adren,icell,icellfp,zfpt,afpt,zvert,zfpb,vfpt,cellw,celln,nfp,&
& nb,iqobs,nobs,ibac,nt,nbvfo,icalib
use vars_calib, only: ibcal,qjcal,nqcal,iqcal
implicit none
integer ic,ic2,ifp,nfp2,k,ib
real aux(nc)

!verificar o que acontece quando o posto encontra-se a jusante do exutorio

nqcal=0 ! indica o numero de postos a montante
allocate(iqcal(nb-1)) !vetor que indica as celulas dos exutorios a montante da bacia a ser calibrada
do ib=1,nb-1 ! javier: verificar o que acontece quando calibra a ultima subbacia
	if(ibac(celjus(iexut(ib))) == ibcal) then ! bacia a jusante do exutorio da sub-bacia ib
        nqcal=nqcal+1
        iqcal(nqcal)=iexut(ib) 
    endif
enddo

if(nqcal > 0) then ! existem cc a montante da bacia a ser calibrada
    allocate(qjcal(nqcal,25,nt))
    ncdom=nc
    do ic=1,nc
        icdom(ic)=ic ! simula todo o dominio na rotina simulacao 
    enddo
    !vetor com a sequencia de celulas com areas drenada crescente usados na resolucao na rotina rede
    call sort(nc,adren,idren)
    icalib = -1 ! flag para indicar que precisa gravar cc de montante 
    ! faz uma primeira rodada em todo o dominio de simulacao para criar as condicoes de montante 
    ! e guarda-las em qcontcal(ib,ntrc,it)
    call modelo
    icalib = 1
endif

ncdom=0 ! vetor que identifica as celulas da bacia a ser calibrada
do ic=1,nc
    if(ibac(ic) == ibcal) then
        ncdom=ncdom+1
        icdom(ncdom) = ic ! indica as celulas dentro do dominio
    endif
enddo

! cria vetor auxiliar de area drenada
aux=1.0e15 ! numero muito grande, para que as celulas fora do dominio fiquem no final

do ic=1,ncdom
    aux(icdom(ic))=adren(icdom(ic)) ! assina os valores verdadeiros de area de drenagem apenas no dominio
enddo

!vetor com a sequencia de celulas do dominio com areas drenada crescente usados na resolucao na rotina rede
call sort(nc,aux,idren)

! elimina as celulas de planicie fora do dominio
if(nfp > 0) then ! existem celulas de planicie
    do ifp = 1,nfp ! loop entre as celulas de planicie
        ! icellfp(ifp) = x indica que a celula ifp do arquivo floodplain.hig corresponde
        ! `a x-esima celula do arquivo cell.hig
        ic = icellfp(ifp) 
        do ic2= 1,ncdom ! loop entre as celulas do dominio
            if(icdom(ic2) == ic) then ! a celula de planicie esta dentro do dominio 
                exit
            endif
        enddo
        if(ic2 > ncdom) then ! a celula de planice esta fora do dominio e eh zerada
            ! se icell(ic) = 0 a celula ic de cell.hig nao pertence a planicie
            icell(icellfp(ifp))=0 ! icell indica a posicao das celulas da planicie do arquivo cell.hig no arquivo floodplain.hig
            icellfp(ifp)=0 
        endif
    enddo
	
    ! arruma o arquivo das celulas de planicie em sequencia com as celulas que sobraram
    nfp2 = 0
    do ifp=1,nfp
        ic = icellfp(ifp)
        if(ic == 0) cycle
        nfp2=nfp2+1
        icellfp(nfp2) = icellfp(ifp) ! modifica o arquivo de planicie deslocando as celulas
	    zvert(nfp2) = zvert(ifp)
        zfpb(nfp2) = zfpb(ifp)
        do k=1,50
            zfpt(nfp2,k) = zfpt(ifp,k)
            afpt(nfp2,k) = afpt(ifp,k)
            vfpt(nfp2,k) = vfpt(ifp,k)
        enddo
        celln(nfp2)=celln(ifp)
        cellw(nfp2)=cellw(ifp)
        icell(icellfp(nfp2))=nfp2
    enddo
    nfp = nfp2 ! arquivo de planicie agora eh menor
endif

return
end