!---------------------------------------------------------------------------------------------
!Esta sub-rotina e utilizada pela sub-rotina previsao para fazer a atualizacao
!de algumas variaveis no instante t0, em que inicia o ciclo de previsao. 
!A base para a atualizacao e a comparacao com algum dado observado de vazao.
!(nao serve para atualizar dados de planicie)
!---------------------------------------------------------------------------------------------
subroutine atualiza
use vars_main,only: dir_dados,nobs,qobs,filprev,nc,nu,ntrmax,qcel2,qm2,pmsup,pmssu,pmsub,qj2, &
                  & vrsup,vrssu,vrsub,qrioini,sssu,srad,ssub,qsup,qssu,qsub,it,qr,iqobs,ibac, &
                  & dtp,vsup,tks,vssu,vsub,tkb,adren,acel
implicit none
real,allocatable:: qcelat(:),qmat(:),pmsupat(:),pmssuat(:),pmsubat(:),qjat(:),qrioiniat(:,:) ! vazoes na rede inicial 
real,allocatable:: vrsupat(:), vrssuat(:), vrsubat(:) ! volumes armazenados no rio no incio do intervalo de tempo
real,allocatable:: sssuat(:,:),sradat(:,:),ssubat(:,:) ! armazenamentos atualizados
real,allocatable:: qsupat(:),qssuat(:),qsubat(:) ! fluxos atualizados
real eps(nobs),f1(nobs),f2(nobs),f3(nobs),dif1(nobs),dif2(nobs) ! tolerancia do erro na vazao, fator de correcao e erros no inicio e final da iteracao
real vt,vaux
integer postorio(nobs),ndias(nobs) !indica o numero de posto que sera usado para corrigir a vazao no rio e quantas observacoes serao usadas
integer ic,ib,iu,iter,nbac ! contadores

! le informacoes basicas do arquivo parprev.hig
open(filprev,file=dir_dados//'entrada/ParPrev.hig',status='old')
read(filprev,*)
read(filprev,*)
read(filprev,*)
read(filprev,*) (postorio(ib),ib=1,nobs) ! indica o numero do posto com vazao observada que serao usados para corrigir vazao em cada sub-bacia (0 nao corrige)
read(filprev,*)
read(filprev,*)
read(filprev,*) (ndias(ib),ib=1,nobs) ! numero de intervalos de tempo para calcular a media, o erro e o fator de correcao
read(filprev,*)
read(filprev,*)
read(filprev,*) (eps(ib),ib=1,nobs) !tolerancia admitida no erro da estimativa da vazao
close(filprev)

!   guarda os valores originais antes de iniciar a iteracao
allocate (qcelat(nc),qmat(nc+1),pmsupat(nc+1),pmssuat(nc+1),pmsubat(nc+1))
allocate (qjat(nc),vrsupat(nc), vrssuat(nc), vrsubat(nc),qrioiniat(nc,ntrmax+1))
allocate (sssuat(nc,nu),sradat(nc,nu),ssubat(nc,nu))
allocate (qsupat(nc),qssuat(nc),qsubat(nc))

qcelat=qcel2 ! ultima vazao na celula
qmat=qm2 ! ultima vazao de montante calculada 
pmsupat=pmsup ! proporcao de fluxo superficial
pmssuat=pmssu ! proporcao de fluxo subseperficial
pmsubat=pmsub ! proporcao de fluxo subterraneo    
qjat=qj2 ! ultima vazao de jusante calculada
vrsupat=vrsup ! volume no rio devido ao fluxo superficial no inicio do intervalo
vrssuat=vrssu ! volume no rio devido ao fluxo subsuperficial no inicio do intervalo
vrsubat=vrsub ! volume no rio devido ao fluxo subterraneo no inicio do intervalo
qrioiniat=qrioini ! vazao da condicao inicial de muskingum cunge em cada subtrecho da celula
sssuat=sssu ! armazenamento subsuperficial
sradat=srad ! armazenamento radicular
ssubat=ssub ! armazenamento subterraneo
qsupat=qsup ! fluxo superficial
qssuat=qssu ! fluxo subsuperficial
qsubat=qsub ! fluxo subterraneo    

!f1=0. ! fator de correcao
!f2=0.
!dif1=0.
dif2=0.
iter=0

do while (iter < 15)
    nbac=0 ! indica o numero de subbacias a serem corrigidas
    f3=1.
    do ib=1,nobs
        if(postorio(ib) == 0) cycle ! postorio indica quais postos serao usados para corrigir vazao na bacia ib, 0 nao corrige
        if(qobs(postorio(ib),it) < 0.) cycle ! nao existe dado observado
        
        qr(postorio(ib),it)=qj2(iqobs(postorio(ib)))
        dif2(ib)=(qr(postorio(ib),it)-qobs(postorio(ib),it))/(qobs(postorio(ib),it)+0.001)
        
        if(abs(dif2(ib)) < eps(ib)) cycle ! nao corrige pois o erro eh menor que a tolerancia
        
        nbac=nbac+1 ! numero de subbacias a serem corrigidas
                
!        if(iter > 2) then ! tem pelo menos 2 pontos calculados, aplica a secante   
!            f3(ib)=f2(ib)-dif2(ib)*(f2(ib)-f1(ib))/(dif2(ib)-dif1(ib))
!        else ! aplica uma aproximacao proporcional
            f3(ib)=qobs(postorio(ib),it)/(qr(postorio(ib),it)+0.001)
!        endif
        
    enddo ! fim loop bacias
    
    do ic=1,nc
        ib=ibac(ic)

        if(f3(ib) /= 1.) then ! aplica a correcao somente nas subbacia ib onde o erro for grande
            ! atualiza armazenamento superficial, subsuperficial e subterranea na celula
            do iu=1,nu
    !	        sssuat(ic,iu) = sssuat(ic,iu)*f3(ib)
    !	        sradat(ic,iu) = sradat(ic,iu)*f3(ib) 
                ssubat(ic,iu) = ssubat(ic,iu)*f3(ib)	            
            enddo
            ! atualiza a vazao superficial, subsuperficial e subterranea na celula
            qsupat(ic)=qsupat(ic)*f3(ib)
            qssuat(ic)=qssuat(ic)*f3(ib)
            qsubat(ic)=qsubat(ic)*f3(ib)
            qmat(ic)=qmat(ic)*f3(ib)
            qjat(ic)=qjat(ic)*f3(ib)
            qrioiniat(ic,1:ntrmax+1)=qrioiniat(ic,1:ntrmax+1)*f3(ib)

            
            if(adren(ic)/acel(ic)>1.5)then
                ! corrige os volumes no canal
                vt = vrsupat(ic)+vrssuat(ic)+vrsubat(ic)
                vaux = (1.-1/f3(ib))*dtp*(pmsupat(ic)*qmat(ic)-vrsupat(ic)*qjat(ic)/vt)
                vrsupat(ic) = max(vrsupat(ic)+vaux,0.)
                vaux = (1.-1/f3(ib))*dtp*(pmssuat(ic)*qmat(ic)-vrssuat(ic)*qjat(ic)/vt)
                vrssuat(ic) = max(vrssuat(ic)+vaux,0.)
                vaux = (1.-1/f3(ib))*dtp*(pmsubat(ic)*qmat(ic)-vrsubat(ic)*qjat(ic)/vt)
                vrsubat(ic) = max(vrsubat(ic)+vaux,0.)
            end if
        endif
    enddo
    
    if(nbac == 0) exit ! nenhuma subbacia precisando correcao

    iter = iter + 1    
    ! no inicio do intervalo de tempo adota os valores atualizados
!   sssu = sssuat
!   srad = sradat
    ssub = ssubat
    qsup = qsupat
    qssu = qssuat
    qsub = qsubat
    ! atualiza os volumes superficial, subsuperficial, e subterraneo
!    vsup(1:nc)=qsup(1:nc)*dtp*tks(1:nc)/(1-tks(1:nc))
!    vssu(1:nc)=qssu(1:nc)*dtp*tks(1:nc)/(1-tks(1:nc))
!    vsub(1:nc)=qsub(1:nc)*dtp*tkb(1:nc)/(1-tkb(1:nc))
    vsup(1:nc)=qsup(1:nc)*dtp/(1-tks(1:nc))
    vssu(1:nc)=qssu(1:nc)*dtp/(1-tks(1:nc))
    vsub(1:nc)=qsub(1:nc)*dtp/(1-tkb(1:nc))    
    
    ! atualiza com os valores de vazoesz do incio do intervalo
    qm2=qmat
    qj2=qjat
    qcel2=qcelat    
    pmsup=pmsupat
    pmssu=pmssuat
    pmsub=pmsubat
    vrsup=vrsupat
    vrssu=vrssuat
    vrsub=vrsubat
    qrioini=qrioiniat
    
    call celula
    call rede ! propaga
    
    ! substitui as vazoes calculadas
    do ib=1,nobs
        qr(ib,it)=qj2(iqobs(ib)) !qr(ib,it) vai ser comparado a qobs(ib,it)
!        f1(ib)=f2(ib)
!        f2(ib)=f3(ib)
!        dif1(ib)=dif2(ib)
    enddo
enddo

deallocate (qcelat,qmat,pmsupat,pmssuat,pmsubat,qjat,qrioiniat) ! vazoes na rede anteriores 
deallocate (vrsupat,vrssuat,vrsubat) ! volumes armazenados no rio anteriores
deallocate (sssuat,sradat,ssubat)
deallocate (qsupat,qssuat,qsubat)
	
return
end