!---------------------------------------------------------------------------------------------
!subrotina de calibracao
!
!a calibracao eh baseada na metodologia
!Shuffled Complex Evolution Method for Global Optimization - SCE-UA version 2.1
!desenvolvida por Qingyun Duan et al.
!Department of Hydrology & Water Resources
!University of Arizona, Tucson, AZ 85721
!(602)621-9360, email: duan@hwr.arizona.edu
!
!declaracao do autor:
!     This general purpose global optimization program is developed at
!     the Department of Hydrology & Water Resources of the University
!     of Arizona. Further information regarding the SCE-UA method can
!     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V. K. Gupta
!     at the address and phone number listed above. We request all
!     users of this program make proper reference to the paper entitled
!     'Effective and Efficient Global Optimization for Conceptual
!     Rainfall-runoff Models' by Duan, Q.; S. Sorooshian and V. K. Gupta,
!     Water Resources Research, vol 28(4), pp.1015-1031, 1992.
!---------------------------------------------------------------------------------------------
subroutine calibra

 	use vars_main,only: dir_dados,filcal,filajuste,filevol,nbvfo,nb,titulo,cb
	use vars_calib,only: maxn,kstop,pcento,ngs,iseed,npg,nps,nspl,mings,iniflg,iprint,ibcal, &
	                   & pos,npar,nopt,a,blow,bupp,vpar,cb_ref,xname,npt
	implicit none
	integer i,j,jseed(10),ib 
    integer ierror,iwarn,ideflt
    integer nrun
    character*10 pcntrl,deflt,usrsp
    character*4 reduc,initl,ysflg,noflg
    character*200 fmt !fromatacao
    data deflt/' default  '/
    data usrsp/'user spec.'/
    data ysflg/'yes '/
    data noflg/'no  '/
    data jseed/2,3,5,7,11,13,17,19,23,29/
    
    !---------------------------------------------------------------------------------------------
    ! lista das variaveis de entrada:
    !     a(.) = conjunto de parametros iniciais
    !     blow(.) = limite inferior dos parametros
    !     bupp(.) = limite superior dos parametros
    !     nopt = numero de parametros a serem otimizados
    !     pos(i) = loc(i) = indicador da ordem dos parametros a serem otimizados
    !
    ! lista dos parâmetros de controle do algoritimo SCE:
    !     ngs = numero de complexos para a populacao inicial
    !     npg = numero de pontos em cada complexo
    !     npt = numero total de pontos na populacao inicial (npt=ngs*npg)
    !     nps = numero de pontos em um sub-complexo
    !     nspl = numero de evolucoes permitidas para cada complexo antes
    !            da combinacao entre complexos 
    !     mings = numero minimo de complexos necessario, se o numero de
    !             complexos permite reduzir o rendimento da otimizacao
    !     iseed = semente aleatoria inicial
    !     iniflg = indicador sobre a possibilidade de inclusao de valores iniciais na populacao
    !         = 0, nao inclui
    !         = 1, inclui
    !     iprint = indicador para o controle de impressao de saidas apos cada ciclo de combinacao 
    !         = 0, imprime informacoes do melhor ponto da populacao
    !         = 1, imprime informacoes de cada ponto da populacao
    !
    !  parametros de verificacao de convergencia:
    !     maxn = numero maximo de tentativas permitidas antes da otimizacao ser encerrada
    !     kstop = numero de ciclos de combinacao nos quais o valor de criterio deva ser alterado,
    !             considerando uma percentagem dada, antes da otimizacao ser encerrada
    !     pcento = percentagem que o valor de criterio deve mudar em um dado numero de 
    !              ciclos de combinacao
    !     ipcnvg = indicador de convergencia, ou seja, verifica se gnrng eh inferior a 0.001
    !         = 0, parametro de convergencia nao satisfeito
    !         = 1, parametro de convergencia satisfeito
    !---------------------------------------------------------------------------------------------

    ! leitura do arquivo com parametros de calibracao
	open(filcal,file=dir_dados//'entrada/Calibra.hig',status='old') 
    read(filcal,*)
    read(filcal,*) maxn,kstop,pcento,ngs,iseed,ideflt
    
    if (iseed==0) iseed = 1969 ! valor padrao

    read(filcal,*)
    read(filcal,*)
    if (ideflt == 1) then ! faz a leitura dos parametros de controle do SCE definidos pelo usuario
        read(filcal,*) npg,nps,nspl,mings,iniflg,iprint
        pcntrl = usrsp ! parametros definidos pelo usuario
    else
        read(filcal,*)
        pcntrl = deflt ! parametros padrao
    end if
    
    read(filcal,*)
    read(filcal,*)
	read(filcal,*) ibcal,nbvfo ! indica qual bacia sera calibrada e a estacao associada cujos dados serao utilizados na calibracao
	
	call dominio ! selecionando apenas o dominio da bacia que se deseja calibrar
	
	read(filcal,*)
	read(filcal,*)(pos(i),i=1,npar) ! verifica quais variaveis serao calibradas (1 = calibra, 0 = nao calibra)
	read(filcal,*)
	
    ! calcula o numero de parametros que serao otimizados
    nopt=0
    do i=1,npar
        if(pos(i) > 0) then
            nopt=nopt+1
            pos(nopt)=i ! pos passa a ser um vetor de posicao
        endif
    enddo

    allocate(a(nopt),blow(nopt),bupp(nopt)) 

    ! leitura dos parametros de contorno
    read(filcal,*) (vpar(i),i=1,npar)
    do i=1,nopt
        blow(i)=vpar(pos(i)) ! valor minimo
    enddo
    read(filcal,*) (vpar(i),i=1,npar)
    do i=1,nopt
        bupp(i)=vpar(pos(i)) ! valor maximo
    enddo
    
	close(filcal)  ! fim leitura arquivo de calibracao
	
	! definindo velores padrao dos parametros de controle do SCE
	if (ideflt == 0) then
        npg = 2*nopt + 1
        nps = nopt + 1
        nspl = npg
        mings = ngs
        iniflg = 0
        iprint = 0
    endif
    
    
    ! leitura dos fatores de ajuste originais do arquivo ParAjuste.hig
	open(filajuste,file=dir_dados//'entrada/ParAjuste.hig',status='old')
    read(filajuste,*)
    do ib=1,nb
        if(ib == ibcal) then
            read(filajuste,'(a10,10g10.0)') titulo,(vpar(i),i=1,npar)
        else
            read(filajuste,*)
        endif
    enddo
	close (filajuste)

    ! calcula o valor de cb de referencia em dias para atualizacao de parametros
    cb_ref=cb(ibcal)/(vpar(npar)*86400.)


    ! abre arquivo evolucao
    write(*,*)filevol,dir_dados//'saida/evolucao.hig'
    open(filevol,file=dir_dados//'saida/evolucao.hig')
    write(filevol,'(10x,a46,/,10x,46(1h=))' ) 'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION'
    write(*,'(a46,/,46(1h=))' ) 'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION'
    write(filevol,'(//,10x,a59,/,10x,59(1h=))') 'read and write the input information for hydrological model'
    write(*,'(//,a59,/,59(1h=))') 'read and write the input information for hydrological model'
    write(filevol,'(//,a25,5x,i10,//,a30,i10)') 'sub-basin to be optimized:',ibcal,'station used in the optimization:',nbvfo
    write(*,'(//,a25,5x,i10,//,a30,i10)') 'sub-basin to be optimized:',ibcal,'station used in the optimization:',nbvfo
    write(fmt,'(2(a,i2.2),a)')"(//,a25,//,",npar,"(4x,a6),/,",npar,"(f10.4))"
    write(filevol,fmt) 'initial parameter values:',(xname(j),j=1,npar),(vpar(j),j=1,npar)
    write(*,fmt) 'initial parameter values:',(xname(j),j=1,npar),(vpar(j),j=1,npar)
    write(fmt,'(a,i2.2,a)')"(//,a31,//,",nopt,"(4x,a6))"
    write(filevol,fmt) 'the parameters to be optimized:',(xname(pos(j)),j = 1, nopt)
    write(*,fmt) 'the parameters to be optimized:',(xname(pos(j)),j = 1, nopt)

    ! verifica se os parametros de controle do SCE sao validos
    ierror = 0
    iwarn = 0
    if (ngs < 1 .or. ngs >= 1320) then
        write(filevol,900) ngs
        write(*,900) ngs
900     format(//,1x,'**error** number of complexes in initial population ',&
     &  i5,' is not a valid choice')
        ierror = ierror + 1
    endif

    if (kstop < 0 .or. kstop >= 20) then
        write(filevol,901) kstop
        write(*,901) kstop
901     format(//,1x,'**warning** the number of shuffling loops in',&
     &  ' which the criterion value must change ',/,13x,'should be',&
     &  ' greater than 0 and less than 10.  ','kstop = ',i2,&
     &  ' was specified.'/,13x,'but kstop = 5 will be used instead.')
        iwarn = iwarn + 1
        kstop=5
    endif

    if (mings < 1 .or. mings > ngs) then
        write(filevol,902) mings
        write(*,902) mings
902     format(//,1x,'**warning** the minimum number of complexes ',&
     &         i2,' is not a valid choice. set it to default')
        iwarn = iwarn + 1
        mings = ngs
    endif

    if (npg < 2 .or. npg > 1320/max(ngs,1)) then
        write(filevol,903) npg
        write(*,903) npg
903     format(//,1x,'**warning** the number of points in a complex ',&
     &         i4,' is not a valid choice, set it to default')
        iwarn = iwarn + 1
        npg = 2*nopt+1
    endif

    if (nps<2 .or. nps>npg .or. nps>50) then
        write(filevol,904) nps
        write(*,904) nps
904     format(//,1x,'**warning** the number of points in a sub-','complex ',&
     &  i4,' is not a valid choice, set it to default')
        iwarn = iwarn + 1
        nps = nopt + 1
    endif

    if (nspl < 1) then
        write(filevol,905) nspl
        write(*,905) nspl
905     format(//,1x,'**warning** the number of evolution steps taken in each complex before shuffling ',i4,/,13x,&
     &         'is not a valid choice, set it to default')
        iwarn = iwarn + 1
        nspl = npg
    endif

    ! calculo do numero total de pontos da populacao inicial
    npt = ngs * npg
    if (npt > 1320) then
        write(filevol,906) npt
        write(*,906) npt
906     format(//,1x,'**warning** the number of points in initial population ',i5,' exceed the population limit,',/,13x,&
     &         'set ngs to 2, and npg, nps and nspl to defaults')
        iwarn = iwarn + 1
        ngs = 2
        npg = 2*nopt + 1
        nps = nopt + 1
        nspl = npg
    endif

    do i=1,nopt
        a(i)=vpar(pos(i)) ! aloca os parametros a serem otimizados no vetor a
        if(a(i) < blow(i) .or. a(i) > bupp(i)) then ! verifica se os parametros estao dentro do intervalo admissivel
            ierror=ierror+1
            write(filevol,917) xname(pos(i))
            write(*,917) xname(pos(i))
917         format(//,8x,'*** parametros ',a6,' fora do intervalo admissivel ****')            
        endif
    enddo

    ! imprimindo o numero total de erros e mensagens de aviso
    if (ierror >= 1) write(filevol,907) ierror
907 format(//,1x,'*** total number of error messages is ',i2)

    if (iwarn >= 1) write(filevol,908) iwarn
908 format(//,1x,'*** total number of warning messages is ',i2)

    if (mings < ngs) then
        reduc = ysflg
    else
        reduc = noflg
    endif

    if (iniflg .ne. 0) then
        initl = ysflg
    else
        initl = noflg
    endif

    ! imprimindo as opcoes de otimizacao da evolucao dos complexos combinados
    write(filevol,910)
910 format(//,2x,'SCE control',5x,'max trials',5x,'required improvement',5x,'random',/,3x,'parameter',8x,&
   &'allowed',6x,'percent',4x,'no. loops',6x,'seed',/,2x,11(1h-),5x,10(1h-),5x,7(1h-),4x,9(1h-),5x,6(1h-))
    write(filevol,912) pcntrl,maxn,pcento*100.,kstop,iseed
912 format(3x,a10,7x,i5,10x,f3.1,9x,i2,9x,i5)
    write(filevol,914) ngs,npg,npt,nps,nspl
914 format(//,18x,'SCE algorithm control parameters',/,18x,32(1h=),//,2x,'number of',5x,'points per',5x, &
   & 'points in',6x,'points per',4x,'evol. steps',/,2x,'complexes',6x,'complex',6x,'ini. popul.',5x, &
   & 'sub-complx',4x,'per complex',/,2x,9(1h-),5x,10(1h-),4x,11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
    write(filevol,915) reduc,mings,initl
915 format(//,15x,'complx no.',5x,'min complex',5x,'ini. point',/,15x,'reduction',6x,'no. allowed',6x,'included',/,&
   &15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
    write(filevol,916)
916 format(//,8x,'initial parameter values and parameter bounds',/,8x,45(1h=),//,2x,'parameter',5x,'initial value',5x,&
   &'lower bound',5x,'upper bound',/,2x,9(1h-),5x,13(1h-),5x,11(1h-),5x,11(1h-))
    do i = 1, nopt
        write(filevol,'(2x,a6,4x,3(6x,f10.4))') xname(pos(i)),a(i),blow(i),bupp(i)
    enddo
    if (ierror >= 1) then
        write(filevol,922)
922     format(//,'*** the optimization search is not conducted because of input data error ***')
        stop
    endif

    ! transformacao logaritmica das variaveis de otimizacao
    a(1:nopt)=log(a(1:nopt))
    blow(1:nopt)=log(max(blow(1:nopt),0.0001))
    bupp(1:nopt)=log(max(bupp(1:nopt),0.0001)) 

   ! otimiza
   if (iseed > 0) then
        nrun = min(iseed,10)
   else
        nrun = 1
   endif
   j=npar
   do i=1, nrun
        if (nrun .ne. 1) iseed = -jseed(i)
        write (*,*) '@ SCE-UA run number',i,' random seed value',iseed
        call sceua(filevol)
    enddo

    deallocate(a,blow,bupp)
    return
end subroutine calibra

!====================================================================
subroutine sceua(filevol)
    use vars_main,only: nash,lnash
    use vars_calib,only: nopt,ngs,npg,iseed,npt,blow,bupp,xname,a,pos,iniflg,maxn, &
                       & iprint,nspl,mings,xnstd,bound,x,cx,cf,nps,kstop,pcento,ibcal
    implicit none
    integer filevol !entrada da rotina
    integer i,j,k,k1,k2
    integer loop,nloop,igs,icall,ipcnvg,lpos,lcs(nps)
    integer npt1,ngs1,ngs2,iseed1
    real*8,dimension(nopt):: unit,xx,bestx,worstx
    real*8 s(nps,nopt),criter(20),xf(ngs*npg),sf(nps)
    real*8 bestf,worstf,rand,ran1,gnrng,fa,functn,denomi,timeou
    character*200 fmt !fromatacao
    
    allocate(xnstd(nopt),bound(nopt),x(ngs*npg,nopt),cx(ngs*npg,nopt),cf(ngs*npg))
    
    !---------------------------------------------------------------------------------------------
    !  lista de variaveis
    !     x(.,.) = coordinates of points in the population
    !     xf(.) = function values of x(.,.)
    !     xx(.) = coordinates of a single point in x
    !     cx(.,.) = coordinates of points in a complex
    !     cf(.) = function values of cx(.,.)
    !     s(.,.) = coordinates of points in the current simplex
    !     sf(.) = function values of s(.,.)
    !     bestx(.) = best point at current shuffling loop
    !     bestf = function value of bestx(.)
    !     worstx(.) = worst point at current shuffling loop
    !     worstf = function value of worstx(.)
    !     xnstd(.) = standard deviation of parameters in the population
    !     gnrng = normalized geometric mean of parameter ranges
    !     lcs(.) = indices locating position of s(.,.) in x(.,.)
    !     bound(.) = bound on ith variable being optimized
    !     ngs1 = number of complexes in current population
    !     ngs2 = number of complexes in last population
    !     iseed1 = current random seed
    !     criter(.) = vector containing the best criterion values of the last
    !         10 shuffling loops
    !---------------------------------------------------------------------------------------------
    
    write (*,*) ' enter the sceua subroutine --- '     
    
    ! inicializando variaveis
    nloop = 0
    loop = 0
    igs = 0
    
    ! inicializando semente aleatoria com um inteiro negativo
    iseed1 = -abs(iseed)
    
    ! calculando o numero total de pontos da populacao inicial
    npt = ngs * npg
    ngs1 = ngs
    npt1 = npt
    
    write(filevol,400)
400 format(//,2x,50(1h=),/,2x,'enter the shuffled complex evolution global search',/,2x,50(1h=))
    write (*,*) ' ***  evolution loop number ',nloop
    
    ! calculo do limite dos parametros que serao otimizados
    bound = bupp - blow
    unit = 1.0
    
    ! calculo da funcao do ponto inicial
    fa = functn(a)
    
    ! imprime o ponto inicial e seu valor de criterio
    write(filevol,500)
500 format(//,'*** print the initial point and its criterion value ***')
    write(fmt,'(a,i2.2,a,i3.3,a)')"(/,' criterion','    nash','   lnash',",nopt,"(4x,a6),/1x,",10*nopt+25,"(1h-))"
    write(filevol,fmt) (xname(pos(j)),j=1,nopt)
    write(fmt,'(a,i2.2,a)')"(g10.3,2f8.2,",nopt,"f10.4)"
    write(filevol,fmt) fa,nash(ibcal),lnash(ibcal),dexp(a(1:nopt))
    
    ! gera um conjunto inicial de npt1 pontos dentro do espaco de parametros
    if (iniflg == 1) then ! utilizar parametros definidos pelo ParAjuste
        x(1,1:nopt) = a(1:nopt)
        xf(1) = fa
    else ! gera valores aleatorios para os parametros
        call getpnt(1,iseed1,xx,unit,blow)
        x(1,1:nopt) = xx(1:nopt)
        xf(1) = functn(xx)
    endif
    
    
    icall = 1
    if (icall >= maxn) go to 900
    
    ! gera npt1-1 pontos aleatorios distribuidos uniformemente dentro do
    ! espaco de parametros e calcula as funcoes correspondentes
    do i = 2, npt1
        call getpnt(1,iseed1,xx,unit,blow)
        x(i,1:nopt) = xx(1:nopt)
        xf(i) = functn(xx)
        icall = icall + 1
        if (icall >= maxn) then
            npt1 = i
            exit
        endif
    enddo
    
    ! ordena os pontos na ordem crescente da funcao
    call sort1(npt,npt1,nopt,x,xf)
    
    ! guarda o melhor e o pior pontos
    bestx(1:nopt) = x(1,1:nopt)
    worstx(1:nopt) = x(npt1,1:nopt)
    bestf = xf(1)
    worstf = xf(npt1)
    
    fa = functn(bestx(1:nopt)) ! recalculando Nash, LNash para o melhor ponto
    
    ! calcula o intervalo de parametros para a populacao inicial
    call parstt(npt1,gnrng,ipcnvg)
    
    ! imprime os resultados para a populacao inicial
    write(filevol,600)
600 format(//,1x,'*** print the results of the SCE search ***')
    write(fmt,'(2a,i2.2,a)')"(/,1x,'loop',1x,'trials',1x,'complxs',2x,' nash   ','lnash   ','best f',3x,'worst f',", &
                         & "2x,'par rng',2x,",nopt,"(4x,a6))" !610
    write(filevol,fmt) (xname(pos(j)),j=1,nopt)
    write(fmt,'(a,i2.2,a)')"(i5,1x,i5,3x,i5,2f8.2,3g10.3,",nopt,"f10.4)" !612
    write(filevol,fmt) nloop,icall,ngs1,nash(ibcal),lnash(ibcal),bestf,worstf,gnrng,dexp(bestx(1:nopt))
    if (iprint == 1) then
        write(filevol,650) nloop
        write(fmt,'(2a,i2.2,a)')"(15x,2f8.2,g10.3,20x,",nopt,"(f10.4))"
        do i = 1, npt1
            write(filevol,fmt) xf(i),nash(ibcal),lnash(ibcal),dexp(x(i,1:nopt))
        enddo
    endif
    
    if (icall >= maxn) go to 900
    if (ipcnvg == 1) go to 920
    
    ! inicio do loop principal
1000 continue
    nloop = nloop + 1
    
    write (*,*) ' ***  evolution loop number ',nloop
    
    ! inicio do loop de complexos
    do igs = 1, ngs1
        
        ! atribui pontos para os complexos
        do k1 = 1, npg
            k2 = (k1-1) * ngs1 + igs
            cx(k1,1:nopt) = x(k2,1:nopt)
            cf(k1) = xf(k2)
        enddo
        
        ! inicio do loop interno - selecao aleatoria de sub-complexos
        do loop = 1,nspl
            
            ! escolhe um sub-complexo (nps pontos) de acordo com uma distribuicao
            ! linear de probabilidades
            if (nps==npg) then
                !lcs(1:nps) = (/1:nps/)
                do i = 1,nps
                    lcs(i) = i
                end do
            else
                rand = ran1(iseed1)
                lcs(1) = 1 + dint(npg + 0.5 - dsqrt( (npg+.5)**2 - npg * (npg+1) * rand ))
                do k = 2, nps
    60              rand = ran1(iseed1)
                    lpos = 1 + dint(npg + 0.5 - dsqrt((npg+.5)**2 - npg * (npg+1) * rand ))
                    do k1 = 1, k-1
                        if (lpos == lcs(k1)) go to 60
                    enddo
                    lcs(k) = lpos
                enddo
                
                ! ordena os sub-complexos em ordem crescente da funcao
                call sort2(nps,lcs)
            endif
            
            ! cria o conjunto de sub-complexos
            do k = 1, nps
                s(k,1:nopt) = cx(lcs(k),1:nopt)
                sf(k) = cf(lcs(k))
            enddo
            
            ! usa o sub-complexo para gerar novos pontos
            call cce(s,sf,icall,iseed1)
            
            ! se o sub-complexo eh aceito, substitui ele no complexo
            do k = 1, nps
                cx(lcs(k),1:nopt) = s(k,1:nopt)
                cf(lcs(k)) = sf(k)
            enddo
            
            ! ordena os pontos
            call sort1(npt,npg,nopt,cx,cf)
            
            ! se o numero maximo de rodadas for excedida, encera-se o loop
            if (icall >= maxn) go to 2222
            
        enddo ! fim do loop interno
2222    continue
        
        ! substitui o novo complexo dentro do conjunto original x(.,.)
        do k1 = 1, npg
            k2 = (k1-1) * ngs1 + igs
            x(k2,1:nopt) = cx(k1,1:nopt)
            xf(k2) = cf(k1)
        end do
        if (icall >= maxn) go to 3333
    
    enddo ! fim do loop de complexos
    
    ! reordena os pontos
3333 call sort1(npt,npt1,nopt,x,xf)
    
    ! guarda o melhor e o pior pontos
    bestx(1:nopt) = x(1,1:nopt)
    worstx(1:nopt) = x(npt1,1:nopt)
    bestf = xf(1)
    worstf = xf(npt1)
    fa = functn(bestx) ! recalculando Nash, LNash para o melhor ponto
    
    ! teste de convergencia dos parametros da populacao
    call parstt(npt1,gnrng,ipcnvg)
    
    ! imprime os resultados da populacao atual
    write(fmt,'(2a,i2.2,a)')"(/,1x,'loop',1x,'trials',1x,'complxs',2x,' nash   ','lnash   ','best f',3x,'worst f',", &
                         & "2x,'par rng',2x,",nopt,"(4x,a6))" !610
    if (mod(nloop,5) == 0) then
        write(filevol,fmt) (xname(pos(j)),j=1,nopt)
    endif
    write(fmt,'(a,i2.2,a)')"(i5,1x,i5,3x,i5,2f8.2,3g10.3,",nopt,"f10.4)" !612
    write(filevol,fmt) nloop,icall,ngs1,nash(ibcal),lnash(ibcal),bestf,worstf,gnrng,dexp(bestx(1:nopt))
    if (iprint == 1) then
        write(filevol,650) nloop
650     format(/,1x,'population at loop ',i3,/,1x,22(1h-))
        write(fmt,'(a,i2.2,a)')"(15x,2f8.2,g10.3,20x,",nopt,"(f10.4))"
        do i = 1, npt1
            write(filevol,fmt) xf(i),nash(ibcal),lnash(ibcal),dexp(x(i,1:nopt))
        enddo
    endif
    
    ! verifica se o numero maximo de avaliacoes da funcao foi excedido
    if (icall >= maxn) go to 900
    
    ! calcula a quantidade de loop sucessivos que nao apresentaram melhorias na funcao
    criter(20) = bestf
    if (nloop >= (kstop+1)) then
        denomi = dabs(criter(20-kstop) + criter(20)) / 2.
        timeou = dabs(criter(20-kstop) - criter(20)) / denomi
        if (timeou < pcento) go to 910
    endif
    criter(1:19) = criter(2:20)
    
    ! se a populacao convergiu dentro de um espaco suficientemente pequeno
    if (ipcnvg == 1) go to 920
    
    ! nenhum dos criterios foi satisfeito, continua a busca
    ! verifica a reducao do numero de complexos
    if (ngs1 > mings) then
        ngs2 = ngs1
        ngs1 = ngs1 - 1
        npt1 = ngs1 * npg        
        call comp(npt1,ngs1,ngs2,xf)
    endif
    
    go to 1000 ! fim do loop principal
    
    ! busca finalizada
900 write(filevol,800) maxn,loop,igs,nloop
800 format(//,1x,'*** optimization search terminated because the',&
    &       ' limit on the maximum',/,5x,'number of trials ',i5,&
    &       ' exceeded.  search was stopped at',/,5x,'sub-complex ',&
    &       i3,' of complex ',i3,' in shuffling loop ',i3,' ***')
    go to 999
910 write(filevol,810) pcento*100.,kstop
    810 format(//,1x,'*** optimization terminated because the criterion',&
    &       ' value has not changed ',/,5x,f5.2,' percent in',i3,&
    &       ' shuffling loops ***')
    go to 999
920 write(filevol,820) gnrng*100.
820 format(//,1x,'*** optimization terminated because the population',&
    &       ' has converged into ',/,4x,f5.2,' percent of the',&
    &       ' feasible space ***')
    
    ! imprime a estimativa final dos parametros e sua funcao
999 write(filevol,830)
830 format(//,'*** print the final parameter estimate and its',&
    &       ' criterion value ***')
    
    write(fmt,'(a,i2.2,a,i3.3,a)')"(/,' criterion','    nash','   lnash',",nopt,"(4x,a6),/1x,",10*nopt+25,"(1h-))"
    write(filevol,fmt) (xname(pos(j)),j=1,nopt)
    write(fmt,'(a,i2.2,a)')"(g10.3,2f8.2,",nopt,"f10.4)"
    write(filevol,fmt) bestf,nash(ibcal),lnash(ibcal),dexp(bestx(1:nopt))

    deallocate (x,xnstd,bound)
    deallocate (cx,cf)
    return
end ! fim da sub-rotina sceua

!====================================================================
subroutine cce(s,sf,icall,iseed)
! algoritmo que gera um novo ponto a partir de um sub-complexo
    use vars_calib, only: nopt,nps,blow,bupp,maxn,xnstd
    implicit none
    real*8,parameter:: c1=0.8,c2=0.4
    integer j,n,m,ibound
    integer icall,iseed
    real*8 s(nps,nopt),sf(nps)
    real*8 sb(nopt),sw(nopt),ce(nopt),snew(nopt)
    real*8 alpha,beta,fw,fnew,functn
    
    !-----------------------------------------------------------------
    ! lista de variaveis locais
    !    sb(.) = o melhor ponto do Simplex
    !    sw(.) = o piot ponto do Simplex
    !    fw = valor da funcao do pior ponto
    !    ce(.) = centroide do Simplex excluindo o pior ponto
    !    snew(.) = novo ponto gerado a aprtir do Simplex
    !-----------------------------------------------------------------
    
    ! inicializando variaveis
    n = nps
    m = nopt
    alpha = 1.0 !esse alpha nao eh o mesmo do vars_main!
    beta = 0.5
    
    ! identificacao do pior ponto do sub-complexo s
    sb(1:m) = s(1,1:m)
    sw(1:m) = s(n,1:m)
    ! calculo do centroide dos demais pontos
    ! calculo do passo, vetor entre o pior ponto e o centroide
    do j = 1, m
        ce(j) = 0.0
        ce(j) = (ce(j) + sum(s(1:n-1,j)))/dble(n-1)
    enddo
    ! identificacao do pior valor da funcao
    fw = sf(n)
    
    ! calculo do novo ponto (snew)
    ! tentando primeiro o reflection step
    snew(1:m) = ce(1:m) + alpha * (ce(1:m) - sw(1:m))
    
    ! verifica se snew satisfaz todas as restricoes
    call chkcst(nopt,snew,blow,bupp,ibound)
    
    ! se snew estiver fora do contorno,
    ! escolher um ponto ao acaso dentro da regiao viavel de acordo com
    ! uma distribuicao normal com o melhor ponto do sub-complexo
    ! como media e desvio padrao da populacao como std
    if (ibound>=1) call getpnt(2,iseed,snew,xnstd,sb)
    
    ! calculo da funcao para snew
    fnew = functn(snew)
    icall = icall + 1
    
    ! compara fnew com o pior valor da funcao (fw)
    if (fnew<=fw) go to 20 ! aceita o novo ponto
    if (icall>=maxn) return
    
    ! se fnew é maior que fw, tenta o contraction step
    snew(1:m) = ce(1:m) - beta * (ce(1:m) - sw(1:m))
    
    ! calculo da funcao para o ponto contraction step
    fnew = functn(snew)
    icall = icall + 1
    
    ! compara fnew com o pior valor da funcao (fw)
    if (fnew <= fw) go to 20 ! aceita o novo ponto
    if (icall >= maxn) return
    
    ! se o reflection step e o contraction step falharem, escolher outro ponto
    ! de acordo com a distribuicao normal com o melhor ponto do sub-complexo
    ! como media e desvio padrao da populacao como std
    call getpnt(2,iseed,snew,xnstd,sb)
    
    ! calculo da funcao para o ponto aleatorio
    fnew = functn(snew)
    icall = icall + 1
    
    ! substitui o pior ponto pelo novo ponto
20  continue
    s(n,1:m) = snew(1:m)
    sf(n) = fnew
    return

end subroutine cce

!===================================================================
subroutine getpnt(idist,iseed,x,std,xi)
! essa subrotina gera um novo ponto dentro da regiao viavel
    use vars_calib, only: nopt,blow,bupp
    implicit none
    integer idist,iseed
    integer j,ibound
    real*8 x(nopt),std(nopt),xi(nopt)
    real*8 rand
    real*8 ran1,gasdev !funcoes

    ! lista de variaveis:
    !     x(.) = novo ponto
    !     xi(.) = ponto foco
    !     blow(.) = limite inferior
    !     bupp(.) = limite superior
    !     std(.) = desvio padrao da distribuicao de probabilidades
    !     idist = indicador de probabilidade probability flag
    !           = 1 - distribuicao uniforme
    !           = 2 - distribuicao gaussiana
    !
!    implicit real*8 (a-h,o-z)
!    dimension x(nopt),std(nopt),xi(nopt)
    
    ibound=1
    do while (ibound == 1)
        do j=1,nopt
            ibound=1
            do while (ibound == 1)
                if (idist == 1) then
                    rand = ran1(iseed)
                elseif (idist == 2) then
                    rand = gasdev(iseed)
                endif
                x(j) = xi(j) + std(j) * rand * (bupp(j) - blow(j))
                
                ! verifica restricoes explicitas        
                call chkcst(1,x(j),blow(j),bupp(j),ibound)              
            enddo
        enddo
        
        ! verifica restricoes implicitas     
        call chkcst(nopt,x,blow,bupp,ibound)
    enddo
    
    return
end subroutine getpnt

!===================================================================
subroutine parstt(npt,gnrng,ipcnvg)
! subrotina que verifica a convergencia dos parametros
    use vars_calib, only: x,xnstd,bound,nopt
    implicit none
    real*8,parameter:: delta=1.0d-20, peps=1.0d-3
    integer ipcnvg,npt
    integer k
    real*8 xmax(nopt),xmin(nopt),xmean(nopt)
    real*8 gsum,gnrng
    
    ! calculo do maximo, minimo e desvio padrao dos parametros
    gsum = 0.d0
    do k = 1, nopt
        xmean(k) = sum(x(1:npt,k)) / dble(npt)
        xnstd(k) = sum(x(1:npt,k)*x(1:npt,k))  / dble(npt) - xmean(k)*xmean(k)
        xmax(k) = maxval(x(1:npt,k))
        xmin(k) = minval(x(1:npt,k))
        if (xnstd(k) <= delta) xnstd(k) = delta 
        xnstd(k) = dsqrt(xnstd(k)) / bound(k)
        gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
    enddo
    gnrng = dexp(gsum/dble(nopt))
    
    ! verifica se o desvio padrao normalizado dos parametros é <= eps
    ipcnvg = 0
    if (gnrng <= peps) then
        ipcnvg = 1
    endif
    return
end subroutine parstt

!====================================================================
subroutine comp(npt1,ngs1,ngs2,af)
! essa subrotina reduz a matriz de entrada a(n,ngs2*npg) para a 
! matriz b(n,ngs1*npg) e o vetor af(ngs2*npg) para o vetor bf(ngs1*npg)
    use vars_calib, only: npg,a => x,n => nopt,npt,b => cx, bf => cf
    implicit none
    integer npt1,ngs1,ngs2,nopt
    integer igs,ipg,j,k1,k2
    real*8 af(npt)
    
    do igs=1, ngs1
        do ipg=1, npg
          k1=(ipg-1)*ngs2 + igs
          k2=(ipg-1)*ngs1 + igs
          b(k2,1:n) = a(k1,1:n)
          bf(k2) = af(k1)
        end do
    end do
    
    do j=1, npt1
        a(j,1:n) = b(j,1:n)
        af(j) = bf(j)
    end do
    
    return
end subroutine comp

!===================================================================
subroutine sort1(nmax,n,m,rb,ra)
! subrotina de ordenacao adapatada de "Numerical Recipes"
! W.H. Press et al., p. 233-234
    implicit none
    integer j,n,m,nmax
    integer iwk(nmax)
    real*8 ra(nmax),rb(nmax,m),wk(nmax,m)
    
    ! lista de variaveis:
    !     ra(.) = array to be sorted
    !     rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
    !     wk(.,.), iwk(.) = local varibles
    
    call indexx(nmax,n, ra, iwk)
    wk(1:n,1) = ra(1:n)
    ra(1:n) = wk(iwk(1:n),1)
    do j = 1, m
        wk(1:n,j) = rb(1:n,j)
    enddo
    do j = 1, m
        rb(1:n,j) = wk(iwk(1:n),j)
    enddo
    return
    
end subroutine sort1

!===========================================================
subroutine sort2(n,ra)
! subrotina de ordenacao adapatada de "Numerical Recipes"
! W.H. Press et al., p. 231
    implicit none
    integer i,j,l,n,ir,rra
    integer ra(n)
    
    !  list of variables
    !     ra(.) = integer array to be sorted
    
    l = (n / 2) + 1
    ir = n
10  continue
    if (l > 1) then
        l = l - 1
        rra = ra(l)
    else
        rra = ra(ir)
        ra(ir) = ra(1)
        ir = ir - 1
        if (ir == 1) then
            ra(1) = rra
            return
        endif
    endif
    i = l
    j = l + l
20  continue
    if (j <= ir) then
        if (j < ir) then
            if (ra(j) < ra(j + 1)) j = j + 1
        endif
        if (rra < ra(j)) then
            ra(i) = ra(j)
            i = j
            j = j + j
        else
            j = ir + 1
        endif
        goto 20
    endif
    ra(i) = rra
    goto 10
    
end subroutine sort2

!=======================================================
subroutine indexx(nmax,n, arrin, indx)
! essa subrotina foi obtida do "Numerical Recipes"
! Press et al.
    implicit none
    integer nmax,n,i,j,l,ir,indxt
    integer indx(nmax)
    real*8 arrin(nmax)
    real*8 q
    
    do j = 1, n
        indx(j) = j
    enddo
    l = (n / 2) + 1
    ir = n
10  continue
    if (l > 1) then
        l = l - 1
        indxt = indx(l)
        q = arrin(indxt)
    else
        indxt = indx(ir)
        q = arrin(indxt)
        indx(ir) = indx(1)
        ir = ir - 1
        if (ir == 1) then
            indx(1) = indxt
            return
        end if
    end if
    i = l
    j = l + l
20  continue
    if (j <= ir) then
        if (j < ir) then
            if (arrin(indx(j)) < arrin(indx(j + 1))) j = j + 1
        endif
        if (q < arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
        else
            j = ir + 1
        endif
        goto 20
    endif
    indx(i) = indxt
    goto 10

end subroutine indexx

!==============================================================
real*8 function ran1(idum)
! essa funcao foi obtida do "Numerical Recipes"
! Press et al.
    implicit none
    integer,parameter:: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32
    integer,parameter:: ndiv=1+(im-1)/ntab
    real*8,parameter:: eps=1.2e-7,am=1./im,rnmx=1.-eps
    integer idum,j,k
    integer iv(ntab),iy
    save iv,iy
    data iv /ntab*0/, iy /0/
    
    if (idum<=0.or.iy==0) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum<0) idum=idum+im
            if (j<=ntab) iv(j)=idum
        enddo
        iy=iv(1)
    endif
    
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    
    if (idum<0) idum=idum+im
    
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=idum
    ran1=min(am*iy,rnmx)
    
    return
end function ran1       

!===============================================================
real*8 function gasdev(idum)
! essa funcao foi obtida do "Numerical Recipes"
! Press et al.
    implicit none
    integer idum,iset
    real*8 v1,v2,r,fac,gset
    real*8 ran1 !funcao
    common /gasblk/ iset
    data iset / 0 /
    
    if (iset == 0) then
1       v1 = (2. * ran1(idum)) - 1.
        v2 = (2. * ran1(idum)) - 1.
        r = (v1 ** 2) + (v2 ** 2)
        if (r >= 1.) goto 1
        fac = sqrt(- ((2. * log(r)) / r))
        gset = v1 * fac
        gasdev = v2 * fac
        iset = 1
        else
        gasdev = gset
        iset = 0
    end if
    
    return
end function gasdev

!==================================================================
subroutine chkcst(nopt,x,blow,bupp,ibound)
! essa subrotina verifica se o conjunto de pontos satisfaz todas
! as restricoes
implicit none
integer ibound,nopt,i
real*8 x(nopt),blow(nopt),bupp(nopt)

! lista de variaveis:
!     ibound = indicador de violacao
!            = -1 valor inicial
!            = 0  nao violado
!            = 1  violado
!     nopt = numero de variaveis otimizadas
!     i = a inesima variavel dos arrays x, bl, and bu    

ibound = -1

! verifica se as restricoes explicitas sao violadas
do i=1, nopt
    if (x(i) < blow(i) .or. x(i) > bupp(i)) then
        ibound = 1
        return  
    endif
enddo
if (nopt == 1) then
    ibound = 0
    return
endif

return
end subroutine chkcst

!==================================================================
real*8 function functn(xpar)
! esta funcao prepara os dados para chamar o modelo e calcular a
! funcao
!use ieee_arithmetic
use vars_main, only: nu,ssmax,srmax,smax,alpha,kss,tsub,mu,csi, &
                   & d1,d2,d3,cs,cb,ncdom,icdom,ibac,lambdam, lambdai, &
                   & ac,tkb,dtp,tks,tcon,fob
use vars_calib, only: vpar,pos,solo,ibcal,ths,thr,alpha_ref,ksat, &
                    & tsub_ref,mu_ref,csi_ref,cs_ref,cb_ref,nopt,npar
implicit none
real*8 xpar(nopt)
real par(npar)
integer  iu,ic,ic2,k,iaux

par=vpar ! vetor original de parametros

do k = 1, nopt
    iaux=dexp(xpar(k))*10000.
    par(pos(k))=float(iaux)/10000. ! limita o valor de para 4 digitos significativos
enddo

! ksat condutividade hidraulica saturada (mm/h)
! b parametro da curva de retencao
! psib pressao de entrada do ar (kpa)
! ths umidade volumetrica na saturacao
! thr umidade volumetrica residual     

do iu=1,nu
    if(solo(iu) /= 14) then
        ! calcula os parametros de cada sub-bacia
        ssmax(ibcal,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*par(1) ! armazenamento maximo camada superior do solo
        srmax(ibcal,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*par(2) ! armazenamento maximo em mm na segunda camada
        smax(ibcal,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*par(3) ! armazenamento maximo em mm na terceira camada
        alpha(ibcal,iu)=alpha_ref*par(7) ! coeficiente de anisotropia
        kss(ibcal,iu)=(ksat(solo(iu))*24./1000.)*par(4) ! condutividade hidraulica saturada camada superior m/dia
        tsub(ibcal,iu)=tsub_ref*par(5) ! transmissividade maxima
        mu(ibcal,iu)=mu_ref*par(6) ! 
        csi(ibcal,iu)=min(csi_ref(iu)*par(8),1.) ! capacidade de campo           
    endif
enddo
d1(ibcal)=par(1)
d2(ibcal)=par(2)
d3(ibcal)=par(3)
cs(ibcal)=cs_ref*par(9)
cb(ibcal)=cb_ref*par(10)*86400. ! cb em dias foi multiplicado por 86400 para cb ficar em seg


do ic2=1,ncdom
    ic=icdom(ic2)
    if(ibac(ic) == ibcal) then ! celula que pertence a bacia sendo otimizada
        do iu=1,nu
	        lambdam(ic,iu)=0.0  ! atualiza o valor de lambdam
		    if(smax(ibcal,iu)>0) then ! pula agua
		        do k=2,50
				    lambdam(ic,iu)=lambdam(ic,iu)+0.5*(lambdai(ic,k)**(1./mu(ibcal,iu))+ &
                                  &lambdai(ic,k-1)**(1./mu(ibcal,iu)))*(ac(ic,k)-ac(ic,k-1))
			    enddo
			    lambdam(ic,iu)=lambdam(ic,iu)**mu(ibcal,iu)
		    endif
	    enddo
	    !atualiza os coeficientes do rls 
   		tkb(ic)=exp(-dtp/(cb(ibcal)))	    !cb definido em segundos em lesolo.f90
    	tks(ic)=exp(-dtp/(cs(ibcal)*tcon(ic)))   !cs definido em lesolo.f90	
	endif
enddo

! chama o modelo
call modelo
call fobj

!if(isnan(fob)) then ! caso a combinacao de parametros gere um nan
if(fob/=fob)then ! caso a combinacao de parametros gere um nan
    fob=1.
    fob=huge(fob) ! gera um numero grande
endif

functn = fob
return

end function functn

