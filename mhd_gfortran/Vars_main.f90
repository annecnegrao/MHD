	!---------------------------------------------------------------------------------------------
    ! modulo de declaracao de variaveis do programa principal
    !---------------------------------------------------------------------------------------------
	module vars_main
    implicit none
	save
	character(*),parameter:: dir_dados = 'C:\Users\proclima\Documents\dados_'
	
	integer,parameter:: ntrmax=20 !numero maximo de subtrechos de propagacao muskingum cunge
	real,parameter:: rugman=0.030 !coeficiente de manning e considerado fixo em 0,03 para todos os trechos e rios
	real,parameter:: pi=3.141592
    
	!numeros 	 
	integer nc,nu,ncdom ! numero de celulas, de usos, e de celulas do dominio simulado/calibrado
	integer nt !numero de intervalos de tempo
	integer ntprev ! numero de intervalos de tempo da previsao
	integer nstep ! numero de intervalos de tempo ate o inicio da simulacao (nstep>0, a simulacao se inicia apos o inicio das observacoes) 	
	integer nb !numero de bacias
    integer nobs !numero de series de vazao observada
    integer,allocatable:: iqobs(:) !numero das celulas com dados observados
	integer nhidg !numero de hidrogramas que se deseja gravar
	integer,allocatable:: ihidg(:) !numero das celulas em que se deseja os hidrogramas
	integer numsubst !numero de series de vazao que serao substituidas
	integer,allocatable:: isubst(:) !numero das celulas onde a vazao calculada deve ser substituida pela lida de um arquivo
	integer nmapas !numero de mapas de uso do solo
	integer nbvfo ! numero da bacia para a qual se deseja imprimir as funcoes objetivo e os valors medios das variaveis

	!contadores
	integer it !intervalo de tempo

	!-------variaveis com numeros dos arquivos------------
	!arquivos de entrada
	integer,parameter:: filfix=10 !arquivo dos parametros fixos
	integer,parameter:: filveg=20 !arquivo dos parametros de vegetacao
	integer,parameter:: filuso=30 !arquivo dos parametros associados ao uso
	integer,parameter:: filhig=40 !arquivo das celulas
	integer,parameter:: filobs=60 !arquivo de vazoes observadas
	integer,parameter:: filsubs=70 !arquivo de hidrogramas lidos para subst. calculados
	integer,parameter:: filtop=90 !arquivo de indices topograficos nas celulas
	integer,parameter:: filajuste=100 !arquivo dos fatores de ponderacao dos parametros de solo e da bacia
	integer,parameter:: filbacia=110 !arquivo dos parametros de escoamento da bacia
	integer,parameter:: filsolo=120 !arquivo parametros de solo
	integer,parameter:: filbloc=130 !arquivo indicando a combinacao de vegetacao e solo do arquivo filuso	
	integer,parameter:: filplu=140 !arquivo dos dados de chuva	
	integer,parameter:: filmet=150 !arquivo de dados meteorologicos
	integer,parameter:: filccini=160 !arquivo de dados com condicoes iniciais do sistema usados na previsao
	integer,parameter:: filprev=170 !arquivo com dados da previsao e da rotina de assimilacao
	integer,parameter:: filpluprev=180 !arquivo de dados de previsao de chuva	
	integer,parameter:: filmetprev=190 !arquivo com dados meteorologico previstos
	integer,parameter:: filprevfin=200 !arquivo com dados meteorologico previstos
	integer,parameter:: filcal=210 !arquivo de dados da bacia a ser calibrada

	!arquivos de saida	
	integer,parameter:: filhid=220 !arquivo dos hidrogramas das sub-bacias
	integer,parameter:: filprp=230 !arquivo de hidrograma com proporcoes
	integer,parameter:: filaju=240 !arquivo de ajuste
	integer,parameter:: filbac=250 !arquivo de saida de dados medios das bacias 
	integer,parameter:: filsol=260 !arquivo do solo (1 bloco)
	integer,parameter:: filevol=290 !arquivo com parametros de evolucao da calibracao
	integer,parameter:: filfp=300 !arquivo contendo informacoes da topografia do floodplain
	integer,parameter:: filpmedia=750 !arquivo com prec media diaria na bacia
	!-----------------------------------------------------

	integer icalib,iniaju ! opcao da rodada, e numero de deltats antes inicio da comparacao calculado x observado (spin-up)
	real fob !valor da funcao objetivo de calibracao
	real,allocatable:: nash(:),lnash(:) ! coef. Nash e Log-Nash
	real,allocatable:: qr(:,:),qobs(:,:) !vazao calculada e observada nos exutorios das bacias
	real,allocatable:: qlido(:,:) !vazao em alguns pontos para substituicao de valores calculados
	integer,allocatable:: kcb(:) !numero de celulas em cada sub-bacia 

	!valores medios por sub-bacia
	real,allocatable:: dsupb(:,:),dssub(:,:),dsubb(:,:)
	real,allocatable:: ecb(:,:),etpb(:,:),ettb(:,:),sssub(:,:),preb(:,:),ssubb(:,:),sradb(:,:),asatb(:,:) !valores medios de evap, transp, arma e prec em cada subbacia	
	real,allocatable:: qbsub(:),qbssu(:),qbsup(:)
	
	!datas
	integer idiaArquivo,imesArquivo,ianoArquivo,ihoraArquivo !data e hora do inicio do arquivo meteorologico
    integer idiafimArquivo,imesfimArquivo,ianofimArquivo,ihorafimArquivo !data e hora do fim do arquivo meteorologico	
	integer idia,imes,iano,ihora,ihoraini,iiano,iihoraini !data e hora do inicio da simulacao  
	integer idiafin,imesfin,ianofin,ihorafin !data e hora do fim da simulacao  
	integer idini ! dia em que se inicia a simulacao
	integer jdia !dia do calendario juliano
	real dtp !intervalo de tempo de calculo principal (dado em segundos no arquivo de entrada)
	character cdata*10 ! data em formato character
	character,allocatable:: cdatag(:)*10 ! arquivo com todas as datas em formato character para impressao
	
	character*120 titulo !variavel de leitura dos titulos dos arquivos de entrada
	
	real,allocatable:: qesp(:) !vazao especifica de base (m3/s/km2)
	real,allocatable:: qcons(:) !consumo de agua na bacia no exutorio (m3/s)
	real,allocatable:: alb(:,:),iaf(:,:),zveg(:,:),rc(:) !albedo, indice de area foliar e altura media, resistencia minima superficial
	real,allocatable:: prad(:,:) !profundidade radicular
	real,allocatable:: scmax(:) ! capacidade maxima de armazenamento no dossel
	real,allocatable:: cover(:,:) ! fracao de cobertura do dossel
	real,allocatable:: dveg(:,:),z0(:,:) ! plano de deslocamento zero e rugosidade			
	real:: omegac,fr ! valor critico do fator de ponderacao de estresse e fator de distribuicao de raizes
	
	!numero de subtrechos e intervalo de tempo de calculo muskingum-cunge
	integer,allocatable:: ntrc(:)
	real,allocatable:: dt(:)
	real,allocatable:: cel(:),tks(:),tkb(:) !celeridade no rio e coeficientes do rls da celula
	
	character (40) arqobs 	!nome do arquivo das vazoes observadas
	character (40) arqsubst !nome do arquivo das vazoes lidas que substituirao as calculadas em algumas celulas

	!parametros relacionados ao solo
	real,allocatable:: smax(:,:),kss(:,:),tsub(:,:),alpha(:,:),csi(:,:),ssmax(:,:),mu(:,:),srmax(:,:)
	real,allocatable:: neta(:),d1(:),d2(:),d3(:)
	real,allocatable:: ssub(:,:) !armazenamento na camada inferior do solo
	real,allocatable:: srad(:,:) !armazenamento na zona radicular profunda
	real,allocatable:: secrbs(:) !parametro de transpiracao de solo nu 
	real,allocatable:: secrjah(:),secrjal(:) !parametros de transpiracao de solo nu 	
	real,allocatable:: tcon(:),cb(:),cs(:) !parametros da propagacao na celula
	
	integer,allocatable:: iexut(:) !indica celula do exutorio da bacia
	real,allocatable:: brio(:) !largura do rio
	real,allocatable:: prec(:) !chuva acumulada na celula
	integer,allocatable:: icell(:),ibac(:),celjus(:),icdom(:) !bacia, celula de jusante, celula do dominio simulado/calibrado
	real,allocatable:: xcel(:),ycel(:) !coordenadas do centro da celula
	real,allocatable:: hcel(:),lcel(:),acel(:),adren(:),srio(:),decl(:) !hmax, hmin, e area da celula, area drenada em celulas, comprimento do rio e declividade do rio
	integer,allocatable:: idren(:) !vetor que ordena os valores de area de drenagem em ordem crescente para a rotina rede
	real,allocatable:: puso(:,:) !proporcao de usos na celula
	real,allocatable:: lambda(:,:),lambdam(:,:),lambdan(:,:),ac(:,:),lambdai(:,:),tanb(:) !indices topograficos e histograma na celula
	real,allocatable:: qref(:) !vazao de referencia
	real qmesp,bc1,bc2,bc3 ! vazao media especifica de toda a bacia e coef de calculo da largura do rio na celula

	real,allocatable:: tontem(:) !temperatura no dia anterior
	real,allocatable:: pre(:) !chuva no intervalo na celula
	!valores em cada celula usados no intervalo de tempo
    real,allocatable:: ta(:),tdew(:),vv(:),roc(:),patm(:)
    !valores das variaveis meteorologicas para todo o intervalo de simulacao
    real,allocatable:: preall(:,:),taall(:,:),tdewall(:,:),vvall(:,:),rocall(:,:),patmall(:,:)    
    
	!uso do solo
    character*25,allocatable:: arqmap(:) ! nome do arquivo de uso da terra
    integer,allocatable:: itmap(:) ! intervalo de tempo onde troca o mapa do uso da terra
    
	!-------------------------------------------------
	
	real,allocatable:: qrg(:,:)	!armazena hidrogramas nas celulas designadas	
	real,allocatable:: qm1(:),qj1(:),qm2(:),qj2(:) !vazoes a montante e a jusante na celula i no tempo 1 e 2
	real,allocatable:: qcel1(:),qcel2(:) !vazoes originadas na celula nos instantes t e t+1 na celula i
	real,allocatable:: pmsub(:),pmssu(:),pmsup(:),pjsub(:),pjssu(:),pjsup(:) !proporcoes de origem das vazoes no rio
	real,allocatable:: vrsub(:),vrssu(:),vrsup(:) !volumes das proporcoes no trecho
	real,allocatable:: ett(:,:) !evapotraspiracao total
	real,allocatable:: etp(:,:) !transpiracao potencial
	real,allocatable:: ec(:,:) ! evaporacao no dossel de agua interceptada
	real,allocatable:: asat(:,:) !percentagem de area saturada
	real,allocatable:: sssu(:,:) !armazenamento na camada superior do solo
	real,allocatable:: qsub(:),qssu(:),qsup(:) !vazao na celula
	real,allocatable:: vsub(:),vssu(:),vsup(:) !volume na celula

	real e0,e0ag ! evaporacao potencial de agua interceptada e de uma superficie de agua livre
	real gamma ! constante psicometrica	em kpa/c (eq. 4.2.28)
	real delta ! declividade da curva de pressao de vapor (eq. 4.2.3)
	real ra ! resistencia aerodicamica da vegetacao e da agua
	real,allocatable:: evq(:) !evaporacao direta das superficies liquidas (bloco agua) da celula em m3/s

	real,allocatable:: qcontorm(:,:),qcontorj(:,:) !vazao condicao de contorno de montante e jusante de cada celula em cada dt de mc 
	real,allocatable:: qrioini(:,:) !vazao da condicao inicial de muskingum cunge em cada subtrecho da celula
	integer,allocatable:: icodmusk(:) !codigo que indica se e linear (0) ou nao linear (1)
	
	real,allocatable:: zfpt(:,:),afpt(:,:),vfpt(:,:) ! tabelas de cota, area alagada e volume para celulas de planicie
	real,allocatable:: zfp(:),afp(:),vfp(:) ! cota da agua na planicie, area alagada e volume para celulas de planicie
	real,allocatable:: zfpb(:),zvert(:) ! cota de fundo da celula de planicie, altura de vertimento
	integer,allocatable:: icellfp(:),celln(:),cellw(:) ! celulas de planicie e celulas adjacentes a uma celula de planicie ao norte e oeste
	real,allocatable:: qrfp(:) ! vazao de troca entre a planicie e o rio
	real cva,bva,cvl,bvl,fch ! coeficiente e largura de vertedor afogado e de vertedor livre, condutancia hidraulica da planicie 
	integer nfp ! numero de celulas que compoem a planicie
	
	real dsup ! escoamento superficial
	real dd   ! drenagem profunda
	real dss  ! drenagem supsuperficial
	real dsub ! escoamento base
	real dtfp ! passo de tempo para o modelo de planicie

    integer:: imesah, fmesah           !mes inicial/final do ano hidrologico
    integer:: ibalanco                 !flag p/ executar ou nao calculo do balanco hidrico

	end module
