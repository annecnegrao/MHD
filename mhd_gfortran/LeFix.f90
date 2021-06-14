!---------------------------------------------------------------------------------------------
!subrotina que leh os parametros gerais fixos, como o numero de celulas, o numero de usos do
!solo, o numero de intervalos de tempo e de postos pluviometricos...
!ParFix.hig
!---------------------------------------------------------------------------------------------
subroutine lefix
use vars_main, only: filfix,dir_dados,idia,imes,iano,ihoraini,dtp,nt,nc,nu,nb,icalib,ntprev,imesah,fmesah, &
                   & nbvfo,nobs,arqobs,iqobs,nhidg,ihidg,nmapas,arqmap,itmap,numsubst,iniaju,iiano,ibalanco, &
                   & arqsubst,isubst,idiafimArquivo,imesfimArquivo,ianofimArquivo,ihorafimArquivo,iihoraini, &
                   & idiaArquivo,imesArquivo,ianoArquivo,ihoraArquivo,idiafin,imesfin,ianofin,ihorafin	
implicit none
integer k,idiapre,imespre,ianopre,ihorapre,jul,julday
character titulo*100

!formatacao
71	format(6i10)
75	format(a100)

open(filfix,file=dir_dados//'entrada/ParFix.hig',status='old')
read(filfix,75) titulo
write(*,*) titulo
read(filfix,*)
read(filfix,*)
read(filfix,*)
read(filfix,75) titulo
write(*,*) titulo

read(filfix,71)idiaArquivo,imesArquivo,ianoArquivo,ihoraArquivo !data e hora do inicio do arquivo de dados meteorologicos
write(*,71)idiaArquivo,imesArquivo,ianoArquivo,ihoraArquivo
read(filfix,*)
read(filfix,75) titulo
write(*,*) titulo	

read(filfix,71)idiafimArquivo,imesfimArquivo,ianofimArquivo,ihorafimArquivo !data e hora do fim do arquivo de dados meteorologicos
write(*,71)idiafimArquivo,imesfimArquivo,ianofimArquivo,ihorafimArquivo
read(filfix,*)
read(filfix,75) titulo
write(*,*) titulo	

read(filfix,71)idia,imes,iano,ihoraini !data e hora do inicio da simulacao
write(*,71)idia,imes,iano,ihoraini
read(filfix,*)
read(filfix,75) titulo
write(*,*) titulo	
read(filfix,71)idiafin,imesfin,ianofin,ihorafin !data e hora do fim da simulacao
write(*,71)idiafin,imesfin,ianofin,ihorafin	
read(filfix,*)

iihoraini=ihoraini
iiano=iano
read(filfix,75) titulo
write(*,*) titulo	
read(filfix,71)iniaju ! spin-up: inicia a comparacao observado x calculado dois ano apos
write(*,*)iniaju
read(filfix,*)

read(filfix,75) titulo	
write(*,*) titulo
read(filfix,*) dtp !intervalo de tempo em segundos
write(*,*)dtp
! calcula o numero de intervalos de tempo
jul=(julday(imesfin,idiafin,ianofin)-julday(imes,idia,iano))*86400+int(dtp)+(ihorafin-ihoraini)*3600
nt=jul/int(dtp)
read(filfix,*)

!verificacoes das datas
if( (julday(imes,idia,iano)*86400+int(dtp)+ihoraini*3600) < &
&   (julday(imesArquivo,idiaArquivo,ianoArquivo)*86400+int(dtp)+ihoraArquivo*3600) ) &
&  stop '**Data Inicial da Simulacao incompativel com arquivos prec e met.bin'

if( (julday(imesfin,idiafin,ianofin)*86400+int(dtp)+ihorafin*3600) > &
&   (julday(imesfimArquivo,idiafimArquivo,ianofimArquivo)*86400+int(dtp)+ihorafimArquivo*3600) ) &
&  stop '**Data Final da Simulacao incompativel com arquivos prec e met.bin'

read(filfix,75) titulo 
write(*,*) titulo
read(filfix,71)nc,nu,nb !numero de celulas, usos, bacias  
write(*,71)nc,nu,nb
read(filfix,*)

read(filfix,75) titulo
write(*,*) titulo
read(filfix,71) icalib !opcao de execucao
write(*,71) icalib
if(icalib == 2) then ! faz previsao
    read(filfix,*)
	read(filfix,75) titulo
    write(*,*) titulo
    read(filfix,71)idiapre,imespre,ianopre,ihorapre !data e hora do fim da previsao
    write(*,71)idiapre,imespre,ianopre,ihorapre
    ! calcula alcance da previsao em intervalos de tempo
    jul=(julday(imespre,idiapre,ianopre)-julday(imesfin,idiafin,ianofin))*86400+(ihorapre-ihorafin)*3600
    ntprev=jul/int(dtp) !numero de intervalos de previsao
    if(ntprev < 0) then
        write(*,*) '*** data de fim de simulacao ou de previsao errada(s)***'
        stop
    endif
else
    ntprev=0
    read(filfix,*)
    read(filfix,*)
    read(filfix,*)
endif
read(filfix,*)

read(filfix,75)titulo
write(*,*)titulo
read(filfix,71)nbvfo
write(*,71)nbvfo
read(filfix,*)

read(filfix,75)titulo
write(*,*)titulo
read(filfix,*)nobs,arqobs !numero de postos com dados observados e nome do arquivo dos dados
write(*,*)nobs,arqobs
read(filfix,*)

allocate(iqobs(nobs))

read(filfix,75)titulo
write(*,*)titulo
read(filfix,*)(iqobs(k),k=1,nobs) !número da celula em que existem dados observados
write(*,*)(iqobs(k),k=1,nobs)
read(filfix,*)

read(filfix,75)titulo
write(*,*)titulo
read(filfix,*)nhidg	!numero de pontos em que se deseja gravar o hidrograma
write(*,*)nhidg
read(filfix,*)

allocate(ihidg(nhidg))

read(filfix,75)titulo
write(*,*)titulo
do k=1,nhidg
	read(filfix,*)ihidg(k) ! celulas correspondentes aos pontos em que se deseja gravar os hidrogramas
	write(*,*)ihidg(k)
enddo
read(filfix,*)

read(filfix,75)titulo
write(*,*)titulo
read(filfix,*) nmapas !numero de mapas de uso da terra
write(*,*) nmapas
read(filfix,*)

allocate(arqmap(nmapas),itmap(nmapas))

read(filfix,75) titulo
write(*,*)titulo
read(filfix,75) titulo
write(*,*)titulo
do k=1,nmapas
	read(filfix,'(2(i2,1x),i4,1x,a25)') idiafin,imesfin,ianofin,arqmap(k) !le data do mapa de uso da terra e seu nome
	! determina o numero do intervalo de tempo onde muda o mapa de uso da terra
	jul = (julday(imesfin,idiafin,ianofin) - julday(imes,idia,iano))*86400 + int(dtp)
	! se a data inicial da simulacao ou previsao for posterior a ultima data do uso da terra, itmap<0 e vale o ultimo mapa lido
	itmap(k) = jul/int(dtp) 
	write(*,'(2(i2,a),i4,2x,a25)') idiafin,'/',imesfin,'/',ianofin,arqmap(k)
enddo
read(filfix,*)

read(filfix,75)titulo
write(*,*)titulo
read(filfix,*) numsubst,arqsubst !numero de pontos que se deseja subsituir a vazao calculada pela lida
write(*,*) numsubst,arqsubst
read(filfix,*)

allocate(isubst(numsubst))

read(filfix,75)titulo
write(*,*)titulo
if(numsubst > 0) then
    write(*,*)'numsubst',numsubst
    do k=1,numsubst
        read(filfix,*)isubst(k) !celulas cuja vazao calculada sera substituida pela lida
        write(*,*)isubst(k)
    enddo
else
    read(filfix,*)
endif

read(filfix,*)
read(filfix,75)titulo
write(*,*)trim(titulo)
read(filfix,*)ibalanco          !calcular balanco hidrico(opcoes: 1=sim, 0=nao)
write(*,*)ibalanco
read(filfix,75)titulo
write(*,*)trim(titulo)
read(filfix,*)imesah,fmesah
write(*,*)imesah,fmesah
close (filfix)

write(*,*)
write(*,*) ' leu parametros gerais fixos'

return
end