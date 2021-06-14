!---------------------------------------------------------------------------------------------
! Subrotina que calcula a simulacao  
! [nao calibra nem faz previsao, opcao "0" do ParFix.hig]  
!---------------------------------------------------------------------------------------------

subroutine simula
use vars_main,only: icdom,ncdom,nc,adren,filpmedia,dir_dados,filsol,prec,preb,ecb,ettb,&
&sssub,sradb,ssubb,dsupb,dssub,etpb,dsubb,qbsup,qbssu,qbsub,asatb,it,cdatag,nb,filbac,filprp,&
&filhid,qrg,nhidg,idren,dtp,ibac,nbvfo,nt

implicit none
integer ic,ib
integer kbanual(nb)
real prebanual(nb)
character subbac*2,cdatal(nt)*16
character*100 fmt !fromatacao
	
! dominio da simulacao
ncdom=nc
do ic=1,nc
    icdom(ic)=ic ! por enquanto, simula todo o dominio na rotina simulacao 
enddo
!vetor com a sequencia de celulas com areas drenada crescente usados na resolucao na rotina rede
call sort(nc,adren,idren)

! arquivo de chuva diaria da bacia
open(filpmedia,file=dir_dados//'saida/chuva_diaria_media.txt',status='unknown') 

!************** nosolo.hig *************************************************************** 
!arquivo de saida com dados do solo de uma celula       
open(filsol,file=dir_dados//'saida/nosolo.hig',status='unknown')
! titulos dos arquivos
write(filsol,'(2a)') 'celula    precip       ett        ec       etp      asat      ', &
                   & 'sssu      srad      ssub      qsub      qssu      qsup'
	
!zerando variaveis do balanco hidrico da bacia
prec=0.0 ! precipitacao acumulada na celula
preb=0.0 !precipitacao na bacia
ecb=0.0 ! lamina interceptada na bacia
ettb=0.0 ! potranspiracao na bacia
etpb=0.0 ! transpiracao potencial na bacia
sssub=0.0 ! lamina armazenada media camada superior na bacia
sradb =0.0 ! lamina armazenada media na camada intermediaria da bacia
ssubb=0.0 !lamina armazenada media camada infeiror na bacia
dsupb=0.0 !escomento superficial da bacia
dssub=0.0 !escomento intermediario da bacia
dsubb=0.0 !escomento subterraneo da bacia
qbsup=0.0 !vazao superficial
qbssu=0.0 !vazao sub-superficial
qbsub=0.0 !vazao de base
asatb=0.0 !area saturada media na bacia
	
call modelo !chama a rotina que comanda o loop do tempo

! cria as datas longas para o excel
do it=1,nt
    cdatal(it)=cdatag(it)(1:2)//'/'//cdatag(it)(3:4)//'/'//cdatag(it)(5:8)//' '//cdatag(it)(9:10)//':00'
enddo

! abre arquivos com resultados da simulacao
do ib=1,nb
    write(subbac,'(i2.2)') ib
    !arquivo de saida com dados das bacias
    open(filbac+ib-1,file=dir_dados//'saida/bacia'//subbac//'.hig',status='unknown') 
    write(filbac+ib-1,'(a23,i5)') 'dados medios da bacia: ',ib
    write(filbac+ib-1,'(2a)') '            data        sssu        srad        ssub        asat        prec', &
                              & '          ec         etp         ett        qsup        qssu        qsub'
    write(filbac+ib-1,'(2a)') '                          mm          mm          mm           %          mm', &
                              & '          mm          mm          mm          mm          mm          mm'

enddo
open(filprp,file=dir_dados//'saida/qprop.hig',status='unknown') !arquivo de saida com vazoes segundo a origem
open(filhid,file=dir_dados//'saida/vazao.hig',status='unknown') !arquivo de saida com hidrogramas calculados

!grava hidrogramas nos nb exutorios desejados
write(fmt,'(a,i2.2,a)')"(a16,",nhidg,"f10.3)"
write(filhid,fmt) (cdatal(it),(qrg(ib,it),ib=1,nhidg),it=1,nt)
!grava valores medios das variaveis das sub-bacias
do ib=1,nb
    write(filbac+ib-1,'(a16,11f12.3)')  (cdatal(it),sssub(ib,it),sradb(ib,it),ssubb(ib,it),asatb(ib,it), &
             & preb(ib,it),ecb(ib,it),etpb(ib,it),ettb(ib,it),dsupb(ib,it),dssub(ib,it),dsubb(ib,it),it=1,nt)
    close(filbac+ib-1)
enddo
!grava arquivo com as vazoes segundo a origem
write(filprp,'(a16,3f10.3)') (cdatal(it),qbsup(it),qbssu(it),qbsub(it),it=1,nt)
close (filhid)
close (filprp)

call fobj !calcula os valores das funcoes objetivo
	
!calcula precipitacao media anual por bacia
prebanual=0.0
kbanual=0
prec=(prec/nt)*(365.*24.*3600./dtp)

do ic=1,nc
	ib=ibac(ic)
	prebanual(ib)=prebanual(ib)+prec(ic)
	kbanual(ib)=kbanual(ib)+1.
enddo
do ib=1,nb
	prebanual(ib)=prebanual(ib)/kbanual(ib)
enddo

!fim do calculo da precipitacao media anual por bacia
write(*,*)'   chuva anual na bacia ',nbvfo,prebanual(nbvfo)

close (filpmedia) ! arquivo de chuva da bacia
return
end