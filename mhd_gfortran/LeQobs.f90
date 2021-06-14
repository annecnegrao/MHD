!---------------------------------------------------------------------------------------------
!subrotina que leh os dados observados de vazao
!---------------------------------------------------------------------------------------------
subroutine leqobs
use vars_main, only: filobs,dir_dados,arqobs,titulo,nstep,qobs,imes,idia,iano,ihoraini,dtp, &
                   & it,nt,nobs
implicit none
integer k,julday,idia1,imes1,iano1,ihora1,iaux
character*100 fmt !fromatacao

71	format(a80)

open(filobs,file=dir_dados//'entrada/'//arqobs,status='old')
read(filobs,*)
! le primeiro registro e extrai a data, que deve estar em formato dd/mm/aaaahh
read(filobs,71) titulo
iaux=index(titulo,'/') ! identifica a primeira barra da data
read(titulo(iaux-2:iaux-1),'(i2)') idia1 
read(titulo(iaux+1:iaux+2),'(i2)') imes1
read(titulo(iaux+4:iaux+7),'(i4)') iano1
if(titulo(iaux+8:iaux+9) == '  ') then ! data nao inclui hora
    ihora1=0
    iaux=iaux+7
else ! data inclui hora
    read(titulo(iaux+8:iaux+9),'(i2)') ihora1
    iaux=iaux+9
endif
backspace(filobs)


! calcula o numero de intervalos de tempo ate o registro correspondente ao primeiro dia da simulacao
nstep=(julday(imes,idia,iano)-julday(imes1,idia1,iano1))*86400+(ihoraini-ihora1)*3600
nstep=nstep/int(dtp)
    
if(nstep < 0) then ! os dados meteorologicos se iniciam apos o periodo de simulacao
    write(*,*) ' ***inicio da simulacao anterior ao inicio dos dados disponiveis*** '
	stop 
else ! se nstep for maior que zero despreza os primeiros nstep valores
    do it=1,nstep
        read(filobs,*)
    enddo
endif

! leitura do restante do arquivo
write(fmt,'(a,i2,a,i2.2,a)')"(a",iaux,",1x,",nobs,"g10.0)"
do it=1,nt
	read(filobs,fmt) titulo(1:iaux),(qobs(k,it),k=1,nobs)
enddo
write(*,*)
write(fmt,'(a,i2,a,i2.2,a)')"(a20,/,a",iaux,",",nobs,"f10.3)"
write(*,fmt) ' ultima vazao lida  ',titulo(1:iaux),(qobs(k,nt),k=1,nobs)
write(*,*)
close (filobs)
!fim da leitura de dados de vazao observada --------------------!

return
end
