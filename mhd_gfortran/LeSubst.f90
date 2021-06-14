!---------------------------------------------------------------------------------------------
!subrotina que leh os dados de vazao que substituem a vazao calculada
!---------------------------------------------------------------------------------------------
subroutine lesubst
use vars_main,only: dir_dados,filsubs,arqsubst,titulo,nstep,it,nt,qobs,numsubst
implicit none
integer k,idia1,imes1,iano1,ihora1,iaux
character*100 fmt !fromatacao

71	format(a80)

open(filsubs,file=dir_dados//'entrada\'//arqsubst,status='old')
read(filsubs,*)
! le primeiro registro e extrai a data, que deve estar em formato dd/mm/aaaahh
read(filsubs,71) titulo
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
backspace(filsubs)

! se a simulacao se inicia apos a data do inicio das observacoes, o nstep, calculado em leqobs.prn, utilizado eh > 0
if(nstep > 0) then ! a simulacao se inicia apos o inicio dos dados
    ! despreza os primeiros nstep valores
    do it=1,nstep
        read(filsubs,*)
    enddo
endif

! le o restante do arquivo
write(fmt,'(a,i2,a,i2.2,a)')"(a",iaux,",1x,",numsubst,"g10.0)"
do it=1,nt
	read(filsubs,fmt) titulo(1:iaux),(qobs(k,it),k=1,numsubst)
enddo
write(*,*)
write(fmt,'(a,i2,a,i2.2,a)')"(a20,/,a",iaux,",",numsubst,"f10.3)"
write(*,fmt) ' ultima vazao lida  ',titulo(1:iaux),(qobs(k,nt),k=1,numsubst)
write(*,*)
close (filsubs)

write(*,*)
write(*,*)'fim da leitura de dados de vazao para a substituicao de valores calculados'

return
end