!---------------------------------------------------------------------------------------------
!subrotina que leh dados meteorologicos dos arquivos(real*4) de dados meteorologicos, sao eles:
!prec.bin
!met.bin
!---------------------------------------------------------------------------------------------
subroutine lemet
use vars_main, only: filplu,filmet,nc,nstep,preall,taall,tdewall,vvall,patmall,rocall, &
                   & dir_dados,nt,pre,ta,tdew,vv,patm,roc,it
implicit none
integer ic,itaux,irec
real lixo

! abre arquivos de precipitacao e meteorologicos de todo o periodo de dados
! o inicio da simulacao pode ser posterior ao inicio do periodo com dados
! eh necessario que a data de incio das observacoes nos arquivos qobs.prn, prec.bin e met.bin sejam as mesmas!!!
open(filplu,file=dir_dados//'entrada/prec.bin', &
                  & status='old',form='unformatted',access='direct',recl=4*nc)  !dados de chuva interpolados	
open(filmet,file=dir_dados//'entrada/met.bin', &
     & status='old',form='unformatted',access='direct',recl=4*nc) !dados dados meteorologicos para calcular evaporacao

irec = 1
! se a simulacao se inicia apos a data do inicio das observacoes, o nstep, calculado em leqobs.prn, utilizado eh > 0 
if(nstep > 0) then ! a simulacao se inicia apos o inicio dos dados
    ! despreza os nstep primeiros dados                  
    do itaux=1,nstep
        read(filplu,rec=itaux)(lixo,ic=1,nc) ! precipitacao
        
        read(filmet,rec=irec)(lixo,ic=1,nc) ! temperatura do ar
        irec = irec+1
        read(filmet,rec=irec)(lixo,ic=1,nc) ! temperatura ponto de orvalho
        irec = irec+1
        read(filmet,rec=irec)(lixo,ic=1,nc) ! velocidade do vento
        irec = irec+1
        read(filmet,rec=irec)(lixo,ic=1,nc) ! pressao atmosferica
        irec = irec+1
        read(filmet,rec=irec)(lixo,ic=1,nc) ! radiacao global incidente
        irec = irec+1
    enddo
endif

! le o restante do arquivo
do itaux=1,nt
    read(filplu,rec=itaux)(pre(ic),ic=1,nc) ! precipitacao
    
    read(filmet,rec=irec)(ta(ic),ic=1,nc) ! temperatura do ar
    irec = irec+1
    read(filmet,rec=irec)(tdew(ic),ic=1,nc) ! temperatura ponto de orvalho
    irec = irec+1
    read(filmet,rec=irec)(vv(ic),ic=1,nc) ! velocidade do vento
    irec = irec+1
    read(filmet,rec=irec)(patm(ic),ic=1,nc) ! pressao atmosferica
    irec = irec+1
    read(filmet,rec=irec)(roc(ic),ic=1,nc) ! radiacao global incidente
    irec = irec+1
    
    it = itaux
    
    call qqmet !verifica a qualidade dos dados meteorologicos
    
    ! guardando nas variaveis all
    preall(:,itaux) = pre
    taall(:,itaux) = ta
    tdewall(:,itaux) = tdew
    vvall(:,itaux) = vv
    patmall(:,itaux) = patm
    rocall(:,itaux) = roc
enddo
close (filmet)
close(filplu)

write(*,*)
write(*,*) ' leu dados meteorologicos'

return
end
