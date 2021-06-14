!---------------------------------------------------------------------------------------------
!subrotina que le os arquivos de previsão de chuva para cada célula do modelo.
!---------------------------------------------------------------------------------------------
subroutine lemetprev
use vars_main, only: filpluprev,filmetprev,dir_dados,cdatag,pre,ta,tdew,vv,patm,roc,nc,it,nt
implicit none
integer ic,irec

! abre arquivo de dados de precipitacao e meteorologicos previstos
! formato dos arquivos: preprev_ddmmyyyyhh.bin e metddmmyyyyhh.bin
open(filpluprev,file=dir_dados//'previsao/'//cdatag(nt+1)//'/preprev_'//cdatag(it)//'.bin', &
   & status='old',form='unformatted',access='direct',recl=4*nc)
open(filmetprev,file=dir_dados//'previsao/'//cdatag(nt+1)//'/metprev_'//cdatag(it)//'.bin', &
   & status='old',form='unformatted',access='direct',recl=4*nc)

irec = 1

read(filpluprev,rec=irec)(pre(ic),ic=1,nc) ! precipitacao

read(filmetprev,rec=irec)(ta(ic),ic=1,nc) ! temperatura do ar
irec = irec+1
read(filmetprev,rec=irec)(tdew(ic),ic=1,nc) ! temperatura ponto de orvalho
irec = irec+1
read(filmetprev,rec=irec)(vv(ic),ic=1,nc) ! velocidade do vento
irec = irec+1
read(filmetprev,rec=irec)(patm(ic),ic=1,nc) ! pressao atmosferica
irec = irec+1
read(filmetprev,rec=irec)(roc(ic),ic=1,nc) ! radiacao global incidente
    
call qqmet !verifica a qualidade dos dados meteorologicos

close(filpluprev)
close(filmetprev)	

return
end
