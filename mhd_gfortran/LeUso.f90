!---------------------------------------------------------------------------------------------
!subrotina que leh dados de uso do solo
!para cada celula eh feita a leitura da porcentagem de cada uso
!---------------------------------------------------------------------------------------------
subroutine leuso(imapas)
use vars_main, only: filuso, dir_dados, arqmap, nc, puso, nu
implicit none
integer ic,iu,imapas,iaux
real aux

open(filuso,file=dir_dados//'entrada/'//arqmap(imapas),status='old')	
do ic =1,nc
	read(filuso,*) iaux,aux,aux,iaux,(puso(ic,iu),iu=1,nu)
enddo
close(filuso)

!write(*,*)
!write(*,*) ' leu mapa(s) de usos da terra'

return
end