!---------------------------------------------------------------------------------------------
!subrotina interpmusk faz uma interpolação rápida para muskingun cunge nao linear
!interpola a largura da secao (bmusknl) e a celeridade (celmusknl)
!---------------------------------------------------------------------------------------------
subroutine interpmusk(qmusknl,celmusknl,bmusknl,ic)
use vars_main, only: nc
implicit none
integer,parameter:: ntmup=20 !numero de pontos da tabela do muskingun-cunge nao linear
integer kh,ic
real qmusknl,celmusknl,bmusknl
real qmup(nc,ntmup),amup(nc,ntmup),bmup(nc,ntmup),cmup(nc,ntmup) !tabela muskingun nao linear

do kh=2,ntmup-1 !percorre tabela
	if(qmusknl<qmup(ic,kh)) then
	 write(*,*)'qmusknl<qmup rutina interpmusk', ic, kh
	 exit
	endif
enddo
!sempre assume que o ponto está entre kh-1 e kh (mesmo nos extremos)
celmusknl=cmup(ic,kh-1)+(cmup(ic,kh)-cmup(ic,kh-1))*(qmusknl-qmup(ic,kh-1))/(qmup(ic,kh)-qmup(ic,kh-1))
bmusknl=bmup(ic,kh-1)+(bmup(ic,kh)-bmup(ic,kh-1))*(qmusknl-qmup(ic,kh-1))/(qmup(ic,kh)-qmup(ic,kh-1))

return
end