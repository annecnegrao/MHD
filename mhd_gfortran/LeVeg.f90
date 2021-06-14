!---------------------------------------------------------------------------------------------
!subrotina que leh dados de vegetacao
!sao dados mensais de 16 tipos de vegetacao
!---------------------------------------------------------------------------------------------
subroutine leveg
use vars_main, only: filbloc,filveg,nu,dir_dados,omegac,fr,alb,iaf,zveg,rc,cover,prad,dveg,z0,scmax
use vars_calib, only: solo,veg
implicit none 
real,allocatable:: aux(:,:) ! variavel auxiliar de leitura
integer iu,im,iveg,nveg,iaux

allocate(veg(nu+1),solo(nu+1))

! leh blocos da bacia
open(filbloc,file=dir_dados//'entrada/Blocos.hig',status='old')
read(filbloc,*)
read(filbloc,*)
read(filbloc,*)(solo(iu),veg(iu),iu=1,nu)
close(filbloc)

!leh parametros de vegetacao
open(filveg,file=dir_dados//'entrada/ParVeg.hig',status='old')
read(filveg,*)
read(filveg,*)
read(filveg,*)omegac ! fator de compensacao de jarvis
read(filveg,*)fr ! fator de distribuicao de raizes
read(filveg,*)nveg ! numero de tipos funcionais de vegetacao

allocate(aux(nveg,12))

veg(nu+1)=nveg ! reserva o uso nu + 1 para evaporacao de area alagada

! compatibiliza o tipo de vegetacao com o tipo de solo para agua
do iu=1,nu+1
    if(veg(iu)==nveg) solo(iu)=14 ! solo correspondente a agua
end do

!albedo
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    alb(iu,im)=aux(veg(iu),im)
	end do
end do

!indice de area foliar
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    iaf(iu,im)=aux(veg(iu),im)
	end do
end do

!zveg (altura media da vegetacao)
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    zveg(iu,im)=aux(veg(iu),im)
	end do
end do

!percentagem de cobertura da vegetacao - cover 
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    cover(iu,im)=aux(veg(iu),im)
	end do
end do

!prad (profundidade radicular)
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    prad(iu,im)=aux(veg(iu),im)
	end do
end do    

!dveg - plano de deslocamento zero (m)
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    dveg(iu,im)=aux(veg(iu),im)
	end do
end do

!rugosidade z0 (m) 
read(filveg,*)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,(aux(iveg,im),im=1,12)
end do
do iu=1,nu+1
    do im=1,12
	    z0(iu,im)=aux(veg(iu),im)
	end do
end do   

!rc (resistencia superficial minima da vegetacao)
read(filveg,*)
do iveg=1,nveg
	read(filveg,'(a3,g8.0)')iaux,aux(iveg,1)
end do
do iu=1,nu+1
	rc(iu)=aux(veg(iu),1)
end do

!scmax (capacidade maxima do dossel)
read(filveg,*)
do iveg=1,nveg
	read(filveg,*)iaux,aux(iveg,1)
end do
do iu=1,nu+1
    scmax(iu)=aux(veg(iu),1)
end do

close(filveg)

deallocate(aux)

write(*,*)
write(*,*) ' leu parametros de vegetacao'

return
end