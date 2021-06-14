!---------------------------------------------------------------------------------------------
!subrotina que faz a leitura dos dados das celulas e parametros topograficos
!tambem analisa as planicies de alagamento caso existam
!cell.hig
!partop.hig
!floodplain.hig
!---------------------------------------------------------------------------------------------
subroutine lecell
use vars_main, only: dir_dados,filhig,filtop,filfp,nfp,icell,nc,xcel,ycel,ibac,adren,acel, &
                   & hcel,lcel,srio,decl,celjus,tanb,ac,lambda,lambdam,lambdan,lambdai, &
                   & iexut,kcb,nb,nu,smax,mu,neta,zfpt,afpt,vfpt,zfp,zvert,zfpb,afp,vfp, &
                   & icellfp,celln,cellw,qrfp,cvl,bvl,cva,bva,fch
implicit none
integer flagfp(nc) !indicador 0 ou 1 das celulas que fazem parte do floodplain
integer ibant,ic,k,ib,iu,kb,ifp
real res,xc,yc,dist
logical existe

open(filhig,file=dir_dados//'entrada/cell.hig',status='old')
open(filtop,file=dir_dados//'entrada/ParTop.hig',status='old')
ibant=1
ib=1
nfp=0
icell=0
do ic=1,nc
	read(filhig,*) k,xcel(ic),ycel(ic),ibac(ic),adren(ic),acel(ic),hcel(ic),lcel(ic),srio(ic),decl(ic),celjus(ic),flagfp(ic)
    nfp=nfp+flagfp(ic)
    !leitura da distribuicao do indice topográfico
	read(filtop,'(45x,g10.0,50g10.0,50g20.0)') tanb(ic),(ac(ic,k),k=1,50),(lambdai(ic,k),k=1,50)
	!define exutorios (verifica qual celula define o fim de uma sub-bacia)
	if(ibac(ic) /= ibant)then
		iexut(ib)=ic-1
		ibant=ibac(ic)
		ib=ib+1
	endif
enddo

iexut(ib)=nc !o exutorio da ultima sub-bacia é a ultima celula
celjus(nc)=nc+1 !a celula a jusante da ultima subbacia corresponde a celula nc+1

!verifica o numero de celulas em cada sub-bacia
kcb(1)=iexut(1)
do kb=2,nb
	kcb(kb)=iexut(kb)-iexut(kb-1)
enddo

!calcula os valores medios de lambda e lambdan
lambdam=0.
lambdan=0.
lambda=0.
do ic=1,nc
	do iu=1,nu
		if(smax(ibac(ic),iu).gt.0) then ! pula agua
			do k=2,50
				lambda(ic,iu)=lambda(ic,iu)+0.5*(lambdai(ic,k)+lambdai(ic,k-1))*(ac(ic,k)-ac(ic,k-1))
				lambdam(ic,iu)=lambdam(ic,iu)+0.5*(lambdai(ic,k)**(1./mu(ibac(ic),iu))+ &
                              &lambdai(ic,k-1)**(1./mu(ibac(ic),iu)))*(ac(ic,k)-ac(ic,k-1))
				lambdan(ic,iu)=lambdan(ic,iu)+0.5*(lambdai(ic,k)**(1./neta(iu))+ &
                              &lambdai(ic,k-1)**(1./neta(iu)))*(ac(ic,k)-ac(ic,k-1))
			enddo
			lambdam(ic,iu)=lambdam(ic,iu)**mu(ibac(ic),iu)				
			lambdan(ic,iu)=lambdan(ic,iu)**neta(iu)
		endif
	enddo
enddo
 
if(nfp > 0) then ! existem celulas com planicie de alagamento
    allocate (zfpt(nfp,50),afpt(nfp,50),vfpt(nfp,50),zfp(nfp),zvert(nfp),zfpb(nfp))
    allocate (afp(nfp),vfp(nfp),icellfp(nfp),celln(nfp),cellw(nfp),qrfp(nfp))
    inquire(file=dir_dados//'entrada/floodplain.hig',exist=existe)
	if(.not. existe)then ! verifica existencia de arquivo com dados da planicie alagada   
        stop 'nao existe arquivo com dados topograficos da planicie'
    else
        open(filfp,file=dir_dados//'entrada/floodplain.hig',status='old')
    endif
    read(filfp,*)
    read(filfp,*)
    read(filfp,'(2g10.0)') cvl,bvl ! coeficiente e largura de vertedor livre
    read(filfp,'(2g10.0)') cva,bva ! coeficiente e largura de vertedor afogado
    read(filfp,'(2g10.0)') fch ! condutancia hidraulica da planicie de alagamento
    read(filfp,*)
    
    ! le dados da altura, da fracao da celula alagada e do volume alagado
    do ifp=1,nfp
        ! icellfp indica quais celulas do arquivo cell.hig pertencem a planicie de inundacao
        ! icellfp(ifp) = x indica que a celula ifp do arquivo floodplain.hig corresponde a x-esima celula do arquivo cell.hig
        read(filfp,'(i10,2g10.0)') icellfp(ifp),zvert(ifp),zfpb(ifp)
        read(filfp,'(10x,50g10.0)') (zfpt(ifp,k),k=1,50)
        read(filfp,'(10x,50g10.0)') (afpt(ifp,k),k=1,50)
        read(filfp,'(10x,50g10.0)') (vfpt(ifp,k),k=1,50)
        ! icell indica a posicao das celulas da planicie do arquivo cell.hig no arquivo floodplain.hig
        ! se icell(ic) = 0 a celula ic de cell.hig nao pertence a planicie
        ! se icell(ic) = x > 0, a celula ic pertence a planicie e corresponde a x-esima celula do arquivo floodplain.hig
        icell(icellfp(ifp))=ifp
    enddo
    close(filfp)

    ! determina a conectividadade da planicie usando as coordenadas da celula
    celln=0 !celula da planicie localizada ao norte da celula analisada
    cellw=0 !celula da planicie localizada ao oeste da celula analisada

    ! determina a resolucao do modelo procurando as minimas distancias entre celulas contiguas no meio do dominio
    res=9999.
    do ic=1,10
        dist=abs(xcel(nc/2+ic)-xcel(nc/2+2*ic))
        if(dist /= 0) res=min(res,dist)
    enddo

    do ifp=1,nfp
        ! testa ao norte
        xc=xcel(icellfp(ifp))
        yc=ycel(icellfp(ifp))+res
        do k=1,nfp
            dist=(xc-xcel(icellfp(k)))**2.+(yc-ycel(icellfp(k)))**2.
            if(dist < 0.01) then
                celln(ifp)=k
                exit
            endif
        enddo
        ! testa ao oeste
        xc=xcel(icellfp(ifp))-res
        yc=ycel(icellfp(ifp))
        do k=1,nfp
            dist=(xc-xcel(icellfp(k)))**2.+(yc-ycel(icellfp(k)))**2.
            if(dist < 0.01) then
                cellw(ifp)=k
                exit
            endif
        enddo
    enddo
endif

close (filhig)
close (filtop)

write(*,*)
write(*,*) ' leu dados das celulas e topograficos'

return
end