!---------------------------------------------------------------------------------------------
!rotina que calcula Muskingum Cunge nao linear
!---------------------------------------------------------------------------------------------
subroutine musk_nl(qjx,ntr,ntc,ntclast,ic)
use vars_main,only: qcontorm,iexut,qcons,evq,dt,qrioini,srio,ntrc,decl,qcontorj, &
                  & ibac,nc
implicit none
integer itc,itr,ic,itclast
real c1,c2,c3 !coeficientes muskingum
real dtc,dtclast1,dtclast2 !variaveis auxiliares
real qjx
integer ntr,ntc,ntclast !numero de trechos, de intervalos de tempo da celula, e da ultima celula
!vazoes internas do trecho	
real,allocatable:: qc(:,:)
real qmusknl,bmusknl,celmusknl ! vazao, largura e celeridade de mc n�o linear
real dtcal,dx,kmc,xmc,den ! delta t e delta x, parametros k e x de mc
real qperdas
allocate (qc(ntr+1,ntc+1))

qc=0.0 !zera matriz

if(ic == iexut(ibac(ic))) then
    qperdas=qcons(ibac(ic)) ! perdas totais devido a irrigacao, cons. humano, etc
else
    qperdas=0.
endif

!condicoes de contorno de montante	
if(ntc == ntclast) then		
    do itc=1,ntc+1
	    qc(1,itc)=max(qcontorm(ic,itc)-evq(ic)-qperdas,0.0)
    enddo
else
    do itc=1,ntc+1
        dtc=float(itc-1)*dt(ic)
        dtclast1=0.
        do itclast=2,ntclast+1
            dtclast2=float(itclast-1)*dt(nc)
            if(dtc <= dtclast2) then
                qc(1,itc)=qcontorm(ic,itclast-1)+(dtc-dtclast1)*(qcontorm(ic,itclast)-qcontorm(ic,itclast-1))/(dtclast2-dtclast1)
                qc(1,itc)=max(qc(1,itc)-evq(ic)-qperdas,0.0)
                exit
            endif
            dtclast1=dtclast2
        enddo
    enddo
endif

!condicoes iniciais	
do itr=1,ntr+1
	qc(itr,1)=qrioini(ic,itr)
enddo

do itc=1,ntc
	do itr=1,ntr
		!aqui est� a diferenca entre linear e nao linear
		dtcal=dt(ic)
		dx=srio(ic)*1000./ntrc(ic)
		qmusknl=qc(itr+1,itc) !considera v�lida a �ltima vazao calculada
		call interpmusk(qmusknl,celmusknl,bmusknl,ic)
		kmc=dx/celmusknl
		xmc=0.5*(1.0-(qmusknl/(bmusknl*decl(ic)*celmusknl*dx)))
		den=2.*kmc*(1.-xmc)+dtcal
		c1=(2.*kmc*xmc+dtcal)/den
		c2=(dtcal-2.*kmc*xmc)/den
		c3=(2.*kmc*(1.-xmc)-dtcal)/den
		qc(itr+1,itc+1)=c1*qc(itr,itc)+c2*qc(itr,itc+1)+c3*qc(itr+1,itc)
		qc(itr+1,itc+1)=max(qc(itr+1,itc+1),0.0) !evita vaz�es negativas
	enddo
	!guarda hidrograma de sa�da completo do ultimo subtrecho em cada intervalo de tempo de mc 
	if(ntc == ntclast) then		
		qcontorj(ic,itc+1)=qc(ntr+1,itc+1)
	else
        dtc=float(itc)*dt(ic)
        dtclast1=0.
        do itclast=2,ntclast+1
            dtclast2=float(itclast-1)*dt(nc)
            if(dtc <= dtclast2) then
                qcontorj(ic,itclast)=qc(ntr+1,itc+1)+(dtc-dtclast1)*(qc(ntr+1,itc+1)-qc(ntr+1,itc))/(dtclast2-dtclast1)
                exit
            endif
            dtclast1=dtclast2
        enddo
    endif
enddo

!guarda valores da condi��o inicial para a proxima vez que vai calcular Muskingum Cunge
qrioini(ic,1)=qcontorm(ic,ntc+1)
do itr=2,ntr+1
	qrioini(ic,itr)=qc(itr,ntc+1)
enddo

!guarda valor do final do per�odo de propaga��o no extremo de jusante do trecho
qjx=qc(ntr+1,ntc+1) !este valor volta para programa principal

deallocate (qc)
return
end