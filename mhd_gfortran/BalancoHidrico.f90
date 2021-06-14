subroutine balanco
use vars_main, only: nb,nc,it,nt,iiano,jdia,ihoraini,dtp,dir_dados,preall,icdom,ncdom,qr,iihoraini, &
& idini,qobs,adren,iexut,ettb,imesah,fmesah,ibac,kcb,preb,taall,tdewall,vvall,patmall,rocall
implicit none
real:: tempo,qobsmedia,qsimmedia,ettsimmedia,areadrenagem,erro,precmedia
real:: tamedia,tdewmedia,vvmedia,patmmedia,rocmedia
integer, dimension(1:12):: diasmes=(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer:: ib, ic, ic2, itiniah, itfimah, julday, ianoaux, ibacaux,ibflag,ntaux,ndiasobs,ihora,incrementoano
integer, parameter:: filbalanco=666

open(filbalanco, file=dir_dados//'saida/balanco_hidrico.hig', status='unknown')

!loop sobre as bacias
do ib=1, nb
    ianoaux = iiano
    ihora=iihoraini
    tamedia = 0.
    tdewmedia = 0.
    vvmedia = 0.
    patmmedia = 0.
    rocmedia = 0.
    precmedia = 0.
    qobsmedia = 0.
    qsimmedia = 0.
    ettsimmedia = 0.
    erro = 0.
    ntaux=0
    ndiasobs=0
    
    areadrenagem = adren(iexut(ib))        
    ibacaux=ib
    call dominioBacia(ibacaux)

    write(filbalanco,*)'sb:',ib
    write(filbalanco,'(11a12)')'ano','qobs','ta','tdew','vv','patm','roc','prec','qsim','ettsim','%erro'
       

    !ano hidrologico          
    incrementoano = ianoaux+1 
    if(imesah==1) incrementoano=ianoaux
    itiniah = julday(imesah,1,ianoaux) !converte o dia de inicio em calendario juliano                   
    itfimah = julday(fmesah,diasmes(fmesah),incrementoano) !converte o dia de inicio em calendario juliano
        
    !loop sobre o tempo        
    do it=1, nt        
        ! arruma as datas
        tempo=ihora/24.+float(it-1)*dtp/86400        
	    jdia=idini+int(tempo) ! dia juliano do calendário 	    	    	    	   
        
        if(jdia >= itiniah)then                      
            !write(*,*)it, jdia, itiniah
            !pause
            !mudou o ano
            if(jdia > itfimah)then
                tamedia = tamedia/real(ncdom)/real(ntaux)
                tdewmedia = tdewmedia/real(ncdom)/real(ntaux)
                vvmedia = vvmedia/real(ncdom)/real(ntaux)
                patmmedia = patmmedia/real(ncdom)/real(ntaux)
                rocmedia = rocmedia/real(ncdom)/real(ntaux)
                precmedia = precmedia/real(ncdom)
                ettsimmedia = ettsimmedia/real(ncdom)/real(ntaux)
                qsimmedia = qsimmedia/real(ntaux)
                erro = (precmedia-qsimmedia*real(ntaux)-ettsimmedia*real(ntaux))/precmedia*100
                
                if(ndiasobs==0) then 
                    qobsmedia=-1.
                else
                    qobsmedia=qobsmedia/real(ndiasobs)
                endif
                
                !write(*,*)ib,ianoaux,qobsmedia, ndiasobs, it
                !pause
                    
                write(filbalanco,'(i12,10f12.2)')ianoaux,qobsmedia,tamedia,tdewmedia,vvmedia,patmmedia,&
                                                & rocmedia,precmedia,qsimmedia,ettsimmedia,erro
                ianoaux = ianoaux + 1
                erro=0.
                tamedia = 0.
                tdewmedia = 0.
                vvmedia = 0.
                patmmedia = 0.
                rocmedia = 0.
                precmedia = 0.            
                qobsmedia = 0.
                qsimmedia = 0.
                ettsimmedia = 0.
                
                !ano hidrologico
                 incrementoano = ianoaux+1 
                if(imesah==1) incrementoano=ianoaux
                itiniah = julday(imesah,1,ianoaux) !converte o dia de inicio em calendario juliano
                itfimah = julday(fmesah,diasmes(fmesah),incrementoano) !converte o dia de inicio em calendario juliano
                ntaux=0
                ndiasobs=0
            endif
            
            if(qobs(ib,it)>=0.)then                            
                qobsmedia = qobs(ib,it)/areadrenagem*86.4 + qobsmedia
                ndiasobs = ndiasobs + 1
            end if
            
            qsimmedia = qr(ib,it)/areadrenagem*86.4 + qsimmedia                        
            
            ibflag=0            
            !loop sobre as celulas do dominio da bacia
            do ic2=1, ncdom    
                ic = icdom(ic2)                    
                
                tamedia = taall(ic,it) + tamedia
                tdewmedia = tdewall(ic,it) + tdewmedia
                vvmedia = vvall(ic,it) + vvmedia
                patmmedia = patmall(ic,it) + patmmedia
                rocmedia = rocall(ic,it) + rocmedia
                
                if(ibflag/=ibac(ic))then                    
                    ibflag=ibac(ic)
                    ettsimmedia = ettb(ibac(ic),it)*kcb(ibac(ic)) + ettsimmedia                        
                    precmedia = preb(ibac(ic),it)*kcb(ibac(ic)) + precmedia
                endif                
            enddo                                    
           ntaux=ntaux+1            
        endif       		                                    
    enddo   !nt     
enddo !nb

return
end

!---------------------------------------------------------------------------------------------
!subrotina que eh chamada na opcao calibra e seleciona como dominio total da bacia
!---------------------------------------------------------------------------------------------
subroutine dominioBacia(ibcal)
use vars_main, only: celjus,iexut,icdom,nc,ncdom
implicit none
integer ic,ic2,ibcal,ilast

icdom=0 ! vetor que identifica as celulas da bacia a ser simulado/calibrada
ilast=iexut(ibcal) ! exutorio da bacia a ser simulada/calibrada 
ncdom = 0 !numero de celulas do dominio a ser simulado/calibrado
do ic=1,nc        
    ic2=ic                     
    do while (celjus(ic2) <= nc)            
        ic2 = celjus(ic2)                        
        if(ic2 == ilast) then
            ncdom=ncdom+1
            icdom(ncdom) = ic ! indica as celulas dentro do dominio
            exit
        endif
    enddo        
enddo
ncdom=ncdom+1   
icdom(ncdom)=ilast ! o exutorio tambem faz parte do dominio
return
end
