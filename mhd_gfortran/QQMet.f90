!---------------------------------------------------------------------------------------------
!subrotina que verifica a qualidade dos dados meteorologicos
!e faz a tranformacao da unidade de pressao
!---------------------------------------------------------------------------------------------
subroutine qqmet
use vars_main,only: it,nt,ihoraini,dtp,ihora,jdia,idini,imes,idia,iano,cdata,nc,ta, &
                  & tdew,vv,roc,patm
implicit none
integer ic
real tempo ! tempo de simulacao em dias


! verifica a data
tempo=ihoraini/24.+float(it-1)*dtp/86400
ihora=nint((tempo-int(tempo))*24.) ! hora da simulacao
jdia=idini+int(tempo) ! dia juliano do calendário 
call caldat(jdia,imes,idia,iano)
write(cdata(1:2),'(i2.2)') idia
write(cdata(3:4),'(i2.2)') imes
write(cdata(5:8),'(i4)') iano
write(cdata(9:10),'(i2.2)') ihora

do ic=1,nc
    !temperatura do ar em graus centigrados
    if(ta(ic)<-25.0 .or. ta(ic)>50.0)then 
        write(*,*)'erro no valor da temperatura na celula:', ic, ' na data: ', & 
        & cdata(1:2)//'/'//cdata(3:4)//'/'//cdata(5:8)//' '//cdata(9:10)//':00',' valor: ',ta(ic)
        stop
    endif
    if(tdew(ic)>ta(ic)) tdew(ic)=ta(ic) !comparando temperaturas

    !ponto de orvalho em graus centigrados
    if(tdew(ic)<-25 .or. tdew(ic)>50.)then 
	    write(*,*)'erro no ponto de orvalho:', ic, ' na data: ', & 
	    & cdata(1:2)//'/'//cdata(3:4)//'/'//cdata(5:8)//' '//cdata(9:10)//':00',' valor: ',tdew(ic)
	    stop
    endif
	    
    !velocidade do vento em m/s
    if(vv(ic)<0.0 .or. vv(ic)>30.)then
	    write(*,*) 'erro na velocidade do vento na celula:', ic,' na data: ', & 
	    & cdata(1:2)//'/'//cdata(3:4)//'/'//cdata(5:8)//' '//cdata(9:10)//':00',' valor: ',vv(ic)
        stop
    endif	
    	
	!radiacao de onda curta global incidente									
    if(roc(ic)<0.0)then
	    write(*,*)'erro na radicao de onda curta na celula ', ic, ' na data: ', & 
	    & cdata(1:2)//'/'//cdata(3:4)//'/'//cdata(5:8)//' '//cdata(9:10)//':00',' valor:',roc(ic)
	    stop
    endif

    ! prmd pressao em milibares
    if(patm(ic)>1300 .or. patm(ic)<400)then
	    write(*,*)'erro na pressao na celula:', ic, ' na data:', & 
	    & cdata(1:2)//'/'//cdata(3:4)//'/'//cdata(5:8)//' '//cdata(9:10)//':00',' valor:',patm(ic)
        stop
    endif
    patm(ic)=0.1*patm(ic) !transforma pressao de mbar para kPa
    
enddo

return
end