!---------------------------------------------------------------------------------------------
! funcao que determina o dia do calendario juliano correspondente ao dia, mes e ano dados
!---------------------------------------------------------------------------------------------
      function julday(mm,id,iyyy)
      implicit none
      integer,parameter:: igreg=15+31*(10+12*1582)
      integer julday,id,iyyy,mm
      integer ja,jm,jy
      
      jy=iyyy
      if (jy.eq.0)then
          write(*,*) 'julday: there is no year zero'
          read(*,*)
      end if
      if (jy.lt.0) jy=jy+1
      
      if (mm.gt.2) then
          jm=mm+1
      else
          jy=jy-1
          jm=mm+13
      endif
      
      julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
      
      if (id+31*(mm+12*iyyy).ge.igreg) then
          ja=int(0.01*jy)
          julday=julday+2-ja+int(0.25*ja)
      endif

      return
      end
