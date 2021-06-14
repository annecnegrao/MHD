    program mhd
	!******************************************************************
	!***************modelo hidrol�gico distribu�do*********************
	!******************************************************************
	!*		    baseado no modelo mgb do                              *
	!*			instituto de pesquisas hidr�ulicas                    *
	!*          desenvolvido pelo inpe                                *
	!******************************************************************
	!use portlib !biblioteca para calcular o tempo de processamento
	use vars_main !m�dulo que cont�m as vari�veis principais
	use vars_calib !m�dulo que cont�m as vari�veis da calibra��o
	implicit none
	integer julday
		
    ! diretorio de trabalho (especificado em vars_main)
    write(*,*) 'diretorio de trabalho'
    write(*,*) dir_dados

	!____________leitura de arquivos principais________________
	call lefix 	!subrotina de leitura de parametros fixos	
	call alloca_vars(0) !alloca vari�veis principais
	idini=julday(imes,idia,iano) !converte o dia de inicio em calendario juliano
	call leveg !subrotina de leitura de parametros da vegetacao (albedo, iaf, rs e z mensais)
	call lesolo !subrotina associada aos parametros de solo
	call lecell !subrotina de leitura do arquivo das celulas
	call leqobs	!leitura de dados de vazao observada		
	call lemet ! le os dados meterologicos de todos o periodo de simulacao	
	if(numsubst > 0) call lesubst	
	!___________fim da leitura dos arquivos principais ______________

	!_______________preparacao de dados_____________________
	call parcel !calcula alguns par�metros das c�lulas e do rio
	call parcunge !calcula parametros para propagacao muskingum-cunge
	
	!____________fim da preparacao de dados_____________________

	!decide se vai simular, calibrar ou fazer previs�o (icalib=0; icalib=1; icalib=2) 
	calibra_case: select case (icalib) !verifica se calibra ou n�o
	case (0) ! n�o calibra, apenas simula
		write(*,*)
		write(*,*)' fazendo simulacao'
		write(*,*)
		call simula	
	    close (filsol) !arquivo de saida dados do solo
	case (1) ! calibra
		write(*,*)
		write(*,*)' fazendo calibracao'	
		write(*,*)
		call calibra
	case (2) ! simula previs�o
	    write(*,*)
		write(*,*)' fazendo previsao'
		write(*,*)
		call previsao
	case default
		stop ' erro: icalib desconhecido no arquivo parfix.hig!!!'
	end select calibra_case
	
	call alloca_vars(1) !dealloca vari�veis principais

	write(*,*)'programa terminou'
	stop
	end program
