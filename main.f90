PROGRAM modelisation_2
    IMPLICIT NONE

    INTEGER :: nbr_constituant
    DOUBLE PRECISION :: charge,temperature,PUISSANCE,P,HLIQ,H_V,debit_vapeur,T_EB
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: x_l,x_g,Teb,Hvap,Cp_vap,Cp_liq,gamma,K
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coef_binaires,Alpha,Coef_Psat,margule_coef,Coef_unquiac,coef_wilson

    INTEGER :: test_choix_teb,test_choix_depart,choix_enthalpie,choix_melange
    INTEGER :: i

    DOUBLE PRECISION :: cp_liq_melange

    !VARIABLES DE DISCO©
    EXTERNAL :: DISCO, FX, DFDX, DFDXDOT, GEX
    INTEGER :: N
    INTEGER :: LIP,LRP,ITASK,LRW,IOPT,ISTATE,LIW,MF,LG,ITOL
    DOUBLE PRECISION :: T,TOUT,RTOL
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: G,X,XDOT,RPAR,F,RWORK,ATOL,BOOL
    INTEGER, DIMENSION(:),ALLOCATABLE:: IPAR,IWORK,JROOT


    CALL execute_command_line("clear")

    PRINT*,'------------------------------------------------------------------------------------------------------------'
    PRINT*,'                                MODELISATION DE LA DISTILLATION DE RAYLEIGH'
    PRINT*,'------------------------------------------------------------------------------------------------------------'
    PRINT*,"                    EQUILIBRE LIQUIDE VAPEUR : CONVENTION SYMETRIQUE, APPROCHE GAMMA_i/PHI"
    PRINT*,"        RESOLUTION DU SYSTEME ALGEBRICO-DIFFERENTIEL PAR LA METHODE DE GEAR DANS L'INTEGRATEUR D.I.S.C.o ® "
    PRINT*,'------------------------------------------------------------------------------------------------------------'
    PRINT*,''
    PRINT*,''
    PRINT*,'Cliquer ENTRER Pour continuer'
    read(*,*)

    OPEN(UNIT=12,FILE="Donnees.txt")
    READ(12,*)nbr_constituant
    CLOSE(12)
    N=2*nbr_constituant+4 !! VALEUR POUR ALLOCATE

    LIP=4  !DIMENSION de IPAR
    ALLOCATE(IPAR(LIP))
    IPAR=0
    !IPAR : vecteur des paramètres entiers Donc ici de nbr_constituant + gestion de la phase
    IPAR(1)=nbr_constituant

    LRP=4*IPAR(1)**2+15*IPAR(1)+3    !DIMENSION de RPAR
    !RPAR : vecteur des paramètres réels qui ne changent pas dans le temps
    ALLOCATE(RPAR(LRP))             !RPAR : vecteur des paramètres réels qui ne changent pas dans le temps
    RPAR=0

    ALLOCATE(x_l(IPAR(1)))
    ALLOCATE(x_g(IPAR(1)))
    ALLOCATE(Teb(IPAR(1)))
    ALLOCATE(Hvap(IPAR(1)))
    ALLOCATE(Cp_vap(IPAR(1)))
    ALLOCATE(Cp_liq(IPAR(1)))
    ALLOCATE(Coef_Psat(5,IPAR(1)))
    ALLOCATE(Coef_binaires(IPAR(1),IPAR(1)))
    ALLOCATE(coef_wilson(IPAR(1),IPAR(1)+1))
    ALLOCATE(Alpha(IPAR(1),IPAR(1)))
    ALLOCATE(margule_coef(IPAR(1),IPAR(1)))
    ALLOCATE(Coef_unquiac(IPAR(1),5))
    ALLOCATE(K(IPAR(1)))
    ALLOCATE(gamma(IPAR(1)))

    x_l=0
    x_g=0
    Teb=0
    Hvap=0
    Cp_vap=0
    Cp_liq=0
    Coef_Psat=0
    Coef_binaires=0
    coef_wilson=0
    Alpha=0
    margule_coef=0
    Coef_unquiac=0
    K=0
    gamma=0

    
    CALL lecture_fichier(nbr_constituant,x_l,x_g,Teb,Hvap,Cp_vap,Cp_liq,Coef_Psat,Coef_binaires,charge,&
    & temperature,PUISSANCE,P,Alpha,margule_coef,Coef_unquiac,coef_wilson,RPAR,LRP,T_EB)


    !CALL execute_command_line("clear")
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                                           Depart de la modelisation ? "
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                      TAPER 1 >>>>>>>>>>>>>>>>> DEPART T=",temperature,'K'
    PRINT*,"                      TAPER 2 >>>>>>>>>>>>>>>>> DEPART T=T_{ebulition,melange}"
73  READ(*,*)test_choix_depart
    IF (test_choix_depart == 1) THEN
        IPAR(2)=1
        PRINT*,""
        PRINT*, "                                               CHOIX VALIDE !"
        PRINT*,""
    ELSEIF (test_choix_depart==2) THEN
        IPAR(2)=2
        PRINT*,""
        PRINT*, "                                               CHOIX VALIDE !"
        PRINT*,""
    ELSE
        PRINT*,""
        PRINT*, "                                           ERREUR ! ENTRER 1 OU 2"
        PRINT*,""
        GO TO 73
    ENDIF

    !----------------------------------------------------------X----------------------------------------------------------
    ALLOCATE(X(N))
    X=0
    X(1)=charge
    debit_vapeur=0
    X(2)=debit_vapeur
    !X(3)=temperature TEMPÉRATURE à t=0 : NON CONNU À CE STADE
    !X(4)=HLIQ TEMPÉRATURE à t=0 : NON CONNU À CE STADE ====>> CALCUL IMPOSSIBLE
    DO i=5,nbr_constituant+4
        X(i)=x_l(i-4)
    ENDDO
    DO i=nbr_constituant+5,(nbr_constituant+4)+nbr_constituant
        X(i)=x_g(i-(nbr_constituant+4))
    ENDDO

    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                                          Choix du modele de melange "
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                      TAPER 1 >>>>>>>>>>>>>>>>> MODELE NRTL "
    PRINT*,"                      TAPER 2 >>>>>>>>>>>>>>>>> MODELE MARGULES (SI NC=2) !! "
    PRINT*,"                      TAPER 3 >>>>>>>>>>>>>>>>> MODELE IDEALE"
    PRINT*,"                      TAPER 4 >>>>>>>>>>>>>>>>> MODELE UNIQUAC (à verifier!)"
    PRINT*,"                      TAPER 5 >>>>>>>>>>>>>>>>> MODELE WILSON"

190 READ(*,*)choix_melange
    IF (choix_melange/=1 .and. choix_melange/=2 .and. choix_melange/=3 &
    &.and. choix_melange/=4 .and. choix_melange/=5) THEN
        PRINT*,""
        PRINT*, "                                       ERREUR ! ENTRER 1 OU 2 OU 3 OU 4 OU 5"
        PRINT*,""
        GO TO 190
    ENDIF
    PRINT*,""
    PRINT*, "                                               CHOIX VALIDE !"
    PRINT*,""
    IPAR(4)=choix_melange

    !CALL execute_command_line("clear")
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                                          Determination de T_{ebulition} "
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                      TAPER 1 >>>>>>>>>>>>>>>>> LECTURE DANS LE FICHIER : 'Donnees.txt' (SELON COMPO !!) "
    PRINT*,"                      TAPER 2 >>>>>>>>>>>>>>>>> CALCUL PAR DICHOTOMIE"
111 READ(*,*)test_choix_teb

    IF (test_choix_teb == 2) THEN
        CALL DICHOTOMIE_T_EBU(RPAR,IPAR,LRP,LIP,X,N,T_EB) !X utlisé pour les x_l
        RPAR(2*IPAR(1)**2+9*IPAR(1)+3)=T_EB
    ELSEIF (test_choix_teb == 1) THEN
        !Temperature = Tfichier
        RTOL=1d-2
    ELSE
        PRINT*,""
        PRINT*, "                                           ERREUR ! ENTRER 1 OU 2"
        PRINT*,""
        GO TO 111
    ENDIF
    PRINT*,""
    PRINT*, "                                               CHOIX VALIDE !"
    PRINT*, "                                     T_{ebulition} =",T_EB,"K"
    PRINT*,""


!--------------------------------REMANIMANT DE X SELON LES CONDITIONS DE DÉPART DE LA DISTILATION--------------------------------
    IF (test_choix_depart == 1) THEN !PHASE 1!
        IPAR(2)=1 ! Phase = 1
        !T départ= T lu dans le DOcument DOnnes.txt
    ELSEIF (test_choix_depart==2) THEN !PHASE 2!
        IPAR(2)=2 !Phase = 2
        temperature=T_EB !T départ= Teb lu dans le Document Donnes.txt ou calculé par dichotomie
    ENDIF
    X(3)=temperature

    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                                          Choix du modele d'enthalpie "
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                      TAPER 1 >>>>>>>>>>>>>>>>> MODELE DE ENTHALPIE IDEAL "
    PRINT*,"                      TAPER 2 >>>>>>>>>>>>>>>>> MODELE DE ENTHALPIE REEL : CALCUL D'ENTHALPIE D'EXCES"
169 READ(*,*)choix_enthalpie
    IF (choix_enthalpie/=1 .and. choix_enthalpie/=2) THEN
        PRINT*,""
        PRINT*, "                                           ERREUR ! ENTRER 1 OU 2"
        PRINT*,""
        GO TO 169
    ENDIF
    PRINT*,""
    PRINT*, "                                               CHOIX VALIDE !"
    PRINT*,""
    IPAR(3)=choix_enthalpie

    !Une fois la température de départ connue on calcul HLIQ
    CALL H_LIQ(nbr_constituant,temperature,Cp_liq,x_l,choix_enthalpie,choix_melange,RPAR,LRP,HLIQ)
    X(4)=HLIQ
    !Y_i

    !CHOIX MODELE MELANGE
    IF (IPAR(4) == 1) THEN
        CALL NRTL_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 2) THEN
        CALL margule_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 3) THEN
        gamma(:)=1
    ELSEIF (IPAR(4) == 4) THEN
        CALL UNIQUAC_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 5) THEN
        CALL wilson_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    END IF

    !Une fois les gammas obtenu à la température de départ connue on calcul K
    CALL K_EQUILIBRE(IPAR(1),gamma,RPAR(2),X(3),RPAR,LRP,K)
    DO i=nbr_constituant+5,(nbr_constituant+4)+nbr_constituant
        X(i)=K(i-(nbr_constituant+4))*X(i-nbr_constituant)
    ENDDO

    !----------------------------------------------------------XDOT----------------------------------------------------------!
    ALLOCATE(XDOT(N))
    XDOT=0
    CALL H_VAP(nbr_constituant,temperature,Cp_vap,Cp_liq,x_g,Teb,Hvap,H_V)   !Utilisation des varibles LUES et non celle de RPAR, car t=0 <=>
    CALL calcul_CP_LIQ_MEL(nbr_constituant,temperature,x_l,cp_liq,&
        &choix_enthalpie,choix_melange,RPAR,LRP,cp_liq_melange)

    XDOT(1)=-X(2)  !dU/dt=-V
    XDOT(2)=0              !dV/dt=0 à t=0 pour les deux phases
    XDOT(3)= (X(2)*(X(4)-H_V)+ RPAR(1))/(X(1)*cp_liq_melange)   !dT/dt=(V.mHv+Q)/(U.cp_liq)
    XDOT(4)= (X(2)*(X(4)-H_V)+ RPAR(1))/X(1)                    !d(U.h)/dt=-V.mHv+Q <=> dh/dt=(-V.mHv+Q)/U
    DO i=5,nbr_constituant+4
        XDOT(i)=(X(2)/X(1))*(X(i)-X(IPAR(1)+i))                   !dxi/dt=(-V.yi-xi.dU/dt)/U=V.(xi-yi)/U
    ENDDO
    DO i=nbr_constituant+5,(nbr_constituant+4)+nbr_constituant
        XDOT(i)=0
    ENDDO
    !----------------------------------------------------------DISCO----------------------------------------------------------
    
    !----------------------------------------------- DP et INT de DISCo -------------------------------------------------------------!

    !LRP fait
    !LIP fait
    LG=1                        !DIMENSION de G soit le nombre d'évenments possibles
    ITASK=1                     !Calcul de X à T=TOUT par dépassement et interpolation
    LRW=278+11*N+5*(N*N)+3*LG   !DIMENSION de RWORK pour une matrice pleine (NE=N*N)
    IOPT=1                      !Existance d'entrées optionnelles car =1
    ISTATE=1                    !Valeur de DX/dT DOnnées par utilisateur
    LIW=251+36*N+5*(N*N)        !DIMENSION de IWORK
    MF=22                       !MF=10*METH+MITER or METH=2 : méthode de GEAR et MITER=2 :mat dense calc num.
    ITOL=2                      !RTOL : scalaire & ATOL : vecteur
    RTOL=1d-8                   !RTOL : paramètre de tolérance sur erreurs relatives

    T=0
    TOUT=10

    !----------------------------------------------- Vecteur et Matrice de DISCo -------------------------------------------------------------!

    ALLOCATE (JROOT(LG),G(LG),BOOL(N),F(N),ATOL(N),IWORK(LIW),RWORK(LRW))
    JROOT=0                     !EVENEMENT
    G=0
    BOOL=0
    F=0
    ATOL=0                      !Vecteurs de paramètres de tolérance sur erreurs absolues
    IWORK=0                     !Vecteurs d'entiers : espace de travail de DISCo
    RWORK=0                     !Vecteurs de réels : espace de travail de DISCo


    ATOL(1)=1
    ATOL(2)=1
    ATOL(3)=1d-2
    ATOL(4)=1
    DO i=1,2*IPAR(1)
        ATOL(i+4)=1d-8
    ENDDO

    IWORK(70)=6                 !IWORK(70) : entrée optionnelle , =6 : affiche erreur ISTATE à l'écran

    JROOT(1)=0

    BOOL=1                      !BOOL(i)=1 alors l'équation i est prise en compte lors du cal. d'erreur

    OPEN(UNIT=17,FILE="resultats_modelo.xls")
    OPEN(UNIT=18,FILE="resultats_modelo.txt")

    WRITE(17,'(*(g15.9,:,"' //ACHAR(9)// '"))') 'Temps (s)','Charge (kmol)','x_{1}',&
    &'x_{2}','x_{3}','y_{1}','y_{2}','y_{3}','Température (K)','Débit Vapeur (kmol/s)'
    WRITE(18,'(*(g15.9,:,"' //ACHAR(9)// '"))') 'Temps (s)','Charge (kmol)','x_{1}',&
    &'x_{2}','x_{3}','y_{1}','y_{2}','y_{3}','Température (K)','Débit Vapeur (kmol/s)'


    OPEN(UNIT=295,FILE="evolution_temps.xls")
    WRITE(295,'(*(g15.9,:,"' //ACHAR(9)// '"))')"Temps(s)","gamma_1","gamma_2","gamma_3","K_1","K_2","K_3"
    !ECRITURE de t=0 SUR EXCEL
    CALL ecriture_excel(N,IPAR,LIP,RPAR,LRP,t,X)

    PRINT*," Appel de DISCo, VEUILLEZ PATIENTEZ..."

    DO WHILE (ISTATE==1 .or. ISTATE==2)
        CALL DISCO(FX, DFDX, DFDXDOT,N, X, XDOT, T, TOUT,ITOL, RTOL, ATOL,ITASK, ISTATE, IOPT,RWORK, LRW, IWORK, &
        & LIW, MF,GEX, LG, JROOT, BOOL,RPAR, LRP, IPAR, LIP)

        !ECRITURE SUR EXCEL
        CALL ecriture_excel(N,IPAR,LIP,RPAR,LRP,t,X)

        !Test de la valeur de ISTATE
        IF (ISTATE<=1) THEN
            PRINT*,""
            PRINT*,'                      ---------------------------------------------------------------------'
            PRINT*,"                                              Erreur d'integration "
            PRINT*,'                      ---------------------------------------------------------------------'
            PRINT*,"                       ISTATE=",ISTATE
            STOP
        ELSE IF (ISTATE==2) THEN
        !Tout s'est bien passé, on passe au temps d'après
            TOUT=TOUT+10  !Incrémentation de TOUT de 10 secondes   CHANGEMENT DU PAS ???
        ELSE IF (ISTATE==3) THEN
            !Occurence d'un évenement
            IF (IPAR(2)==1) THEN
                IPAR(2)=2  !Passage à la phase diphasique
                ISTATE=1 !Valeur de istate pour la prochaine entrée de DISCo =>> Suite des calculs non basé sur les anciens
            ELSE
                !Fin de le distillation, on sort les valeurs finales
            ENDIF
        ENDIF
    ENDDO
    
    CLOSE(17)
    CLOSE(18)
    CLOSE(295)

    PRINT*,""
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"                                                   Distillation FINIE !"
    PRINT*,'                      ---------------------------------------------------------------------'
    PRINT*,"RESULTATS DANS :"
    PRINT*,"      -resultats_modelo.txt"
    PRINT*,"      -resultats_modelo.xls"
    !CALL execute_command_line("open resultats_modelo.xls")

    DEALLOCATE(x_l)
    DEALLOCATE(x_g)
    DEALLOCATE(Teb)
    DEALLOCATE(Hvap)
    DEALLOCATE(Cp_vap)
    DEALLOCATE(Cp_liq)
    DEALLOCATE(Coef_Psat)
    DEALLOCATE(Coef_binaires)
    DEALLOCATE(coef_wilson)
    DEALLOCATE(margule_coef)
    DEALLOCATE(Coef_unquiac)
    DEALLOCATE(Alpha)
    DEALLOCATE(K)
    DEALLOCATE(gamma)
END PROGRAM modelisation_2

SUBROUTINE ecriture_excel(N,IPAR,LIP,RPAR,LRP,t,X)
    IMPLICIT NONE

    INTEGER :: LIP,N,LRP
    INTEGER, DIMENSION(LIP) :: IPAR
    DOUBLE PRECISION,DIMENSION(LRP)::RPAR
    DOUBLE PRECISION, DIMENSION(N) :: X
    DOUBLE PRECISION :: t

    INTEGER :: i
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: gamma,K

    ALLOCATE(gamma(IPAR(1)),K(IPAR(1)))


    IF (IPAR(2)==1)THEN !PHASE 1 Yi = 0
        WRITE(17,'(*(g15.9,:,"'// ACHAR(9)// '"))')t,X(1),(X(I),I=5,IPAR(1)+4),(0.0,I=IPAR(1)+5,2*IPAR(1)+4),x(3),X(2),X(4)
        WRITE(18,'(*(g15.9,:,"'// ACHAR(9)// '"))')t,X(1),(X(I),I=5,IPAR(1)+4),(0.0,I=IPAR(1)+5,2*IPAR(1)+4),x(3),X(2),X(4)
    ELSE
        WRITE(17,'(*(g15.9,:,"'// ACHAR(9)// '"))')t,X(1),(X(I),I=5,IPAR(1)+4),(X(I),I=IPAR(1)+5,2*IPAR(1)+4),x(3),X(2),X(4)
        WRITE(18,'(*(g15.9,:,"'// ACHAR(9)// '"))')t,X(1),(X(I),I=5,IPAR(1)+4),(X(I),I=IPAR(1)+5,2*IPAR(1)+4),x(3),X(2),X(4)
    END IF

    IF (IPAR(4) == 1) THEN
        CALL NRTL_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 2) THEN
        CALL margule_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 3) THEN
        gamma(:)=1
    ELSEIF (IPAR(4) == 4) THEN
        CALL UNIQUAC_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 5) THEN
        CALL wilson_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    END IF

    CALL K_EQUILIBRE(IPAR(1),gamma,RPAR(2),X(3),RPAR,LRP,K)

    WRITE(295,'(*(g15.9,:,"'// ACHAR(9)// '"))')t,(gamma(I),I=1,IPAR(1)),(K(I),I=1,IPAR(1))

    DEALLOCATE(gamma)
END SUBROUTINE ecriture_excel


SUBROUTINE lecture_fichier(nbr_constituant,x_l,x_g,Teb,Hvap,Cp_vap,Cp_liq,Coef_Psat,Coef_binaires,charge,&
    & temperature,PUISSANCE,P,Alpha,margule_coef,Coef_unquiac,coef_wilson,RPAR,LRP,T_EB)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbr_constituant,LRP

    DOUBLE PRECISION, INTENT(OUT),DIMENSION(nbr_constituant) :: x_l,x_g,Teb,Hvap,Cp_vap,Cp_liq
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(5,nbr_constituant) :: Coef_Psat
    DOUBLE PRECISION, INTENT(OUT) :: charge,temperature,PUISSANCE,P,T_EB
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant,nbr_constituant):: Coef_binaires,Alpha,margule_coef
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant,nbr_constituant+1):: coef_wilson
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant,5):: Coef_unquiac 
    DOUBLE PRECISION,DIMENSION(LRP), INTENT(OUT)::RPAR

    INTEGER :: i,j,k

    OPEN(UNIT=12,FILE="Donnees.txt")

    READ(12,*)

    RPAR=0

    READ(12,*)(x_l(j), j=1,nbr_constituant)
    READ(12,*)(x_g(j), j=1,nbr_constituant)
    READ(12,*)charge
    READ(12,*)temperature
    READ(12,*)PUISSANCE
    RPAR(1)=PUISSANCE
    READ(12,*)P
    RPAR(2)=P

    DO i=1,nbr_constituant
        READ(12,*)Teb(i)
        RPAR(i+2)=Teb(i)
    ENDDO

    DO i=1,nbr_constituant
        READ(12,*)Hvap(i)
        RPAR(2+nbr_constituant+i)=Hvap(i)
    ENDDO

    DO i=1,nbr_constituant
        READ(12,*)Cp_vap(i)
        RPAR(2+2*nbr_constituant+i)=Cp_vap(i)
    ENDDO

    DO i=1,nbr_constituant
        READ(12,*)Cp_liq(i)
        RPAR(2+3*nbr_constituant+i)=Cp_liq(i)
    ENDDO
    READ(12,*)
    READ(12,*)
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            READ(12,*)Coef_binaires(i,j)
            RPAR(2+4*nbr_constituant+ k)=Coef_binaires(i,j)
        ENDDO
    ENDDO

    READ(12,*)
    READ(12,*)
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            READ(12,*)Alpha(i,j)
            RPAR(2+4*nbr_constituant+nbr_constituant*nbr_constituant+ k )=Alpha(i,j)
        ENDDO
    ENDDO

    READ(12,*)
    READ(12,*)
    READ(12,*)
    k=0
    DO j=1,nbr_constituant
        DO i=1,5
            k=k+1
            READ(12,*)Coef_Psat(i,j)
            RPAR(2+4*nbr_constituant+2*nbr_constituant*nbr_constituant + k )=Coef_Psat(i,j)
        ENDDO
    ENDDO

    READ(12,*)
    READ(12,*)T_EB
    RPAR(2*nbr_constituant**2+9*nbr_constituant+3)=T_EB

    READ(12,*)
    READ(12,*)
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            READ(12,*)margule_coef(i,j)
            RPAR(2*nbr_constituant**2+9*nbr_constituant+3 + k )=margule_coef(i,j)
        ENDDO
    ENDDO

    READ(12,*)
    READ(12,*)
    k=0
    DO j=1,5
        DO i=1,nbr_constituant
            k=k+1
            READ(12,*)Coef_unquiac(i,j)
            RPAR(3*nbr_constituant**2+9*nbr_constituant+3 + k )=Coef_unquiac(i,j)
        ENDDO
    ENDDO

    READ(12,*)
    READ(12,*)
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            READ(12,*)coef_wilson(i,j)
            RPAR(3*nbr_constituant**2+14*nbr_constituant+3 + k )=coef_wilson(i,j)
        ENDDO
    ENDDO
    k=0
    DO i=1,nbr_constituant
        k=k+1
        READ(12,*)coef_wilson(i,nbr_constituant+1)
        RPAR(4*nbr_constituant**2+14*nbr_constituant+3 + k )=coef_wilson(i,nbr_constituant+1)
    ENDDO
    
    CLOSE(12)
END SUBROUTINE lecture_fichier

FUNCTION p_sat(nbr_constituant,constituant,temperature,RPAR,LRP) result(PSAT)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: constituant,nbr_constituant,LRP
    DOUBLE PRECISION, INTENT(IN) :: temperature
    DOUBLE PRECISION, DIMENSION(LRP), INTENT(IN) :: RPAR

    DOUBLE PRECISION :: PSAT

    INTEGER :: i,j,k_variable
    DOUBLE PRECISION :: a,b,c,d,e
    DOUBLE PRECISION, DIMENSION(:,:),allocatable :: Coef_Psat

    ALLOCATE(Coef_Psat(5,nbr_constituant))
    Coef_Psat=0

    !--- Réécriture de Coef_Psat ---!
    k_variable=0
    DO j=1,nbr_constituant
        DO i=1,5
            k_variable=k_variable+1
            Coef_Psat(i,j)=RPAR(2+4*nbr_constituant+2*nbr_constituant*nbr_constituant + k_variable )
        ENDDO
    ENDDO

    a=Coef_Psat(1,constituant)
    b=Coef_Psat(2,constituant)
    c=Coef_Psat(3,constituant)
    d=Coef_Psat(4,constituant)
    e=Coef_Psat(5,constituant)

    PSAT=exp(A+(B/(temperature))+C*log((temperature))+D*(temperature)**(E))

    DEALLOCATE(Coef_Psat)
END FUNCTION p_sat

SUBROUTINE margule_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gamma)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbr_constituant
    INTEGER, INTENT(IN) :: LRP
    DOUBLE PRECISION, INTENT(IN) :: temperature
    DOUBLE PRECISION ,INTENT(IN),DIMENSION(nbr_constituant):: x_l
    DOUBLE PRECISION, INTENT(IN),DIMENSION(LRP) :: RPAR

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant) :: gamma

    INTEGER :: i,j,k
    DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: margule_coef

    ALLOCATE(margule_coef(nbr_constituant,nbr_constituant))
    margule_coef=0

    gamma=0

    !ECRITURE DES COEF DE MARGULES BIEN!
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            margule_coef(i,j)=RPAR(2*nbr_constituant**2+9*nbr_constituant+3 + k )
        ENDDO
    ENDDO

    gamma(1)=exp( (1/(8.314*temperature))* ( (x_l(1)*x_l(2))*(margule_coef(2,1)*x_l(1) + margule_coef(1,2)*x_l(2)) ) )
    gamma(2)=exp( (1/(8.314*temperature))* ( (x_l(2)*x_l(1))*(margule_coef(1,2)*x_l(2) + margule_coef(2,1)*x_l(1)) ) )

    DEALLOCATE(margule_coef)
END SUBROUTINE margule_gamma

SUBROUTINE UNIQUAC_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gamma)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbr_constituant
    INTEGER, INTENT(IN) :: LRP
    DOUBLE PRECISION, INTENT(IN),DIMENSION(LRP) :: RPAR
    DOUBLE PRECISION, INTENT(IN), DIMENSION(nbr_constituant) :: x_l
    DOUBLE PRECISION ,INTENT(IN) :: temperature

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant) :: gamma

    INTEGER :: i,j,k
    DOUBLE PRECISION :: somme,somme2,somme3,somme4,somme5,somme6
    DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: Coef_unquiac,tau_i
    DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: phi_i,teta_i,lngamma_c,lngamma_r

    ALLOCATE(Coef_unquiac(nbr_constituant,5))
    ALLOCATE(phi_i(nbr_constituant),teta_i(nbr_constituant),tau_i(nbr_constituant,nbr_constituant))
    ALLOCATE(lngamma_c(nbr_constituant),lngamma_r(nbr_constituant))
    Coef_unquiac=0
    phi_i=0
    teta_i=0
    tau_i=0
    lngamma_c=0
    lngamma_r=0

    somme=0
    somme2=0
    somme3=0
    somme4=0
    somme5=0
    somme6=0

    gamma=0

    !VALEUR DE COEF "BIEN RANGÉ"!
    k=0
    DO j=1,5
        DO i=1,nbr_constituant
            k=k+1
            Coef_unquiac(i,j)=RPAR(3*nbr_constituant**2+9*nbr_constituant+3 + k )
        ENDDO
    ENDDO

    !VALEUR DE PHI_i,Teta_i!
    DO j=1,nbr_constituant
        DO i=1,nbr_constituant
            somme=somme+x_l(i)*Coef_unquiac(i,2)
            somme2=somme2+x_l(i)*Coef_unquiac(i,1)
        ENDDO
        phi_i(j)=(x_l(j)*Coef_unquiac(j,2))/somme
        teta_i(j)=(x_l(j)*Coef_unquiac(j,1))/somme2
        somme=0
        somme2=0
    ENDDO

    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            tau_i(i,j)=exp(-(Coef_unquiac(i,j))/(8.314*temperature))
        ENDDO
    ENDDO

    !lngamma_c
    DO i=1,nbr_constituant
        lngamma_c(i)=log(phi_i(i)/x_l(i)) + 1-(phi_i(i)/x_l(i)) - (10/2)*Coef_unquiac(i,1)*&
        &( log(phi_i(i)/Coef_unquiac(i,1)) + 1-(phi_i(i)/Coef_unquiac(i,1)))
    ENDDO

    !lngamma_r
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            somme3=0
            DO k=1,nbr_constituant
                somme3=somme3+Coef_unquiac(k,1)*x_l(k)*tau_i(k,j)
            ENDDO
            somme4=somme4+Coef_unquiac(j,1)*x_l(j)*tau_i(i,j)/somme3
            somme5=somme5+Coef_unquiac(j,1)*x_l(j)*tau_i(j,i)
            somme6=somme6+Coef_unquiac(j,1)*x_l(j)
        ENDDO
        lngamma_r(i)=Coef_unquiac(i,1)*( 1 - log(somme5/somme6) -  somme4)
        somme5=0
        somme6=0
        somme4=0
    ENDDO

    DO i=1,nbr_constituant
        gamma(i)=exp(lngamma_c(i) + lngamma_r(i))
    ENDDO

    DEALLOCATE(Coef_unquiac,lngamma_c,lngamma_r,phi_i,teta_i,tau_i)
END SUBROUTINE UNIQUAC_gamma

SUBROUTINE wilson_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gamma)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbr_constituant
    INTEGER, INTENT(IN) :: LRP
    DOUBLE PRECISION, INTENT(IN),DIMENSION(LRP) :: RPAR
    DOUBLE PRECISION, INTENT(IN), DIMENSION(nbr_constituant) :: x_l
    DOUBLE PRECISION ,INTENT(IN) :: temperature

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant) :: gamma

    INTEGER :: i,j,k
    DOUBLE PRECISION :: somme,somme1,somme2,somme3
    DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE :: coef_wilson,Lambda_wilson

    ALLOCATE(coef_wilson(nbr_constituant,nbr_constituant+1))
    ALLOCATE(Lambda_wilson(nbr_constituant,nbr_constituant))
    coef_wilson=0
    Lambda_wilson=0

    somme=0
    somme1=0
    somme2=0
    somme3=0

    gamma=0

    !COEF WILSON BIEN RANGÉ!
    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            coef_wilson(i,j)=RPAR(3*nbr_constituant**2+14*nbr_constituant+3 + k )
        ENDDO
    ENDDO
    k=0
    DO i=1,nbr_constituant
        k=k+1
        coef_wilson(i,nbr_constituant+1)=RPAR(4*nbr_constituant**2+14*nbr_constituant+3 + k )
    ENDDO

    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            Lambda_wilson(i,j)=(coef_wilson(j,nbr_constituant+1)/coef_wilson(i,nbr_constituant+1))*&
            &exp(-(coef_wilson(i,j))/(8.314*temperature))
        ENDDO
    ENDDO

    DO k=1,nbr_constituant
        DO i=1,nbr_constituant
            DO j=1,nbr_constituant
                somme=somme + x_l(j)*Lambda_wilson(j,i)
            ENDDO
            somme2=somme2+ x_l(i)*Lambda_wilson(i,k)/somme
            somme=0
        ENDDO
        DO j=1,nbr_constituant
            somme3=somme3+x_l(j)*Lambda_wilson(k,j)
        ENDDO
        somme3=log(somme3)
        gamma(k)=exp(1 - somme3 - somme2)
        somme2=0
        somme3=0
    ENDDO

    DEALLOCATE(coef_wilson,Lambda_wilson)
END SUBROUTINE wilson_gamma

SUBROUTINE NRTL_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gamma)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbr_constituant
    INTEGER, INTENT(IN) :: LRP
    DOUBLE PRECISION ,INTENT(IN):: temperature
    DOUBLE PRECISION, INTENT(IN),DIMENSION(LRP) :: RPAR
    DOUBLE PRECISION, INTENT(IN),DIMENSION(nbr_constituant) :: x_l

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbr_constituant) :: gamma

    INTEGER :: i,j,k,l
    DOUBLE PRECISION :: somme_1, somme_2, somme_3, somme_4, somme_5, somme_6
    DOUBLE PRECISION, DIMENSION(:,:), allocatable :: tau,g_i,Coef_binaires,Alpha


    ALLOCATE(tau(nbr_constituant,nbr_constituant),g_i(nbr_constituant,nbr_constituant))
    ALLOCATE(Coef_binaires(nbr_constituant,nbr_constituant),Alpha(nbr_constituant,nbr_constituant))
    tau=0
    g_i=0
    Coef_binaires=0
    Alpha=0

    !--- Réécriture de Coef Binaire ---!

    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            Coef_binaires(i,j)=RPAR(2+4*nbr_constituant+ k)
        ENDDO
    ENDDO

    !--- Réécriture de Alpha ---!

    k=0
    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            k=k+1
            Alpha(i,j)=RPAR(2+4*nbr_constituant+nbr_constituant*nbr_constituant+ k)
        ENDDO
    ENDDO

    !--- Gamma par NRTL ---!

    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            tau(i,j)=Coef_binaires(i,j)/(temperature*(1.98751)) !1.98751 = R (en cal/mol)
        ENDDO
    ENDDO

    DO i=1,nbr_constituant
        DO j=1,nbr_constituant
            g_i(i,j) = exp(-Alpha(i,j)*tau(i,j))
        ENDDO
    ENDDO


    DO k=1,nbr_constituant
        somme_1=0
        somme_2=0
        DO l=1,nbr_constituant
            somme_1=somme_1+tau(l,k)*g_i(l,k)*x_l(l)
            somme_2=somme_2+g_i(l,k)*x_l(l)
        ENDDO
        somme_6=0
        DO j=1,nbr_constituant
            somme_3=0
            somme_4=0
            somme_5=0
            DO l=1,nbr_constituant
                somme_3 = somme_3 + g_i(l,j)*x_l(l)
                somme_4= somme_4 + tau(l,j)*g_i(l,j)*x_l(l)
                somme_5 = somme_5 + g_i(l,j)*x_l(l)
            ENDDO
            somme_6=somme_6+((x_l(j)*g_i(k,j))/somme_3)*(tau(k,j)-(somme_4/somme_5))
        ENDDO
        gamma(k) = exp((somme_1/somme_2)+somme_6)
    ENDDO

    DEALLOCATE(tau,g_i,Coef_binaires,Alpha)
END SUBROUTINE NRTL_gamma

SUBROUTINE K_EQUILIBRE(nbr_constituant,gamma,P,temperature,RPAR,LRP,K)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbr_constituant,LRP
    DOUBLE PRECISION, INTENT(IN) :: temperature,P
    DOUBLE PRECISION, DIMENSION(nbr_constituant), INTENT(IN) :: gamma
    DOUBLE PRECISION, DIMENSION(LRP), INTENT(IN) :: RPAR

    DOUBLE PRECISION, DIMENSION(nbr_constituant),INTENT(OUT) :: K

    INTEGER :: constituant
    DOUBLE PRECISION :: PSAT,p_sat

    !--- K en approche GAMMA/PHI ---!
    DO constituant=1,nbr_constituant
        PSAT=p_sat(nbr_constituant,constituant,temperature,RPAR,LRP)
        K(constituant)=gamma(constituant)*PSAT/P  !approche gamma/phi => Ptot * K_i = gamma_i * Psat_i
    ENDDO


END SUBROUTINE K_EQUILIBRE

SUBROUTINE DICHOTOMIE_T_EBU(RPAR,IPAR,LRP,LIP,X,N,T_EB)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LIP,LRP,N
    DOUBLE PRECISION, DIMENSION(LRP), INTENT(IN) :: RPAR
    INTEGER, DIMENSION(LIP), INTENT(IN) :: IPAR
    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: X

    DOUBLE PRECISION, INTENT(OUT) :: T_EB

    INTEGER :: i
    DOUBLE PRECISION :: a,b, eps,somme,f_A,f_z
    DOUBLE PRECISION :: p_sat
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: gamma

    ALLOCATE(gamma((N-4)/2))

    eps=1D-8
    a=100
    b=1000

    DO WHILE (abs(b-a)>eps)
        T_EB=(b+a)/2
        somme=0
        DO i=1,IPAR(1)
            !CHOIX MODELE MELANGE
            IF (IPAR(4) == 1) THEN
                CALL NRTL_gamma(IPAR(1),a,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 2) THEN
                CALL margule_gamma(IPAR(1),a,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 3) THEN
                gamma(:)=1
            ELSEIF (IPAR(4) == 4) THEN
                CALL UNIQUAC_gamma(IPAR(1),a,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 5) THEN
                CALL wilson_gamma(IPAR(1),a,X(5),RPAR,LRP,gamma)
            END IF
            somme=somme+gamma(i)* p_sat(IPAR(1),i,a,RPAR,LRP)*X(4+i)
        ENDDO
        f_A=somme-RPAR(2)

        somme=0
        DO i=1,IPAR(1)
            !CHOIX MODELE MELANGE
            IF (IPAR(4) == 1) THEN
                CALL NRTL_gamma(IPAR(1),T_EB,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 2) THEN
                CALL margule_gamma(IPAR(1),T_EB,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 3) THEN
                gamma(:)=1
            ELSEIF (IPAR(4) == 4) THEN
                CALL UNIQUAC_gamma(IPAR(1),T_EB,X(5),RPAR,LRP,gamma)
            ELSEIF (IPAR(4) == 5) THEN
                CALL wilson_gamma(IPAR(1),T_EB,X(5),RPAR,LRP,gamma)
            END IF
            somme=somme+gamma(i)* p_sat(IPAR(1),i,T_EB,RPAR,LRP) *X(4+i)
        ENDDO
        f_z=somme-RPAR(2)

        IF (f_A*f_z > 0) THEN
            a=T_EB
        ELSE
            b=T_EB
        ENDIF
    END DO

    DEALLOCATE(gamma)
END SUBROUTINE DICHOTOMIE_T_EBU

SUBROUTINE H_LIQ(nbr_constituant,temperature,Cp_liq,x_l,choix_enthalpie,choix_melange,RPAR,LRP,HLIQ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbr_constituant,choix_enthalpie,LRP,choix_melange
    DOUBLE PRECISION, DIMENSION(nbr_constituant),INTENT(IN) :: Cp_liq,x_l
    DOUBLE PRECISION, INTENT(IN):: temperature
    DOUBLE PRECISION, DIMENSION(LRP),INTENT(IN) :: RPAR

    DOUBLE PRECISION, INTENT(OUT) :: HLIQ

    INTEGER :: i
    DOUBLE PRECISION :: T_ref,H_REF,temperaturedt
    DOUBLE PRECISION, DIMENSION(:) ,allocatable :: gammaT,gammaTDT

    ALLOCATE(gammaT(nbr_constituant))
    ALLOCATE(gammaTDT(nbr_constituant))
    gammaT=0
    gammaTDT=0

    H_REF=0
    T_ref=25+273  ! 298K
    HLIQ=0

    IF (choix_enthalpie==1) THEN        !mélange idéal
        DO i=1,nbr_constituant
            HLIQ=HLIQ+ H_REF +x_l(i)*Cp_liq(i)*(temperature-T_ref)
        ENDDO
    ELSEIF (choix_enthalpie==2) THEN    !mélange réel avec enthalpie d'excès
        temperaturedt=temperature+(1/temperature)
        IF (choix_melange == 1) THEN
            CALL NRTL_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gammaT)
            CALL NRTL_gamma(nbr_constituant,temperaturedt,x_l,RPAR,LRP,gammaTDT)
        ELSEIF (choix_melange == 2) THEN
            CALL margule_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gammaT)
            CALL margule_gamma(nbr_constituant,temperaturedt,x_l,RPAR,LRP,gammaTDT)
        ELSEIF (choix_melange == 3) THEN
            gammaT(:)=1
            gammaTDT(:)=1
        ELSEIF (choix_melange == 4) THEN
            CALL UNIQUAC_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gammaT)
            CALL UNIQUAC_gamma(nbr_constituant,temperaturedt,x_l,RPAR,LRP,gammaTDT)
        ELSEIF (choix_melange == 5) THEN
            CALL wilson_gamma(nbr_constituant,temperature,x_l,RPAR,LRP,gammaT)
            CALL wilson_gamma(nbr_constituant,temperaturedt,x_l,RPAR,LRP,gammaTDT)
        END IF

        DO i=1,nbr_constituant
            HLIQ=HLIQ+ H_REF+x_l(i)*Cp_liq(i)*(temperature-T_ref)-8.1314d3*(temperature**2)&   !R en J/kmol
            &*((log(gammaTDT(i))-log(gammaT(i)))/(temperaturedt-temperature))     !calcul par perturbation numerique
        ENDDO
    ENDIF

    DEALLOCATE(gammaTDT,gammaT)
END SUBROUTINE H_LIQ

SUBROUTINE H_VAP(nbr_constituant,temperature,Cp_vap,Cp_liq,x_g,Teb,Hvap,H_V)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbr_constituant
    DOUBLE PRECISION, INTENT(IN) :: temperature
    DOUBLE PRECISION, DIMENSION(nbr_constituant), INTENT(IN)::Cp_vap,Teb,Hvap,Cp_liq,x_g

    DOUBLE PRECISION,INTENT(OUT) :: H_V

    DOUBLE PRECISION :: Tref,HREF
    INTEGER :: i

    HREF=0
    Tref=25+273
    H_V=0

    DO i=1,nbr_constituant
        H_V=H_V+x_g(i)*(HREF+Cp_liq(i)*(Teb(i)-Tref))+Hvap(i)*x_g(i)+x_g(i)*(Cp_vap(i)*(temperature-Teb(i)))
    ENDDO

END SUBROUTINE H_VAP

SUBROUTINE calcul_CP_LIQ_MEL(nbr_constituant,temperature,x_l,cp_liq,choix_enthalpie,choix_melange,RPAR,LRP,cp_liq_melange)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbr_constituant,choix_enthalpie,LRP,choix_melange
    DOUBLE PRECISION, INTENT(IN) :: temperature
    DOUBLE PRECISION, DIMENSION(nbr_constituant),INTENT(IN) :: cp_liq,x_l
    DOUBLE PRECISION, DIMENSION(LRP),INTENT(IN) :: RPAR

    DOUBLE PRECISION, INTENT(OUT):: cp_liq_melange

    DOUBLE PRECISION :: T,TdT,HLIQ_T,HLIQ_TdT
    INTEGER :: i

    cp_liq_melange=0
    DO i=1,nbr_constituant
        cp_liq_melange = cp_liq_melange + x_l(i)*cp_liq(i)
    ENDDO

    T=temperature
    TdT=temperature+(1/temperature)
    CALL H_LIQ(nbr_constituant,T,Cp_liq,x_l,choix_enthalpie,choix_melange,RPAR,LRP,HLIQ_T)
    CALL H_LIQ(nbr_constituant,TdT,Cp_liq,x_l,choix_enthalpie,choix_melange,RPAR,LRP,HLIQ_TdT)
    cp_liq_melange=(HLIQ_TdT-HLIQ_T)/(TdT-T)

END SUBROUTINE calcul_CP_LIQ_MEL


!SUBROUTINES DE DISCO - EXTERNAL!

SUBROUTINE FX(N,T,X,XDOT,F,RPAR,LRP,IPAR,LIP)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N,LRP,LIP,T    !T est le temps
    DOUBLE PRECISION, DIMENSION(N),INTENT(IN)  :: X,XDOT
    DOUBLE PRECISION, DIMENSION(LRP),INTENT(IN) :: RPAR
    INTEGER, DIMENSION(LIP),INTENT(IN) :: IPAR

    DOUBLE PRECISION, DIMENSION(N),INTENT(OUT) :: F

    INTEGER :: i
    DOUBLE PRECISION :: HLIQ,h_v
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: gamma,K


    ALLOCATE(gamma((N-4)/2),K((N-4)/2))
    gamma=0
    K=0

    !IPAR(1)=nbr_constituant
    !X(3)=temperature
    !X(5)=x_l           adresse du premier element du vecteur x_l
    !X(IPAR(1)+5)=x_g   adresse du premier element du vecteur x_g

    !---- APPEL DES SUBROUTINES À LA BONNE TEMPÉRATURE ----!
    CALL H_VAP(IPAR(1),X(3),RPAR(3+2*IPAR(1)),RPAR(3+3*IPAR(1)),X(IPAR(1)+5),RPAR(3),RPAR(3+IPAR(1)),H_V)
    IF (IPAR(4) == 1) THEN
        CALL NRTL_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 2) THEN
        CALL margule_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 3) THEN
        gamma(:)=1
    ELSEIF (IPAR(4) == 4) THEN
        CALL UNIQUAC_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    ELSEIF (IPAR(4) == 5) THEN
        CALL wilson_gamma(IPAR(1),X(3),X(5),RPAR,LRP,gamma)
    END IF
    CALL K_EQUILIBRE(IPAR(1),gamma,RPAR(2),X(3),RPAR,LRP,K)
    CALL H_LIQ(IPAR(1),X(3),RPAR(3*IPAR(1)+2 +1),X(5),IPAR(3),IPAR(3),RPAR,LRP,HLIQ)

    !----------- ÉCRITURE DU VECTEUR DES RESIDUS --------
    F(1)=XDOT(1)+X(2)   !Bilan de matière total : dU/dt+V=0
    F(IPAR(1)+2)=(XDOT(1)*X(4)+X(1)*XDOT(4))+X(2)*H_V-RPAR(1)     !Bilan énergétique : d(U.h)/dt+V.H_V-Q=0
    F(N-1)=X(4)-HLIQ !Modèle d'enthalpie liquide

    F(N)=0  !Initialisation de la somme à 0
    DO i=1,IPAR(1)
        F(i+1)= (XDOT(1)*X(i+4)+X(1)*XDOT(i+4)) + X(2)*X(i+IPAR(1)+4)  !Bilan de matière partiel : d(U.x_i)/dt+V.yi=0
        F(i+IPAR(1)+2)=X(IPAR(1)+4+i)-K(i)*X(i+4) !Equations d'équilibre liquide/vapeur : y_i-K(i).x_i=0
        IF ( (IPAR(2)==2) ) THEN !Phase 2 : Diphasique, T=teb
            F(N)=F(N)+X(i+4)-X(i+IPAR(1)+4)   !Σxi-Σyi=0
        ENDIF
    ENDDO
    IF (IPAR(2)==1) THEN  ! Phase 1 : Monophasique
        F(N)=X(2)   !V=0
    ENDIF

    DEALLOCATE(gamma,K)
END SUBROUTINE FX

SUBROUTINE DFDX(N,T,X,XDOT,P,LROWP,ML,MU,RPAR,LRP,IPAR,LIP)
    IMPLICIT NONE
    INTEGER :: N,ML,MU
    INTEGER :: LIP,LROWP,LRP
    DOUBLE PRECISION:: T
    DOUBLE PRECISION, DIMENSION(N)::X,XDOT
    DOUBLE PRECISION, DIMENSION(LRP) :: RPAR
    INTEGER, DIMENSION(LIP) :: IPAR
    DOUBLE PRECISION, DIMENSION(LROWP,N) :: P


END SUBROUTINE DFDX

SUBROUTINE DFDXDOT(N,T,X,XDOT,Q,LROWQ,ML,MU,RPAR,LRP,IPAR,LIP)
    IMPLICIT NONE
    INTEGER :: N,ML,MU
    INTEGER :: LIP,LROWQ,LRP
    DOUBLE PRECISION:: T
    DOUBLE PRECISION, DIMENSION(N)::X,XDOT
    DOUBLE PRECISION, DIMENSION(LRP) :: RPAR
    INTEGER, DIMENSION(LIP) :: IPAR
    DOUBLE PRECISION, DIMENSION(LROWQ,N) :: Q


END SUBROUTINE DFDXDOT

SUBROUTINE GEX(N,T,X,XDOT,G,LG,RPAR,LRP,IPAR,LIP)
    IMPLICIT NONE

    INTEGER :: N,LG,LRP,LIP,gestion_phase
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(N):: X,XDOT
    DOUBLE PRECISION, DIMENSION(LG)::G
    DOUBLE PRECISION, DIMENSION(LRP)::RPAR
    INTEGER, DIMENSION (LIP):: IPAR

    !Fonction : U=0
    !Fonction : T-Teb=0

    gestion_phase=IPAR(2)
    IF (gestion_phase==1) THEN          !G(1)=T-TEB=0 DOnc on est dans la phase 2
        G(1)=X(3)-RPAR(2*IPAR(1)**2+9*IPAR(1)+3)
    ELSEIF (gestion_phase==2) THEN      !G(1)=U=0 : plus rien dans le bouilleur,simulation finie
        G(1)=X(1)
    ENDIF

END SUBROUTINE GEX
!Codé par : Elodie LE GUEN et Erwan BENHENOU
