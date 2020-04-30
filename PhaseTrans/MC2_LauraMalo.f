c	Laura Malo Roset

c -------------- Declaració de variables ----------------

       implicit none
       integer*4 nrand,irand,L,N                    !dimensions 
       integer*4 imc,ipas,ivec,itemp,illav,i,j      !contadors bucles
       integer*4 llav0,mcd,mcini,nllav,mctot        !parametres simulacio
       real*8 temp,temp0,a,delta,suma
       real*8 mag,magne,de,ene,energ
       real*8 sum,sume,sume2,summ,summ2,sumam       !promitjos
       real*8 vare,varm                             !variances
       real*4 rrand(1:196632)                       !vector nombres aleatoris
       integer*4 PBC(0:257)                         !condicions contorn
       integer*2 S(1:256,1:256)                     !matriu configuracuó  
       character*17 nom                             !fitxer sortida
       real*4 time1,time2                           !temps CPU


c ------------- Codi principal --------------------------

       call cpu_time(time1)

       !dades simulacio 
       L=32                     !tamany
       temp0=1.4d0              !temperatura inicial              
       mctot=40000              !número d'iteracions de montecarlo
       mcini=2000               !valor imc a partir del que calculem els promitjos
       mcd=20                   !pas de imc per fer els promitjos
       nllav=200                !número de llavors
       llav0=117654             !llavor inicial
       nom= "MC40000-L032-L200" !fitxer de sortida

       !parametres calculats a partir de les dades de la simulació
       N=L*L
       nrand=3*N

       !calculem les condicions de contorn (primers veins)
       PBC(0)=L
       PBC(L+1)=1
       do i=1,L
         PBC(i)=i
       enddo
      

       !obrim el fitxer on guardarem les dades
       open(unit=13,file=nom//".res")

       
       !bucle de temperatudes
       temp=temp0
       do itemp=0,200
        temp=temp0+0.01*itemp

       !inicialitcem tots els contadors
        sum=0.d0     
        sume=0.0d0           !contador energia
        sume2=0.0d0          !contador energia al quadrat
        summ=0.0d0           !contador magnetització
        summ2=0.0d0          !contador magnetització al quadrat
        sumam=0.0d0          !contador valor absolut magnetització
      
      !comença el bucle de llavors
        do illav=llav0,llav0+nllav-1,1
c         write(*,*) 'llavor=', illav
          call rcarin(illav,rrand,nrand)
          call rcarry(rrand,nrand)      !creem el vector rrand 

       !creem la configuració inicial d'spin      
          irand=1
          do i=1,L
            do j=1,L
              if (rrand(irand).lt.0.5e0) then
                S(i,j)=+1               !assignem S=1 per els valors més petits a 0.5
              else 
                S(i,j)=-1               !assginem S=-1 per els valors més grans a 0.5
              endif
              irand=irand+1
            enddo
          enddo
c       write(*,*) S                  !escrivim per pantalla la matriu S
c       call writeconfig(S,L)         !cridem la subrutina que escriu la configuració

       !calculem l'energia i la magnetització inicial                    
          mag=magne(S,L)
          ene=energ(S,L,PBC)

       !comencem el bucle de les iteracions de MC (algoritme de Metrpolis)
          do imc=1,mctot    
            call rcarry(rrand,nrand)
            ivec=1
       	    do ipas=1,N
              i=int(L*rrand(ivec))+1    !trobem un spin aleatori amb la posció (i,j)
              j=int(L*rrand(ivec+1))+1
              delta=rrand(ivec+2)
              ivec=ivec+3
              suma=S(i,PBC(j+1))+S(i,PBC(j-1))+
     +            S(PBC(i+1),j)+S(PBC(i-1),j)
              de=2*S(i,j)*suma
              if (de.lt.0.d0) then      !acceptem o no el canvi d'spí amb una cerca probabilitat
                  S(i,j)=-S(i,j)
              else 
                  a=dexp(-de/temp)
                  if (delta.lt.a) then
                    S(i,j)=-S(i,j)
                  endif
              endif
            enddo    !acaba el bucle de ipas
            
            !calculem els valors mig només a uns determiats pasos de MC
            if ((imc.gt.mcini).and.(mcd*(imc/mcd).eq.imc)) then
              ene=energ(S,L,PBC)
              mag=magne(S,L)
              sum=sum+1.d0          
              sume=sume+ene
              sume2=sume2+ene*ene
              summ=summ+mag
              sumam=sumam+abs(mag)
              summ2=summ2+mag*mag
            endif

          enddo      !tanquem bucle MC
        enddo        !tanquem bucle llavors

 
        !normalitcem els promitjos
        sume=sume/(sum*N)               
        sume2=sume2/(sum*(N**2))
        summ=summ/(sum*N)
        sumam=sumam/(sum*N)
        summ2=summ2/(sum*(N**2))
        !calculem els estimadors de les variances
        vare= (sume2-sume*sume)      
        varm=(summ2-summ*summ)

       !calculem els erros
c       ee= (1/N)*(sqrt(vare)/sqrt(sum))
c       em= (1/N)*(sqrt(varm)/sqrt(sum))

       !capacitat calorifica i suscebilitat
c       cv= N*(vare/(temp**2))
c       sus=N*(varm/temp)

        !escrivim per pantalla els resultats obtinguts
c       write(*,*) 'N=',N
c       write(*,*) 'Temperatura=',temp
c       write(*,*) '<e>=',sume
c       write(*,*) '<e^2>=',sume2
c       write(*,*) '<m>=', summ
c       write(*,*) '<|m|>=',sumam
c       write(*,*) '<m^2>=', summ2 
        
        !guardem al fitxer les dades obtingudes
        write(13,*) N,temp,sum,sume,sume2,vare,summ,sumam,summ2,varm   
       enddo        !Tanquem bucle de temperatures      
       
       call cpu_time(time2)

       !calculem el temps que triga a executar
       write(*,*) 'cputime= ', time2-time1, 'segons'

       close(13)     !tanquem el fitxer
       
       end


c ------------- Final del codi principal ----------------


c ------------- Definició de subrutines -----------------

c
c     ******************************************************************
c     *                     SUBROUTINE RCARIN                          *
c     ******************************************************************
c

      SUBROUTINE RCARIN(IJKL,RVEC,LENV)

C----------------------------------------------------------------------
C Inicializa valores antes de llamar a la subrutina RCARRY.
C IJKL debe estar en el rango 0<IJKL<900 000 000.
C Para conseguir los valores standar usados por Marsaglia y Zaman en su
C articulo poner IJKL = 54217137 (I=12, J=34, K=56, L=78)
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
      COMMON /RAN1/ CARRY
      DIMENSION RVEC(LENV+24)

      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177,177) + 2
      J = MOD(IJ,177)     + 2
      K = MOD(KL/169,178) + 1
      L = MOD(KL,169)

      DO 2 II=24,1,-1
        S = 0.0
        T = 0.5
        DO 3 JJ=1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF (MOD(L*M,64).GE.32) S = S+T
          T = 0.5*T
3       CONTINUE
        RVEC(II) = S
2     CONTINUE

      CARRY = 0.0

      RETURN
      END
c
c     ******************************************************************
c     *                     SUBROUTINE RCARRY                          *
c     ******************************************************************
c
      SUBROUTINE RCARRY(RVEC,LENV)
C----------------------------------------------------------------------
C Generador de numeros pseudo-aleatorios. Algoritmo de G. Marsaglia y
C A. Zaman. Genera numeros reales de 32-bits con mantisas de 24 bits,
C comprendidos entre 0 y 1 (1, explicitamente excluido).
C Periodo aproximado : 10**171.
C Admite la generacion de subsecuencias disjuntas.
C                   F. James, 1989
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
      DIMENSION RVEC(LENV+24)
      COMMON /RAN1/ CARRY
      PARAMETER (TWOM24=1.0/16777216.0)
C
      DO 100 IVEC=25,LENV+24
        UNI = RVEC(IVEC-24) - RVEC(IVEC-10) - CARRY
        IF (UNI.LT.0.) THEN
          UNI = UNI + 1.0
          CARRY = TWOM24
        ELSE
          CARRY = 0.0
        ENDIF

        IF(UNI.EQ.0.)THEN
          UNI=RVEC(IVEC-24)*TWOM24
            in48=-48
          IF(UNI.EQ.0.)UNI=2**(in48)
        ENDIF

        RVEC(IVEC) = UNI
100   CONTINUE

      DO 200 I=1,24
200   RVEC(I)=RVEC(LENV+I)

      RETURN
      END
  
c
c     ******************************************************************
c     *                     SUBROUTINE WRITECONFIG                     *
c     ******************************************************************
c
c     Subrutina que guarda en un fitxer la configuració d'spins de la matriu S

      subroutine writeconfig(S,L)
       implicit none
       integer i,j,L
       integer*2 S(1:L,1:L)

       open(14,FILE="P1-configuration.conf")     !fitxer on guardem les dades 'P1-configuration.conf'
       do i=1,L
          do j=1,L
            if (S(i,j).eq.1) then
              write(13,*) i,j                    !si S=1 guardem la posció i,j de la matriu
            else 
              continue
            endif
          enddo
        enddo
       close(14)

      end subroutine

c
c     ******************************************************************
c     *                     FUNCTION MAGNE                             *
c     ******************************************************************
c
c     Funció que calcula la magnetització (mag) del sistema donada la matriu S i la mida L 

      real*8 function magne(S,L)
      integer*2 S(1:256,1:256)
      integer*4 i,j,L
      real*8 mag
      mag=0.0d0
      do i=1,L
        do j=1,L
          mag=mag+S(i,j)                      
        enddo
      enddo
      magne=mag
      return
      end

c
c     ******************************************************************
c     *                     FUNCTION ENERGIA                           *
c     ******************************************************************
c
c     Funció que calcula l'energia del sistema 

      real*8 function energ(S,L,PBC)
      integer*2 S(1:256,1:256)
      integer*4 i,j,L
      integer*4 PBC(0:257)                                          
      real*8 ene
      ene=0.0d0
      do i=1,L
        do j=1,L
          ene=ene-S(i,j)*S(PBC(i+1),J)-S(i,j)*S(i,PBC(j+1))       
        enddo
      enddo
      energ=ene
      return
      end

