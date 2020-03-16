!==================================================================================================
!==================================================================================================
!==================================================================================================
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                             I. User Defined variables                                           !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------! 
      Module uArrays
!     Here describe allocatable arrays to be shared with subroutines containing use uArrays

      Integer :: Itercompt            
      Integer , parameter :: Itmax    = 1     ! Max number of internal iteration 
!     The inital pressure is the pressure on the gauss point of the element
      double precision, allocatable :: eVolPrevious(:),eVolumetric(:),b(:), meanStress(:) 
!     vf is the fluid content variation on the gauss point of the element
      double precision, allocatable :: CM(:), vf(:), InitalPressure(:)
      double precision, allocatable :: nodalInitialPressure(:)
      double precision, allocatable :: AUt(:), AUn(:)
      double precision, allocatable :: ATau(:), ASn(:)
      double precision, allocatable :: ADamage(:),totalLength(:)


      character(len=90) :: Fich101   ! Element mechanical Results of rock matrix at each time step
      character(len=90) :: Fich102   ! Nodal mechanical Results of rock matrix at each time step
      character(len=90) :: Fich103   ! Element hydraulic Results of rock matrix at each time step
      character(len=90) :: Fich104   ! Nodal hydraulic Results of rock matrix at each time step


      character(len=90) :: Fich201   ! Mechanical Results of joint elements at each time step
      character(len=90) :: Fich202   ! Average mechanical res of joint elements at each time step
      character(len=90) :: Fich203   ! Hydraulic Results of joint elements at each time step

      character(len=90) :: Fich301   ! Matrix element geometrical information
      character(len=90) :: Fich302   ! Matrix nodal geometrical information
      character(len=90) :: Fich303   ! Joint element geometrical information

!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                             I.1 variables for HM coupling                                       !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      ! variable for controlling hydromechanical coupling in the rock matrix
      Integer , parameter :: matrixHMCoupling   = 0  
      ! variable for controlling hydromechanical coupling in the fracture
      Integer , parameter :: fractureHMCoupling = 0  

      EndModule uArrays
!==================================================================================================      
!==================================================================================================
!==================================================================================================


!==================================================================================================
!==================================================================================================
!==================================================================================================
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                              User Defined module                                                !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Module User
      use Tools
      Contains

!==================================================================================================      
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!             II. Change some options of variables not available in input files                   !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Subroutine uRead
      use Global; implicit none

      !NB(18) = 1 ;      ! At least one material includes plasticity 
      !NB(21) = 1 ;      ! At least one material includes viscous strain
      !NB(105)= 1 ;      ! Exist User-defined Free Strain (in uDefLib)
      !NB(107)= 1 ;      ! Exist User-defined Internal Stress (in uSigInt) 
      !NB(47)=5          ! number of internal variable for Mechanics 
      !NB(77)=1          ! Activate Chemical calculation function
      !NB(67)=1          ! exists sigma interne
      NB(110) = 0        ! Couplage M -> H pris en compte (source SrcH pour variation volumique)
      NB(111) = 0        ! Internal HM iteration process is activated
      print*, 'Calculation Start'
      End Subroutine uRead
!-------------------------------------------------------------------------------------------------! 
!=================================================================================================!



!=================================================================================================! 
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                             User Defined Process (Couplings)                                    !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------! 
!=================================================================================================!


!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                             Before starting Hydro-Thermo-Mechanical loop                        !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                        III. Initialization of variables                                         !
!                   Describe here the complementary initialization commands                       !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uInit
      use Global;use uArrays; implicit none

      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !              III.1. Post treatment variables for joint elements                           !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      integer :: n
      
      if (NB(8).ne.0) then  ! if containing joint elements
            allocate(AUt(1), AUn(1),  ATau(1), ASn(1),            &
            & ADamage(1), totalLength(1) )
            If ((NB(12).eq.1).or.(NB(12).eq.4)) then   ! contains mechanical problem    
                  AUt(1) = 0
                  AUn(1) = 0
                  ATau(1) = 0
                  ASn(1) = 0 
                  ADamage(1) = 0
                  totalLength(1) = 0
            endif
            If ((NB(12).eq.2).or.(NB(12).eq.4)) then   ! contains hydraulic problem    
                  continue
            endif 
      endif

      ! var(18) Total length of joint elements Error: equals 0, total length needs to be calculated
      DO n = 1,NB(2)
            If (ntyp(n).eq.5) then ! joint eleemtns 
                  totalLength(1) = totalLength(1) + s(n)
            endif
      enddo

      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !            III.2.Post treatment variables for rock matrix elements                        !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      ! Aways containing matrix elemetns
      If ((NB(12).eq.2).or.(NB(12).eq.4)) then   ! Contains hydraulic problem
            allocate(InitalPressure(NB(2)),nodalInitialPressure(NB(1)),vf(NB(2)))
            Do n=1,NB(2)
                  InitalPressure(n) = RPg(n)   ! when using resumption
            ENDDO
            Do n=1,NB(1)   ! NB(1)  total number of nodes
                  nodalInitialPressure(n) = RP(n)
            ENDDO
      ENDIF
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !            III.3.Variables for coupling calculation in the joint elements                 !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!



      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !            III.4.Variables for coupling calculation in the rock matrix                    !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      
      allocate(eVolPrevious(NB(2)),eVolumetric(NB(2)),b(NB(2)),meanStress(NB(2)),vf(NB(2)))
      eVolPrevious(:)  = 0  ! Volumetric deformation of the precedent time step
      Itercompt = 0         ! Tempory variable for controlling internal iteration. 
      Do n=1,NB(2)
            ! For the rock matrix in triangular or quatrilateral elements
            If ((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then 
                  eVolumetric(n) = EXY(n,1) + EXY(n,2) + EXY(n,3)
                  meanStress(n) = (SXY(n,1) + SXY(n,2) + SXY(n,3) )/3
                  vf(n) = 0
            Endif
      enddo
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  III.5. File names for joint elements resluts                             !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      if (NB(8).ne.0) then ! if  containing joint elements
 

            Fich201 = trim(foldername)//"\jointMecha.dat"
            Fich202 = trim(foldername)//"\jointMechaAverage.dat"
            Fich203 = trim(foldername)//"\jointHydro.dat"
            Fich303 = trim(foldername)//'jointElemGeo.dat'
            
            open (unit=201, file=Fich201, status='replace')
            open (unit=202, file=Fich202, status='replace')
            open (unit=203, file=Fich203, status='replace')

            open (unit=303, file=Fich303, status='replace')

            ! 201: joint results at each time step
            if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! mechanical problem

                  write(201,*) 'Temps du calcul', char(9), 'NoElem', char(9), 'Ut', char(9),'Un',  &
                  & char(9),'Tau', char(9), 'Sn', char(9), 'Utp', char(9),'Unp', char(9),          &
                  & 'Damage variable'

            else if ((NB(12).eq.2).or.(NB(12).eq.4)) then ! Hydraulic problem

                  write(203,*) 'Temps du calcul', char(9), 'NoElem', char(9), 'Gauss Pressure',    &
                  & char(9), 'Initial Pressure'

            endif
            ! 202: joint average mechanical results
            write(202,*) 'Temps du calcul', char(9), 'NoElem', char(9), 'AUt',char(9), 'AUn',     &
            & char(9), 'ATau',char(9), 'ASn', char(9), 'Average Damage'
            
            ! 303: Joint elem geo data
            write(303,*) 'JointElemNo', char(9), 'Length', char(9), 'Node1', char(9), 'Node2',    &
            & char(9), 'Node3', char(9), 'Node4'

            Do n =1,NB(2)
                  if (ntyp(n).eq.5) then ! joint element
                        write(303,*) n, s(n), konec(n,1),konec(n,2),konec(n,3),konec(n,4)
                  endif 
            enddo
      endif

      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  III.6. File names for rock matrix resluts                                !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!

      
      Fich101 = trim(foldername)//"\matrixMechaElem.dat"
      Fich102 = trim(foldername)//"\matrixMechaNodal.dat"
      Fich103 = trim(foldername)//"\matrixHydroElem.dat"
      Fich104 = trim(foldername)//"\matrixHydroNodal.dat"

      Fich301 = trim(foldername)//"\matrixElemGeo.dat"
      Fich302 = trim(foldername)//"\nodeCoord.dat"

      open (unit=101, file=Fich101, status='replace')
      open (unit=102, file=Fich102, status='replace')
      open (unit=103, file=Fich103, status='replace')
      open (unit=104, file=Fich104, status='replace')
      
      open (unit=301, file=Fich301, status='replace')
      open (unit=302, file=Fich302, status='replace')

      write(101,*)  'Temps du calcul', char(9), 'NoElem', char(9), 'EXX', char(9), 'EYY', char(9),&
      &'EZZ', char(9), '2EXY', char(9), 'Evol', char(9), 'SXX', char(9), 'SYY', char(9),'SZZ',    &
      &char(9),'SXZ', char(9), 'meanStress'
      
      write(102,*)  'Temps du calcul', char(9), 'NoNode', char(9), 'Ux', char(9), 'Uy'
      
      write(103,*)  'Temps du calcul', char(9), 'NoElem', char(9),'Gausse pressure', char(9), 'vf'

      write(104,*)  'Temps du calcul', char(9), 'NoNode', char(9), 'Nodal pressure'

      write(301,*) 'MatrixElemNo', char(9), 'Surface', char(9), 'Node1', char(9), 'Node2', char(9)&
      &, 'Node3'
      write(302,*) 'NodeNo', char(9), 'XNode', char(9), 'YNode'

      Do n =1,NB(2)
            if (ntyp(n).eq.3) then ! triangle element
                  write(301,*) n, s(n), konec(n,1),konec(n,2),konec(n,3)
            elseif (ntyp(n).eq.4) then ! quadrilateral element
                  write(301,*) n, s(n), konec(n,1),konec(n,2),konec(n,3), konec(n,4)
            endif
      enddo

      Do n=1,NB(1)
            write(302,*) n, x(n,1), x(n,2)
      enddo
      End Subroutine uInit
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

       
!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                        IV. Initialization of variables                                          !
! Describe  what should be down at each time increment before calculation of the next increment   !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uCycleEntry
      use Global; use uArrays;implicit none
      

      integer :: n
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                        IV.1  write results at each Time Step for rock matrix                    !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      
      if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! Mecha prob
            Do n = 1,NB(2)
                  if((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                        ! matrix element results mecha
                        write(101,*), Temps, n, EXY(n,1),EXY(n,2),EXY(n,3),EXY(n,4),               &
                        eVolumetric(n), SXY(n,1), SXY(n,2),SXY(n,3),SXY(n,4),meanStress(n)
                  endif
            enddo
            Do n=1,NB(1)
                        ! matrix nodal results mecha
                        write(102,*) Temps, n, RU(n,1),  RU(n,2) 
            enddo
      elseif ((NB(12).eq.2).or.(NB(12).eq.4)) then  ! Hydro prob
            Do n = 1,NB(2)
                  if((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                        ! matrix element results hydro
                        write(103,*), Temps, n, RPg(n), vf(n)
                  endif
            enddo
            Do n=1,NB(1)
                        ! matrix nodal results hydro
                        write(104,*) Temps, n, RP(n) 
            enddo
      end if



!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                        IV.2  write results at each Time Step for joint elements                 !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      if (NB(8).ne.0) then ! Containing joint elemants
            if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! Mecha prob
                  Do n = 1,NB(2)
                        if(ntyp(n).eq.5) then
                              ! joint element results mecha
                              write(201,*), Temps, n, EXY(n,1),EXY(n,2), SXY(n,1), SXY(n,2), Ep(n,1),   &
                              & Ep(n,2), Vinm(n,1)  
                              write(202,*) Temps, n, Aut, AUn, ATau, ASn, ADamage
                        endif
                  enddo

            elseif ((NB(12).eq.2).or.(NB(12).eq.4)) then  ! Hydro prob
                  Do n = 1,NB(2)
                        if(ntyp(n).eq.5) then
                              ! matrix element results hydro
                              write(203,*), Temps, n, RPg(n)
                        endif
                  enddo
            end if
      endif


      End Subroutine uCycleEntry
!
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!


!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 V.Coupling from Hydraulic to Thermal and Mechanics                              !
!      Describe here the effect of Hydraulic variables on the Thermal and Mechical processes      !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Subroutine uHydroTM
      use Global; use uArrays;implicit none

      End Subroutine uHydroTM
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 VI. Coupling from Thermal to Mechanics and Hydraulic                            !
!      Describe here the effect of Thermal variables on the Mechical and Hydraulic processes      !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Subroutine uThermoMH
      use Global; implicit none
!      
!
      End Subroutine uThermoMH
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 VII. Coupling from Mechanics to Hydraulic and Thermal                           !
!      Describe here the effect of Mechanical variables on the Hydraulic and Thermal processes    !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uMechaHT
      use Global;use uArrays; implicit none
!      
!
      integer :: n
      allocate(CM(NB(2)))

      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  VII.1. Coupling for joint elements                                       !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      if (NB(12).eq.4) then
            if (fractureHMCoupling.eq.1) then
                  continue
            endif
      endif
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  VII.2. Coupling for rock matrix                                          !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      if (NB(12).eq.4) then
            if (matrixHMCoupling.eq.1) then
                  DO n=1,NB(2)
                        ! rock matrix could be triangular or quatrilateral elements
                        If ((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then 
                              b(n) = vmat(mat(n),4,2)
                              SrcH(n) = -b(n)*(EXY(n,1) + EXY(n,2) + EXY(n,3)-eVolumetric(n))/dt
                              eVolumetric(n) = EXY(n,1) + EXY(n,2) + EXY(n,3)
                              ! SrcH(n) = -b(n)*(EXY(n,1) + EXY(n,2) + EXY(n,3)-2*eVolumetric(n)+eVolPrevious(n))/dt
                              ! eVolPrevious(n)= eVolumetric(n);eVolumetric(n) = EXY(n,1) + EXY(n,2) + EXY(n,3)
                              ! InitalPressure: Intial pressure on gauss points
                        endif
                  ENDDO
                  ! Mise en place iteration
                  IF (NB(111).eq.1) THEN  ! Procedure iteration HM active
                        if (Itercompt.lt.ItMax) then
                              NB(103)   = 0 ! on n'ajoute plus les charges dues aux conditions limites
                              NB(112)   = 1 ! on retourne au caclul hydraulic au lieu de passer a l'instant suivant
                              Itercompt = Itercompt + 1
                        else
                              Itercompt=0; NB(103)=1; NB(112)=0
                        endif
                  ENDIF
            endif
      endif
 
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  VII.3. Calculate post treatment variables: Joint elements                !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      if (NB(8).ne.0) then ! if containing joint elements
            if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! containing mechanical problem
                  DO n=1,NB(2)
                        if (ntyp(n).eq.5) then ! for joint elements only
                              ! S(n) element's surface
                              AUt(1) = AUt(1) + EXY(n,1)*S(n)/totalLength(1)
                              AUn(1) = AUn(1) + EXY(n,2)*S(n)/totalLength(1)
                              ATau(1) = ATau(1) + SXY(n,1)*S(n)/totalLength(1)
                              ASn(1) = ASn(1) + SXY(n,2)*S(n)/totalLength(1)
                              ADamage(1) = ADamage(1) + Vinm(n,1)*S(n)/totalLength(1)
                        endif
                  ENDDO
            endif
      endif
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  VII.4. Calculate post treatment variables: Rock matrix                   !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!
      Do n=1,NB(2)
            If ((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                  CM(n) = vmat(mat(n),2,2)! the storage coefficient
                  meanStress(n) = (sxy(n,1)+sxy(n,2)+sxy(n,3))/3.
                  if (NB(12).eq.2) then  ! Hyrdo
                        vf(n) = (RPg(n) - InitalPressure(n))*CM(n)
                  elseif (NB(12).eq.4) then  ! Hydromecha
                        vf(n) = b(n)*eVolumetric(n) + (RPg(n) - InitalPressure(n))*CM(n)
                  endif
            Endif
      enddo
      

      End Subroutine uMechaHT
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

      
!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 VIII. Coupling from Chemical to Hydraulic and Thermal and Mechanics             !
!           Describe the effect of Chemical variables changes on the THM processes                !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Subroutine uChemoHTM
      use Global; implicit none
      
      End Subroutine uChemoHTM
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                            After Hydro-Thermo-Mechanical is finished                            !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 IX. Wrie the results                                                            !
!           Describe here the complementary commands before leaving the process                   !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!
      Subroutine uWriteF
      use Global; use uArrays;implicit none

      integer :: n
      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  IX.1. Results of joint elements                                          !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!

      if (NB(8).ne.0) then ! Containing joint elemants
            if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! Mecha prob
                  Do n = 1,NB(2)
                        if(ntyp(n).eq.5) then
                              ! joint element results mecha
                              write(201,*), Temps, n, EXY(n,1),EXY(n,2), SXY(n,1), SXY(n,2), Ep(n,1),   &
                              & Ep(n,2), Vinm(n,1)  
                              write(202,*) Temps, n, Aut, AUn, ATau, ASn, ADamage
                        endif
                  enddo

            elseif ((NB(12).eq.2).or.(NB(12).eq.4)) then  ! Hydro prob
                  Do n = 1,NB(2)
                        if(ntyp(n).eq.5) then
                              ! matrix element results hydro
                              write(203,*), Temps, n, RPg(n)
                        endif
                  enddo
            end if
      endif
      
      

      !-------------------------------------------------------------------------------------------!
      !                                                                                           !
      !                  IX.2. Results of rock matrix                                             !
      !                                                                                           !
      !-------------------------------------------------------------------------------------------!

      if ((NB(12).eq.1).or.(NB(12).eq.4)) then ! Mecha prob
            Do n = 1,NB(2)
                  if((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                        ! matrix element results mecha
                        write(101,*), Temps, n, EXY(n,1),EXY(n,2),EXY(n,3),EXY(n,4),               &
                        eVolumetric(n), SXY(n,1), SXY(n,2),SXY(n,3),SXY(n,4),meanStress(n)
                  endif
            enddo
            Do n=1,NB(1)
                  if ((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                        ! matrix nodal results mecha
                        write(102,*) Temps, n, RU(n,1),  RU(n,2) 
                  endif
            enddo
      elseif ((NB(12).eq.2).or.(NB(12).eq.4)) then  ! Hydro prob
            Do n = 1,NB(2)
                  if((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then
                        ! matrix element results hydro
                        write(103,*), Temps, n, RPg(n), vf(n)
                  endif
            enddo
            Do n=1,NB(1)
                        ! matrix nodal results hydro
                        write(104,*) Temps, n, RP(n) 
            enddo
      end if



      close(101)
      close(102)
      close(103)
      close(104)
      close(201)
      close(202)
      close(203)
      close(301)
      close(302)
      close(303)
      
      End Subroutine uWriteF
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!      


!=================================================================================================! 
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                             User Defined Models                                                 !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------! 
!=================================================================================================!



!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                        Subroutine for defining material's damage criterion                      !
!                   Describe here the complementary initialization commands                       !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uCalCrit(n,matn,modmatn,sxyl,f0)
            use Global; implicit none
            integer :: i, j, n, matn, modmatn
            double precision :: sxyl(NB(28)), f0
      !-----------------------------------------------  
      f0=0.D0; goto 1000
      !      if (natmatn.eq.91130) then
      !       f0=0.D0; Vinm(n,0)=f0
      !       goto 1000 
      !       endif
1000  Continue
      end Subroutine uCalCrit
!-------------------------------------------------------------------------------------------------! 
!=================================================================================================!

      !==============================================================
      !     Subroutine for defining material's elastic tensor
      !--------------------------------------------------------------
            Subroutine ucalCmat(n,matn,natmatn,C)
            use Global; implicit none
            character ch
            integer :: i, j, n, matn, natmatn
            double precision, dimension(2,2) :: kjr
            double precision C(NB(28),NB(28))
            double precision :: xfact,EY, XNU
            double precision :: ALPHA,D, XD, beta, sigR, sigPc, tau0, sigman0,tanfi,Coh,f0
      !-----------------------------------------------  
            goto 1000 
             if (natmatn.eq.91130) then
              C(1,1)=vmat(matn,1,1); C(2,2)=vmat(matn,1,2)
              C(1,2)=vmat(matn,1,3); C(2,1)=C(1,2)
              goto 1000 
             endif
             goto 1000
      !       print*, 'Error1 materiau; element, mat, natmat: ', n,matn, natmatn
      !       stop
      !510   continue ; ! Linear and isotropic elastic material	   
      !       EY=vmat(mat(n),1,1); XNU=vmat(mat(n),1,2); xfact = EY/((1+XNU)*(1-2*XNU))
      !          C(1,1)=1-XNU;  C(1,2)=XNU;  C(1,3)=XNU; C(1,4) =0.D0
      !          C(2,1)= XNU;  C(2,2)=1-XNU;  C(2,3)=XNU; C(2,4) =0.D0
      !          C(3,1)= XNU;  C(3,2)= XNU;  C(3,3)=1-XNU; C(3,4) =0.D0
      !          C(4,1)= 0.D0;  C(4,2)= 0.D0;  C(4,3)=0.D0; C(4,4) =(1-2*XNU)/2.D0
      !          Do i=1, 4; do j=1, 4; C(i,j) = xfact*C(i,j); enddo; Enddo
      !      GOTO 1000
      !
            goto 1000 
      !
      1000    Continue
              end Subroutine ucalCmat
      !==============================================================
      !     Subroutine for defining plastic strain
      !--------------------------------------------------------------
            Subroutine udefplas(n,matn,natmatn,CC,SXYL,iter,fNmax)
      !	
            use Global; implicit none
            double precision :: CC(NB(28),NB(28)), SXYL(NB(28)), dFdsig(NB(28))
          ! general variables 
            integer :: i, j, n, ng, n1, n2, n3, n4, matn, natmatn, ia, ib, ndimC, iter
            double precision gama, Gkapa, rI1, rJ2, f0, h, xlam, fNmax, tanfi, sig, epsp
      !
             Goto 10000
            epsp=var(3)
      !-----------
             if (natmatn.eq.91110) then
              gama = vmat(matn,1,7); Gkapa = vmat(matn,1, 8); goto 520
             endif
             goto 10000
      520    Continue ; ! Elastoplasti Drucker-Prager Material
      !       given for example
             ndimC = 4
             rI1 = sxyl(1)+sxyl(2)+sxyl(3)
             rJ2 = sqrt((sxyl(1)**2+sxyl(2)**2+sxyl(3)**2+2*sxyl(4)**2-rI1*rI1/3)/2)
             f0 = gama*rI1 +rJ2-Gkapa; if (f0/Gkapa.gt.fNmax) fNmax = f0/Gkapa
             if (f0.gt.0.) then
              do i = 1, 3; dFdsig(i) = gama +(sxyl(i)-rI1/3)/2/rJ2; enddo
              dFdsig(4) = sxyl(4)/rJ2 ;! multiple par deux car epsil4=2*epsilonxy
              H = 0.
              Do ia=1,ndimC; do ib = 1, ndimC; H = H + dFdsig(ia)*CC(ia, ib)*dFdsig(ib); enddo; Enddo
              xlam = f0/H
              Do ia = 1, ndimC; dEp(n,ia) =xlam*dFdsig(ia); Enddo
             endif
             Goto 10000
      !
      10000  Continue
             end Subroutine udefplas
      !==============================================================
      !     Subroutine for defining viscous strain
      !--------------------------------------------------------------
            Subroutine udefvis(n,matn,natmatn)
            use Global; implicit none
            character(len=80) :: fich
            integer :: n, j, matn, natmatn
      !      double precision bt,bn,qn,alfa,xsitp,xsinp,xalt,xaln,s1,s2,deltat,xsimin,epsmin
      !
            goto 1000
      !---------------------------------------
      !
      1000  Continue
            end Subroutine udefvis
      !==============================================================
      !     Subroutine for defining Free Strain
      !--------------------------------------------------------------
             Subroutine uDefLib
            use Global; implicit none
            character(len=80) :: fich
            integer :: n, j, matn, natmatn, ndimC
            double precision const  
      !
      !     if (natmatn.eq.943000) then
      !  Calcul de defrmation de gonflement
      !      ndimC = nodgaus(ntyp(n),3)
      !      do j=1, ndimC; EpL(n,j)=0; enddo
      !  Deformation thermique, Axisymetrie
      !      E = vamt(matn,1) 
      !     acier :
      !!       const = 5.e-5     
      !!      EpL(n,1) = const*(-50.); EpL(n,2) = 0.; EpL(n,3) = EpL(n,1); EpL(n,4) = 0.
      !      EpL(n,2) = const*pe(n)*(1.2)
      !!      EpL(n,3) = 0.D0
      !!      Endif
      !      Enddo
      !
      1000  Continue
            end Subroutine uDefLib
      !==============================================================
      !     Subroutine for defining volume forces
      !--------------------------------------------------------------
            Subroutine uForVol(n,matn,natmatn)  
            use Global; implicit none
            character(len=80) :: fich
            integer :: n, ia, j, ng, matn, natmatn
            double precision coefBiot, p, yg
      !
      !--- lecture connectivite
      !
             goto 1000
             print*, 'Calculation of User volume forces'  
      1000  Continue
            end Subroutine uForVol
      !
      !==============================================================
      !     Subroutine for defining Initial Stresses
      !--------------------------------------------------------------
            Subroutine uSigInt
            use Global; implicit none
      !
      !      SXYi(:,:) = 0.D0
      !
      !5000  Continue
            end Subroutine uSigInt
      !
      !================================================================================
      !     Subroutine for defining Hydraulic diffusion parameters
      !--------------------------------------------------------------------------------
            Subroutine ucalvkc(n,matn,natmatn,VVK,vcinfn)
            use Global; implicit none 
            double precision VVK(2,2), vcinfn      
            integer :: n,ii, matn,natmatn
      !
      
      !     if (natmatn.eq.95050) then
      !      VVK(1,1)=1.;  VVK(2,2)=1.;  VVK(1,2)=0.;  VVK(2,1)=0.;
      !      Vcinfn=0.
      !     endif
            end Subroutine ucalvkc
      !
      !================================================================================
      !     Subroutine for defining Thermal diffusion parameters
      !--------------------------------------------------------------------------------
            Subroutine ucalqkc1(n,matn,natmatn,qKTh,qcinfn)
            use Global; implicit none 
            double precision qKTh(2,2), qcinfn      
            integer :: n,ii, matn,natmatn
      !
      !     if (natmatn.eq.95050) then
      !      qKTh(1,1)=1.;  qKTh(2,2)=1.;  qKTh(1,2)=0.;  qKTh(2,1)=0.;
      !      qcinfn=0.
      !     endif
            end Subroutine ucalqkc1
      !
      !================================================================================
      !     Subroutine for defining chemical process parameters
      !--------------------------------------------------------------------------------
            Subroutine uDphyC
            use Global; implicit none
            integer n, matn
      !
      !   Debut intervention--------------
             Do n=1,NB(2); matn=mat(n)   
      !       if (natmat(mtn,1,1).eq.931100) goto 93100
      !       if (natmat(mtn,1,1).eq.921300) goto 92100
             goto 1000
      !       
      92100 Continue; ! Element joint
             Difc(n,1,1)=1.0;Difc(n,2,2)=Difc(n,1,1);Difc(n,1,2)=0.0;Difc(n,2,1)=Difc(n,1,2);srcC(n)=1.0; asrcC(n)= 0.1
             PhyC(n) = 1; !  PhyC(n) = vmat(mtn,3,11)
             goto 1000
      93100 Continue; ! Element massif
      !       Difc(n,1,1)=1.0;Difc(n,2,2)=2.0;Difc(n,1,2)=0.0;Difc(n,2,1)=Difc(n,1,2);srcC(n)=1.0; asrcC(n)= 0.1
      !        PhyC(n) = vmat(mtn,6)
            goto 1000
      !     
      1000   Continue
            Enddo
            end Subroutine uDphyC
      !================================================================================

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 X. Subroutine to create time-evolution plots                                    !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uCurve(m)
      use Global; implicit none
      character(len=90) :: Fich1
      double precision :: ux, uy
      integer m ,i, nbdl, n1, n2
!
!    ----------------   m=1    Opening the file et writing title line      
      if (m.gt.1) goto 105  
      Fich1 = trim(foldername)//"\Ucurve.dat"
      open (unit=19, file=Fich1, status='replace')
      write(19,*) 't, Ux(n), Uy(n) ,   n=',nb(35)
!      call ulit_front
105   continue
!    ----------------  m=2   Writing for the time t   
       n1 = NB(35)
       ux= RU(n1,1);uy=RU(n1,2)
     write(19,305) Temps, ux, uy
!    ----------------  m=3   Closing the file
      if (m.eq.3) close(19)
305   Format (15(1X,E14.7))
1000  Continue
      end Subroutine uCurve
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
!                 XI. Subroutine for completing the post-process output file                      !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!

      Subroutine uWriteGID
      use Global; use uArrays; implicit none
!      double precision sxynod(NB(1),4)
!      double precision EPnod(NB(1)),EPval
!      double precision nx, ny, sigL
      integer :: i, n, nj, j , n1, n2, kj, njb
!
!      Goto 1000
!------------
!      if ((NB(8)+NB(9)+NB(65)).eq.0) goto 150
!      write(17,*) 'Result "Joint Stress (E)" "analysis name" ', Temps, '  Vector OnGausspoints "joints"'
!      write(17,*) 'ComponentNames "Tau","Sn", "|Tau|"'
!      write(17,*) 'Values'
!      DO n=1,NB(2)
!      If ((ntyp(n).eq.5).or.(ntyp(n).eq.6).or.(ntyp(n).eq.8)) then
!       write(17,*) n, SXY(n,1), SXY(n,2), abs(SXY(n,1))
!      Endif
!      ENDDO
!      write(17,*) 'end values'
!150   Continue
!-----------------------
      if ((NB(12).eq.1).or.(NB(12).eq.4).or.(NB(12).eq.6)) then ! contains mechanical problem 
         if ((NB(8)+NB(9)+NB(65)).eq.0) goto 150   ! if not containing joint/cable/bolt elements
            write(17,*) 'Result "Joint displacement (E)" "analysis name" ', Temps, '  Vector OnGausspoints "joints"'
            write(17,*) 'ComponentNames "Ut","Un", "|Ut|"'
            write(17,*) 'Values'
            DO n=1,NB(2)
               If ((ntyp(n).eq.5).or.(ntyp(n).eq.6).or.(ntyp(n).eq.8)) then ! joint/cable/bolt elements
                  write(17,*) n, EXY(n,1), EXY(n,2), abs(EXY(n,1))
               Endif
            ENDDO
          write(17,*) 'end values'
150      Continue
      endif
!------------
!        if ((NB(8)+NB(9)+NB(65)).eq.0) goto 150
!        write(17,*) 'Result "Joint Displacement (E)" "analysis name" ', Temps, '  Vector OnGausspoints "joints"'
!        write(17,*) 'ComponentNames "Ut","Un","|Ut|"'
!        write(17,*) 'Values'
!        DO n=1,NB(2)
!        If ((ntyp(n).eq.5).or.(ntyp(n).eq.6).or.(ntyp(n).eq.8)) then
!         write(17,*) n, EXY(n,1), EXY(n,2), abs(EXY(n,1))
!        Endif
!        ENDDO
!        write(17,*) 'end values'
! 150    Continue
! !------------
!      If ((NB(5)+NB(9)+NB(61)+NB(65)).ne.0) then      
!      write(17,*) 'Result "Bar Force Vector (E)" "analysis name" ', Temps, '  Vector OnGaussPoints "Bars"'
!      write(17,*) 'ComponentNames "SLx","SLy"'
!      write(17,*) 'Values'
!      do n=1,NB(2)    
!       n1 = konec(n,1); n2 = konec(n,2)
!       nx = -(x(n2,2)-x(n1,2))/S(n); ny = (x(n2,1)-x(n1,1))/S(n)
!      if ((ntyp(n).eq.2).or.(ntyp(n).eq.7)) write(17,*) n, SXY(n,1)*nx,SXY(n,1)*ny
!      if ((ntyp(n).eq.6).or.(ntyp(n).eq.8)) write(17,*) n, SXY(n,3)*nx,SXY(n,3)*ny
!      end do
!      write(17,*) 'end values'
!      Endif       
!------------
      If ((NB(12).eq.2).or.(NB(12).eq.4)) then   ! Contains hydraulic problem
         write(17,*) 'Result "Fluid Pressures Variation(N)" "analysis name" ',Temps,' Scalar OnNodes'
         write(17,*) 'ComponentNames "Fluid Pressure"'; write(17,*) 'Values'
         do i=1,NB(1); write(17,*) i, RP(i)-nodalInitialPressure(i); end do
         write(17,*) 'end values'
      endif

!------------
!      write(17,*) 'Result "Fluid Velocity (N)" "analysis name" ',Temps,' Vector OnNodes'
!      write(17,*) 'ComponentNames "Vx","Vy"'; write(17,*) 'Values'
!      do i=1,NB(1); write(17,*) i, VXY(i,1), VXY(i,2); end do
!      write(17,*) 'end values'
!------------
!      write(17,*) 'Result "Fluid Velocities (E)" "analysis name" ', Temps, '  Vector OnGaussPoints "Qtriang"'
!      write(17,*) 'ComponentNames "Vx", "Vy"'; write(17,*) 'Values'
!      do n=1,NB(2)
!       if ((ntyp(n).eq.3).or.(ntyp(n).eq.4)) then 
!        write(17,*) n,VXYL(n,1),VXYL(n,2)
!       endif 
!      end do
!      write(17,*) 'end values'
!------------
!      IF (NB(8).NE.0) THEN
!       write(17,*) 'Result "Joint Flux (E)" "analysis name" ', Temps, '  Vector OnGaussPoints "joints"'
!       write(17,*) 'ComponentNames "Qx", "Qy"'; write(17,*) 'Values'
!       do n=1, NB(2)
!        if (ntyp(n).eq.5) then
!         write(17,*) n, VXYL(n,1),VXYL(n,2)
!        endif
!       end do
!       write(17,*) 'end values'
!      ENDIF
!------------
!1000   Continue
      end Subroutine uWriteGID
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

!=================================================================================================!
!-------------------------------------------------------------------------------------------------!
!                                                                                                 !
! XII. Subroutine to show how works the Subroutine Vsort for sorting a table in increasing order  !
!                                                                                                 !
!-------------------------------------------------------------------------------------------------!


!    Minax(n) provides the initial position of the n-th rank if sorted in increasing order
      Subroutine ExampleSort
      implicit none
      integer, allocatable :: minax(:)
      double precision, allocatable :: sc(:)
      integer :: nbv,i,j,k,n
!
      nbv=5   ! Array dimension
      ! if (allocated(sc)) deallocate (sc, minax)
      allocate (sc(nbv), minax(nbv))
      sc(1)=1.3; sc(2)=1.5; sc(3)=0.3; sc(4)=2.3; sc(5)=2.0
      print*, (sc(i), i=1,nbv)    !  -->  1.3, 1.5, 0.3, 2.3, 2.0
      call VSort(nbv,sc,minax)
      print*, (sc(minax(i)), i=1,nbv)    !  -->  0.3, 1.3, 1.5, 2.0, 2.3
      END subroutine ExampleSort
!-------------------------------------------------------------------------------------------------!
!=================================================================================================!

      End Module User
!==================================================================================================
!==================================================================================================
!==================================================================================================