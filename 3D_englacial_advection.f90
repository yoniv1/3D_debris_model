subroutine  debris_cover_main_englacial(H, HT, mask, G, G_obs, h_debris, margin, VX, VY, VZAGE, meltout_debris, inputlocation_debris, ui, vi, wi, za, ti, ti_new,actual_mass_hist_eng, expected_mass_hist_eng, expected_debris_mass_eng,c_debris_out,DHDT_yearly,ti_add,t_s,pert_status,H_old,TIME,sum_mass_meltout,mass_in_old)

      use        PARAMETER_LIST

      implicit none

!--   Parameter List in                                                                                                                                                                  
      real, dimension(NX,NY)               :: HT,DHDT_yearly, H_old
      double precision                     :: TIME
      double precision, dimension(NX,NY,3) :: H
      real,dimension(NX,NY)                :: mask
      real,dimension(NX,NY)                :: margin 
      real,dimension(NX,NY)                :: h_debris, c_debris_out,count_debris,cells_debris
      real,dimension(NX,NY)                :: meltout_debris, inoutdebris,mass_in,mass_out,loop_number,mass_in_old       
      real,dimension(NX,NY)                :: inputlocation_debris, t_s,sum_C,ti_corr,th_per_piece,meltout_pieces,sum_c_meltout,sum_mass_meltout
      double precision,dimension(NX,NY,NZ) :: VX,VY,VZAGE, ui, vi, wi, ti, ti_new, tinit, advx_eng, advy_eng, advz_eng, advel_x, advel_y, advel_z, ti_diff, vel_diff, actual_mass_debris_eng, sink_debris, new_debris_mass_eng, expected_debris_mass_eng, old_debris_mass_eng, term1_debris_eng_x, term1_debris_eng_y, term1_debris_eng_z, new_C,ti_add,ti_add2,ti_add3,c_meltout
      real,dimension(NZ)                   :: za, za_eng, dzeta
      real,dimension(1)                    :: actual_mass_hist_eng, expected_mass_hist_eng, mass_ratio_hist_eng, expected_mass_before_eng

!--   Parameter List out                                                                                                                                                                 
      real,dimension(NX,NY) :: G,G_obs

!--   Variables                                                                                                                                                                          
      integer               :: i,j,k,yr_deng, it_deng,it_smolar_eng,c,pert_status,W
      real                  :: startTime, stopTime
      real, allocatable     :: ti_corr2(:,:)

!-----------------------------------------------------------------------                              
! Initialization of velocities                                                                          
!-----------------------------------------------------------------------                    

! Velocities                                                      
! ui = x-component
! vi = y-component
! wi = z-component

! Vertical grid spacing                                                                                                    
! za_eng

!-----------------------------------------------------------------------        
! RUN THE ENGLACIAL DEBRIS CONCENTRATION MODEL
!-----------------------------------------------------------------------    

write(*,*),'starting the loop for the englacial debris concentration...'

yr_deng = 1

do while (yr_deng.le.numberofyears_eng)

  it_deng = 1

  do while (it_deng.le.(nyears_deng/dt_eng))

!-----------------------------------------------------------------------                        
! DEBRIS INPUT AND INITIALIZATION                         
!-----------------------------------------------------------------------     

  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! ccccccccccc INITIALIZE MATRICES cccccccccccccc  
  ! cccccccccccccccccccccccccccccccccccccccccccccc         

  ! Initialize

      do J=1,NY 
         do I=1,NX      
            inoutdebris(I,J) = 0.
         end do
       end do

  ! Initialize                                                                                                                                                                 

   if(it_deng.eq.1)then

      do J=1,NY
         do I=1,NX
            sum_c_meltout(I,J) = 0.
         end do
       end do

  ! Initialize      

     do J=1,NY
        do I=1,NX
          do K=1,NZ
             c_meltout(I,J,K) = 0
          end do
         end do
       end do
       
   endif

  ! Initialize                                                                                                   

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              ti_add2(I,J,K) = 0
           end do
         end do
      end do

  ! Zeta coordinate   
  
      do K=1,NZ-1
          dzeta(K)=za_eng(K+1)-za_eng(K)
      enddo
      dzeta(NZ) = dzeta(NZ-1)

  ! cccccccccccccccccccccccccccccccccccccccccccccc
  ! ccccccccccccc DEFINE THE INPUT ccccccccccccccc                       
  ! cccccccccccccccccccccccccccccccccccccccccccccc               

  ! Debris input

      do J=1,NY
         do I=1,NX

            if (I.eq.14.and.J.ge.83.and.J.le.95)then
                inputlocation_debris(I,J) = 1
           end if
            if (I.eq.15.and.J.ge.83.and.J.le.95)then
                inputlocation_debris(I,J) = 1
           end if
            if (I.eq.16.and.J.ge.83.and.J.le.95)then
                inputlocation_debris(I,J) = 1
           end if
            if (I.eq.17.and.J.ge.83.and.J.le.95)then
                inputlocation_debris(I,J) = 1
           end if

           if (inputlocation_debris(I,J).gt.0)then
              inoutdebris(I,J) = depositionrate_debris_acc
           end if 

         end do
       end do

   ! Recalculate to concentrations      

      do J=1,NY
         do I=1,NX
           if(H(I,J,3).gt.0)then
                tinit(I,J,NZ) = (inoutdebris(I,J)*rho_debris*nz*dt_eng)/(H(I,J,3))
           endif
         enddo
      enddo

  ! Set top boundary (input) condition

    if(TIME.gt.tdebris.and.TIME.lt.tdebris+1)then
      do J=1,NY
         do I=1,NX
                t_s(I,J) = tinit(I,J,NZ)
                mass_in_old(I,J) = tinit(I,J,NZ)*dx*dx*(H(I,J,3)/NZ)
         end do
      end do
   endif

  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! cccccc APPLY INPUT BOUNDARY CONDITION cccccccc  
  ! cccccccccccccccccccccccccccccccccccccccccccccc    

  ! Initialization of ti and ti_new

      do J=1,NY
         do I=1,NX
                ti(I,J,NZ) = t_s(I,J)
	        ti_new(I,J,NZ) = t_s(I,J)
         end do
      end do

   ! Mass input   

      do J=1,NY
         do I=1,NX
              mass_in(I,J) = t_s(I,J)*dx*dx*(H(I,J,3)/NZ)
         end do
      end do

  ! cccccccccccccccccccccccccccccccccccccccccccccc 
  ! ccccccccccc CORRECTION FOR DHDT cccccccccccccc     
  ! cccccccccccccccccccccccccccccccccccccccccccccc 

   if(it_deng.eq.1.and.TIME.gt.tdebris+1)then

  ! Initialize

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              ti_add2(I,J,K) = 0
           end do
         end do
      end do

   ! Take into account dHdt

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              if(ti(I,J,K).gt.0)then
                ti_add2(I,J,K) = -ti(I,J,K) * (NZ/H(I,J,3)) * (-DHDT_yearly(I,J)/NZ)
              else
                ti_add2(I,J,K) = 0
           end if
         end do
       end do
     end do

  ! Delete false pixels

      do J=1,NY
         do I=1,NX
           do K=1,NZ
              if(H(I,J,3).eq.0)then
                ti_add2(I,J,K) = 0
              else
                ti_add2(I,J,K) = ti_add2(I,J,K)
           end if
         end do
       end do
     end do

   ! Adjust debris concentration accordingly

      ti = ti + ti_add2
      ti_new = ti_new + ti_add2

       do J=1,NY
         do I=1,NX
           if(wi(I,J,NZ).lt.0)then
               t_s(I,J) = ti(I,J,NZ)
           endif
         end do
       end do

   ! Mass input     

      do J=1,NY
         do I=1,NX
              mass_in(I,J) = t_s(I,J)*dx*dx*(H(I,J,3)/NZ)
         end do
      end do

   ! cccccccccccccccccccccccccccccccccccccccccccccc 
   ! cccccc MELT OUT BOUNDARY CONDITION  cccccccccc  
   ! cccccccccccccccccccccccccccccccccccccccccccccc    

    do J=1,NY
      do I=1,NX

        th_per_piece(I,J) = H(I,J,3)/NZ
        if((H(I,J,3).gt.0))then
           meltout_pieces(I,J) = -G(I,J) / th_per_piece(I,J)
        endif
        loop_number(I,J) = ceiling(meltout_pieces(I,J))
        if((loop_number(I,J).gt.NZ-1))then
            loop_number(I,J) = NZ-1
        endif
        if((H(I,J,3).gt.0).and.(G(I,J).lt.0).and.(ti(I,J,NZ-1).gt.0))then
         do W=1,loop_number(I,J)
           if((meltout_pieces(I,J)-W).gt.0)then
              c_meltout(I,J,W+1) = ti(I,J,NZ-W)
           elseif((meltout_pieces(I,J)-W).lt.0.and.(loop_number(I,J)).eq.1)then
              c_meltout(I,J,W+1) = ti(I,J,NZ-W)*abs(meltout_pieces(I,J))
           elseif((meltout_pieces(I,J)-W).lt.0.and.(loop_number(I,J)).gt.1)then
              c_meltout(I,J,W+1) = ti(I,J,NZ-W)*abs(1-meltout_pieces(I,J))
           endif
         end do
        endif

      end do
    end do

    do J=1,NY
     do I=1,NX
        if((H(I,J,3).gt.0).and.(G(I,J).lt.0).and.(ti(I,J,NZ-1).gt.0))then
           do K=1,NZ
             c_meltout(I,J,K) = c_meltout(I,J,K)
           end do
        else
             c_meltout(I,J,K) = 0
        end if
     end do
    end do

  ! Summation

     do J=1,NY
        do I=1,NX
          if((H(I,J,3).gt.0).and.(G(I,J).lt.0).and.(ti(I,J,NZ-1).gt.0))then
            do K=1,NZ
              if (K.eq.1)then
                 sum_c_meltout(I,J) = c_meltout(I,J,1)
              elseif (K.gt.1)then
                 sum_c_meltout(I,J) = sum_c_meltout(I,J) + c_meltout(I,J,K)
              endif
            enddo
          else
              sum_c_meltout(I,J) = 0
          endif
       enddo
    enddo

   endif !endif it_deng = 1

      !-----------------------------------------------------------------------    
      ! STEP 1: DEBRIS CONCENTRATION CALCULATION
      !-----------------------------------------------------------------------                                                                
 
      do J=1,NY 
         do I=1,NX
	   do K=1,NZ                                                                                                 
    
	   ti_new(I,J,K) = ti(I,J,K)

	   IF((H(I,J,3).ge.0).and.(I.gt.1).and.(J.gt.1).and.(K.gt.1).and.(I.lt.NX).and.(J.lt.NY).and.(K.lt.NZ))then

            ! cccccccccccccccccccccccccccccccccccccccccccccc
            ! cccccccccc START ADVECTION SCHEME cccccccccccc
            ! cccccccccccccccccccccccccccccccccccccccccccccc

            ! cccccccccccccccccccccccccccccccccccccccccccccc 
            ! cccccccccccccccc Advection in X cccccccccccccc   
            ! cccccccccccccccccccccccccccccccccccccccccccccc                                                             
 
            if((I.eq.2).or.(K.eq.2).or.(I.eq.NX-1).or.(K.eq.NZ-1))then
             if (ui(I,J,K).lt.0)then
                advx_eng(I,J,K) = -ui(I,J,K) * (ti(I+1,J,K)-ti(I,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) * (1 - za_eng(K)) * ui(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * (ui(I+1,J,K)-ui(I,J,K)) /deltax_d - (1/H(I,J,3)) * (ui(I,J,K+1)-ui(I,J,K))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d)
             else if (ui(I,J,K).ge.0)then
                advx_eng(I,J,K) = -ui(I,J,K) * (ti(I,J,K)-ti(I-1,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * ui(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d) - ti(I,J,K) * (ui(I,J,K)-ui(I-1,J,K)) /deltax_d - (1/H(I,J,3)) * (ui(I,J,K)-ui(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d)
             endif
            else
             if (ui(I,J,K).ge.0)then
                advx_eng(I,J,K) = -ui(I,J,K) * (((3*ti(I,J,K))-(4*ti(I-1,J,K))+ti(I-2,J,K))/(2*deltax_d)) - (1/H(I,J,3)) * (((3*ti(I,J,K))-(4*ti(I,J,K-1))+ti(I,J,K-2))/(2*dzeta(k))) * (1 - za_eng(k)) * ui(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d) - ti(I,J,K) * (((3*ui(I,J,K))-(4*ui(I-1,J,K))+ui(I-2,J,K))/(2*deltax_d)) - (1/H(I,J,3)) * (((3*ui(I,J,K))-(4*ui(i,J,K-1))+ui(I,J,K-2))/(2*dzeta(k))) * (1 - za_eng(k)) * ti(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d)
             else if (ui(I,J,K).lt.0)then
                advx_eng(I,J,K) = -ui(I,J,K) * ((-ti(I+2,J,K)+(4*ti(I+1,J,K))-(3*ti(I,J,K)))/(2*deltax_d)) - (1/H(I,J,3)) * (-ti(I,J,K+2)+(4*ti(I,J,K+1))-(3*ti(I,J,K)))/(2*dzeta(k)) * (1 - za_eng(k)) * ui(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * ((-ui(I+2,J,K)+(4*ui(I+1,J,K))-(3*ui(I,J,K)))/(2*deltax_d)) - (1/H(I,J,3)) * (-ui(I,J,K+2)+(4*ui(I,J,K+1))-(3*ui(I,J,K)))/(2*dzeta(k)) * (1 - za_eng(k)) * ti(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d)
             endif
	          endif

            ! cccccccccccccccccccccccccccccccccccccccccccccc               
            ! cccccccccccccccc Advection in Y cccccccccccccc                    
            ! cccccccccccccccccccccccccccccccccccccccccccccc  

            if((J.eq.2).or.(K.eq.2).or.(J.eq.NY-1).or.(K.eq.NZ-1))then                                                                                    
             if (vi(I,J,K).lt.0)then
                advy_eng(I,J,K) = -vi(I,J,K) * (ti(I,J+1,K)-ti(I,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) * (1 - za_eng(K)) * vi(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * (vi(I,J+1,K)-vi(I,J,K)) /deltax_d - (1/H(I,J,3)) * (vi(I,J,K+1)-vi(I,J,K))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d)
             else if (vi(I,J,K).ge.0)then
                advy_eng(I,J,K) = -vi(I,J,K) * (ti(I,J,K)-ti(I,J-1,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * vi(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d) - ti(I,J,K) * (vi(I,J,K)-vi(I,J-1,K)) /deltax_d - (1/H(I,J,3)) * (vi(I,J,K)-vi(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d)
            end if
           else
             if (vi(I,J,K).ge.0)then
                advy_eng(I,J,K) = -vi(I,J,K) * (((3*ti(I,J,K))-(4*ti(I,J-1,K))+ti(I,J-2,K))/(2*deltax_d)) - (1/H(I,J,3)) * (((3*ti(I,J,K))-(4*ti(I,J,K-1))+ti(I,J,K-2))/(2*dzeta(k))) * (1 - za_eng(k)) * vi(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d) - ti(I,J,K) * (((3*vi(I,J,K))-(4*vi(I,J-1,K))+vi(I,J-2,K))/(2*deltax_d)) - (1/H(I,J,3)) * (((3*vi(I,J,K))-(4*vi(i,J,K-1))+vi(I,J,K-2))/(2*dzeta(k))) * (1 - za_eng(k)) * ti(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d)
             else if (vi(I,J,K).lt.0)then
                advy_eng(I,J,K) = -vi(I,J,K) * ((-ti(I,J+2,K)+(4*ti(I,J+1,K))-(3*ti(I,J,K)))/(2*deltax_d)) - (1/H(I,J,3)) * (-ti(I,J,K+2)+(4*ti(I,J,K+1))-(3*ti(I,J,K)))/(2*dzeta(k)) * (1 - za_eng(k)) * vi(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * ((-vi(I,J+2,K)+(4*vi(I,J+1,K))-(3*vi(I,J,K)))/(2*deltax_d)) - (1/H(I,J,3)) * (-vi(I,J,K+2)+(4*vi(I,J,K+1))-(3*vi(I,J,K)))/(2*dzeta(k)) * (1 - za_eng(k)) * ti(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d)
             endif
           endif

            ! cccccccccccccccccccccccccccccccccccccccccccccc  
            ! cccccccccccccccc Advection in Z cccccccccccccc   
            ! cccccccccccccccccccccccccccccccccccccccccccccc              

            if((K.eq.2).or.(K.eq.NZ-1))then                                                                                    
             if (wi(I,J,K).lt.0)then
                advz_eng(I,J,K) = -wi(I,J,K) * (-1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) - ti(I,J,K) * (-1/H(I,J,3)) * (wi(I,J,K+1)-wi(I,J,K))/dzeta(K)
             else if (wi(I,J,K).ge.0)then
                advz_eng(I,J,K) = -wi(I,J,K) * (-1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) - ti(I,J,K) * (-1/H(I,J,3)) * (wi(I,J,K)-wi(I,J,K-1))/dzeta(K)
             end if
	          else
             if (wi(I,J,K).ge.0)then
                advz_eng(I,J,K) = -wi(I,J,K) * (-1/H(I,J,3)) * (((3*ti(I,J,K))-(4*ti(I,J,K-1))+ti(I,J,K-2))/(2*dzeta(k))) - ti(I,J,K) * (-1/H(I,J,3)) * (((3*wi(I,J,K))-(4*wi(I,J,K-1))+wi(I,J,K-2))/(2*dzeta(k)))
             else if (wi(I,J,K).lt.0)then
                advz_eng(I,J,K) = -wi(I,J,K) * (-1/H(I,J,3)) * (-ti(I,J,K+2)+(4*ti(I,J,K+1))-(3*ti(I,J,K)))/(2*dzeta(k)) - ti(I,J,K) * (-1/H(I,J,3)) * ((-wi(I,J,K+2)+(4*wi(I,J,K+1))-(3*wi(I,J,K)))/(2*dzeta(k)))
             end if
            endif

           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! cccccccccccccccccccc END OF ADVECTION SCHEME ccccccccccccccccccc
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
            ! Calculate new debris concentration
  
            ti(I,J,K) = ti_new(I,J,K) + dt_eng * (advx_eng(I,J,K) + advy_eng(I,J,K) + advz_eng(I,J,K))

            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
            ! cccccccccccccccc UPDATE BOUNDARY CONDITIONS cccccccccccccccccccc 
            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                      
              ! Boundary conditions on top                                                                                          
               if (wi(I,J,K).le.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
               endif
               if (wi(I,J,K).gt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               end if
              ! Boundary conditions on sides                                                                                                       
                ti_new(NX,J,K) = 0.0
                ti(NX,J,K) = 0.0
                ti_new(I,NY,K) = 0.0
                ti(I,NY,K) = 0.0
                ti_new(1,J,K) = 0.0
                ti(1,J,K) = 0.0
                ti_new(I,1,K) = 0.0
                ti(I,1,K) = 0.0
              ! Boundary conditions at bottom                                                             
                ti(I,J,1) = 0.0
                ti_new(I,J,1) = 0.0

            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! cccccccccccccccc ADJUST FALSE CONCENTRATION cccccccccccccccccccc  
            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                     
               if (ti(I,J,K).lt.1e-3)then
                    ti(I,J,K) = 0.
	                  ti_new(I,J,K) = 0.
               end if

               if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(ti(I,J,K).gt.0))then
                    ti(I,J,K) = 0.
		                ti_new(I,J,K) = 0.
               end if

               if (H(I,J,3).eq.0)then
                    ti(I,J,K) = 0.
		                ti_new(I,J,K) = 0.
               end if
               
               if (isnan(ti(I,J,K)))then
                  ti(I,J,K) = 0
                  ti_new(I,J,K) = 0
               end if
           end if
         end do
      end do
    end do

! Update boundary condition at the top             
 
      do J=1,NY
         do I=1,NX
           do K=1,NZ
               if (wi(I,J,NZ).ge.0)then
                ti(I,J,NZ) = ti(I,J,NZ-1)
                ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               else
                ti(I,J,NZ) = ti(I,J,NZ)
                ti_new(I,J,NZ) = ti_new(I,J,NZ)
               endif
            end do
        end do
     end do

 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! ccccccccccccccccccccccc SAVE FOR MASS CONSERVATION cccccccccccccccccccccc
 ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do J=2,NY-1
         do I=2,NX-1
           do K=2,NZ-1

            actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
            if (wi(I,J,NZ-1).gt.0)then
               sink_debris(I,J,1)=sum_c_meltout(I,J)
            else
               sink_debris(I,J,1)=0.
            endif
            new_debris_mass_eng(I,J,NZ-1) = (H(I,J,3)/nz)*(deltax_d*deltax_d*t_s(I,J))*dt_eng - (H(I,J,3)/nz)*(deltax_d*deltax_d*sink_debris(I,J,1))*dt_eng
	          expected_debris_mass_eng(I,J,K) = expected_debris_mass_eng(I,J,K) + new_debris_mass_eng(I,J,K)

            end do
        end do
     end do

      !-----------------------------------------------------------------------    
      ! STEP 2: ANTI-DIFFUSION VELOCITY CALCULATION
      !-----------------------------------------------------------------------                                                                

        ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! cccccccccccccc CALCULATE ANTI-DIFFUSION VELOCITIES AFTER SMOLARKIEWICZ cccccccccccc
        ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do J=1,NY 
         do I=1,NX
	       do K=1,NZ                                                                                                 
	           IF((H(I,J,3).ge.0).and.(I.gt.1).and.(J.gt.1).and.(K.gt.1).and.(I.lt.NX).and.(J.lt.NY).and.(K.le.NZ))then
                     if((ti(I,J,K).gt.0))then

                     ! X-DIRECTION 
                       advel_x(I,J,K) = -(0.5*((abs(ui(I,J,K))*deltax_d - (dt_eng*((ui(I,J,K))**2))) / ti(I,J,K))) * (((ti(I,J,K)+ti(I+1,J,K))/2)-((ti(I,J,K)+ti(I-1,J,K))/2))/deltax_d - (((((dt_eng*(wi(I,J,K))/((-H(I,J,3)*dzeta(k))))*(dt_eng*(ui(I,J,K))/deltax_d)*deltax_d*(-H(I,J,3)*dzeta(k))))/(2*dt_eng)) / (ti(I,J,K)) * ((((ti(I,J,K)+ti(I,J,K+1))/2)-((ti(I,J,K)+ti(I,J,K-1))/2))/(-H(I,J,3)*dzeta(k)))) - (((vi(I,J,K)*ui(I,J,K)*dt_eng)/2) / ti(I,J,K) * ((((ti(I,J,K)+ti(I,J+1,K))/2)-((ti(I,J,K)+ti(I,J-1,K))/2))/deltax_d))
                     ! Y-DIRECTION 
                       advel_y(I,J,K) = -(0.5*((abs(vi(I,J,K))*deltax_d - (dt_eng*((vi(I,J,K))**2))) / ti(I,J,K))) * (((ti(I,J,K)+ti(I,J+1,K))/2)-((ti(I,J,K)+ti(I,J-1,K))/2))/deltax_d - (((((dt_eng*(wi(I,J,K))/((-H(I,J,3)*dzeta(k))))*(dt_eng*(vi(I,J,K))/deltax_d)*deltax_d*(-H(I,J,3)*dzeta(k))))/(2*dt_eng)) / (ti(I,J,K)) * ((((ti(I,J,K)+ti(I,J,K+1))/2)-((ti(I,J,K)+ti(I,J,K-1))/2))/(-H(I,J,3)*dzeta(k)))) - (((ui(I,J,K)*vi(I,J,K)*dt_eng)/2) / ti(I,J,K) * ((((ti(I,J,K)+ti(I+1,J,K))/2)-((ti(I,J,K)+ti(I-1,J,K))/2))/deltax_d))
                     ! Z-DIRECTION
                       advel_z(I,J,K) = -(0.5*((abs(wi(I,J,K))*(-H(I,J,3)*dzeta(k)) - (dt_eng*((wi(I,J,K))**2))) / ti(I,J,K))) * (((ti(I,J,K)+ti(I,J,K+1))/2)-((ti(I,J,K)+ti(I,J,K-1))/2))/(-H(I,J,3)*dzeta(k)) - ((((((dt_eng*(wi(I,J,K))/((-H(I,J,3)*dzeta(k))))*(dt_eng*(ui(I,J,K))/deltax_d)*deltax_d*(-H(I,J,3)*dzeta(k))))/(2*dt_eng)) / (ti(I,J,K)) * (((ti(I,J,K)+ti(I+1,J,K))/2)-((ti(I,J,K)+ti(I-1,J,K))/2))/deltax_d)) - ((((((dt_eng*(wi(I,J,K)))/(-H(I,J,3)*dzeta(k)))*(dt_eng*(vi(I,J,K))/deltax_d))*deltax_d*(-H(I,J,3)*dzeta(k)))/(2*dt_eng)) / (ti(I,J,K)) * ((((ti(I,J,K)+ti(I,J+1,K))/2)-((ti(I,J,K)+ti(I,J-1,K))/2))/(deltax_d)))
                    
                     ! ANTI-DIFFUSION VELOCITY
                     vel_diff(I,J,K) = sqrt((advel_x(I,J,K)**2)+(advel_y(I,J,K)**2)+(advel_z(I,J,K)**2))

                     ! BOUNDARY CONDITIONS
                     if (wi(I,J,K).ge.0)then
                        advel_x(I,J,NZ) = advel_x(I,J,NZ-1)
                        advel_y(I,J,NZ) = advel_y(I,J,NZ-1)
                        advel_z(I,J,NZ) = advel_z(I,J,NZ-1)
                        vel_diff(I,J,NZ) = vel_diff(I,J,NZ-1)
                     endif
                     
                     ! DELETE HIGH VALUES
                     if (vel_diff(I,J,K).gt.10000)then
                        vel_diff(I,J,K) = 0
                     endif

                end if
             end if
         end do
      end do
    end do

	!-----------------------------------------------------------------------    
	! STEP 3: RECALCULATE ADVECTION EQUATION WITH THE SMOLARKIEWICZ SCHEME
	!-----------------------------------------------------------------------                                                                

        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ccccccc CALCULATE DEBRIS CONCENTRATION WITH NUMERICAL ANTI-DIFFUSION CORRECTION cccccccc
        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(it_deng.eq.(nyears_deng/dt_eng))then

        it_smolar_eng = 1	

         do while (it_smolar_eng.le.(it_smolar_engl_max))
                                                                                                                                                                      
           do J=1,NY 
            do I=1,NX
	            do K=1,NZ

	              ti_new(I,J,K) = ti(I,J,K)                                                                                                 

	              IF((H(I,J,3).ge.0).and.(I.gt.1).and.(J.gt.1).and.(K.ge.1).and.(I.lt.NX).and.(J.lt.NY).and.(K.le.NZ))then

                ! cccccccccccccccccccccccccccccccccccccccccccccc
                ! cccccccccc START ADVECTION SCHEME cccccccccccc
                ! cccccccccccccccccccccccccccccccccccccccccccccc

		! cccccccccccccccccccccccccccccccc
                ! cccccccccc X-DIRECTION ccccccccc
		! cccccccccccccccccccccccccccccccc

                  if (advel_x(I,J,K).lt.0)then
                     term1_debris_eng_x(I,J,K) = -advel_x(I,J,K) * (ti(I+1,J,K)-ti(I,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) * (1 - za_eng(K)) * advel_x(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * (advel_x(I+1,J,K)-advel_x(I,J,K)) /deltax_d - (1/H(I,J,3)) * (advel_x(I,J,K+1)-advel_x(I,J,K))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I+1,J,3)-H(I,J,3)) /deltax_d)
                  else if (advel_x(I,J,K).ge.0)then
                     term1_debris_eng_x(I,J,K) = -advel_x(I,J,K) * (ti(I,J,K)-ti(I-1,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * advel_x(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d) - ti(I,J,K) * (advel_x(I,J,K)-advel_x(I-1,J,K)) /deltax_d - (1/H(I,J,3)) * (advel_x(I,J,K)-advel_x(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J,3)-H(I-1,J,3)) /deltax_d)
                  endif

		! cccccccccccccccccccccccccccccccc
                ! cccccccccc Y-DIRECTION ccccccccc
	        ! cccccccccccccccccccccccccccccccc

                  if (advel_y(I,J,K).lt.0)then
                     term1_debris_eng_y(I,J,K) = -advel_y(I,J,K) * (ti(I,J+1,K)-ti(I,J,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) * (1 - za_eng(K)) * advel_y(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d) - ti(I,J,K) * (advel_y(I,J+1,K)-advel_y(I,J,K)) /deltax_d - (1/H(I,J,3)) * (advel_y(I,J,K+1)-advel_y(I,J,K))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J+1,3)-H(I,J,3)) /deltax_d)
                  else if (advel_y(I,J,K).ge.0)then
                     term1_debris_eng_y(I,J,K) = -advel_y(I,J,K) * (ti(I,J,K)-ti(I,J-1,K)) /deltax_d - (1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * advel_y(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d) - ti(I,J,K) * (advel_y(I,J,K)-advel_y(I,J-1,K)) /deltax_d - (1/H(I,J,3)) * (advel_y(I,J,K)-advel_y(I,J,K-1))/dzeta(K) * (1 - za_eng(K)) * ti(I,J,K) * ((H(I,J,3)-H(I,J-1,3)) /deltax_d)
                 end if

	        ! cccccccccccccccccccccccccccccccc
                ! cccccccccc Z-DIRECTION ccccccccc
		! cccccccccccccccccccccccccccccccc

                  if (advel_z(I,J,K).lt.0)then
                     term1_debris_eng_z(I,J,K) = -advel_z(I,J,K) * (-1/H(I,J,3)) * (ti(I,J,K+1)-ti(I,J,K))/dzeta(K) - ti(I,J,K) * (-1/H(I,J,3)) * (advel_z(I,J,K+1)-advel_z(I,J,K))/dzeta(K)
                  else if (advel_z(I,J,K).ge.0)then
                     term1_debris_eng_z(I,J,K) = -advel_z(I,J,K) * (-1/H(I,J,3)) * (ti(I,J,K)-ti(I,J,K-1))/dzeta(K) - ti(I,J,K) * (-1/H(I,J,3)) * (advel_z(I,J,K)-advel_z(I,J,K-1))/dzeta(K)
		              end if
                                                                         
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
               ! cccccccccccccccccccc BOUNDAARY CONDITIONS cccccccccccccccccccccc
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

                  if (wi(I,J,K).ge.0)then
                     term1_debris_eng_x(I,J,NZ) = term1_debris_eng_x(I,J,NZ-1)
                     term1_debris_eng_y(I,J,NZ) = term1_debris_eng_y(I,J,NZ-1)
                     term1_debris_eng_z(I,J,NZ) = term1_debris_eng_z(I,J,NZ-1)
                  endif

               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ! cccccccccccccccccccc END OF ADVECTION SCHEME ccccccccccccccccccc
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
                ! Calculate new debris concentration

                 if (wi(I,J,K).le.0)then
                    ti_diff(I,J,K) = deltat_smolar_engl*(term1_debris_eng_x(I,J,K) + term1_debris_eng_y(I,J,K) + term1_debris_eng_z(I,J,K))
                 elseif (wi(I,J,K).gt.0)then
                    ti_diff(I,J,K) = deltat_smolar_engl*(term1_debris_eng_x(I,J,K) + term1_debris_eng_y(I,J,K) + term1_debris_eng_z(I,J,K))
                 endif

                 if (ti_diff(I,J,K).lt.0)then
                    ti_diff(I,J,K)=0
                 end if
                 ti(I,J,K) = ti_new(I,J,K) - (ti_diff(I,J,K))

            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
            ! cccccccccccccccc UPDATE BOUNDARY CONDITIONS cccccccccccccccccccc   
            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
 
              ! Boundary conditions on top
               if (wi(I,J,K).le.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
               endif
               if (wi(I,J,K).gt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               end if
              ! Boundary conditions on sides   
                ti_new(NX,J,K) = 0.0
                ti(NX,J,K) = 0.0
                ti_new(I,NY,K) = 0.0
                ti(I,NY,K) = 0.0
                ti_new(1,J,K) = 0.0
                ti(1,J,K) = 0.0
                ti_new(I,1,K) = 0.0
                ti(I,1,K) = 0.0
              ! Boundary conditions at bottom  
                ti(I,J,1) = 0.0
                ti_new(I,J,1) = 0.0

              ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ! cccccccccccccccc ADJUST FALSE CONCENTRATION cccccccccccccccccccc  
              ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                     
                if (ti(I,J,K).lt.1e-3)then
                    ti(I,J,K) = 0.
                end if

                if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(ti(I,J,K).gt.0))then
                    ti(I,J,K) = 0.
                end if

                if (H(I,J,3).eq.0)then
                    ti(I,J,K) = 0.
                end if
               
                if (isnan(ti(I,J,K)))then
                  ti(I,J,K) = 0
                end if

              end if

             end do
            end do
           end do

          it_smolar_eng = it_smolar_eng + 1

    	 end do
	 
    end if

    it_deng = it_deng + 1 

end do

! Update boundary condition at the top  
  
      do J=1,NY
         do I=1,NX
           do K=1,NZ
               if (wi(I,J,K).le.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
               endif
               if (wi(I,J,K).gt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               end if
            end do
        end do
     end do

     !-----------------------------------------------------------------------    
     ! STEP 4: ENSURE MASS CONSERVATION  
     !-----------------------------------------------------------------------   

      ! Update mass in glacier after diffusion correction

      do J=2,NY-1
         do I=2,NX-1
           do K=2,NZ-1
              actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
            end do
         end do
      end do

  ! Ensure mass conservation

      actual_mass_hist_eng = sum(actual_mass_debris_eng/1e6)
      expected_mass_hist_eng = sum(expected_debris_mass_eng(:,:,NZ-1)/1e6)
      mass_ratio_hist_eng = (expected_mass_hist_eng/actual_mass_hist_eng)

      do J=2,NY-1
        do I=2,NX-1
          do K=2,NZ-1
            ti(I,J,K) = ti(I,J,K)*mass_ratio_hist_eng(1)
          end do
        end do
      end do

  ! Check if it worked

      do J=2,NY-1
        do I=2,NX-1
          do K=2,NZ-1
           actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
          end do
        end do
      end do

      actual_mass_hist_eng = sum(actual_mass_debris_eng/1e6)
      mass_ratio_hist_eng = (expected_mass_hist_eng/actual_mass_hist_eng)

  ! Update boundary condition at the top                

      do J=1,NY
         do I=1,NX
           do K=1,NZ
               if (wi(I,J,K).le.0)then
                ti(I,J,NZ) = t_s(I,J)
                ti_new(I,J,NZ) = t_s(I,J)
               endif
               if (wi(I,J,K).gt.0)then
                 ti(I,J,NZ) = ti(I,J,NZ-1)
                 ti_new(I,J,NZ) = ti_new(I,J,NZ-1)
               end if
            end do
        end do
     end do

  ! Update sources and sinks

      do J=2,NY-1
         do I=2,NX-1
           do K=2,NZ-1
               actual_mass_debris_eng(I,J,K) = (H(I,J,3)/nz)*(deltax_d*deltax_d*ti(I,J,K))
              if (wi(I,J,NZ-1).gt.0)then
                 sink_debris(I,J,1)=sum_c_meltout(I,J)
              else
                 sink_debris(I,J,1)=0.
              endif
           end do
        end do
     end do

! cccccccccccccccccccccccccccccccccccccccccc      
! ccccccccccc MASS CONSERVED? cccccccccccccc   
! cccccccccccccccccccccccccccccccccccccccccc        

   ! Mass output  
 
      do J=1,NY
         do I=1,NX
             mass_out(I,J) = sum_c_meltout(I,J)*dx*dx*(H(I,J,3)/NZ)
         end do
      end do

    write(*,*),'mass in engl = ',sum(mass_in)
    write(*,*),'mass out engl = ',sum(mass_out)
    write(*,*),'mass in/out ratio engl = ',sum(mass_out)/sum(mass_in)

! cccccccccccccccccccccccccccccccccccccccccc
! cccccccccc END OF THE LOOP ccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccc

yr_deng = yr_deng + 1

end do

! cccccccccccccccccccccccccccccccccccccccccc   
! cccccccc PREPARE FOR OUTPUT cccccccccccccc    
! cccccccccccccccccccccccccccccccccccccccccc    

sum_mass_meltout = mass_out

end subroutine debris_cover_main_englacial
