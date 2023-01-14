      subroutine  debris_cover_main (H, HT, mask, G, G_obs, h_debris, inputlocation_debris, inoutdebris, meltout_debris, h_debrisini, term1_debris_x, term1_debris_y, term2_debris, term3_debris, fdebris, margin, VX, VY, area_debris, actual_mass_debris,expected_mass_debris,mass_ratio_hist,c_debris_out,h_debris_margin,sum_mass_meltout)

      use        PARAMETER_LIST

      implicit none

!--   Parameter List in
      real, dimension(NX,NY)               :: HT
      double precision, dimension(NX,NY,3) :: H
      real,dimension(NX,NY)                :: mask
      real,dimension(NX,NY)                :: margin
      real,dimension(NX,NY)                :: inoutdebris,c_debris_out,sum_mass_meltout
      real,dimension(NX,NY)                :: h_debris, h_debris_margin, h_debris_diff	
      real,dimension(NX,NY)                :: meltout_debris	
      real,dimension(NX,NY)                :: h_debrisini
      real,dimension(NX,NY)                :: term1_debris_x, fl_debris_x, term1a_debris_x, term1b_debris_x	
      real,dimension(NX,NY)                :: term1_debris_y, fl_debris_y, term1a_debris_y, term1b_debris_y
      real,dimension(NX,NY)                :: term2_debris, advel_x, advel_y, vel_diff
      real,dimension(NX,NY)                :: term3_debris
      real,dimension(NX,NY)                :: fdebris, area_debris, actual_mass_debris, new_debris_mass, expected_mass_debris, old_debris_mass
      real,dimension(1)                    :: actual_mass_hist, expected_mass_hist, mass_ratio_hist
      real,dimension(NX,NY)                :: inputlocation_debris
      double precision,dimension(NX,NY,NZ) :: VX,VY

!--   Parameter List out
      real,dimension(NX,NY) :: G,G_obs

!--   Variables
      integer               :: i,j,it_d,yr_d,it_smolar

!-----------------------------------------------------------------------                               
! Run the debris cover model
!-----------------------------------------------------------------------    

write(*,*),'starting the loop for the debris cover...'

yr_d = 1

do while (yr_d.le.numberofyears)

  it_d = 1

  do while (it_d.le.(nyears_d/deltat_d))
                                                                                    
!-----------------------------------------------------------------------    
! DEBRIS INPUT
!-----------------------------------------------------------------------                                                                                                                                                                   
  ! Initialize

      do J=1,NY 
         do I=1,NX                                                              
            inoutdebris(I,J) = 0.
         end do
       end do

!-----------------------------------------------------------------------    
! DEBRIS MELTOUT INPUT
!-----------------------------------------------------------------------  

      do J=1,NY
         do I=1,NX
            if (G(I,J).lt.0.and.sum_mass_meltout(I,J).gt.0)then
               meltout_debris(I,J) = -((sum_mass_meltout(I,J)) / (rho_debris)) / (deltax_d*deltax_d*deltat_d)
            elseif (G(I,J).ge.0)then
               meltout_debris(I,J) = 0.
            endif
         end do
       end do

    if(sum(meltout_debris).eq.0)then
      do J=1,NY
         do I=1,NX
            if (G(I,J).lt.0.and.sum_mass_meltout(I,J).gt.0)then
!              meltout_debris(I,J) = (c_debris_out(I,J)*G(I,J)) / ((1-phi_debris)*(rho_debris))
               meltout_debris(I,J) = -((c_debris_out(I,J)*deltax_d*deltax_d*(H(I,J,3)/nz)) / (rho_debris)) / (deltax_d*deltax_d*deltat_d)
            elseif (G(I,J).ge.0)then
               meltout_debris(I,J) = 0.
            endif
         end do
       end do
     endif

       do J=1,NY
          do I=1,NX
             if (H(I,J,3).eq.0)then
                meltout_debris(I,J)=0
             endif
          end do
       end do

   if (it_d.eq.(nyears_d/deltat_d))then
    write(*,*),'mass in supra = ',sum(sum_mass_meltout)
   endif
!-----------------------------------------------------------------------
! DEBRIS FLUX CALCULATION
!-----------------------------------------------------------------------

      do J=1,NY
         do I=1,NX                                                                                                                                                                                    
             fl_debris_x(I,J) = vx(I,J,1)*h_debris(I,J)
             fl_debris_y(I,J) = vy(I,J,1)*h_debris(I,J)
         end do
      end do

      !-----------------------------------------------------------------------    
      ! STEP 1: DEBRIS THICKNESS CALCULATION
      !-----------------------------------------------------------------------                                                                
 
      ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! cccc CALCULATE INITIAL DEBRIS THICKNESS WITHOUT NUMERICAL DIFFUSION CORRECTION cccc
      ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do J=1,NY 
         do I=1,NX                                                                                                 
    
	   h_debrisini(I,J) = h_debris(I,J)

	   if((I.gt.1).and.(J.gt.1).and.(I.lt.NX).and.(J.lt.NY))then

            ! cccccccccccccccccccccccccccccccccccccccccccccc
            ! cccccccccc START ADVECTION SCHEME cccccccccccc
            ! cccccccccccccccccccccccccccccccccccccccccccccc

            ! cccccccccc FIRST ORDER UPWIND SCHEME ccccccccc

             if (vx(I,J,1).lt.0)then                            
                term1_debris_x(I,J) = -vx(I,J,1)*((h_debris(I+1,J)-h_debris(I,J))/(2*deltax_d)) - h_debris(I,J)*((vx(I+1,J,1)-vx(I,J,1))/(2*deltax_d))                      !      
             else if (vx(I,J,1).ge.0)then                                                                         
                term1_debris_x(I,J) = -vx(I,J,1)*((h_debris(I,J)-h_debris(I-1,J))/(2*deltax_d)) - h_debris(I,J)*((vx(I,J,1)-vx(I-1,J,1))/(2*deltax_d))            
             endif                                                                                                                                                      
             if (vy(I,J,1).lt.0)then  
                term1_debris_y(I,J) = -vy(I,J,1)*((h_debris(I,J+1)-h_debris(I,J))/(2*deltax_d)) - h_debris(I,J)*((vy(I,J+1,1)-vy(I,J,1))/(2*deltax_d)) 
             else if (vy(I,J,1).ge.0)then                 
                term1_debris_y(I,J) = -vy(I,J,1)*((h_debris(I,J)-h_debris(I,J-1))/(2*deltax_d)) - h_debris(I,J)*((vy(I,J,1)-vy(I,J-1,1))/(2*deltax_d))
             endif                                                                       
                                                                         
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! cccccccccccccccccccc END OF ADVECTION SCHEME ccccccccccccccccccc
           ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            ! Meltout of debris material from the ice
  
             term2_debris(I,J) = abs(meltout_debris(I,J));
       
            ! In and output
     
             term3_debris(I,J) = inoutdebris(I,J);
  
            ! Calculate new debris thickness
  
            h_debris(I,J) = h_debrisini(I,J) + deltat_d*(term1_debris_x(I,J) + term1_debris_y(I,J) + term2_debris(I,J) + term3_debris(I,J));

            ! Adjust false debris thickness

            if (h_debris(I,J).lt.0)then
                 h_debris(I,J) = 0.
            end if

            if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(h_debris(I,J).gt.0))then
                 h_debris(I,J) = 0.
            end if

         ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! cccc REMOVAL OF DEBRIS AT THE GLACIER MARGINS (UNLOADING PARAMETERIZATION) cccc
         ! ccCCCCCcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            h_debris_margin(I,J) = 0

            if (margin(I,J).eq.1)then
                h_debris_margin(I,J) = h_debris(I,J)
            else if (margin(I,J).eq.0)then
                h_debris_margin(I,J) = 0
            end if

         ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! cccccccccccccccccccccccc SAVE TERMS FOR MASS CONSERVATION ccccccccccccccccccccc
         ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           
           actual_mass_debris(I,J) = rho_debris*(deltax_d*deltax_d*h_debris(I,J))
           new_debris_mass(I,J) = rho_debris*(deltax_d*deltax_d*(inoutdebris(I,J))*deltat_d) - rho_debris*(deltax_d*deltax_d*(meltout_debris(I,J))*deltat_d) 
           old_debris_mass(I,J) = rho_debris*(deltax_d*deltax_d*(h_debris_margin(I,J)))
           expected_mass_debris(I,J) = expected_mass_debris(I,J) + new_debris_mass(I,J) - old_debris_mass(I,J)

          ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! cccccccccccccccccccccccc REMOVE DEBRIS AT MARGINS ccccccccccccccccccccccccccccc  
          ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

             if (margin(I,J).eq.1)then
                 h_debris(I,J) = 0
             endif 

           end if
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
	     if((I.gt.1).and.(J.gt.1).and.(I.lt.NX).and.(J.lt.NY))then
                if((h_debris(I,J).gt.0).and.(inoutdebris(I,J).eq.0).and.(margin(I,J).eq.0))then
                     ! X-DIRECTION 
                     advel_x(I,J) = -(0.5*((abs(vx(I,J,1))*deltax_d - (deltat_d*((vx(I,J,1))**2))) / h_debris(I,J))) * (((h_debris(I,J)+h_debris(I+1,J))/2)-((h_debris(I,J)+h_debris(I-1,J))/2))/deltax_d
                     ! Y-DIRECTION 
                     advel_y(I,J) = -(0.5*((abs(vy(I,J,1))*deltax_d - (deltat_d*((vy(I,J,1))**2))) / h_debris(I,J))) * (((h_debris(I,J)+h_debris(I,J+1))/2)-((h_debris(I,J)+h_debris(I,J-1))/2))/deltax_d
                     ! VELOCITY
                     vel_diff(I,J) = sqrt(advel_x(I,J)**2 + advel_y(I,J)**2)
                end if
             end if
          end do
        end do

	!-----------------------------------------------------------------------    
	! STEP 3: RECALCULATE ADVECTION EQUATION WITH ANTI-DIFFUSION VELOCITY
	!-----------------------------------------------------------------------                                                                

        ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! cccccccccc CALCULATE DEBRIS THICKNESS WITH NUMERICAL DIFFUSION CORRECTION ccccccccc
        ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if(it_d.eq.(nyears_d/deltat_d))then

        it_smolar = 1	

         do while (it_smolar.le.(0))

          do J=1,NY 
             do I=1,NX                                                                                                 
    
	       h_debrisini(I,J) = h_debris(I,J)

	       if((I.gt.1).and.(J.gt.1).and.(I.lt.NX).and.(J.lt.NY))then

                ! cccccccccccccccccccccccccccccccccccccccccccccc
                ! cccccccccc START ADVECTION SCHEME cccccccccccc
                ! cccccccccccccccccccccccccccccccccccccccccccccc

                ! cccccccccc FIRST ORDER UPWIND SCHEME ccccccccc

                 if (vx(I,J,1).lt.0)then                            
                    term1_debris_x(I,J) = -advel_x(I,J)*((h_debris(I+1,J)-h_debris(I,J))/(2*deltax_d)) - h_debris(I,J)*((advel_x(I+1,J)-advel_x(I,J))/(2*deltax_d))                                        !      
                 else if (vx(I,J,1).ge.0)then                                                                         
                    term1_debris_x(I,J) = -advel_x(I,J)*((h_debris(I,J)-h_debris(I-1,J))/(2*deltax_d)) - h_debris(I,J)*((advel_x(I,J)-advel_x(I-1,J))/(2*deltax_d))                                        !      
                 endif                                                                                                                                                      
                 if (vy(I,J,1).lt.0)then  
                    term1_debris_y(I,J) = -advel_y(I,J)*((h_debris(I,J+1)-h_debris(I,J))/(2*deltax_d)) - h_debris(I,J)*((advel_y(I,J+1)-advel_y(I,J))/(2*deltax_d)) 
                 else if (vy(I,J,1).ge.0)then                 
                    term1_debris_y(I,J) = -advel_y(I,J)*((h_debris(I,J)-h_debris(I,J-1))/(2*deltax_d)) - h_debris(I,J)*((advel_y(I,J)-advel_y(I,J-1))/(2*deltax_d)) 
                 endif                                                                       
                                                                         
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ! cccccccccccccccccccc END OF ADVECTION SCHEME ccccccccccccccccccc
               ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
                ! Calculate new debris thickness

                 h_debris_diff(I,J) = deltat_smolar*(term1_debris_x(I,J) + term1_debris_y(I,J))
                 if (h_debris_diff(I,J).lt.0)then
                    h_debris_diff(I,J)=0
                 end if
                 h_debris(I,J) = h_debrisini(I,J) - (h_debris_diff(I,J))

                ! Adjust false debris thickness

                if (h_debris(I,J).lt.0)then
                     h_debris(I,J) = 0.
                end if

                if ((sqrt(vx(I,J,1) * vx(I,J,1) + vy(I,J,1) * vy(I,J,1)).eq.0).and.(h_debris(I,J).gt.0))then
                     h_debris(I,J) = 0.
                end if

              end if

             end do
          end do

          it_smolar = it_smolar + 1

	      end do
	 
       end if

      it_d = it_d + 1 

  end do

     !-----------------------------------------------------------------------    
     ! STEP 4: ENSURE MASS CONSERVATION  
     !-----------------------------------------------------------------------   

      write(*,*),'mass out supra = ',rho_debris*(deltax_d*deltax_d*(sum(h_debris_margin)))
      write(*,*),'mass in/out ratio supra = ',(sum(sum_mass_meltout))/(rho_debris*(deltax_d*deltax_d*(sum(h_debris_margin))))

      ! Remove debris at the margins

      do J=1,NY
         do I=1,NX
           if (margin(I,J).eq.1)then
               h_debris(I,J) = 0
           endif
         end do
      end do

      ! Update mass on glacier after diffusion correction

      do J=1,NY
         do I=1,NX
            actual_mass_debris(I,J) = rho_debris*(deltax_d*deltax_d*h_debris(I,J)) 
         end do 
      end do

      ! Calculate sum of masses

      actual_mass_hist = (sum(actual_mass_debris(:,:))/1e9)
      expected_mass_hist = (sum(expected_mass_debris(:,:))/1e9)
      mass_ratio_hist = (expected_mass_hist/actual_mass_hist)

      ! Ensure mass conservation

      do J=1,NY
         do I=1,NX
            h_debris(I,J) = h_debris(I,J)*mass_ratio_hist(1)
         end do
      end do

      ! Check if it worked

      do J=1,NY
         do I=1,NX
            actual_mass_debris(I,J) = rho_debris*(deltax_d*deltax_d*h_debris(I,J))
         end do
      end do

      ! Calculate sum of masses                                                                                                                                                                                             
      actual_mass_hist = (sum(actual_mass_debris(:,:))/1e9)
      expected_mass_hist = (sum(expected_mass_debris(:,:))/1e9)
      mass_ratio_hist = (expected_mass_hist/actual_mass_hist)

  !-----------------------------------------------------------------------        
  ! STEP 5: DETERMINE DEBRIS-RELATED MELT-MODIFICATION FACTOR        
  !-----------------------------------------------------------------------              

  ! Calculate the debris-related melt-modification factor fdebris                                                                                 

   do J=1,NY
      do I=1,NX
         if(h_debris(I,J).le.0.015)then
            fdebris(I,J) = (26.667*h_debris(I,J)) + 1
         else if(h_debris(I,J).gt.0.015.and.h_debris(I,J).le.0.04)then
            fdebris(I,J) = (-16*h_debris(I,J)) + 1.64
         else if(h_debris(I,J).gt.0.04)then
            fdebris(I,J) = 0.1061*(h_debris(I,J)**(-0.7205))-0.07925
         end if
      end do
   end do

   do J=1,NY
      do I=1,NX
         if(fdebris(I,J).lt.0)then
            fdebris(I,J)=0
         end if
      end do
   end do

  !-----------------------------------------------------------------------  
  ! STEP 6: DETERMINE FRACTIONAL DEBRIS-COVERED AREA
  !-----------------------------------------------------------------------  

  ! Calculate fractional debris-covered area

   do J=1,NY
      do I=1,NX
         area_debris(I,J)=(1-exp(-20*h_debris(I,J)))
      end do
  end do

!area_debris = meltout_debris !!!!!!!!!

yr_d = yr_d + 1

end do

      end subroutine debris_cover_main
