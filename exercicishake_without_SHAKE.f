*****************************************************************************
*     Programa de Dinamica Molecular per colectivitat NVT (algoritme
*     de Berendsen) per un sistema de N molecules triatomiques dins una 
*     capsa cubica. Els atoms son de Ar (per exemple) i interactuen 
*     de la seguent forma: 
*
*          1) Els atoms d'una mateixa molecula no interactuen 
*             ja que fan SHAKE i formen un triangle equilater. 
*          2) Si pertanyen a molecules diferents interactuen amb un 
*             potencial de Lennard-Jones truncat a 2.5 sigmes
*
*           Exercici SHAKE, en que els alumnes han de fer una
*           subrutina que implementi la funcio SHAKE a fi de que 
*           els atoms formin molecules triatomiques amb forma 
*           de triangle equilater 
*****************************************************************************

      implicit double precision(a-h,o-z)
      double precision massa,lambda
      real*8,dimension(:,:), allocatable::gr_mat_atm, r_sum
      include 'exercicishake.dim'

c     1. Dimensionament de magnituds.

      dimension r(3,nmax,nmaxmol),rpro(3,nmax,nmaxmol)
      dimension vinf(3,nmax,nmaxmol)
      dimension accel(3,nmax,nmaxmol)


c     2. Lectura de dades i calcul de quantitats relacionades

      open(1,file='exercicishake.dades',status='old')
         read(1,*) nconf
         read(1,*) deltat,taut 
         read(1,*) nmolecules,tempref
         read(1,*) sigma,epsil
         read(1,*) massa
         read(1,*) r0
         read(1,*) num_bins
      close(1)

      natoms = 3
      nf = 6*nmolecules-3 !nombre de graus de 
c           llibertat nmolecules*(9-3)-3 
      rc = 2.5d0 !abast del potencial en unitats reduides

c     3. Lectura de la configuracio anterior en A i A/ps

      open(2,file='conf.data',status='old')
         do ic = 1,nmolecules
            do is = 1,natoms
               read(2,*) (r(l,is,ic),l=1,3)
               read(2,*) (vinf(l,is,ic),l=1,3)
            end do         
         end do         
         read(2,*) costat
      close(2)
c     4. Expressa les quantitats en unitats reduides

      call reduides(nmolecules,natoms,r,vinf,costat,deltat,
     &taut,tempref,epsil,sigma,massa,r0,uvel,utemps)
      tol=0.0000001d0
c     5. Comen�a el bucle de la generacio de configuracions 
      open(99, file='thermdata_withoutSHAKE.out', status='replace')
      write(99,*) ("Time, ekin, epot, etot, temp, lambda")

      allocate(gr_mat_atm(2,num_bins), r_sum(natoms, nmolecules))
      call g_r_atm(nmolecules,natoms,gr_mat_atm,r,num_bins,costat,1) 
      call dist_atomic(natoms, nmolecules,r, r_sum, 1)
      
      do i = 1,nconf
         call forces(nmolecules,natoms,r,costat,accel,rc,epot)
         call factlambda(nmolecules,natoms,vinf,nf,tempref,
     &deltat,taut,lambda)
         call velpospro(nmolecules,natoms,vinf,accel,deltat,lambda,
     &r,rpro)

         !call shake(nmolecules, r, rpro, natoms, deltat, r0, tol)
         
         call velocitat(nmolecules,natoms,r,rpro,deltat,vinf,
     &temperatura,nf,ecin)
         
         sim_temps = i*deltat

         write(99,*) sim_temps,ecin,epot,ecin+epot,temperatura,lambda
         
         call g_r_atm(nmolecules,natoms,gr_mat_atm,r,num_bins,costat,2)
         call dist_atomic(natoms, nmolecules,r,r_sum, 2)
     	 
      end do
      close(99)

c     Escriptura g(r)
      open(22, file='radial_func_withoutSHAKE.out', status='replace')
      write(22,*) ("dr, g(r)")
      call g_r_atm(nmolecules,natoms,gr_mat_atm,r,num_bins,costat,3)
      do  j=1,num_bins 
         write(22,*)  gr_mat_atm(1,j),gr_mat_atm(2,j)
      enddo
      close(22)
      call dist_atomic(natoms, nmolecules,r,r_sum, 3)
      
      
c     5. Escriptura de la darrera configuracio en A i A/ps 

      open(11,file='confnova_withoutSHAKE.data',status='unknown')
      do ic = 1,nmolecules
         do is = 1,natoms
            write(11,*) (r(l,is,ic)*sigma,l=1,3)
            write(11,*) (vinf(l,is,ic)*uvel,l=1,3)
         end do         
      end do         
      write(11,*) costat*sigma
      close(11)


      open(33,file='confnova_vmd_withoutSHAKE.data',status='unknown')
      write(33,*) natoms*nmolecules
      write(33,*) " "
      do ic = 1,nmolecules
         do is = 1,natoms
            write(33,*) 'A ', (r(l,is,ic)*sigma,l=1,3)
         end do         
      end do         
      close(33)
      stop
      end

*********************************************************
*********************************************************
c              subrutina reduides
*********************************************************
*********************************************************

      subroutine reduides(nmolecules,natoms,r,vinf,costat,deltat,
     &taut,tempref,epsil,sigma,massa,r0,uvel,utemps)
      implicit double precision(a-h,o-z)
      double precision massa
      include 'exercicishake.dim'
      dimension r(3,nmax,nmaxmol),vinf(3,nmax,nmaxmol)

c              unitat de temps expressada en ps 

      rgas = 8.314472673d0 !J/(mol*K) 
      utemps = sigma*dsqrt(massa/epsil)*dsqrt(10.d0/rgas)
      uvel = sigma/utemps !unitat de velocitat expressada en A/ps

      costat = costat/sigma
      r0 = r0/sigma
      deltat = deltat/utemps
      taut = taut/utemps
      tempref = tempref/epsil
      do ic = 1,nmolecules
         do is = 1,natoms
            do l = 1,3
               r(l,is,ic) = r(l,is,ic)/sigma
               vinf(l,is,ic) = vinf(l,is,ic)/uvel 
            end do
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subrutina forces
*********************************************************
*********************************************************

      subroutine forces(nmolecules,natoms,r,costat,accel,rc,epot)
      implicit double precision(a-h,o-z)
      include 'exercicishake.dim'
      dimension r(3,nmax,nmaxmol),accel(3,nmax,nmaxmol)

      do ic = 1,nmolecules
         do is = 1,natoms
            do l = 1,3
               accel(l,is,ic) = 0.d0 !fa 0 les acceleracions
            end do 
         end do 
      end do 
      epot = 0.d0 

c        interaccio entre atoms de molecules diferents

      do ic = 1,nmolecules-1
         do is = 1,natoms
            do jc = ic+1,nmolecules
               do js = 1,natoms
                  call lj(is,ic,js,jc,r,costat,accel,rc,pot)
                  epot = epot + pot
               end do
            end do
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subrutina Lennard-Jones
*********************************************************
*********************************************************

      subroutine lj(is,ic,js,jc,r,costat,accel,rc,pot)
      implicit double precision(a-h,o-z)
      include 'exercicishake.dim'
      dimension r(3,nmax,nmaxmol),accel(3,nmax,nmaxmol)
      dimension rij(3)

      rr2 = 0.d0
      do l = 1,3
         rijl = r(l,js,jc) - r(l,is,ic)
         rij(l) = rijl - costat*dnint(rijl/costat)
         rr2 = rr2 + rij(l)*rij(l)
      end do

      rr = dsqrt(rr2)
      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forsadist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         pot = 4.d0*(ynvrr12-ynvrr6)  
         do l = 1,3
            accel(l,is,ic) = accel(l,is,ic) - forsadist*rij(l)
            accel(l,js,jc) = accel(l,js,jc) + forsadist*rij(l)
         end do
      end if

      return
      end

*********************************************************
*********************************************************
c              subrutina factlambda
*********************************************************
*********************************************************

c     calcul de l'energia cinetica de la configuracio a fi
c     de determinar la temperatura, i per tant retocar el
c     sistema per adaptar-ho a la temperatura dessitjada

      subroutine factlambda(nmolecules,natoms,vinf,nf,tempref,
     &deltat,taut,lambda)
      implicit double precision(a-h,o-z)
      double precision lambda
      include 'exercicishake.dim'
      dimension vinf(3,nmax,nmaxmol)

      ecin = 0.d0
      do ic = 1,nmolecules
         do is = 1,natoms
            v2 = 0.d0
            do l = 1,3
               v2 = v2 + vinf(l,is,ic)*vinf(l,is,ic)
            end do
            ecin = ecin + 0.5d0*v2
         end do
      end do

      temperatura = 2.d0*ecin/dfloat(nf) 
      lambda = dsqrt(1.d0+deltat/taut*(tempref/temperatura-1.d0))

      return
      end

*********************************************************
*********************************************************
c              subrutina velpospro
*********************************************************
*********************************************************

c       Calcul de la velocitat als instants t+delta/2 i t i
c       la posicio provisional l'instant t + deltat.

      subroutine velpospro(nmolecules,natoms,vinf,accel,deltat,
     &lambda,r,rpro)
      implicit double precision(a-h,o-z)
      double precision lambda
      include 'exercicishake.dim'
      dimension vinf(3,nmax,nmaxmol),accel(3,nmax,nmaxmol),
     &r(3,nmax,nmaxmol),rpro(3,nmax,nmaxmol)

      do ic = 1,nmolecules
         do is = 1,natoms
            do l = 1,3
               vsup = vinf(l,is,ic) + accel(l,is,ic)*deltat
               vsup = vsup*lambda
               rpro(l,is,ic) = r(l,is,ic) + vsup*deltat
            end do
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subrutina velocitat
*********************************************************
*********************************************************

c       Calcul de la velocitat als instants t+delta/2 i t i
c       la posicio l'instant t + deltat.

      subroutine velocitat(nmolecules,natoms,r,rpro,deltat,vinf,
     &temperatura,nf,ecin)
      implicit double precision(a-h,o-z)
      include 'exercicishake.dim'
      dimension v(3,nmax,nmaxmol),vinf(3,nmax,nmaxmol),
     &rpro(3,nmax,nmaxmol),r(3,nmax,nmaxmol)

      ecin = 0.d0
      do ic = 1,nmolecules
         do is = 1,natoms
            v2 = 0.d0
            do l = 1,3
               vsup = (rpro(l,is,ic) - r(l,is,ic))/deltat
               v(l,is,ic) = (vinf(l,is,ic)+vsup)/2.d0 !v(t)
               vinf(l,is,ic) = vsup
               r(l,is,ic) = rpro(l,is,ic)
               v2 = v2 + v(l,is,ic)*v(l,is,ic)
            end do
            ecin = ecin + 0.5d0*v2
         end do
      end do
      temperatura = 2.d0*ecin/dfloat(nf)

      return
      end
      
      
      
*********************************************************
*********************************************************
c              subrutina radial distribution
*********************************************************
*********************************************************
            
      subroutine g_r_atm(nmolecules, natoms,gr_mat_atm, r, num_bins, 
     &   box_size, switch_case)
            
            implicit double precision(a-h,o-z)
            integer, intent(in)      :: switch_case
            real*8, parameter        :: pi = 4.d0 * datan(1.d0)
            integer, save            :: index_mat,total_part, n_gdr
            real*8, save             :: dr, dist, dv, ndg, dens
            include 'exercicishake.dim'
            dimension r(3,nmax,nmaxmol), gr_mat_atm(2,num_bins), rij(3)
            
            
            select case (switch_case)
               case (1)
                  ! SWITCH 1 => Memmory initialization
                  n_gdr = 0
                  total_part=natoms*nmolecules
               
                  dr = box_size / (2.d0*dble(num_bins))
                  dens = dble(total_part) / (box_size ** 3)

                  gr_mat_atm(1,:) = [(dble(i)*dr, i=1, num_bins)]
                  gr_mat_atm(2,:) = 0.d0
               
               case (2)
                  ! SWITCH 2 => Positions binning
                  n_gdr = n_gdr + 1
                  
                  ! Cae of atoms with other molecules
                  do ic = 1,nmolecules-1
                     do is = 1,natoms
                        do jc = ic+1,nmolecules
                           do js = 1,natoms
                              rij(:) = r(:,js,jc)-r(:,is,ic)
                              rij = rij - box_size*dnint(rij/box_size)
               
                              dist=dsqrt(sum(rij**2))
                              ! Apliquem el cutoff de maxima distancia
                              if (dist .lt. box_size/2.d0) then
                                 index_mat = int(dist/dr) + 1
                                 gr_mat_atm(2,index_mat) = 
     &                           gr_mat_atm(2,index_mat) + 2.d0
                              end if
                           end do
                        end do
                     end do
                  end do

                  ! Case with atoms of the same molecule
                  do ic = 1, nmolecules
                     do is = 1, natoms-1
                        do js = is + 1, natoms
                           rij(:) = r(:,js,ic) - r(:,is,ic)
                           rij = rij - box_size*dnint(rij/box_size)
               
                           dist=dsqrt(sum(rij**2))
                           ! Apliquem el cutoff de maxima distancia
                           if (dist .lt. box_size/2.d0) then
                              index_mat = int(dist/dr) + 1
                              gr_mat_atm(2,index_mat) = 
     &                        gr_mat_atm(2,index_mat) + 2.d0
                           endif
                        enddo
                     enddo
                  enddo

               case (3)
                  ! SWITCH = 3 => Normalization of g(r)
                     
                  do i = 1, num_bins
                     associate(gdr => gr_mat_atm(2,i))
                     dv = (((dble(i) + 1.d0)**3)-(dble(i)**3))*(dr**3)
                     ndg = (4.d0/ 3.d0) * pi * dv * dens
                     
                     gdr = gdr / (dble(total_part) * ndg * dble(n_gdr))
                     end associate
                  end do
                     
               end select
                  
            end subroutine g_r_atm
    

      subroutine shake(nmolecules, r, rpro, natoms, deltat, r0, tol)
         implicit double precision(a-h,o-z)
         real*8 :: lambda_shake 
         include 'exercicishake.dim'
         dimension r(3,nmax,nmaxmol), rpro(3,nmax,nmaxmol)
         dimension rnova(3,3,nmaxmol)
         dimension rpij(3)
         dimension rij(3)
         
         dimension rpro_p(3, 3), r_p(3, 3)

         do im = 1, nmolecules

            rpro_p = rpro(:, :, im)
            r_p = r(:, :, im)

            r12_m = dsqrt(sum((rpro_p(:, 1) - rpro_p(:, 2))**2))
            r13_m = dsqrt(sum((rpro_p(:, 1) - rpro_p(:, 3))**2))
            r23_m = dsqrt(sum((rpro_p(:, 2) - rpro_p(:, 3))**2))

            do while ((abs(r12_m**2-r0**2).gt. tol) .and. 
     &      (abs(r13_m**2-r0**2).gt. tol) .and. (abs(r23_m**2-r0**2)
     &      .gt. tol))
            
               i_contador_pij = 1
               do iai = 1, natoms - 1
                  do iaj = iai + 1, natoms
                  rpij(:) = rpro_p(:, iaj) - rpro_p(:, iai)
                  rij(:) = r_p(:, iaj) - r_p(:, iai)
                  
                  rpij_norm_sq = (sum(rpij(:)**2))
                  
                  lambda_shake= (rpij_norm_sq - r0**2)/
     &            (8.0d0 * deltat**2 * sum(rpij(:)*rij(:)))
                  
                  rnova(:, iai, im) = rpro_p(:, iai)+
     &            (2.0d0 * lambda_shake* deltat**2 * rij(:))
                  rnova(:, iaj, im) = rpro_p(:, iaj)-
     &            (2.0d0 * lambda_shake* deltat**2 * rij(:))
                  rpro_p(:,iai) = rnova(:, iai, im)
                  rpro_p(:,iaj) = rnova(:, iaj, im)
                  enddo
               enddo

               
               r12_m = dsqrt(sum((rpro_p(:, 1) - rpro_p(:, 2))**2))
               r13_m = dsqrt(sum((rpro_p(:, 1) - rpro_p(:, 3))**2))
               r23_m = dsqrt(sum((rpro_p(:, 2) - rpro_p(:, 3))**2))


            enddo
            
         enddo

         rpro(:, :, :) = rnova(:, :, :)
            

      end subroutine shake

      subroutine dist_atomic(natoms, nmolecules,r,r_sum, switch_case)
         implicit double precision(a-h,o-z)
         integer, intent(in)      :: switch_case
         integer, save            :: ntimes
         real*8                   :: r_sum

         include 'exercicishake.dim'
         dimension r(3,nmax,nmaxmol),r_sum(natoms, nmolecules)
         dimension r_aux(natoms, nmolecules), std(nmolecules)
         dimension dist_mean(nmolecules)

         select case (switch_case)
            case (1)
               ntimes = 0
               r_sum = 0.d0
               std = 0.d0
               mean= 0.d0

            case (2)
               ntimes = ntimes+1
               do imol = 1, nmolecules
                  index=1
                  do i= 1,natoms-1
                     do j= i+1, natoms
                        rij = dsqrt(sum(((r(:,j, imol)-
     &                  r(:,i, imol))**2)))
                        
                        r_sum(index, imol)= r_sum(index, imol)+rij
                        index = index + 1
                     enddo
                  enddo
               enddo

            case (3)
               r_aux = r_sum(:, :)/dble(ntimes)
               dist_mean = sum(r_aux(:,:), dim=2)/3.d0

               do imol = 1,nmolecules
                  do j=1, natoms
                     std(imol)= std(imol)+
     &              (dist_mean(imol)-r_aux(j, imol))**2
                  enddo
               enddo
               std = dsqrt(1.d0/dble(natoms)*std)

               open(44, file= 'distance_atoms_withoutSHAKE.data'
     &         , status='replace',action='write')
               do imol=1, nmolecules
                  write(44, *) imol, dist_mean(imol), std(imol)
               enddo

               close(44)

            end select



      end subroutine dist_atomic