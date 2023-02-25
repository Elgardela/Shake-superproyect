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
      real*8,dimension(:,:), allocatable::gr_mat
      include 'exercicishake.dim'

c     1. Dimensionament de magnituds.

      dimension r(3,nmax,nmaxmol),rpro(3,nmax,nmaxmol),
     &rnova(3,nmax,nmaxmol)
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

c     5. Comença el bucle de la generacio de configuracions 
      open(99, file='thermdata.out', status='replace')
      write(99,*) ("Iter, ekin, epot, etot, temp, lambda")

      allocate(gr_mat(num_bins))
      do i = 1,nconf
         call forces(nmolecules,natoms,r,costat,accel,rc,epot)
         call factlambda(nmolecules,natoms,vinf,nf,tempref,
     &deltat,taut,lambda)
         call velpospro(nmolecules,natoms,vinf,accel,deltat,lambda,
     &r,rpro)
c
c         IDEA: afegiu una rutina "shake" aqui...
c         call shake(nmolecules,r,rpro,rnova,r0)
c
         call velocitat(nmolecules,natoms,r,rpro,deltat,vinf,
     &temperatura,nf,ecin)
         write(99,*) (i,ecin,epot,ecin+epot,temperatura,lambda,l=1,1)
         
         call g_r(gr_mat, r, num_bins, costat, switch_case)
      end do
      close(99)
c     5. Escriptura de la darrera configuracio en A i A/ps 

      open(11,file='confnova.data',status='unknown')
      do ic = 1,nmolecules
         do is = 1,natoms
            write(11,*) (r(l,is,ic)*sigma,l=1,3)
            write(11,*) (vinf(l,is,ic)*uvel,l=1,3)
         end do         
      end do         
      write(11,*) costat*sigma
      close(11)

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
      
      subroutine g_r(gr_mat, r, num_bins, box_size, switch_case)
      implicit double precision(a-h,o-z)
      integer(kind=i64), intent(in)      :: switch_case
      
      integer(kind=i64), save                      :: i, index_mat, j, total_part, n_gdr
      real(kind=dp), save                          :: dr, dist, dv, ndg, dens
      include 'exercicishake.dim'
      dimension r(3,nmax,nmaxmol), gr_mat(2,num_bins), rij(3)
      
      
       select case (switch_case)
            case (1_i64)
            	n_gdr = 0_i64
      		total_part=nmax*nmaxmol
      
      		dr = box_size / real(2*num_bins, kind=dp)
     		dens = real(total_part, kind=dp) / (box_size ** 3)

       		gr_mat(1,:) = [(real(i, kind=dp)*dr, i=1, num_bins)]
       		gr_mat(2,:) = 0.0_dp
       		
       	    case (2_i64)
       	    	n_gdr = n_gdr + 1_i64
       	    	
	        do ic = 1,nmolecules-1
		  do is = 1,nmax
		    do jc = ic+1,nmolecules
		       do js = 1,natoms
			  rij(:)=r(:,js,jc)-r(:,is,ic)
			  rij=rij - box_size*dnint(rij/box_size)
			  
			  dist=dsqrt(sum(rij**2))
			  
			  ! Apliquem el cutoff de maxima distancia
                          if (dist .lt. box_size/2) then
                              index_mat = int(dist/dr, kind=i64) + 1_i64
                              gr_mat(2, index_mat) = gr_mat(2, index_mat) + 2.0_dp
                          end if
		       end do
		    end do
		  end do
	        end do
        
              
             case (3_i64)
        
                ! SWITCH = 3 => Calculem g(r)
                open(11, file='radial_distribution.txt',status='unknown',access='append', action="write")
 
                do i = 1, num_bins
                    associate(gdr => gr_mat(2,i))
                        dv = (((real(i, kind=dp) + 1.0_dp) ** 3) - (real(i, kind=dp) ** 3)) * (dr ** 3)
                        ndg = (4.0_dp / 3.0_dp) * pi * dv * dens
                        gdr = gdr / (real(total_part, kind=dp) * ndg * real(n_gdr, kind=dp))
                    end associate
                    write(11,*)   gr_mat(1,i), gr_mat(2,i)
                end do
                write(1,*) " "
                close(1)
            end select
            
      end subroutine g_r
