      program writeWW

      double complex ww,wwsq
      double precision wreal,wimag,Q
      double precision pi,tot
 
      open(10,file='ww.dat.ascii')
      open(20,file='omega.dat')
      rewind(10)
      rewind(20)

ca      write(6,*) 'Give dimension'
ca      read(5,*) nsum
ca      write(6,*) 'Give rad (1) or mHz (2)'
ca      read(5,*) ifreq
      ifreq=2
      pi=acos(-1.0)
      fac=1.0/(pi*2.0)

      ione=1
      itwo=2

      wrealsum=0.0
      wimagsum=0.0
      wrsum=0.0
      wisum=0.0
      tot=0.0
      Qsum=0.0

      do i=1,5000
         read(10,*,END=100) wreal,wimag
         wrsum=wrsum+wreal
         wisum=wisum+wimag
         ww=dcmplx(wreal,wimag)
         wwsq=zsqrt(ww)
         if(ifreq.eq.ione) then
           write(20,*) dreal(wwsq),dimag(wwsq)
         else if(ifreq.eq.itwo) then
           wimag=dimag(wwsq)
           wreal=dble(wwsq)
           Q=wreal/(2.d0*wimag)
           wreal=1000.0*fac*dble(wwsq)
c          wimag=1000.0*0.5*wreal/wimag
           write(20,*) wreal,Q
           Qsum=Qsum+Q
           wrealsum=wrealsum+wreal
           wimagsum=wimagsum+wimag
           tot=tot+1.0
         endif
      enddo

100   continue

      close(10)
      close(20)

      wrealav=wrealsum/tot
      wimagav=wimagsum/tot
      Qav=wrealav/(2.d0*wimagav*fac*1000.0)
      Qsumav=Qsum/tot

      write(6,*) 'Average values are:'
      write(6,*) 'av freq = ',wrealav
ca      write(6,*) 'av Q    = ',wimagav,Qav,Qsumav
      write(6,*) 'av Q    = ',Qav
      write(6,*) 'diagsum is',wrsum,wisum

      end
