!  Harmsynthloadeffect.f90 
!
!  FUNCTIONS:
!  Harmsynthloadeffect - Entry point of console application.
!
!****************************************************************************

      program Harmsynthloadeffect
      implicit none
	character*80000::line
	character*4 fh
	integer i,j,k,n,m,nn,maxn,sn,kln,astat(5),kk
      real*8 tdn(24),rec(8000),flv(40000,3),BLH(3),rln(3),pi,RAD
	real*8 GRS(6),dcin(40),dout(40),rw,re,fk,t,NFD(5),gr,rlat0
	integer::status=0
	real*8,allocatable::cnm(:),snm(:)
	real*8,allocatable::pnm(:),dpt1(:),dpt2(:)
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
	pi=datan(1.d0)*4.d0;RAD=pi/180.d0; rw=1.025d3; re=5.517d3
	maxn=360!Input maximum trucated degree 输入负荷球谐系数最大计算阶数
      GRS(6)=1.d0!=1.d0 land water or sea level variation load, =-1.d0 surface atmosphere load
      !打开地面等效水高球谐系数文件
      !read the load EWH spherical harmonic coefficient model
 	allocate(cnm((maxn+2)**2), stat=astat(1))
 	allocate(snm((maxn+2)**2), stat=astat(2))
	if (sum(astat(1:2)) /= 0) goto 902
      cnm=0.d0;snm=0.d0
      open(unit=8,file="landwater2018011012cs.dat",status="old",iostat=status)
      if(status/=0)goto 904
      read(8,'(a)') line
	do while(.not.eof(8))
         read(8,'(a)') line
         call PickReclong(line,kln,rec,sn)
         if(sn<4)goto 1011
         if(rec(2)>dble(maxn)+0.1)goto 1011
	   n=nint(rec(1));m=nint(rec(2))
         if(n>maxn)goto 1011
         cnm(n*(n+1)/2+m)=rec(3);snm(n*(n+1)/2+m)=rec(4)
1011     continue
      enddo
901   close(8)
      if(n<maxn)maxn=n
      !将等效水高球谐系数转换为位系数直接影响
      !convert the EWH spherical harmonic coefficient to the direct influence of geopotential (Stokes) coefficients
      fk=3.d0*rw/re
      do n=1,maxn
        do m=0,n
            j=n*(n+1)/2+m
            cnm(j)=cnm(j)*fk/dble(2*n+1);snm(j)=snm(j)*fk/dble(2*n+1)
         enddo
	enddo
	!Read load love numbers 读负荷勒夫数
      flv=0.d0
      open(unit=8,file="Love_load_cm.dat",status="old",iostat=status)
      if(status/=0) goto 902 
      do i=1,6
        read(8,'(a)') line
      enddo
      n=0
	do while(.not.eof(8).and.n<3600)
        n=n+1
	  read(8,*,end=903)i,(flv(n,j),j=1,3)
      enddo
903   close(8)
 	allocate(pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2))
      open(unit=8,file="calcpnt.txt",status="old",iostat=status)
      if(status/=0)goto 904
      open(unit=10,file="reslt.txt",status="replace")!Output file 输出文件
      read(8,'(a)') line
      write(10,'(a)')trim(line)
      rlat0=999.d0;kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickRecord(line,kln,rec,sn)
         if(sn<4)goto 906
         BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4);kk=kk+1
         call BLH_RLAT(GRS,BLH,rln);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
!The Legendre function is calculated in advance, which can be speeded up.一次性计算勒让德函数，提速
         if(dabs(rln(2)-rlat0)>1.d-6)then
           t=dsin(rln(2)*RAD); call BelPnmdt(pnm,dpt1,dpt2,maxn,t);rlat0=rln(2)
         endif
         !calculate all-element load effects 计算全要素负荷效应
         call Loadeffectpnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr,BLH(3))
         tdn(11)=tdn(10)-tdn(1)  !orthometric height 正（常）高变化（mm）
         tdn(12:14)=tdn(12:14)*1.d2 !10μE
         write(10,'(a,40F12.4)')trim(line),(tdn(i),i=1,14)
         if(kk/500*500==kk)write(*, '(a,i9)'), '    Calculated point number: ',kk
906      continue
	enddo
905	close(8)
      close(10)
      deallocate(pnm,dpt1,dpt2)
904   deallocate(cnm,snm)
902	continue
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
