#include "fintrf.h"
C======================================================================   
#if 0
C  fine1.F - Gateway function for fine1.f on vidar
C
c  Magnetics function to determine the fields from an 
c  axisymmetric rectangular uniform current region 
c  (ri,ro,zl,zu,cur) to a point (rp,zp).
c  routine operates on matrix of data: ri,ro,zl,zu,cur,rp,zp
c  Does not use %VAL() calling convention to call fine1.f which calls fine.f
c  (for some reason this doesn't work and isn't necessary on vidar)
c
c  SYNTAX: [hr,hz,fl] =fine1g(ri,ro,zl,zu,cur,rp,zp,(idebug))
c
c  INPUTS:
c           ri,ro,zl,zu = inner,outer,lower,upper dimensions of coil [matrix]
c           cur         = current [matrix]
c           rp,zp       = receiving point [matrix]
c           idebug      = debug parameter, 0=no, >0 = debug print out
c  OUTPUTS:
c           hr, hz = fields at points rp,zp [Amp/m]
c           fl     = flux at point [Vs/mu_o]
c                    mu_o= 4*pi*1e-7
c
C     
C     Compile with mex by calling link_fine1.m
C     
#endif
C     
C     fine1.f
c     routine to call fine.f many times for vectors of coils and field points
c     assumes the correct address is passed in arrays using %VAL()
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2011 The MathWorks, Inc.
C     $Revision: 1.13.2.5 $
C
c GA Matthew Lanctot
C VERSION @(#)fine1.F	1.2 05/26/14
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN
 
C     Pointers to input/output mxArrays:
      mwPointer ri_8pr, ro_8pr, zl_8pr, zu_8pr, cu_8pr
      mwPointer rp_8pr, zp_8pr
      mwPointer hr_8pr, hz_8pr, fl_8pr
 
C     Array information:
      mwPointer m, n
      mwSize size
 
C     Arguments for communication with computational routine:
C     Maximum size = numel
      integer*4 numel, num
      parameter(numel = 2000000)
      real*8 ri_8(numel),ro_8(numel)
      real*8 zl_8(numel), zu_8(numel), cu_8(numel)
      real*8 rp_8(numel), zp_8(numel)
      real*8 hr_8(numel), hz_8(numel), fl_8(numel)
      integer idebug_pt

c     debug parameter
      real*8 idebug_8
      integer idebug
      
C     Debugging variables
C     Uncomment lines below and relocate write and k lines
C
      character(len=240) :: line
      integer*4 k
      integer*4, external :: mexPrintf
C      write(line,*) 'size = ',size
C      k=mexPrintf(line//achar(13)) 

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .lt. 7) then
         call mexErrMsgIdAndTxt ('MATLAB:fine1:nInput',
     +                           'At least seven inputs required.')
      elseif(nlhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:fine1:nOutput',
     +                           'Three outputs required.')
      endif

C     Get the size of the input arrays.
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m*n

C     Column * row must be smaller than num
      if(size .gt. numel) then
         write(line,*) 'size = ',size
         k=mexPrintf(line//achar(13))
         call mexErrMsgIdAndTxt ('MATLAB:fine1:mSize',
     +                           'Row * column must be <= 500000.')
      endif
      
C     Check that the array is numeric (not strings).
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:NonNumeric',
     +                           'Input must be a numeric array.')
      endif

C     Create matrices for the return arguments.
      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      plhs(2) = mxCreateDoubleMatrix(m,n,0)
      plhs(3) = mxCreateDoubleMatrix(m,n,0)
      hr_8pr = mxGetPr(plhs(1))
      hz_8pr = mxGetPr(plhs(2))
      fl_8pr = mxGetPr(plhs(3))
      
C     Create pointers for input arguments and copy to inputs
      ri_8pr = mxGetPr(prhs(1))
      ro_8pr = mxGetPr(prhs(2))
      zl_8pr = mxGetPr(prhs(3))
      zu_8pr = mxGetPr(prhs(4))
      cu_8pr = mxGetPr(prhs(5))
      rp_8pr = mxGetPr(prhs(6))
      zp_8pr = mxGetPr(prhs(7))
                
      call mxCopyPtrToReal8(ri_8pr ,ri_8,size)
      call mxCopyPtrToReal8(ro_8pr ,ro_8,size)
      call mxCopyPtrToReal8(zl_8pr ,zl_8,size)
      call mxCopyPtrToReal8(zu_8pr ,zu_8,size)
      call mxCopyPtrToReal8(cu_8pr ,cu_8,size)
      call mxCopyPtrToReal8(rp_8pr ,rp_8,size)
      call mxCopyPtrToReal8(zp_8pr ,zp_8,size)

C     Call the computational subroutine.
      call fine1 (size,ri_8, ro_8,
     +                 zl_8, zu_8, cu_8,
     +                 rp_8, zp_8,
     +                 hr_8, hz_8, fl_8)

C     Load the data into output pointers for communication to MATLAB.
      call mxCopyReal8ToPtr(hr_8,hr_8pr,size) 
      call mxCopyReal8ToPtr(hz_8,hz_8pr,size) 
      call mxCopyReal8ToPtr(fl_8,fl_8pr,size)	  

      return
      end


C-----------------------------------------------------------------------
C     Computational routine

      SUBROUTINE fine1 (num,ri_8,ro_8,zl_8,zu_8,cu_8,
     +                      rp_8,zp_8,hr_8,hz_8,fl_8)
c
c     routine to call fine many times for vectors of coils and field points
c     assumes the correct address is passed in arrays using %VAL()
c
c  INPUTS:
c      num                 = size of each vector
c      ri_8,ro_8,zl_8,zu_8 = inner,outer,lower,upper dimensions of coil [vector]
c      cu_8               = current [vector]
c      rp_8,zp_8           = receiving point [vector]
c
c  OUTPUTS:
c           hr_8, hz_8     = fields at points rp,zp [Amp/m]
c           fl_8           = flux at point [Vs/amu_o]
c                            amu_o= 4*pi*1e-7
c
c  CALLED BY:              fine1g.f = gateway rotuine to produce MEX file
c  CALLS:                  fine.f = group of routines to do magnetics
c
c GA Jim Leuer 4-95
C ---------------------------------------------------------------------
c
c     MATLAB real*8 (num) elements storage provided by calling routine

      real*8 ri_8(num), ro_8(num), zl_8(num), zu_8(num)
      real*8 cu_8(num), rp_8(num), zp_8(num)
      real*8 hr_8(num), hz_8(num), fl_8(num)
c
c     FORTRAN input temporary storage variables 

      real*4 ri_4, ro_4, zl_4, zu_4, cu_4
      real*4 rp_4, zp_4
      real*4 hr_4, hz_4, fl_4

c     some integer variables
      integer num,i

c      integer ithru

C     Debugging variables
C     Uncomment lines below and relocate write and k lines
C
C      character(len=240) :: line
C      integer*4 k
C      integer*4, external :: mexPrintf
C      write(line,*) 'num = ',num
C      k=mexPrintf(line//achar(13)) 

C ---------------------------------------------------------------------
c
c        data ithru/0/
c        ithru= ithru+1
c        print*,' '
c        print*,' ithru=',ithru
C
C ---------------------------------------------------------------------
C Main DO loop around FINE
C
c       imhere= 4

       do i=1,num
c       convert from real*8 to real*4
        ri_4=  ri_8(i)
        ro_4=  ro_8(i)
        zl_4=  zl_8(i)
        zu_4=  zu_8(i)
        cu_4=  cu_8(i)
        rp_4=  rp_8(i)
        zp_4=  zp_8(i)
     
        CALL FINE(m1,ri_4,ro_4,zl_4,zu_4,cu_4,rp_4,zp_4,
     1               hr_4, hz_4, fl_4)

C      write(line,*) 'ri_4 = ',ri_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'ro_4 = ',ro_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'zl_4 = ',zl_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'zu_4 = ',zu_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'cu_4 = ',cu_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'rp_4 = ',rp_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'zp_4 = ',zp_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'hr_4 = ',hr_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'hz_4 = ',hz_4
C      k=mexPrintf(line//achar(13))
C      write(line,*) 'fl_4 = ',fl_4
C      k=mexPrintf(line//achar(13))

c
c      use existing arrays to store new data (?could dump into ri,ro,zl?)

        hr_8(i)= hr_4
        hz_8(i)= hz_4
        fl_8(i)= fl_4

       end do
c

      RETURN
      END

