C***********************************************************************
      subroutine istfin(ntot,subj,m,ist,ifin)
C creates vectors of starting and finishing positions for subjects
      integer ntot,subj(ntot),m,ist(m),ifin(m),scur,icur,i
      scur=-999
      icur=0
      do 100 i=1,ntot
         if(subj(i).ne.scur) then
            icur=icur+1
            ist(icur)=i
            scur=subj(i)
         endif
 100  continue
      do 200 i=1,m-1
         ifin(i)=ist(i+1)-1
 200  continue
      ifin(m)=ntot
      return
      end
C**************************************************************************
      subroutine chfc(p,pw,s)
C overwrites s (upper tri) with cholesky factor
      integer p,pw,i,j,k
      double precision s(p,p),sum
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i)**2
 10      continue
         s(i,i)=dsqrt(s(i,i)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i)*s(k,j)
 30         continue
            s(i,j)=(s(i,j)-sum)/s(i,i)
 40      continue
 50   continue
      return
      end
C***********************************************************************
      subroutine chl(p,pw,m,s,l)
C overwrites lth layer of s (upper tri) with cholesky factor
      integer p,pw,m,l,i,j,k
      double precision s(p,p,m),sum
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i,l)**2
 10      continue
         s(i,i,l)=dsqrt(s(i,i,l)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i,l)*s(k,j,l)
 30         continue
            s(i,j,l)=(s(i,j,l)-sum)/s(i,i,l)
 40      continue
 50   continue
      return
      end
C***********************************************************************
      subroutine bkslv(p,pw,s)
C inverts an upper triangular matrix
      integer p,pw,i,j,k
      double precision s(p,p),sum
      s(1,1)=dble(1.)/s(1,1)
      do 10 k=2,pw
         s(k,k)=dble(1.)/s(k,k)
         do 5 j=1,k-1
            sum=dble(0.)
            do 3 i=j,k-1
               sum=sum+s(j,i)*s(i,k)
 3          continue
            s(j,k)=-sum*s(k,k)
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine mm(p,pw,wm,cm)
C calculates upper tri part of cm = wm%*%t(wm), where wm is upper tri
      integer p,pw,i,j,k
      double precision wm(p,p),cm(p,p),sum
      do 10 i=1,pw
         do 5 j=i,pw
            sum=dble(0.)
            do 2 k=max(i,j),pw
               sum=sum+wm(i,k)*wm(j,k)
 2          continue
            cm(i,j)=sum
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine bkslvl(p,pw,m,s,l)
C inverts an upper triangular matrix in layer l of s
      integer p,pw,m,l,i,j,k
      double precision s(p,p,m),sum
      s(1,1,l)=dble(1.)/s(1,1,l)
      do 10 k=2,pw
         s(k,k,l)=dble(1.)/s(k,k,l)
         do 5 j=1,k-1
            sum=dble(0.)
            do 3 i=j,k-1
               sum=sum+s(j,i,l)*s(i,k,l)
 3          continue
            s(j,k,l)=-sum*s(k,k,l)
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine mmul(p,pw,m,wm,l,cm)
C calculates upper tri part of cm = wm%*%t(wm), where wm is upper tri
C(works in layer l)
      integer p,pw,m,l,i,j,k
      double precision wm(p,p,m),cm(p,p),sum
      do 10 i=1,pw
         do 5 j=i,pw
            sum=dble(0.)
            do 2 k=max(i,j),pw
               sum=sum+wm(i,k,l)*wm(j,k,l)
 2          continue
            cm(i,j)=sum
 5       continue
 10   continue
      return
      end
C***********************************************************************
      subroutine chle(p,pw,m,s,l,err)
C overwrites lth layer of s (upper tri) with cholesky factor
C If it fails, err is set to one.
      integer p,pw,m,l,err,i,j,k
      double precision s(p,p,m),sum
      err=0
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i,l)**2
 10      continue
         if(sum.ge.s(i,i,l)) then
            err=1
            goto 999
         endif
         s(i,i,l)=dsqrt(s(i,i,l)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i,l)*s(k,j,l)
 30         continue
            s(i,j,l)=(s(i,j,l)-sum)/s(i,i,l)
 40      continue
 50   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine chfce(p,pw,s,err)
C Overwrites s (upper tri) with cholesky factor
c If s is not positive definite, err is set to one.
      integer p,pw,err,i,j,k
      double precision s(p,p),sum
      err=0
      do 50 i=1,pw
         sum=dble(0.)
         do 10 k=1,i-1
            sum=sum+s(k,i)**2
 10      continue
         if(sum.ge.s(i,i)) then
            err=1
            goto 999
         endif
         s(i,i)=dsqrt(s(i,i)-sum)
         do 40 j=i+1,pw
            sum=dble(0.)
            do 30 k=1,i-1
               sum=sum+s(k,i)*s(k,j)
 30         continue
            s(i,j)=(s(i,j)-sum)/s(i,i)
 40      continue
 50   continue
 999  continue
      return
      end
C***********************************************************************
      subroutine bdiag(a,m,sig)
C fills in elements of each layer of sig below the diagonal
      integer a,m,s,i,j
      double precision sig(a,a,m)
      do 10 s=1,m
         do 5 i=1,a
            do 1 j=1,i-1
               sig(i,j,s)=sig(j,i,s)
 1          continue
 5       continue
 10   continue
      return
      end
C***********************************************************************
C The following are used only by mlmem() and mlmembd().
C***********************************************************************
      subroutine prefem(ntot,subj,m,ist,ifin,pcol,pred,p,q,
     /     xcol,zcol,ztz,patt,nstar,nstari,wkpp,xtxinv,err)
C Preliminary manipulations for pan. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     ztz  = t(z_i)%*%z, i=1,...,m
C     nstar = total number of rows in y containing data
C     nstari = total number of rows in y_i containing data
C     wkpp = t(pred)%*%pred
C     xtxinv = inv(t(pred)%*%pred)
      implicit none
      integer ntot,subj(ntot),m,ist(m),ifin(m),pcol,
     /     p,q,xcol(p),zcol(q),patt(ntot),nstar,nstari(m),err,
     /     st,fin,s,i,j,k
      double precision pred(ntot,pcol),ztz(q,q,m),wkpp(p,p),
     /     xtxinv(p,p),sum
      call istfin(ntot,subj,m,ist,ifin)
      do 2 s=1,m
         st=ist(s)
         fin=ifin(s)
         nstari(s)=0
         do 1 i=st,fin
            if(patt(i).ne.0) nstari(s)=nstari(s)+1
 1       continue
 2    continue
      nstar=0
      do 10 i=1,ntot
         if(patt(i).ne.0) nstar=nstar+1
 10   continue
      do 100 s=1,m
         do 90 i=1,q
            do 80 j=i,q
               sum=dble(0.)
               do 60 k=ist(s),ifin(s)
                  if(patt(k).ne.0) then
                     sum=sum+pred(k,zcol(i))*pred(k,zcol(j))
                  endif
 60            continue
               ztz(i,j,s)=sum
               if(i.ne.j) ztz(j,i,s)=sum
 80         continue
 90      continue
 100  continue
      do 150 i=1,p
         do 140 j=i,p
            sum=dble(0.)
            do 130 k=1,ntot
               if(patt(k).ne.0) then
                  sum=sum+pred(k,xcol(i))*pred(k,xcol(j))
               endif
 130        continue
            wkpp(i,j)=sum
 140     continue
 150  continue
      call chfce(p,p,wkpp,err)
      if(err.eq.1) then
         goto 999
      endif
      call bkslv(p,p,wkpp)
      call mm(p,p,wkpp,xtxinv)
      do 160 i=1,p
         do 155 j=i,p
            xtxinv(j,i)=xtxinv(i,j)
 155     continue
 160  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mimpy(ntot,r,y,patt,npatt,rmat)
C unconditional mean imputation for variables in y
      implicit none
      integer ntot,r,patt(ntot),npatt,rmat(npatt,r),denom,rij,j,i
      double precision y(ntot,r),sum,mean
      do 200 j=1,r
         sum=dble(0.)
         denom=0
         do 100 i=1,ntot
            if(patt(i).ne.0) then
               rij=rmat(patt(i),j)
               sum=sum+dfloat(rij)*y(i,j)
               denom=denom+rij
            endif
 100     continue
         mean=sum/dfloat(denom)
         do 150 i=1,ntot
            if(patt(i).ne.0) then
               if(rmat(patt(i),j).eq.0) y(i,j)=mean
            endif
 150     continue
 200  continue
      return
      end
C***********************************************************************
      subroutine mkbeta(p,r,xtxinv,xty,beta)
C calculates betahat
      integer p,r,i,j,k
      double precision xtxinv(p,p),xty(p,r),beta(p,r),sum
      do 100 i=1,p
         do 90 j=1,r
            sum=dble(0.)
            do 50 k=1,p
               sum=sum+xtxinv(i,k)*xty(k,j)
 50         continue
            beta(i,j)=sum
 90      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine mkxty(ntot,r,y,pcol,pred,p,xcol,patt,xty)
      implicit none
      integer ntot,r,pcol,p,xcol(p),patt(ntot),i,j,k
      double precision y(ntot,r),pred(ntot,pcol),xty(p,r),sum
      do 100 i=1,p
         do 80 j=1,r
            sum=dble(0.)
            do 50 k=1,ntot
               if(patt(k).ne.0) then
                  sum=sum+pred(k,xcol(i))*y(k,j)
               endif
 50         continue
            xty(i,j)=sum
 80      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,epsi,patt)
C calculates eps = y - X%*% beta
      implicit none
      integer ntot,r,patt(ntot),pcol,p,xcol(p),i,j,k
      double precision y(ntot,r),pred(ntot,pcol),beta(p,r),epsi(ntot,r),
     /     sum
      do 100 i=1,ntot
         if(patt(i).ne.0) then
            do 90 j=1,r
               sum=dble(0.)
               do 80 k=1,p
                  sum=sum+pred(i,xcol(k))*beta(k,j)
 80            continue
               epsi(i,j)=y(i,j)-sum   
 90         continue
         endif
 100  continue
      return
      end
C***********************************************************************
      subroutine mksigma(ntot,r,epsi,nstar,sigma,patt)
C calculates (1/nstar)*sum of t(eps_i)%*%eps_i
      implicit none
      integer ntot,r,patt(ntot),nstar,i,j,k
      double precision epsi(ntot,r),sigma(r,r)
      do 2 i=1,r
         do 1 j=i,r
            sigma(i,j)=dble(0.)
 1       continue
 2    continue
      do 100 i=1,ntot
         if(patt(i).ne.0) then
            do 90 j=1,r
               do 80 k=j,r
                  sigma(j,k)=sigma(j,k)+epsi(i,j)*epsi(i,k)
 80            continue
 90         continue
         endif
 100  continue
      do 110 i=1,r
         do 108 j=i,r
            sigma(i,j)=sigma(i,j)/dfloat(nstar)
            if(i.ne.j) sigma(j,i)=sigma(i,j)
 108     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mkpsi0(r,q,m,psi,sig,wkrqrq1)
C calculates initial estimate of psi given sig.
      integer r,q,m,s,i,j
      double precision psi(r*q,r*q),sig(r*q,r*q,m),wkrqrq1(r*q,r*q)
      do 5 i=1,r*q
         do 4 j=i,r*q
            psi(i,j)=dble(0.)
 4       continue
 5    continue
      do 100 s=1,m
         call mmul(r*q,r*q,m,sig,s,wkrqrq1)
         do 20 i=1,r*q
            do 10 j=i,r*q
               psi(i,j)=psi(i,j)+wkrqrq1(i,j)
 10         continue
 20      continue
 100  continue
      do 110 i=1,r*q
         do 105 j=i,r*q
            psi(i,j)=psi(i,j)/dfloat(m)
            if(i.ne.j) psi(j,i)=psi(i,j)
 105     continue
 110  continue
      return
      end
C***********************************************************************
      subroutine mksigbd(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,wkqq1bd,wkqq2bd)
C Version of mksig for the block-diagonal pan.
C calculates the square root of
C       sig = inv( inv(psi) + ( inv(sigma) Otimes ztz ) ), i=1,..m
C Use zflag = 1 to set psi equal to its prior guess.
C When finished, wkrr2 contains inv(sigma).
      implicit none
      integer r,q,m,s,l,i,j,ii,jj,ia,ja
      double precision psi(q,q,r),sigma(r,r),ztz(q,q,m),
     /     sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q),wkqq1bd(q,q,r),wkqq2bd(q,q),sum
      do 2 l=1,r
         sum=dble(0.)
         do 1 i=1,q
            sum=sum+dble(1.)
            psi(i,i,l)=sum
 1       continue
 2    continue
C invert psi and put into wkrqrq2, using the fact that it's block-diagonal
      do 102 i=1,r*q
         do 101 j=i,r*q
            wkrqrq2(i,j)=dble(0.)
 101     continue
 102  continue
      do 110 l=1,r
         do 105 i=1,q
            do 104 j=i,q
               wkqq1bd(i,j,l)=psi(i,j,l)
 104        continue
 105     continue
         call chl(q,q,r,wkqq1bd,l)
         call bkslvl(q,q,r,wkqq1bd,l)
         call mmul(q,q,r,wkqq1bd,l,wkqq2bd)
         do 108 i=1,q
            do 107 j=i,q
               ii=(l-1)*q + i
               jj=(l-1)*q + j
               wkrqrq2(ii,jj)=wkqq2bd(i,j)
 107        continue
 108     continue
 110  continue
C invert sigma and put into wkrr2
      do 7 i=1,r
         do 6 j=i,r
            wkrr1(i,j)=sigma(i,j)
 6       continue
 7    continue
      call chfc(r,r,wkrr1)
      call bkslv(r,r,wkrr1)
      call mm(r,r,wkrr1,wkrr2)
      do 9 i=1,r
         do 8 j=i+1,r
            wkrr2(j,i)=wkrr2(i,j)
 8       continue
 9    continue
      do 100 s=1,m
C initialize sig(,,s) to inv(sigma) Otimes t(Z_i)%*%Z_i
         do 30 i=1,r
            do 20 j=i,r
               do 15 ii=1,q
                  do 10 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     sig(ia,ja,s)=wkrr2(i,j)*ztz(ii,jj,s)
 10               continue
 15            continue
 20         continue
 30      continue
C add in inv(psi) and take inverse square-root
         do 40 i=1,r*q
            do 35 j=i,r*q
               sig(i,j,s)=sig(i,j,s)+wkrqrq2(i,j)
 35         continue
 40      continue
         call chl(r*q,r*q,m,sig,s)
         call bkslvl(r*q,r*q,m,sig,s)
 100  continue
      return
      end
C***********************************************************************
      subroutine mkpsi0bd(r,q,m,psi,sig,wkrqrq1)
C calculates initial estimate of psi given sig.
      implicit none
      integer r,q,m,s,l,i,j,ii,jj
      double precision psi(q,q,r),sig(r*q,r*q,m),wkrqrq1(r*q,r*q)
      do 6 l=1,r
         do 5 i=1,q
            do 4 j=i,q
               psi(i,j,l)=dble(0.)
 4          continue
 5       continue
 6    continue
      do 100 s=1,m
         call mmul(r*q,r*q,m,sig,s,wkrqrq1)
         do 30 l=1,r
            do 20 i=1,q
               do 10 j=i,q
                  ii=(l-1)*q + i
                  jj=(l-1)*q + j
                  psi(i,j,l)=psi(i,j,l)+wkrqrq1(ii,jj)
 10            continue
 20         continue
 30      continue
 100  continue
      do 120 l=1,r
         do 110 i=1,q
            do 105 j=i,q
               psi(i,j,l)=psi(i,j,l)/dfloat(m)
               if(i.ne.j) psi(j,i,l)=psi(i,j,l)
 105        continue
 110     continue
 120  continue
      return
      end
C***********************************************************************
      subroutine mksig(r,q,m,psi,sigma,ztz,sig,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2)
C calculates the square root of
C       sig = inv( inv(psi) + ( inv(sigma) Otimes ztz ) ), i=1,..m
C Use zflag = 1 to set psi equal to its prior guess.
C When finished, wkrr2 contains inv(sigma).
      implicit none
      integer r,q,m,s,i,j,ii,jj,ia,ja
      double precision psi(r*q,r*q),sigma(r,r),ztz(q,q,m),
     /     sig(r*q,r*q,m),wkrr1(r,r),wkrr2(r,r),wkrqrq1(r*q,r*q),
     /     wkrqrq2(r*q,r*q)
C invert psi and put into wkrqrq2
      do 2 i=1,r*q
         psi(i,i)=dble(1.)
         do 1 j=i,r*q
            wkrqrq1(i,j)=psi(i,j)
 1       continue
 2    continue
      call chfc(r*q,r*q,wkrqrq1)
      call bkslv(r*q,r*q,wkrqrq1)
      call mm(r*q,r*q,wkrqrq1,wkrqrq2)
C invert sigma and put into wkrr2
      do 6 i=1,r
         do 5 j=i,r
            wkrr1(i,j)=sigma(i,j)
 5       continue
 6    continue
      call chfc(r,r,wkrr1)
      call bkslv(r,r,wkrr1)
      call mm(r,r,wkrr1,wkrr2)
      do 9 i=1,r
         do 8 j=i+1,r
            wkrr2(j,i)=wkrr2(i,j)
 8       continue
 9    continue
      do 100 s=1,m
C initialize sig(,,s) to inv(sigma) Otimes t(Z_i)%*%Z_i
         do 30 i=1,r
            do 20 j=i,r
               do 15 ii=1,q
                  do 10 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     sig(ia,ja,s)=wkrr2(i,j)*ztz(ii,jj,s)
 10               continue
 15            continue
 20         continue
 30      continue
C add in inv(psi) and take inverse square-root
         do 40 i=1,r*q
            do 35 j=i,r*q
               sig(i,j,s)=sig(i,j,s)+wkrqrq2(i,j)
 35         continue
 40      continue
         call chl(r*q,r*q,m,sig,s)
         call bkslvl(r*q,r*q,m,sig,s)
 100  continue
      return
      end
C*********************************************************************
      subroutine mlmmembd(ntot,m,r,p,q,subj,ist,ifin,nmax,iposn,npatt,
     /     pstfin,patt,nstar,nstari,rmat,pcol,xcol,pred,zcol,w,wkqb2,
     /     vdel,wo,wo1,wm,wom,wkwmm1,wkwmm2,uszxb,usotzo,usotzm,wxbw,
     /     wxbwo,wxbwm,wkeb2,eb,wxbeta,wmusotzm,eybt,wxbetazeb,varb,
     /     wkrrpt,wkrrb21,y,eystar,ey,eyyt,eyxyxt,u,iter,msg,sigma,
     /     beta,ztz,xtw,xtwx,xtwy,xtwxinv,wkrqrq1,wkrqrq2,wkqq1bd,
     /     wkqq2bd,wkqq3,wkrr1,wkrr2,wksigtz,wkqqu,psi,wkqnm,wkeyxyxt,
     /     wkqnm1,cvgd,obeta,osigma,opsi,maxits,llvec,llovec,eps,ggs,
     /     reject,wkg,wkgg,sflag,epsi,wkpr,wkpp,xtxinv)
C Block-diagonal version of EM in mglmm
      implicit none
      integer ntot,m,r,p,q,subj(ntot),ist(m),ifin(m),
     /     nmax,iposn(ntot),
     /     npatt,pstfin(npatt,2),patt(ntot),nstar,nstari(m),
     /     rmat(npatt,r),pcol,zcol(q),xcol(p),err,iter,msg,
     /     cvgd,maxits,ggs,reject(iter),sflag,c1,c2,c3,i,j,l,
     /     nor(100),ormat(100,1000),nmr(100),mrmat(100,1000),
     /     oc(100),oc2(100),mc(100),mc1(100)
      double precision pred(ntot,pcol),w(r*nmax,r*nmax,m),wkqb2(nmax,r),
     /     vdel(r*nmax),wo(r*nmax,r*nmax),wo1(r*nmax,r*nmax),
     /     wm(r*nmax,r*nmax),wom(r*nmax,r*nmax),wkwmm1(r*nmax,r*nmax),
     /     wkwmm2(r*nmax,r*nmax),uszxb(r*q),usotzo(r*q,r*nmax),
     /     usotzm(r*q,r*nmax),wxbw(r*nmax),wxbwo(r*nmax),wxbwm(r*nmax),
     /     wkeb2(r*q,r*nmax),eb(r*q,m),wxbeta(ntot,r),
     /     wmusotzm(r*nmax,r*q),eybt(r*nmax,r*q),wxbetazeb(ntot,r),
     /     varb(r*q,r*q,m),wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),
     /     y(ntot,r),eystar(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),u(r*q,r*q,m),sigma(r,r),beta(p,r),
     /     ztz(q,q,m),xtw(p*r,nmax*r),xtwx(p*r,p*r),xtwy(p*r),
     /     xtwxinv(p*r,p*r),wkrqrq1(r*q,r*q),wkrqrq2(r*q,r*q),
     /     wkqq1bd(q,q,r),wkqq2bd(q,q),wkqq3(r*q,r*q),psi(q,q,r),
     /     wkrr1(r,r),wkrr2(r,r),
     /     wksigtz(r*q,r*nmax,m),wkqnm(r*q,r*nmax,m),
     /     wkqnm1(r*nmax,r*nmax),wkeyxyxt(r*nmax,r*nmax),obeta(p,r),
     /     opsi(q,q,r),osigma(r,r),ldsig,ldu,ldpsi,trdel,trdet,eps,
     /     ll,llo,llvec(maxits),wkg(ggs),wkgg(ggs,ggs),
     /     wkqqu(r*q,r*q,m),
     /     trwex,llovec(maxits),
     /     epsi(ntot,r),wkpr(p,r),wkpp(p,p),xtxinv(p,p)
      msg=0
      iter=0
      ll=dble(0.)
      llo=dble(0.)
      call prefem(ntot,subj,m,ist,ifin,pcol,pred,p,q,
     /     xcol,zcol,ztz,patt,nstar,nstari,wkpp,xtxinv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call mimpy(ntot,r,y,patt,npatt,rmat)
         call mkxty(ntot,r,y,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,epsi,patt)
         call mksigma(ntot,r,epsi,nstar,sigma,patt)
         call mksigbd(r,q,m,psi,sigma,ztz,wkqqu,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,wkqq1bd,wkqq2bd)
         call mkpsi0bd(r,q,m,psi,wkqqu,wkrqrq1)
      endif
      cvgd=0
C****************START OF MAIN ITERATION ****************************
 1    continue
      iter=iter+1
      reject(iter)=0
 3    continue
C     *** save the current parameters ***
      do 5 i=1,p
         do 4 j=1,r
            obeta(i,j)=beta(i,j)
 4       continue
 5    continue
      do 7 i=1,r
         do 6 j=i,r
            osigma(i,j)=sigma(i,j)
 6       continue
 7    continue
      do 10 l=1,r
         do 9 i=1,q
            do 8 j=i,q
               opsi(i,j,l)=psi(i,j,l)
 8          continue
 9       continue
 10   continue
C     *** calculate U_i, U_i%*%(inv(sigma) Otimes t(Z_i)), and W_i ***
      call  mkubd(r,q,m,psi,sigma,ztz,u,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,wkqq1bd,wkqq2bd,wkqqu,ldpsi,ldsig,ldu,err)
      if(err.eq.1) then
         msg=2
         goto 999
      endif
C     *** calculate W_i, inv of the cov of y_i ********
      call mkwkqnm(m,r,q,nmax,ntot,ist,ifin,pcol,zcol,patt,
     /     nstari,pred,wkrr1,u,wksigtz,wkqnm)
      call mkw(m,r,q,nmax,ntot,ist,ifin,
     /     nstari,patt,wkrr1,wksigtz,wkqnm,w)
C  *************************************************************************
      call mkeb(ntot,m,r,p,pcol,q,nmax,xcol,npatt,patt,
     /     ist,ifin,rmat,nor,ormat,nmr,mrmat,nstari,
     /     pred,beta,wkqnm,u,w,wkqb2,wo,wo1,wm,wom,wkwmm1,wkwmm2,
     /     uszxb,usotzo,usotzm,wxbw,wxbwo,wxbwm,wxbeta,wkeb2,y,
     /     err,msg,vdel,trdet,trdel,eb,varb)
      if(err.eq.1) then
         msg=3
         goto 999
      endif
C **** find the current observed loglikelihood ****
      llo=trdet-trdel/dble(2.)
      llovec(iter)=llo
      call mkey(100,100,oc,mc,m,r,ntot,nstari,iposn,npatt,pstfin,
     /     rmat,patt,p,q,xcol,zcol,pcol,ist,ifin,pred,
     /     y,beta,eb,sigma,wkrr1,wkrrpt,wkrrb21,eystar,ey)
      if(iter.gt.1) then
C        *** if observed likelihood has gone down, replace the scoring
C        *** use one step of EM 
         if(reject(iter-1).eq.0) then
            if(llovec(iter).lt.llovec(iter-1)) then
               call psiembd(m,r,q,varb,eb,psi)
               call mkxbeta(ntot,m,ist,ifin,p,r,pcol,xcol,
     /              patt,pred,beta,wxbeta)
               call sigmaem(ntot,nmax,m,r,pcol,q,zcol,ist,ifin,nstari,
     /              100,100,mc,mc1,oc,oc2,
     /              nstar,npatt,patt,rmat,nmr,mrmat,pred,wxbeta,y,ey,
     /              eyyt,wkqnm,usotzm,w,wm,wkwmm1,wkwmm2,wmusotzm,eb,
     /              varb,eybt,wkrrpt,wkrrb21,sigma,err)
               if(err.eq.1) then
                  msg=3
                  goto 999
               endif
               call gls(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /              nstari,patt,pred,w,ey,beta,xtw,xtwx,xtwy,
     /              xtwxinv,err)
               if(err.eq.1) then
                  msg=4
                  goto 999
               endif
               reject(iter-1)=1
               goto 3
            endif
         endif
      endif
C************************************************************************
C************ calculate the fixed effects *******************
      call gls(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /     nstari,patt,pred,w,ey,beta,xtw,xtwx,xtwy,
     /     xtwxinv,err)
      if(err.eq.1) then
         msg=4
         goto 999
      endif
CC*********************************************************************
C ********** PERFORM ONE STEP OF FISHER'S SCORING ON THE **********
C **********          VARIANCE COMPONENETS               **********
C*********************************************************************
      call preyxyxt(ntot,m,ist,ifin,p,q,r,pcol,xcol,zcol,
     /     patt,pred,beta,eb,wxbeta,wxbetazeb)
      call fscovbd(ntot,nmax,npatt,m,r,pcol,q,zcol,ist,ifin,nstari,
     /     patt,rmat,100,100,mc,mc1,oc,oc2,pred,u,ztz,psi,opsi,
     /     wkrr1,wkrr2,wkqnm,wkqnm1,wkeyxyxt,wkrqrq1,wkrqrq2,wkqq1bd,
     /     wkqq2bd,wkqq3,varb,wkrrpt,wkrrb21,wxbeta,ggs,wkg,wkgg,w,
     /     sigma,osigma,y,ey,eyyt,eyxyxt,trwex,msg)
C **** find the current expected loglikelihood ****
      ll = dfloat(nstar)*ldsig + dfloat(m)*ldpsi +ldu - trwex/dble(2.)
      llvec(iter)=ll
      if(msg.ne.0) goto 999
C********CHECK CONVERGENCE  *************************
      c1=0
      do 30 i=1,p
         do 27 j=1,r
            if(dabs(beta(i,j)-obeta(i,j)).gt.(eps*
     /           dabs(obeta(i,j)))) c1=1
 27      continue
 30   continue
      c2=0
      do 50 l=1,r
         do 45 i=1,q
            do 40 j=i,q
               if(dabs(psi(i,j,l)-opsi(i,j,l)).gt.(eps*
     /              dabs(opsi(i,j,l)))) c2=1
 40         continue
 45      continue
 50   continue
      c3=0
      do 60 i=1,r
         do 55 j=i,r
            if(dabs(sigma(i,j)-osigma(i,j)).gt.(eps*
     /           dabs(osigma(i,j)))) c3=1
 55      continue
 60   continue
      if((c1.eq.0).and.(c2.eq.0).and.(c3.eq.0)) cvgd=1
      if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
C********* end of main iteration *****************
      call bdiag(r*q,m,u)
      do 70 i=1,p*r
         do 65 j=i+1,p*r
            xtwxinv(j,i)=xtwxinv(i,j)
 65      continue
 70   continue
 999  continue
      iter=iter
      msg=msg
      return
      end
C*********************************************************************
       subroutine mlmmem(ntot,m,r,p,q,subj,ist,ifin,nmax,iposn,npatt,
     /     pstfin,patt,nstar,nstari,rmat,pcol,xcol,pred,zcol,w,wkqb2,
     /     vdel,wo,wo1,wm,wom,wkwmm1,wkwmm2,uszxb,usotzo,usotzm,wxbw,
     /     wxbwo,wxbwm,wkeb2,eb,wxbeta,wmusotzm,eybt,wxbetazeb,varb,
     /     wkrrpt,wkrrb21,y,eystar,ey,eyyt,eyxyxt,u,iter,msg,sigma,
     /     beta,ztz,xtw,xtwx,xtwy,xtwxinv,wkqq1,wkqq2,wkqq3,wkrr1,wkrr2,
     /     wksigtz,wkqqu,psi,wkqnm,wkeyxyxt,wkqnm1,cvgd,obeta,osigma,
     /     opsi,maxits,llvec,llovec,eps,ggs,reject,wkg,wkgg,sflag,epsi,
     /     wkpr,wkpp,xtxinv)
C  SET sflag=1 if starting values supplied 
C  Returned error codes in msg:
C  0 = no errors
C  1 = t(X)%*%X is not positive definite ( for calculating the starting
C  value of beta)
C  2 = value of psi or sigma or inv(u_i)became 
C      non-positive definite during iterations
C  3 = value of Var(y_i(obs)) became non-positive
C      during the iterations
C  4 = t(X)%*%W%*%X became non-pos.def. during iterations
C  5 = non-positive psi at input scoring step               
C  6 = non-positive sigma at input scoring step                          
C  7 = non-positive wkgg in the scoring step
      implicit none
      integer ntot,m,r,p,q,subj(ntot),ist(m),ifin(m),nmax,iposn(ntot),
     /     npatt,pstfin(npatt,2),patt(ntot),nstar,nstari(m),
     /     rmat(npatt,r),pcol,zcol(q),xcol(p),err,iter,msg,cvgd,maxits,
     /     ggs,reject(iter),c1,c2,c3,i,j,nor(100),ormat(100,1000),
     /     nmr(100),mrmat(100,1000),oc(100),oc2(100),mc(100),mc1(100),
     /     sflag
      double precision pred(ntot,pcol),w(r*nmax,r*nmax,m),wkqb2(nmax,r),
     /     vdel(r*nmax),wo(r*nmax,r*nmax),wo1(r*nmax,r*nmax),
     /     wm(r*nmax,r*nmax),wom(r*nmax,r*nmax),wkwmm1(r*nmax,r*nmax),
     /     wkwmm2(r*nmax,r*nmax),uszxb(r*q),usotzo(r*q,r*nmax),
     /     usotzm(r*q,r*nmax),wxbw(r*nmax),wxbwo(r*nmax),wxbwm(r*nmax),
     /     wkeb2(r*q,r*nmax),eb(r*q,m),wxbeta(ntot,r),
     /     wmusotzm(r*nmax,r*q),eybt(r*nmax,r*q),wxbetazeb(ntot,r),
     /     varb(r*q,r*q,m),wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),
     /     y(ntot,r),eystar(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),u(r*q,r*q,m),sigma(r,r),
     /     beta(p,r),ztz(q,q,m),xtw(p*r,nmax*r),xtwx(p*r,p*r),
     /     xtwy(p*r),xtwxinv(p*r,p*r),wkqq1(r*q,r*q),wkqq2(r*q,r*q),
     /     wkqq3(r*q,r*q),psi(r*q,r*q),wkrr1(r,r),wkrr2(r,r),
     /     wksigtz(r*q,r*nmax,m),wkqnm(r*q,r*nmax,m),
     /     wkqnm1(r*nmax,r*nmax),wkeyxyxt(r*nmax,r*nmax),obeta(p,r),
     /     opsi(r*q,r*q),osigma(r,r),ldsig,ldu,ldpsi,ll,eps,
     /     llvec(maxits),llovec(maxits),wkg(ggs),wkgg(ggs,ggs),
     /     wkqqu(r*q,r*q,m),
     /     epsi(ntot,r),wkpr(p,r),wkpp(p,p),xtxinv(p,p),
     /     trwex,ldsi,ldps,lduu,llo,trdet,trdel
      msg=0
      iter=0
      ll=dble(0.)
      llo=dble(0.)
      call prefem(ntot,subj,m,ist,ifin,pcol,pred,p,q,
     /     xcol,zcol,ztz,patt,nstar,nstari,wkpp,xtxinv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call mimpy(ntot,r,y,patt,npatt,rmat)
         call mkxty(ntot,r,y,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,epsi,patt)
         call mksigma(ntot,r,epsi,nstar,sigma,patt)
         call mksig(r,q,m,psi,sigma,ztz,wkqqu,wkrr1,wkrr2,wkqq1,
     /        wkqq1)
         call mkpsi0(r,q,m,psi,wkqqu,wkqq1)
      endif
      cvgd=0
C****************START OF MAIN ITERATION ****************************
 1    continue
      iter=iter+1
      reject(iter)=0
 3    continue
C     *** save the current parameters ***
      do 5 i=1,p
         do 4 j=1,r
            obeta(i,j)=beta(i,j)
 4       continue
 5    continue
      do 7 i=1,r
         do 6 j=i,r
            osigma(i,j)=sigma(i,j)
 6       continue
 7    continue
      do 9 i=1,r*q
         do 8 j=i,r*q
            opsi(i,j)=psi(i,j)
 8       continue
 9    continue
C     *** calculate U_i, U_i%*%(inv(sigma) Otimes t(Z_i)), and W_i ***
      call  mku(r,q,m,psi,sigma,ztz,u,wkrr1,wkrr2,wkqq1,
     /     wkqq2,wkqqu,ldpsi,ldsig,ldu,err)
      if(err.eq.1) then
         msg=2
         goto 999
      endif
C     *** calculate W_i, inv of the cov of y_i ********
      call mkwkqnm(m,r,q,nmax,ntot,ist,ifin,pcol,zcol,patt,
     /     nstari,pred,wkrr1,u,wksigtz,wkqnm)
      call mkw(m,r,q,nmax,ntot,ist,ifin,
     /     nstari,patt,wkrr1,wksigtz,wkqnm,w)
C  *************************************************************************
      call mkeb(ntot,m,r,p,pcol,q,nmax,xcol,npatt,patt,
     /     ist,ifin,rmat,nor,ormat,nmr,mrmat,nstari,
     /     pred,beta,wkqnm,u,w,wkqb2,wo,wo1,wm,wom,wkwmm1,wkwmm2,
     /     uszxb,usotzo,usotzm,wxbw,wxbwo,wxbwm,wxbeta,wkeb2,y,
     /     err,msg,vdel,trdet,trdel,eb,varb)
      if(err.eq.1) then
         msg=3
         goto 999
      endif
C **** find the current observed loglikelihood ****
      llo=trdet-trdel/dble(2.)
      llovec(iter)=llo
      call mkey(100,100,oc,mc,m,r,ntot,nstari,iposn,npatt,pstfin,
     /     rmat,patt,p,q,xcol,zcol,pcol,ist,ifin,pred,
     /     y,beta,eb,sigma,wkrr1,wkrrpt,wkrrb21,eystar,ey)
      if(iter.gt.1) then
C        *** if observed likelihood has gone down, replace the scoring
C        *** use one step of EM 
         if(reject(iter-1).eq.0) then
            if(llovec(iter).lt.llovec(iter-1)) then
               call psiem(m,r,q,varb,eb,psi)
               call mkxbeta(ntot,m,ist,ifin,p,r,pcol,xcol,
     /              patt,pred,beta,wxbeta)
               call sigmaem(ntot,nmax,m,r,pcol,q,zcol,ist,ifin,nstari,
     /              100,100,mc,mc1,oc,oc2,
     /              nstar,npatt,patt,rmat,nmr,mrmat,pred,wxbeta,y,ey,
     /              eyyt,wkqnm,usotzm,w,wm,wkwmm1,wkwmm2,wmusotzm,eb,
     /              varb,eybt,wkrrpt,wkrrb21,sigma,err)
               if(err.eq.1) then
                  msg=3
                  goto 999
               endif
               call gls(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /              nstari,patt,pred,w,ey,beta,xtw,xtwx,xtwy,
     /              xtwxinv,err)
               if(err.eq.1) then
                  msg=4
                  goto 999
               endif
               reject(iter-1)=1
               goto 3
            endif
         endif
      endif
C************************************************************************
C************ calculate the fixed effects *******************
      call gls(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /     nstari,patt,pred,w,ey,beta,xtw,xtwx,xtwy,
     /     xtwxinv,err)
      if(err.eq.1) then
         msg=4
         goto 999
      endif
CC*********************************************************************
C ********** PERFORM ONE STEP OF FISHER'S SCORING ON THE **********
C **********          VARIANCE COMPONENETS               **********
C*********************************************************************
      call preyxyxt(ntot,m,ist,ifin,p,q,r,pcol,xcol,zcol,
     /     patt,pred,beta,eb,wxbeta,wxbetazeb)
      call fscov(ntot,nmax,npatt,m,r,pcol,q,zcol,ist,ifin,nstari,
     /     patt,rmat,100,100,mc,mc1,oc,oc2,pred,u,ztz,psi,opsi,
     /     wkrr1,wkrr2,wkqnm,wkqnm1,wkeyxyxt,wkqqu,wkqq1,wkqq2,wkqq3,
     /     varb,wkrrpt,wkrrb21,wxbeta,ggs,wkg,wkgg,w,sigma,osigma,
     /     y,ey,eyyt,eyxyxt,trwex,msg,ldsi,ldps,lduu)
C **** find the current expected loglikelihood ****
      ll = dfloat(nstar)*ldsig + dfloat(m)*ldpsi +ldu - trwex/dble(2.)
      llvec(iter)=ll
      if(msg.ne.0) goto 999
C********  CHECK CONVERGENCE  *************************
      c1=0
      do 30 i=1,p
         do 27 j=1,r
            if(dabs(beta(i,j)-obeta(i,j)).gt.(eps*
     /           dabs(obeta(i,j)))) c1=1
 27      continue
 30   continue
      c2=0
      do 40 i=1,r*q
         do 35 j=i,r*q
            if(dabs(psi(i,j)-opsi(i,j)).gt.(eps*
     /           dabs(opsi(i,j)))) c2=1
 35      continue
 40   continue
      c3=0
      do 50 i=1,r
         do 45 j=i,r
            if(dabs(sigma(i,j)-osigma(i,j)).gt.(eps*
     /           dabs(osigma(i,j)))) c3=1
 45      continue
 50   continue
      if((c1.eq.0).and.(c2.eq.0).and.(c3.eq.0)) cvgd=1
      if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
C********* end of main iteration *****************
      call bdiag(r*q,m,u)
      do 70 i=1,p*r
         do 60 j=i+1,p*r
            xtwxinv(j,i)=xtwxinv(i,j)
 60      continue
 70   continue
 999  continue
      iter=iter
      msg=msg
      return
      end
C*********************************************************************
      subroutine psiem(m,r,q,varb,eb,psi)
C replaces psi with its EM estimate
      implicit none
      integer m,r,q,i,j,s
      double precision varb(r*q,r*q,m),eb(r*q,m),psi(r*q,r*q)
      do 5 i=1,r*q
         do 4 j=i,r*q
            psi(i,j)=dble(0.)
 4       continue
 5    continue
      do 100 s=1,m
         do 50 i=1,r*q
            do 40 j=i,r*q
               psi(i,j)=psi(i,j)+varb(i,j,s)+eb(i,s)*eb(j,s)
 40         continue
 50      continue
 100  continue
      do 110 i=1,r*q
         do 105 j=i,r*q
            psi(i,j)=psi(i,j)/dfloat(m)
            if(i.ne.j) psi(j,i)=psi(i,j)
 105     continue
 110  continue
      return
      end
C*********************************************************************
      subroutine psiembd(m,r,q,varb,eb,psi)
C replaces psi with its EM estimate
      implicit none
      integer m,r,q,l,i,j,s,ii,jj
      double precision varb(r*q,r*q,m),eb(r*q,m),psi(q,q,r)
      do 5 l=1,r
         do 4 i=1,q
            do 3 j=i,q
               psi(i,j,l)=dble(0.)
 3          continue
 4       continue
 5    continue
      do 100 s=1,m
         do 90 l=1,r
            do 50 i=1,q
               do 40 j=i,q
                  ii=(l-1)*q + i
                  jj=(l-1)*q + j
                  psi(i,j,l)=psi(i,j,l)+varb(ii,jj,s)+eb(ii,s)*eb(jj,s)
 40            continue
 50         continue
 90      continue
 100  continue
      do 120 l=1,r
         do 110 i=1,q
            do 105 j=i,q
               psi(i,j,l)=psi(i,j,l)/dfloat(m)
               if(i.ne.j) psi(j,i,l)=psi(i,j,l)
 105        continue
 110     continue
 120  continue
      return
      end
C*********************************************************************
      subroutine sigmaem(ntot,nmax,m,r,pcol,q,zcol,ist,ifin,nstari,
     /     lmc,loc,mc,mc1,oc,oc2,
     /     nstar,npatt,patt,rmat,nmr,mrmat,pred,wxbeta,y,ey,eyyt,
     /     wkqnm,usotzm,w,wm,wkwmm1,wkwmm2,wmusotzm,eb,varb,eybt,
     /     wkrrpt,wkrrb21,emsigma,err)
C replaces sigma with its EM estimate
      implicit none
      integer ntot,nmax,m,r,pcol,q,zcol(q),ist(m),ifin(m),nstari(m),
     /     nstar,npatt,patt(ntot),rmat(npatt,r),nmr(r),mrmat(r,nmax),
     /     err,lmc,loc,mc(lmc),mc1(lmc),oc(loc),oc2(loc),
     /     st,fin,ni,i,j,s,ii,gi,kk,jj,ia,ja
      double precision pred(ntot,pcol),wxbeta(ntot,r),y(ntot,r),
     /     ey(ntot,r),eyyt(r*nmax,r*nmax),wkqnm(r*q,r*nmax,m),
     /     usotzm(r*q,r*nmax),w(r*nmax,r*nmax,m),wm(r*nmax,r*nmax),
     /     wkwmm1(r*nmax,r*nmax),wkwmm2(r*nmax,r*nmax),
     /     wmusotzm(r*nmax,r*q),eb(r*q,m),varb(r*q,r*q,m),
     /     eybt(r*nmax,r*q),wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),
     /     emsigma(r,r),sum,sum1,sum2,sum22,sum3,sum4
      do 500 i=1,r
         do 450 j=i,r
            sum=dble(0.)
            sum1=dble(0.)
            sum2=dble(0.)
            sum22=dble(0.)
            sum3=dble(0.)
            sum4=dble(0.)
            do 400 s=1,m
               st=ist(s)
               fin=ifin(s)
               ni=nstari(s)
               call mkyyt(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,
     /              y,eyyt)
               call mkeyyt(ntot,nmax,npatt,m,r,st,fin,ni,patt,
     /              rmat,s,pcol,q,zcol,lmc,loc,mc,mc1,oc,oc2,
     /              pred,varb,wkrrpt,wkrrb21,ey,eyyt)
               call mkeybt(ntot,m,r,q,ni,nmax,npatt,patt,rmat,s,st,fin,
     /              nmr,mrmat,wkqnm,usotzm,w,wm,wkwmm1,wkwmm2,wmusotzm,
     /              y,ey,eb,eybt,err)
               gi=0
               do 390 ii=st,fin
                  if(patt(ii).ne.0) then
                     gi=gi+1
                     sum=sum+eyyt((i-1)*ni+gi,(j-1)*ni+gi)
                  endif
 390           continue
               gi=0
               do 385 kk=st,fin
                  if(patt(kk).ne.0) then
                     gi=gi+1
                     do 380 jj=1,q
                        sum1=sum1+pred(kk,zcol(jj))*
     /                       eybt((i-1)*ni+gi,(j-1)*q+jj)
 380                 continue
                  endif
 385           continue
               gi=0
               do 375 kk=st,fin
                  if(patt(kk).ne.0) then
                     gi=gi+1
                     do 370 jj=1,q
                        sum1=sum1+pred(kk,zcol(jj))*
     /                       eybt((j-1)*ni+gi,(i-1)*q+jj)
 370                 continue
                  endif
 375           continue
               do 365 jj=st,fin
                  if(patt(jj).ne.0) then
                     sum22=dble(0.)
                     do 360 kk=1,q
                        sum22=sum22+pred(jj,zcol(kk))*eb((i-1)*q+kk,s)
 360                 continue
                     sum2=sum2+wxbeta(jj,j)*(ey(jj,i)-sum22)
                  endif
 365           continue
               do 355 jj=st,fin
                  if(patt(jj).ne.0) then
                     sum22=dble(0.)
                     do 350 kk=1,q
                        sum22=sum22+pred(jj,zcol(kk))*eb((j-1)*q+kk,s)
 350                 continue
                     sum2=sum2+wxbeta(jj,i)*(ey(jj,j)-sum22)
                  endif
 355           continue
               do 340 ia=st,fin
                  if(patt(ia).ne.0) then
                     do 330 ja=1,q
                        do 320 jj=1,q
                           sum3=sum3+pred(ia,zcol(ja))*
     /                          pred(ia,zcol(jj))*(
     /                          varb((i-1)*q+ja,(j-1)*q+jj,s)+
     /                          eb((i-1)*q+ja,s)*eb((j-1)*q+jj,s))
 320                    continue
 330                 continue
                  endif
 340           continue
               do 300 ia=st,fin
                  if(patt(ia).ne.0) then
                     sum4=sum4+wxbeta(ia,i)*wxbeta(ia,j)
                  endif
 300           continue
 400        continue
            emsigma(i,j)=(sum-sum1-sum2+sum3+sum4)/dfloat(nstar)
            if(i.ne.j) emsigma(j,i)=emsigma(i,j)
 450     continue
 500  continue
      return
      end
C*********************************************************************
      subroutine mkeybt(ntot,m,r,q,ni,nmax,npatt,patt,rmat,s,st,fin,nmr,
     /     mrmat,wkqnm,usotzm,w,wm,wkwmm1,wkwmm2,wmusotzm,y,ey,eb,eybt,
     /     err)
C calculates E( y^v_i (b^v_i)^T | y_i(obs), theta ) for subject i
      implicit none
      integer ntot,m,r,q,ni,nmax,npatt,patt(ntot),rmat(npatt,r),s,st,
     /     fin,nmr(r),mrmat(r,nmax),err,pdwm,i,j,k,posn,indi
      double precision wkqnm(r*q,r*nmax,m),usotzm(r*q,r*nmax),
     /     w(r*nmax,r*nmax,m),wm(r*nmax,r*nmax),wkwmm1(r*nmax,r*nmax),
     /     wkwmm2(r*nmax,r*nmax),wmusotzm(r*nmax,r*q),y(ntot,r),
     /     ey(ntot,r),
     /     eb(r*q,m),eybt(r*nmax,r*q),sum
      do 3 i=1,r*nmax
         do 2 j=i,r*q
            eybt(i,j)=dble(0.)
 2       continue
         do 1 j=i,r*nmax
            wkwmm1(i,j)=dble(0.)
            wkwmm2(i,j)=dble(0.)
 1       continue
 3    continue
      call getmrmat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nmr,mrmat)
      call mkusotzm(r,m,q,nmax,s,mrmat,nmr,ni,wkqnm,usotzm)
      call mkwm(m,r,nmax,s,mrmat,nmr,ni,pdwm,w,wm)
      do 5 i=1,pdwm
         do 4 j=i,pdwm
            wkwmm1(i,j)=wm(i,j)
 4       continue
 5    continue
      call chfce(r*nmax,pdwm,wkwmm1,err)
      if(err.eq.1) then 
         goto 999
      endif
      call bkslv(r*nmax,pdwm,wkwmm1)
      call mm(r*nmax,pdwm,wkwmm1,wkwmm2)
      do 11 i=1,pdwm
         do 10 j=1,r*q
            sum=dble(0.)
            do 9 k=1,i
               sum=sum+wkwmm2(k,i)*usotzm(j,k)
 9          continue
            do 8 k=i+1,pdwm
               sum=sum+wkwmm2(i,k)*usotzm(j,k)
 8          continue
            wmusotzm(i,j)=sum
 10      continue
 11   continue
      do 200 k=1,r*q
         posn=0
         do 150 i=1,r
            indi=0
            do 100 j=st,fin
               if(patt(j).ne.0) then
                  indi=indi+1
                  if(rmat(patt(j),i).eq.0) then
                     posn=posn+1
                     eybt((i-1)*ni+indi,k)=wmusotzm(posn,k)+
     /                    ey(j,i)*eb(k,s)
                  else
                     eybt((i-1)*ni+indi,k)=y(j,i)*eb(k,s)
                  endif
               endif
 100        continue
 150     continue
 200  continue
 999  continue
      return
      end
C*********************************************************************
      subroutine mkeybt1(ntot,m,r,q,ni,nmax,npatt,patt,rmat,s,st,fin,
     /     nor,ormat,nmr,pcol,zcol,pred,psi,
     /     mrmat,wkqnm,usotzm,w,wom,wm,rozpsi,rozpsio,rozpsim,
     /     y,ey,eb,wmcov,eybt,err)
C calculates E( y^v_i (b^v_i)^T | y_i(obs), theta ) for subject i
      implicit none
      integer ntot,m,r,q,ni,nmax,npatt,patt(ntot),rmat(npatt,r),s,st,
     /     fin,nor(r),ormat(r,nmax),nmr(r),pcol,zcol(q),mrmat(r,nmax),
     /     err,pdwm,pdwo,i,j,k,indi,kk
      double precision pred(ntot,pcol),psi(r*q,r*q),wkqnm(r*q,r*nmax,m),
     /     usotzm(r*q,r*nmax),
     /     w(r*nmax,r*nmax,m),wm(r*nmax,r*nmax),wom(r*nmax,r*nmax),
     /     rozpsi(r*nmax,r*q),
     /     rozpsio(r*nmax,r*q),rozpsim(r*nmax,r*q),
     /     y(ntot,r),
     /     ey(ntot,r),
     /     eb(r*q,m),wmcov(r*nmax,r*q),eybt(r*nmax,r*q),sum
      do 3 i=1,r*nmax
         do 2 j=i,r*q
            eybt(i,j)=dble(0.)
            wmcov(i,j)=dble(0.)
 2       continue
 3    continue
      call getmrmat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nmr,mrmat)
      call getormat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nor,ormat)
      call mkrozpsi(ntot,nmax,r,q,pcol,zcol,ni,st,fin,
     /     patt,pred,psi,rozpsi)
      call mkrozpsim(r,q,nmax,mrmat,nmr,ni,rozpsi,rozpsim)
      call mkrozpsio(r,q,nmax,ormat,nor,ni,rozpsi,rozpsio)
      call mkwm(m,r,nmax,s,mrmat,nmr,ni,pdwm,w,wm)
      call mkwom(m,r,nmax,s,mrmat,ormat,nmr,nor,ni,pdwo,pdwm,
     /     w,wom)
      do 11 i=1,pdwm
         do 10 j=1,r*q
            sum=dble(0.)
            do 9 kk=1,pdwo
               do 8 k=1,i
                  sum=sum+rozpsio(kk,j)*wm(k,i)*wom(kk,k)
 8             continue
               do 7 k=i+1,pdwm
                  sum=sum+rozpsio(kk,j)*wm(i,k)*wom(kk,k)
 7             continue
 9          continue
            wmcov(i,j)=rozpsim(i,j)-sum
 10      continue
 11   continue
      do 100 i=1,r
         indi=0
         do 90 j=st,fin
            if(patt(j).ne.0) then
               indi=indi+1
               do 80 k=1,r*q
                  if(rmat(patt(j),i).eq.0) then
                     eybt((i-1)*ni+indi,k)=wmcov((i-1)*ni+indi,k)+
     /                    ey(j,i)*eb(k,s)
                  else
                     eybt((i-1)*ni+indi,k)=y(j,i)*eb(k,s)
                  endif
 80            continue
            endif
 90      continue
 100  continue
 999  continue
      return
      end
C*********************************************************************
      subroutine mkrozpsi(ntot,nmax,r,q,pcol,zcol,ni,st,fin,
     /     patt,pred,psi,rozpsi)
      implicit none
      integer ntot,nmax,r,q,pcol,zcol(q),ni,st,fin,patt(ntot),i,j,ii,
     /     indi,ia,ja,jj,k
      double precision pred(ntot,pcol),psi(r*q,r*q),rozpsi(r*nmax,r*q),
     /     sum
      do 10 i=1,r*nmax
         do 5 j=1,r*q
            rozpsi(i,j)=dble(0.)
 5       continue
 10   continue
      do 200 i=1,r
         do 150 ii=1,r
            indi=0
            do 100 j=st,fin
               if(patt(j).ne.0) then
                  indi=indi+1
                  ia=(i-1)*ni+indi
                  do 90 jj=1,q
                     ja=(ii-1)*q+jj
                     sum=dble(0.)
                     do 80 k=1,q
                        sum=sum+pred(j,zcol(k))*psi((i-1)*q+k,
     /                       (ii-1)*q+jj)
 80                  continue
                     rozpsi(ia,ja)=sum
 90               continue
               endif
 100        continue
 150     continue
 200  continue
      return
      end
C*********************************************************************
C The following subroutines are required for the E-step
C*********************************************************************
      subroutine getormat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nor,ormat)
C creates a matrix(ormat) of r by nmax that indicates which rows are observed 
C for given column, and a vector of r indicating how many rows are 
C observed for each column
      implicit none
      integer ntot,r,nmax,npatt,patt(ntot),rmat(npatt,r),
     /     st,fin,posn,ormat(r,nmax),nor(r),i,j,obsr
      do 2 i=1,r
         do 1 j=1,nmax
            ormat(i,j)=0
 1       continue
 2    continue
      do 10 i=1,r
         nor(i)=0
         posn=0
         obsr=0
         do 5 j=st,fin
            if(patt(j).ne.0) then
               obsr=obsr+1
               if(rmat(patt(j),i).eq.1)then
                  nor(i)=nor(i)+1
                  posn=posn+1
                  ormat(i,posn)=obsr
               endif
            endif
 5       continue
 10   continue
      return
      end
C*********************************************************************
      subroutine getmrmat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nmr,mrmat)
C does the same the for missing rows
      implicit none
      integer ntot,r,nmax,npatt,patt(ntot),rmat(npatt,r),
     /     st,fin,posn,mrmat(r,nmax),nmr(r),i,j,misr
      do 2 i=1,r
         do 1 j=1,nmax
            mrmat(i,j)=0
 1       continue
 2    continue
      do 10 i=1,r
         nmr(i)=0
         posn=0
         misr=0
         do 5 j=st,fin
            if(patt(j).ne.0) then
               misr=misr+1
               if(rmat(patt(j),i).eq.0)then
                  nmr(i)=nmr(i)+1
                  posn=posn+1
                  mrmat(i,posn)=misr
               endif
            endif
 5       continue
 10   continue
      return
      end
C*********************************************************************
      subroutine mkwo(m,r,nmax,s,ormat,nor,ni,pdwo,w,wo)
C takes the elements of W_i corresponding to observed variables and puts
C them into wo. Note : pdwo is the dimension of wo
      implicit none
      integer m,r,nmax,s,ormat(r,nmax),nor(r),ni,
     /     pdwo,posni,i,j,posnj,ii,jj
      double precision w(r*nmax,r*nmax,m),wo(r*nmax,r*nmax)
      do 10 i=1,r*nmax
         do 5 j=i,r*nmax
            wo(i,j)=dble(0.)
 5       continue
 10   continue
      pdwo=0
      posni=0
      do 50 i=1,r
         do 40 j=1,nor(i)
            posni=posni+1
            posnj=0
            do 30 ii=1,r
               do 20 jj=1,nor(ii)
                  posnj=posnj+1
                  wo(posni,posnj)=w((i-1)*ni+
     /                 ormat(i,j),(ii-1)*ni+ormat(ii,jj),s)
 20            continue
 30         continue
 40      continue
 50   continue
      pdwo=posni
      return
      end
C*********************************************************************
      subroutine mkwm(m,r,nmax,s,mrmat,nmr,ni,pdwm,w,wm)
C takes the elements of W_i corresponding to missing variables and puts
C them into wm. Note : pdwm is the dimension of wm
      implicit none
      integer m,r,nmax,s,posni,posnj,nmr(r),ni,mrmat(r,nmax),pdwm,
     /     i,j,ii,jj
      double precision w(r*nmax,r*nmax,m),wm(r*nmax,r*nmax) 
      do 10 i=1,r*nmax
         do 5 j=i,r*nmax
            wm(i,j)=dble(0.)
 5       continue
 10   continue
      pdwm=0
      posni=0
      do 50 i=1,r
         do 40 j=1,nmr(i)
            posni=posni+1
            posnj=0
            do 30 ii=1,r
               do 20 jj=1,nmr(ii)
                  posnj=posnj+1
                  wm(posni,posnj)=w((i-1)*ni+
     /                 mrmat(i,j),(ii-1)*ni+mrmat(ii,jj),s)
 20            continue
 30         continue
 40      continue
 50   continue
      pdwm=posni
      return
      end
C*********************************************************************
      subroutine mkwom(m,r,nmax,s,mrmat,ormat,nmr,nor,ni,pdwo,pdwm,
     /     w,wom)
C forms the partition of W_i that is the cov matrix of observed and
C missing variables storing in wom.
      implicit none
      integer m,r,nmax,s,posni,posnj,nmr(r),nor(r),
     /     mrmat(r,nmax),ormat(r,nmax),ni,pdwo,pdwm,
     /     i,j,ii,jj,io,im
      double precision w(r*nmax,r*nmax,m),wom(r*nmax,r*nmax)
      do 10 i=1,r*nmax
         do 5 j=i,r*nmax
            wom(i,j)=dble(0.)
 5       continue
 10   continue
      posni=0
       do 50 i=1,r
         do 40 j=1,nor(i)
            posni=posni+1
            posnj=0
            do 30 ii=1,r
               do 20 jj=1,nmr(ii)
                  posnj=posnj+1
                  io=(i-1)*ni+ormat(i,j)
                  im=(ii-1)*ni+mrmat(ii,jj)
                  if(io.le.im) then
                     wom(posni,posnj)=w(io,im,s)
                  else
                     wom(posni,posnj)=w(im,io,s)
                  endif
 20            continue
 30         continue
 40      continue
 50   continue
      pdwo=posni
      pdwm=posnj
      return
      end
C*********************************************************************
      subroutine mkrozpsio(r,q,nmax,ormat,nor,ni,rozpsi,rozpsio)
      implicit none
      integer r,q,nmax,ormat(r,nmax),nor(r),ni,i,ii,j,posni
      double precision rozpsi(r*nmax,r*q),rozpsio(r*nmax,r*q)
      do 10 i=1,r*nmax
         do 5 j=1,r*q
            rozpsio(i,j)=dble(0.)
 5       continue
 10   continue
      posni=0
      do 50 ii=1,r
         do 40 j=1,nor(ii)
            posni=posni+1
            do 30 i=1,r*q
               rozpsio(posni,i)=rozpsi((ii-1)*ni+ormat(ii,j),i)
 30         continue
 40      continue
 50   continue
      return
      end
C*********************************************************************
      subroutine mkrozpsim(r,q,nmax,mrmat,nmr,ni,rozpsi,rozpsim)
      implicit none
      integer r,q,nmax,mrmat(r,nmax),nmr(r),ni,i,ii,j,posni
      double precision rozpsi(r*nmax,r*q),rozpsim(r*nmax,r*q)
      do 10 i=1,r*nmax
         do 5 j=1,r*q
            rozpsim(i,j)=dble(0.)
 5       continue
 10   continue
      posni=0
      do 50 ii=1,r
         do 40 j=1,nmr(ii)
            posni=posni+1
            do 30 i=1,r*q
               rozpsim(posni,i)=rozpsi((ii-1)*ni+mrmat(ii,j),i)
 30         continue
 40      continue
 50   continue
      return
      end
C*********************************************************************
      subroutine mkusotzo(r,m,q,nmax,s,ormat,nor,ni,wkqnm,usotzo)
C takes the observed part of U_i %*% (inv(sigma) Otimes t(Z_i)), storing
C in the upper right side of usotz 
C note that in the scoring part we already calculate 
C U_i %*% (inv(sigma) Otimes t(Z_i))
      implicit none
      integer r,m,q,nmax,s,ormat(r,nmax),nor(r),ni,i,ii,j,posni
      double precision wkqnm(r*q,r*nmax,m),usotzo(r*q,r*nmax)
      do 10 i=1,r*q
         do 5 j=1,r*nmax
            usotzo(i,j)=dble(0.)
 5       continue
 10   continue
      do 60 i=1,r*q
         posni=0
         do 50 ii=1,r
            do 40 j=1,nor(ii)
               posni=posni+1
               usotzo(i,posni)=wkqnm(i,(ii-1)*ni+ormat(ii,j),s)
 40         continue
 50      continue
 60   continue
      return
      end
C*********************************************************************
      subroutine mkusotzm(r,m,q,nmax,s,mrmat,nmr,ni,wkqnm,usotzm)
      implicit none
      integer r,m,q,nmax,s,mrmat(r,nmax),nmr(r),ni,i,ii,j,posni
      double precision wkqnm(r*q,r*nmax,m),usotzm(r*q,r*nmax)
      do 10 i=1,r*q
         do 5 j=1,r*nmax
            usotzm(i,j)=dble(0.)
 5       continue
 10   continue
      do 60 i=1,r*q
         posni=0
         do 50 ii=1,r
            do 40 j=1,nmr(ii)
               posni=posni+1
               usotzm(i,posni)=wkqnm(i,(ii-1)*ni+mrmat(ii,j),s)
 40         continue
 50      continue
 60   continue
      return
      end
C*********************************************************************
      subroutine mkxbw(ntot,m,r,q,p,pcol,nmax,xcol,patt,s,st,
     /     fin,ni,pred,beta,wkqnm,w,wkqb2,wxbw,uszxb)
C calculates (vec(X_i%*%beta)^T %*% W_i)^T= W_i %*% (I_r * X_i)%*%vec(beta)
C and also calculates U_i %*% (inv(sigma) Otimes t(z_i)) %*% vec(X_i%*%beta)
C for individual i
      implicit none
      integer ntot,m,r,q,p,pcol,nmax,xcol(p),
     /     patt(ntot),s,st,fin,ni,i,j,gi,
     /     k,ii,jj,kk,l,indi,indj
      double precision pred(ntot,pcol),beta(p,r),
     /     wkqnm(r*q,r*nmax,m),w(r*nmax,r*nmax,m),
     /     wkqb2(nmax,r),wxbw(r*nmax),uszxb(r*q),sum 
C put X_i%*%(j^th column of beta) into j^th column of wkqb2
      do 5 i=1,r*q
         uszxb(i)=dble(0.)
 5    continue
      gi=0
      do 10 i=1,nmax
         do 7 j=1,r
            gi=gi+1
            wkqb2(i,j)=dble(0.)
            wxbw(gi)=dble(0.)
 7       continue
 10   continue
      gi=0
      do 100 i=st,fin
         if(patt(i).ne.0) then
            gi=gi+1
            do 90 j=1,r
               sum=dble(0.)
               do 80 k=1,p
                  sum=sum+pred(i,xcol(k))*beta(k,j)
 80            continue
               wkqb2(gi,j)=sum
 90         continue
         endif
 100  continue
      do 200 i=1,r
         indi=0
         do 150 j=st,fin
            if(patt(j).ne.0) then
               indi=indi+1
               jj=(i-1)*ni+indi
               sum=dble(0.)
               do 140 k=1,r
                  indj=0
                  do 130 l=st,fin
                     if(patt(l).ne.0) then
                        indj=indj+1
                        kk=(k-1)*ni+indj
                        if(kk.le.jj) then
                           sum=sum+wkqb2(indj,k)*w(kk,jj,s)
                        else
                           sum=sum+wkqb2(indj,k)*w(jj,kk,s)
                        endif
                     endif
 130              continue
 140           continue
               wxbw(jj)=sum
            endif
 150     continue
 200  continue
      do 400 i=1,r
         do 300 j=1,q
            ii=(i-1)*q+j
            sum=dble(0.)
            do 250 k=1,r
               indj=0
               do 220 l=st,fin
                  if(patt(l).ne.0) then
                     indj=indj+1
                     jj=(k-1)*ni+indj
                     sum=sum+wkqnm(ii,jj,s)*wkqb2(indj,k)
                  endif
 220           continue
 250        continue
            uszxb(ii)=sum
 300     continue
 400  continue
      return
      end
C*********************************************************************
      subroutine mkwxbwo(r,nmax,ormat,nor,ni,wxbw,wxbwo)
C takes the part that corresponds to the observed variables in
C W_i %*% (I_r * X_i)%*%vec(beta)
      implicit none
      integer r,nmax,ormat(r,nmax),nor(r),ni,i,ii,j,posni
      double precision wxbw(r*nmax),wxbwo(r*nmax)
      do 10 i=1,r*nmax
         wxbwo(i)=dble(0.)
 10   continue
      posni=0
      do 50 ii=1,r
         do 40 j=1,nor(ii)
            posni=posni+1
            wxbwo(posni)=wxbw((ii-1)*ni+ormat(ii,j))
 40      continue
 50   continue
      return
      end
C*********************************************************************
      subroutine mkwxbwm(r,nmax,ni,mrmat,nmr,wxbw,wxbwm)
C takes the part that corresponds to the missing variables in
C W_i %*% (I_r * X_i)%*%vec(beta)
      implicit none
      integer r,nmax,ni,mrmat(r,nmax),nmr(r),i,ii,j,posni
      double precision wxbw(r*nmax),wxbwm(r*nmax)
      do 10 i=1,r*nmax
         wxbwm(i)=dble(0.)
 10   continue
      posni=0
      do 50 ii=1,r
         do 40 j=1,nmr(ii)
            posni=posni+1
            wxbwm(posni)=wxbw((ii-1)*ni+mrmat(ii,j))
 40      continue
 50   continue
      return
      end
C*********************************************************************
      subroutine mkeb(ntot,m,r,p,pcol,q,nmax,xcol,npatt,patt,
     /     ist,ifin,rmat,nor,ormat,nmr,mrmat,nstari,
     /     pred,beta,wkqnm,u,w,wkqb2,wo,wo1,wm,wom,wkwmm1,wkwmm2,
     /     uszxb,usotzo,usotzm,wxbw,wxbwo,wxbwm,wxbeta,wkeb2,y,
     /     err,msg,vdel,trdet,trdel,eb,varb)
C calculates the expectation of b_i given observed data of i^th subject
C and theta, storing in eb[,s], and the cov matrix of b_i fiven obs data and 
C theta, storing in varb[,,s].
C It also calculates the necessary pieces for evaluating the observed
C log-likelihood :
C -0.5*log(det(cov(y_i(obs)))), and
C trace( (y_i - X_ibeta)^t inv(cov(y_i(obs))) (y_i - X_ibeta) )
      implicit none
      integer ntot,m,r,p,pcol,q,nmax,xcol(p),npatt,
     /     patt(ntot),ist(m),ifin(m),rmat(npatt,r),nor(r),
     /     ormat(r,nmax),nmr(r),mrmat(r,nmax),pdwm,pdwo,nstari(m),
     /     err,msg,s,st,fin,i,j,l,k,posn,ni
      double precision pred(ntot,pcol),beta(p,r),wkqnm(r*q,r*nmax,m),
     /     u(r*q,r*q,m),w(r*nmax,r*nmax,m),wkqb2(nmax,r),
     /     wm(r*nmax,r*nmax),wo(r*nmax,r*nmax),
     /     wo1(r*nmax,r*nmax),wom(r*nmax,r*nmax),
     /     wkwmm1(r*nmax,r*nmax),wkwmm2(r*nmax,r*nmax),
     /     uszxb(r*q),usotzo(r*q,r*nmax),usotzm(r*q,r*nmax),
     /     wxbw(r*nmax),wxbwo(r*nmax),wxbwm(r*nmax),wxbeta(ntot,r),
     /     wkeb2(r*q,r*nmax),y(ntot,r),eb(r*q,m),varb(r*q,r*q,m),
     /     vdel(r*nmax),trdet,trdel,sum
      pdwm=0
      pdwm=0
      trdel=dble(0.)
      trdet=dble(0.)
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
         do 2 i=1,r*q
            do 1 j=1,r*nmax
               wkeb2(i,j)=dble(0.)
 1          continue
 2       continue
         do 4 i=1,r*nmax
            do 3 j=i,r*nmax
               wkwmm1(i,j)=dble(0.)
               wkwmm2(i,j)=dble(0.)
               wo1(i,j)=dble(0.)
 3          continue
 4       continue
         call getormat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nor,ormat)
         call getmrmat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nmr,mrmat)
         call mkxbw(ntot,m,r,q,p,pcol,nmax,xcol,patt,s,st,
     /     fin,ni,pred,beta,wkqnm,w,wkqb2,wxbw,uszxb)
C Now calculate the inverse of W_i(22), whose dim is pdwm, use the workspaces
C wkmm1 and wkmm2 for this inversion, at the end of inversion wkmm2 contains
C the inverse of W_i(22)
         call mkwom(m,r,nmax,s,mrmat,ormat,nmr,nor,
     /        ni,pdwo,pdwm,w,wom)
         call mkwm(m,r,nmax,s,mrmat,nmr,ni,pdwm,w,wm)
         call mkwo(m,r,nmax,s,ormat,nor,ni,pdwo,w,wo)
         do 6 i=1,pdwm
            do 5 j=i,pdwm
               wkwmm1(i,j)=wm(i,j)
 5          continue
 6       continue
         call chfce(r*nmax,pdwm,wkwmm1,err)
         if(err.eq.1) then 
            goto 999
         endif
         call bkslv(r*nmax,pdwm,wkwmm1)
         call mm(r*nmax,pdwm,wkwmm1,wkwmm2)
         do 20 k=1,pdwo
            do 19 l=k,pdwo
              sum=dble(0.)
              do 18 j=1,pdwm
                 do 17 i=1,j
                    sum=sum+wom(k,i)*wkwmm2(i,j)*wom(l,j)
 17              continue
                 do 16 i=j+1,pdwm
                    sum=sum+wom(k,i)*wkwmm2(j,i)*wom(l,j)
 16              continue
 18           continue
              wo1(k,l)=wo(k,l)-sum
 19        continue
 20     continue
        call trdelwdel(nmax,r,ntot,st,fin,npatt,patt,rmat,
     /     p,xcol,pcol,pdwo,pred,beta,y,wxbeta,vdel,wo1,trdel)
        call chfce(r*nmax,pdwo,wo1,err)
        if(err.eq.1) then
           msg=90
           goto 999
        endif
        do 21 i=1,pdwo
           trdet=trdet+dlog(wo1(i,i))
 21     continue
         do 25 i=1,r*q
            eb(i,s)=-uszxb(i)
            do 24 j=i,r*q
               varb(i,j,s)=u(i,j,s)
               if(i.ne.j) varb(j,i,s)=varb(i,j,s)
 24         continue
 25      continue
         call mkusotzo(r,m,q,nmax,s,ormat,nor,ni,wkqnm,usotzo)
         call mkusotzm(r,m,q,nmax,s,mrmat,nmr,ni,wkqnm,usotzm)
         call mkwxbwo(r,nmax,ormat,nor,ni,wxbw,wxbwo)
         call mkwxbwm(r,nmax,ni,mrmat,nmr,wxbw,wxbwm)
C calculate the first part of eb[,s] first: 
         do 50 l=1,r*q
            sum=dble(0.)
            do 40 j=1,pdwm
               do 35 i=1,j
                  sum=sum+wxbwm(j)*usotzm(l,i)*wkwmm2(i,j)
 35            continue
               do 30 i=j+1,pdwm
                  sum=sum+wxbwm(j)*usotzm(l,i)*wkwmm2(j,i)
 30           continue
 40        continue
           eb(l,s)=eb(l,s)+sum
 50     continue
C now add the second part of eb[,s]
        do 100 k=1,r*q
           do 90 l=1,pdwo
              sum=dble(0.)
              do 80 j=1,pdwm
                 do 70 i=1,j
                    sum=sum+usotzm(k,i)*wkwmm2(i,j)*wom(l,j)
 70              continue
                 do 75 i=j+1,pdwm
                    sum=sum+usotzm(k,i)*wkwmm2(j,i)*wom(l,j)
 75              continue
 80           continue
C use workspace wkeb2(r*q,r*nmax) to store usoztm%*%wkwmm2%*%t(wom)
C note we are only using rq*(#obs(i)) part of wkeb2
              wkeb2(k,l)=usotzo(k,l)-sum
 90        continue
 100    continue
C now multiply wkeb2(rq,#obs(i)) by v(y_i(obs))              
        do 200 i=1,r*q
           posn=0
           sum=dble(0.)
           do 150 j=1,r
              do 140 l=st,fin
                 if(patt(l).ne.0) then
                    if(rmat(patt(l),j).eq.1) then
                       posn=posn+1
                       sum=sum+wkeb2(i,posn)*y(l,j)
                    endif
                 endif
 140          continue
 150       continue
           eb(i,s)=eb(i,s)+sum
 200    continue
C now calculate the varb[,s]=cov of b_i given (y_i(obs),theta),
C which is already already equal to u_i
        do 400 k=1,r*q
           do 350 l=1,r*q
              sum=dble(0.)
              do 340 j=1,pdwm
                 do 330 i=1,j
                    sum=sum+usotzm(k,i)*wkwmm2(i,j)*usotzm(l,j)
 330             continue
                 do 335 i=j+1,pdwm
                    sum=sum+usotzm(k,i)*wkwmm2(j,i)*usotzm(l,j)
 335             continue
 340          continue
              varb(k,l,s)=varb(k,l,s)+sum
 350       continue
 400    continue
 500  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkey(loc,lmc,oc,mc,m,r,ntot,nstari,iposn,npatt,pstfin,
     /     rmat,patt,p,q,xcol,zcol,pcol,ist,ifin,pred,
     /     y,beta,eb,sigma,wkrr1,wkrrpt,wkrrb21,eystar,ey)
C computes E( vec(ystar_i) | y_i(obs),theta ), and E(y_i/y_i(obs),theta) where
C ystar=y_i-(X_i%*%beta + Z_i%*%E(b_i/y_i(obs),theta)), and also computes
C wkrrpt containing the matrices of residual covariances from the regression
C of y_i(mis) on y_i(obs) in pattern "pt"
C Work with only partially and/or complete cases.
C Finally at the end it transforms back to y(which is here ey)
      implicit none
      integer loc,lmc,oc(loc),mc(lmc),m,r,ntot,nstari(m),iposn(ntot),
     /     npatt,pstfin(npatt,2),rmat(npatt,r),patt(ntot),
     /     p,q,xcol(p),zcol(q),pcol,ist(m),ifin(m),
     /     i,j,pt,nmc,noc,k,l,s,st,fin,ni,gi
      double precision pred(ntot,pcol),y(ntot,r),beta(p,r),
     /     eb(r*q,m),sigma(r,r),wkrr1(r,r),wkrrpt(r,r,npatt),
     /     wkrrb21(r,r,npatt),eystar(ntot,r),ey(ntot,r),sum
C first transform the observed variables i.e calculate y*_ijk 
      do 60 s=1,m
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
         do 50 i=st,fin
            gi=0
            if(patt(i).ne.0) then
               do 40 j=1,r
                  sum=dble(0.)
                  do 35 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 35               continue
                  do 30 k=1,q
                     gi=gi+1
                     sum=sum+pred(i,zcol(k))*eb(gi,s)
 30               continue
                  if(rmat(patt(i),j).eq.1) then
                     eystar(i,j)=y(i,j)-sum
                  else
                     eystar(i,j)=y(i,j)
                  endif
 40            continue
            endif
 50      continue
 60   continue
C now start filling in the missing portions of y*_i at the same time
C store residual covariance matrix in pt'th layer of wkrrpt
      do 70 i=1,r
         do 75 j=i,r
            wkrr1(i,j)=sigma(i,j)
 75      continue
 70    continue
      do 110 pt=1,npatt
         call swpobs(r,wkrr1,npatt,rmat,pt)
         do 79 i=1,r
            do 78 j=i+1,r
               wkrr1(j,i)=wkrr1(i,j)
 78         continue
 79      continue
         call getmc(r,npatt,rmat,pt,lmc,mc,nmc)
         call getoc(r,npatt,rmat,pt,loc,oc,noc)
         do 85 i=1,nmc
            do 80 j=1,nmc
               wkrrpt(mc(i),mc(j),pt)=wkrr1(mc(i),mc(j))
 80         continue
            do 83 j=1,noc
               wkrrb21(mc(i),oc(j),pt)=wkrr1(mc(i),oc(j))
 83         continue
 85      continue
         do 100 i=pstfin(pt,1),pstfin(pt,2)
            do 95 k=1,nmc
               sum=dble(0.)
               do 90 l=1,noc
                  sum=sum+wkrr1(oc(l),mc(k))*eystar(iposn(i),oc(l))
 90            continue
               eystar(iposn(i),mc(k))=sum
 95         continue
 100     continue
 110  continue
C now calculate E(y_i/y_i(obs),theta) 
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
C         gi=0
         do 450 i=st,fin
            gi=0
            if(patt(i).ne.0) then
               do 400 j=1,r
                  sum=dble(0.)
                  do 350 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 350              continue
                  do 300 k=1,q
                     gi=gi+1
                     sum=sum+pred(i,zcol(k))*eb(gi,s)
 300              continue
                  ey(i,j)=sum+eystar(i,j)
 400           continue
            endif
 450     continue
 500  continue
      return
      end
C*********************************************************************
      subroutine mkyyt(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,
     /     y,eyyt)
C calculates vec(y)%*%t(vec(y)) for only observed positions,
C for subject s
      integer ntot,nmax,r,st,fin,ni,patt(ntot),npatt,rmat(npatt,r),
     /     i,indi,j,ri,k,indj,l,rj
      double precision y(ntot,r),eyyt(r*nmax,r*nmax)
      do 10 i=1,r*nmax
         do 5 j=1,r*nmax
            eyyt(i,j)=dble(0.)
 5       continue
 10   continue
      do 400 i=1,r
         indi=0
         do 350 j=st,fin
            if(patt(j).ne.0) then
               indi=indi+1
               ri=(i-1)*ni+indi
               do 300 k=1,r
                  indj=0
                  do 250 l=st,fin
                     if(patt(l).ne.0) then
                        indj=indj+1
                        rj=(k-1)*ni+indj
                        if(rmat(patt(l),k).eq.1) then
                           if(rmat(patt(j),i).eq.1) then
                              eyyt(ri,rj)=y(j,i)*y(l,k)
                           endif
                        endif
                     endif
 250              continue
 300           continue
            endif
 350     continue
 400  continue
      return
      end
C*********************************************************************
      subroutine mkxbeta(ntot,m,ist,ifin,p,r,pcol,xcol,patt,pred,beta,
     /     wxbeta)
C calculates X_i%*%beta required for the calculations in the
C conventional EM
      integer ntot,m,ist(m),ifin(m),p,r,pcol,xcol(p),
     /     patt(ntot),s,st,fin,i,j,gi,k
      double precision pred(ntot,pcol),beta(p,r),wxbeta(ntot,r),sum
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 50 i=st,fin
            gi=0
            if(patt(i).ne.0) then
               do 45 j=1,r
                  sum=dble(0.)
                  do 30 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 30               continue
                  wxbeta(i,j)=sum
 45            continue
            endif
 50      continue
 100  continue
      return
      end
C*********************************************************************
      subroutine preyxyxt(ntot,m,ist,ifin,p,q,r,pcol,xcol,zcol,
     /     patt,pred,beta,eb,wxbeta,wxbetazeb)
C calculates X_i %*% beta and X_i%*%beta + Z_i%*%E(b_i | y_obs, theta)
      integer ntot,m,ist(m),ifin(m),p,q,r,pcol,xcol(p),zcol(q),
     /     patt(ntot),s,st,fin,i,j,gi,k
      double precision pred(ntot,pcol),beta(p,r),eb(r*q,m),
     /     wxbeta(ntot,r),wxbetazeb(ntot,r),sum
      do 100 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 50 i=st,fin
            gi=0
            if(patt(i).ne.0) then
               do 45 j=1,r
                  sum=dble(0.)
                  do 30 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 30               continue
                  wxbeta(i,j)=sum
                  do 35 k=1,q
                     gi=gi+1
                     sum=sum+pred(i,zcol(k))*eb(gi,s)
 35               continue
                  wxbetazeb(i,j)=sum
 45            continue
            endif
 50      continue
 100  continue
      return
      end
C*********************************************************************
      subroutine mkeyyt(ntot,nmax,npatt,m,r,st,fin,ni,patt,
     /     rmat,s,pcol,q,zcol,lmc,loc,mc,mc1,oc,oc2,
     /     pred,varb,wkrrpt,wkrrb21,ey,eyyt)
C this subroutine fills only the missing portions of the 
C E(vec(y_i) %*%vec(y_i)^T  | y_i(obs),theta )
      implicit none
      integer ntot,nmax,npatt,m,r,st,fin,ni,patt(ntot),
     /     rmat(npatt,r),s,pcol,q,zcol(q),lmc,loc,
     /     mc(lmc),mc1(lmc),oc(loc),nmc,nmc1,noc,i,j,ri,rj,k,
     /     jj,ii,pt,pt1,pt2,gi,gj,noc2,oc2(loc)
      double precision pred(ntot,pcol),varb(r*q,r*q,m),
     /     wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),ey(ntot,r),
     /     eyyt(r*nmax,r*nmax),var1,var2,var3,var4,var22,var32,tmp
      tmp=dble(0.)
      gi=0
      do 40 i=st,fin
         if(patt(i).ne.0) then
            gi=gi+1
            pt=patt(i)
            call getmc(r,npatt,rmat,pt,lmc,mc,nmc)
            call getoc(r,npatt,rmat,pt,loc,oc,noc)
C     both missing in the same row 
            do 30 j=1,nmc
               do 20 k=1,nmc
                  ri=(mc(j)-1)*ni+gi
                  rj=(mc(k)-1)*ni+gi
                  eyyt(ri,rj)=wkrrpt(mc(j),mc(k),pt)+
     /                 ey(i,mc(j))*ey(i,mc(k))+
     /                 var1(ntot,r,q,m,s,i,i,mc(j),mc(k),
     /                 zcol,pcol,pred,varb)-
     /                 var2(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mc(j),
     /                 mc(k),i,zcol,pcol,wkrrb21,varb,pred)-
     /                 var3(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mc(j),
     /                 mc(k),i,zcol,pcol,wkrrb21,varb,pred)+
     /                 var4(ntot,r,q,m,s,pt,pt,noc,noc,loc,oc,oc,
     /                 npatt,mc(j),mc(k),i,i,zcol,pcol,wkrrb21,varb,
     /                 pred)
 20            continue
 30         continue
         endif
 40   continue
C  both missing in the same columns but maybe different rows
      gi=0
      do 70 i=st,fin
         if(patt(i).ne.0) then
            pt1=patt(i)
            gi=gi+1
            call getmc(r,npatt,rmat,patt(i),lmc,mc1,nmc)
            call getoc(r,npatt,rmat,patt(i),loc,oc,noc)
            do 65 ii=1,nmc
               ri=(mc1(ii)-1)*ni+gi
               gj=0
               do 60 j=st,fin
                  if(patt(j).ne.0) then
                     pt2=patt(j)
                     gj=gj+1
                     call getmc(r,npatt,rmat,patt(j),lmc,mc,nmc1)
                     call getoc(r,npatt,rmat,patt(j),loc,oc2,noc2)
                     do 57 k=1,nmc1
                        rj=(mc(k)-1)*ni+gj
                        if(i.ne.j) then
                           eyyt(ri,rj)=ey(i,mc1(ii))*ey(j,mc(k))+
     /                          var1(ntot,r,q,m,s,i,j,mc1(ii),mc(k),
     /                          zcol,pcol,pred,varb)-
     /                          var22(ntot,r,q,m,pt2,s,noc2,oc2,loc,
     /                          npatt,mc1(ii),mc(k),i,j,zcol,pcol,
     /                          wkrrb21,varb,pred)-
     /                          var32(ntot,r,q,m,pt1,s,noc,oc,loc,
     /                          npatt,mc1(ii),mc(k),i,j,zcol,pcol,
     /                          wkrrb21,varb,pred)+
     /                          var4(ntot,r,q,m,s,pt1,pt2,
     /                          noc,noc2,loc,oc,oc2,
     /                          npatt,mc1(ii),mc(k),i,j,zcol,pcol,
     /                          wkrrb21,varb,pred)
                        endif
 57                  continue
                  endif
 60            continue
 65         continue
         endif
 70   continue
C     take care of the observations (one missing one not)
      gi=0
      do 150 i=st,fin
         if(patt(i).ne.0) then
            gi=gi+1
            pt=patt(i)
            call getoc(r,npatt,rmat,pt,loc,oc,noc)
            gj=0
            do 140 j=st,fin
               if(patt(j).ne.0) then
                  gj=gj+1
                  call getmc(r,npatt,rmat,patt(j),lmc,mc,nmc)
                  do 130 ii=1,noc
                     ri=(oc(ii)-1)*ni+gi
                     do 120 jj=1,nmc
                        rj=(mc(jj)-1)*ni+gj
                        eyyt(ri,rj)=ey(i,oc(ii))*ey(j,mc(jj))
 120                 continue
 130              continue
               endif
 140        continue
            call getmc(r,npatt,rmat,pt,lmc,mc,nmc)
            gj=0
            do 145 j=st,fin
               if(patt(j).ne.0) then
                  gj=gj+1
                  call getoc(r,npatt,rmat,patt(j),loc,oc,noc)
                  do 144 ii=1,nmc
                     ri=(mc(ii)-1)*ni+gi
                     do 143 jj=1,noc
                        rj=(oc(jj)-1)*ni+gj
                        eyyt(ri,rj)=ey(i,mc(ii))*ey(j,oc(jj))
 143                 continue
 144              continue
               endif
 145        continue
         endif
 150  continue
      return
      end
C**********************************************************************
      subroutine mkeyxyxt(ntot,nmax,r,st,fin,ni,patt,wxbeta,ey,eyyt,
     /     eyxyxt)
C  calculates 
C  E( vec(y_i - X_i beta)%*%t(vec(y_i - X_i beta)) | y_i(obs),theta )
      implicit none
      integer ntot,nmax,r,st,fin,ni,patt(ntot),i,j,k,gi,gj,di,dj,l
      double precision wxbeta(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax)
C   note that we already calculated E(y_i/y_i(obs)) in the previous routine
      do 10 i=1,r*nmax
         do 5 k=1,r*nmax
            eyxyxt(i,k)=dble(0.)
 5       continue
 10   continue
      di=0
      dj=0
      do 250 j=1,r
         gi=0
         do 200 i=st,fin
            if(patt(i).ne.0) then
               gi=gi+1
               di=(j-1)*ni+gi
               do 190 l=1,r
                  gj=0
                  do 180 k=st,fin
                     if(patt(k).ne.0) then
                        gj=gj+1
                        dj=(l-1)*ni+gj
                        eyxyxt(di,dj)=eyyt(di,dj)-ey(i,j)*wxbeta(k,l)-
     /                       ey(k,l)*wxbeta(i,j)+
     /                       wxbeta(i,j)*wxbeta(k,l)
                     endif
 180              continue
 190           continue
            endif
 200     continue
 250  continue
      return
      end
C***********************************************************************
      subroutine chsub(r,sigma,lmc,mc,nmc,wkrr)
C takes cholesky factor of submatrix of sigma corresponding to mc,
C storing it in the upper-left (nmc x nmc) corner of wkrr
      integer r,lmc,mc(lmc),nmc,i,j
      double precision sigma(r,r),wkrr(r,r)
      do 10 i=1,nmc
         do 5 j=i,nmc
            wkrr(i,j)=sigma(mc(i),mc(j))
 5       continue
 10   continue
      call chfc(r,nmc,wkrr)
      return
      end
C***********************************************************************
      subroutine getmc(r,npatt,rmat,pt,lmc,mc,nmc)
      integer r,npatt,rmat(npatt,r),pt,lmc,mc(lmc),nmc,posn,j
      nmc=0
      posn=0
      do 2 j=1,r
         mc(j)=0
 2    continue
      do 10 j=1,r
         if(rmat(pt,j).eq.0) then
            nmc=nmc+1
            posn=posn+1
            mc(posn)=j
         endif
 10   continue
      return
      end
C***********************************************************************
      subroutine getoc(r,npatt,rmat,pt,loc,oc,noc)
      integer r,npatt,rmat(npatt,r),pt,loc,oc(loc),noc,posn,j
      noc=0
      posn=0
      do 2 j=1,r
         oc(j)=0
 2    continue
      do 10 j=1,r
         if(rmat(pt,j).eq.1) then
            noc=noc+1
            posn=posn+1
            oc(posn)=j
         endif
 10   continue
      return
      end
C***********************************************************************
      subroutine swpobs(r,sigma,npatt,rmat,pt)
C sweeps sigma on positions of rmat(pt,*) containing 1's
      integer r,npatt,rmat(npatt,r),pt,col
      double precision sigma(r,r)
      do 100 col=1,r
         if((rmat(pt,col).eq.1).and.(sigma(col,col).gt.dble(0.))) 
     /        call swp(r,sigma,col)
         if((rmat(pt,col).eq.0).and.(sigma(col,col).le.dble(0.)))
     /        call rsw(r,sigma,col)
 100  continue
      return
      end
C***********************************************************************
      subroutine swp(n,mat,p)
C sweep upper-triangular portion of mat on position p
      integer n,p,i,j
      double precision mat(n,n)
      mat(p,p)=-dble(1.)/mat(p,p)
      do 1 i=1,p-1
         mat(i,p)=-mat(i,p)*mat(p,p)
 1    continue
      do 2 j=p+1,n
         mat(p,j)=-mat(p,j)*mat(p,p)
 2    continue
      do 6 i=1,p-1
         do 4 j=i,p-1
            mat(i,j)=mat(i,j)+mat(i,p)*mat(j,p)/mat(p,p)
 4       continue
         do 5 j=p+1,n
            mat(i,j)=mat(i,j)+mat(i,p)*mat(p,j)/mat(p,p)
 5       continue
 6    continue
      do 8 i=p+1,n
         do 7 j=i,n
            mat(i,j)=mat(i,j)+mat(p,i)*mat(p,j)/mat(p,p)
 7       continue
 8    continue
      return
      end
C***********************************************************************
      subroutine rsw(n,mat,p)
C reverse-sweep upper-triangular portion of mat on position p
      integer n,p,i,j
      double precision mat(n,n)
      mat(p,p)=-dble(1.)/mat(p,p)
      do 1 i=1,p-1
         mat(i,p)=mat(i,p)*mat(p,p)
 1    continue
      do 2 j=p+1,n
         mat(p,j)=mat(p,j)*mat(p,p)
 2    continue
      do 6 i=1,p-1
         do 4 j=i,p-1
            mat(i,j)=mat(i,j)+mat(i,p)*mat(j,p)/mat(p,p)
 4       continue
         do 5 j=p+1,n
            mat(i,j)=mat(i,j)+mat(i,p)*mat(p,j)/mat(p,p)
 5       continue
 6    continue
      do 8 i=p+1,n
         do 7 j=i,n
            mat(i,j)=mat(i,j)+mat(p,i)*mat(p,j)/mat(p,p)
 7       continue
 8    continue
      return
      end
C*********************************************************************
C ********the following subroutines are for the M-step **********
C*********************************************************************
      subroutine mku(r,q,m,psi,sigma,ztz,u,wkrr1,wkrr2,wkqq1,
     /     wkqq2,wkqqu,ldpsi,ldsig,ldu,err)
C calculates U_i, i=1,2,...,m from Sigma and Psi.
C after execution, wkqq1 and wkrr1 contain the inverse of the 
C original psi and sigma, respectively(upper tri only).
C Also calculates ldsig = -0.5*logdet(sigma)
C                 ldpsi = -0.5*logdet(psi) and  
C                 ldu = 0.5*sum(logdet(U_i))
C sets err=1 if the input value of psi and/or sigma is not positive
C definite, or if the value of U_i is not positive definite
      implicit none
      integer r,q,m,err,s,i,j,ii,jj,ia,ja
      double precision psi(r*q,r*q),sigma(r,r),ztz(q,q,m),u(r*q,r*q,m),
     /     wkrr1(r,r),wkrr2(r,r),wkqq1(r*q,r*q),wkqq2(r*q,r*q),
     /     wkqqu(r*q,r*q,m),ldpsi,ldsig,ldu   
      err=0
C invert psi and sigma and put them into wkqq1 and wkrr1
      do 2 i=1,r
         do 1 j=i,r
            wkrr2(i,j)=sigma(i,j)
 1       continue
 2    continue
      do 4 i=1,r*q
         do 3 j=i,r*q
            wkqq2(i,j)=psi(i,j)
 3       continue
 4    continue
      call chfce(r*q,r*q,wkqq2,err)
      if(err.eq.1) goto 999
      call bkslv(r*q,r*q,wkqq2)
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) goto 999
      call bkslv(r,r,wkrr2)
      ldsig=dble(0.)
      ldpsi=dble(0.)
      do 5 i=1,r*q
         ldpsi=ldpsi+dlog(wkqq2(i,i))
 5    continue
      do 10 i=1,r
         ldsig=ldsig+dlog(wkrr2(i,i))
 10   continue
      call mm(r*q,r*q,wkqq2,wkqq1)
      call mm(r,r,wkrr2,wkrr1)
      ldu=dble(0.)
      do 12 i=1,r
         do 11 j=i+1,r
            wkrr1(j,i)=wkrr1(i,j)
 11      continue
 12   continue
      do 500 s=1,m
C initialize inv(u(,,s))=wkqqu(,,s) to inv(sigma)  Otimes t(z_i)%*%z_i
         do 100 i=1,r
            do 50 j=i,r
               do 30 ii=1,q
                  do 15 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     wkqqu(ia,ja,s)=wkrr1(i,j)*ztz(ii,jj,s)
 15               continue
 30            continue
 50         continue
 100     continue
C add inv(psi) to wkqqu(,,s) and take the inverse to get u(,,s)
         do 200 i=1,r*q
            do 150 j=i,r*q
               u(i,j,s)=wkqqu(i,j,s)+wkqq1(i,j)
 150        continue
 200     continue
C at this point wkqqu(,,s) contains inv(u)
         call chle(r*q,r*q,m,u,s,err)
         call bkslvl(r*q,r*q,m,u,s)
         do 300 i=1,r*q
            ldu=ldu+dlog(u(i,i,s))
 300     continue
         call mmul(r*q,r*q,m,u,s,wkqq2)
         do 400 i=1,r*q
            do 350 j=i,r*q
               u(i,j,s)=wkqq2(i,j)
 350        continue
 400     continue
 500  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkubd(r,q,m,psi,sigma,ztz,u,wkrr1,wkrr2,wkrqrq1,
     /     wkrqrq2,wkqq1bd,wkqq2bd,wkqqu,ldpsi,ldsig,ldu,err)
C Version of mku for the block-diagonal mglmmembd
C invert psi and put into wkrqrq2, using the fact that it's block-diagonal
      implicit none
      integer r,q,m,err,s,i,j,l,ii,jj,ia,ja
      double precision psi(q,q,r),sigma(r,r),ztz(q,q,m),u(r*q,r*q,m),
     /     wkrr1(r,r),wkrr2(r,r),wkrqrq1(r*q,r*q),wkrqrq2(r*q,r*q),
     /     wkqq1bd(q,q,r),wkqq2bd(q,q),wkqqu(r*q,r*q,m),ldpsi,ldsig,ldu   
      err=0
      ldsig=dble(0.)
      ldpsi=dble(0.)
      do 102 i=1,r*q
         do 101 j=i,r*q
            wkrqrq2(i,j)=dble(0.)
            wkrqrq1(i,j)=dble(0.)
 101     continue
 102  continue
      do 110 l=1,r
         do 105 i=1,q
            do 104 j=i,q
               wkqq1bd(i,j,l)=psi(i,j,l)
 104        continue
 105     continue
         call chle(q,q,r,wkqq1bd,l,err)
         if(err.eq.1) goto 999
         call bkslvl(q,q,r,wkqq1bd,l)
         do 106 i=1,q
            ldpsi=ldpsi+dlog(wkqq1bd(i,i,l))
 106     continue
         call mmul(q,q,r,wkqq1bd,l,wkqq2bd)
         do 108 i=1,q
            do 107 j=i,q
               ii=(l-1)*q + i
               jj=(l-1)*q + j
               wkrqrq2(ii,jj)=wkqq2bd(i,j)
 107        continue
 108     continue
 110  continue
C invert sigma and put into wkrr1
      do 7 i=1,r
         do 6 j=i,r
            wkrr2(i,j)=sigma(i,j)
 6       continue
 7    continue
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) goto 999
      call bkslv(r,r,wkrr2)
      do 8 i=1,r
         ldsig=ldsig+dlog(wkrr2(i,i))
 8    continue
      call mm(r,r,wkrr2,wkrr1)
      do 10 i=1,r
         do 9 j=i+1,r
            wkrr1(j,i)=wkrr1(i,j)
 9       continue
 10   continue
      ldu=dble(0.)
      do 100 s=1,m
C     initialize inv(u(,,s))=wkqqu(,,s) to inv(sigma)  Otimes t(z_i)%*%z_i
         do 30 i=1,r
            do 20 j=i,r
               do 15 ii=1,q
                  do 14 jj=1,q
                     ia=(i-1)*q+ii
                     ja=(j-1)*q+jj
                     wkqqu(ia,ja,s)=wkrr1(i,j)*ztz(ii,jj,s)
 14               continue
 15            continue
 20         continue
 30      continue
C     add inv(psi) to wkqqu(,,s) and take the inverse to get u(,,s)
         do 40 i=1,r*q
            do 35 j=i,r*q
               u(i,j,s)=wkqqu(i,j,s)+wkrqrq2(i,j)
 35         continue
 40      continue
         call chle(r*q,r*q,m,u,s,err)
         call bkslvl(r*q,r*q,m,u,s)
         do 45 i=1,r*q
            ldu=ldu+dlog(u(i,i,s))
 45      continue
         call mmul(r*q,r*q,m,u,s,wkrqrq1)
         do 60 i=1,r*q
            do 50 j=i,r*q
               u(i,j,s)=wkrqrq1(i,j)
 50         continue
 60      continue
 100  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkwkqnm(m,r,q,nmax,ntot,ist,ifin,pcol,zcol,patt,
     /     nstari,pred,sigmainv,u,wksigtz,wkqnm)
C given U_i, i=1,2,...,m, calculates U_i%*%(inv(sigma) Otimes t(Z_i))
C use wksigtz to store inv(sigma) Otimes t(z_i) for i=1,2,...,m
C note that this subroutine is called after mku, which will have the 
C inverse cholesky factor of the original sigma in wkrr2 and wkqq2 will 
C contain the inverse of the original sigma. 
C In the main subroutine mlmem ( and mlmembd ), sigmainv will be
C replaced by wkrr1(inv(sigma)).
      implicit none
      integer m,r,q,nmax,ntot,ist(m),ifin(m),pcol,
     /     zcol(q),patt(ntot),nstari(m),st,fin,s,i,j,k,l,ll,kk,ind
      double precision pred(ntot,pcol),sigmainv(r,r),u(r*q,r*q,m),
     /     wksigtz(r*q,r*nmax,m),wkqnm(r*q,r*nmax,m),sum     
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
C put inv(sigma) Otimes t(z_i) into wksigtz(,,s)
         do 200 i=1,r
            do 100 j=1,r
               do 50 k=1,q
                  ind=0
                  do 30 l=st,fin
                     if(patt(l).ne.0) then
                        kk=(i-1)*q+k
                        ind=ind+1
                        ll=(j-1)*nstari(s)+ind
                        if(i.le.j) then
                           wksigtz(kk,ll,s)=sigmainv(i,j)*
     /                          pred(l,zcol(k))
                        else
                           wksigtz(kk,ll,s)=sigmainv(j,i)*
     /                          pred(l,zcol(k))
                        endif
                     endif
 30               continue
 50            continue
 100        continue
 200     continue
         do 400 i=1,r*nstari(s)
            do 350 j=1,r*q
               sum=dble(0.)
               do 310 k=1,j-1
                  sum=sum+u(k,j,s)*wksigtz(k,i,s)
 310           continue
               do 330 k=j,r*q
                  sum=sum+u(j,k,s)*wksigtz(k,i,s)
 330           continue
               wkqnm(j,i,s)=sum
 350        continue
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkw(m,r,q,nmax,ntot,ist,ifin,
     /     nstari,patt,sigmainv,wksigtz,wkqnm,w)
C given inv(sigma), wksigtz and wkqnm, calculates (upper-tri part of)
C weight matrices W_i
      implicit none
      integer ntot,m,r,q,nmax,ist(m),ifin(m),nstari(m),patt(ntot),
     /     s,i,j,k,l,kk,st,fin,ni,indi,indj,ll
      double precision sigmainv(r,r),wksigtz(r*q,r*nmax,m),
     /    wkqnm(r*q,r*nmax,m),w(r*nmax,r*nmax,m),sum
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
C put (inv(sigma) Otimes z_i)%*%U_i(inv(sigma) Otimes t(z_i)) into w_i
         do 100 i=1,r*ni
            do 50 j=i,r*ni
               sum=dble(0.)
               do 10 k=1,r*q
                  sum=sum+wksigtz(k,i,s)*wkqnm(k,j,s)
 10            continue
               w(i,j,s)=sum
 50         continue
 100     continue
         do 300 i=1,r
            indi=0
            do 200 k=st,fin
               if(patt(k).ne.0) then
                  indi=indi+1
                  kk=(i-1)*ni+indi
                  do 190 j=i,r
                     indj=0
                     do 180 l=st,fin
                        if(patt(l).ne.0) then
                           indj=indj+1
                           ll=(j-1)*ni+indj
                           if(abs(kk-ll).eq.abs((j-i)*ni)) then
                              w(kk,ll,s)=sigmainv(i,j)-w(kk,ll,s)
                           else
                              w(kk,ll,s)=-w(kk,ll,s) 
                           endif
                        endif
 180                 continue
 190              continue
               endif
 200        continue
 300     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine gls(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /     nstari,patt,pred,w,ey,beta,xtw,xtwx,xtwey,
     /     xtwxinv,err)
C calculates gls estimate of beta, using weights in w
      implicit none
      integer ntot,m,r,ist(m),ifin(m),nmax,pcol,p,
     /     xcol(p),nstari(m),patt(ntot),ni,st,fin,s,
     /     i,j,err,ia,ja
      double precision pred(ntot,pcol),w(r*nmax,r*nmax,m),
     /     ey(ntot,r),beta(p,r),xtw(p*r,nmax*r),xtwx(p*r,p*r),
     /     xtwey(p*r),xtwxinv(p*r,p*r),sum
      err=0
C initialize t(I_r Otimes X_i)%*%W_i%*%vec(y_i) and
C t(I_r Otimes X_i)%*%W_i%*%(I_r Otimes X_i)
      do 10 i=1,r*p
         xtwey(i)=dble(0.)
         do 5 j=i,r*p
            xtwx(i,j)=dble(0.)
 5       continue
 10   continue
      do 100 s=1,m
         ni=nstari(s)
         st=ist(s)
         fin=ifin(s)
         call mkxtw(ntot,r,p,m,pcol,xcol,patt,ni,st,fin,nmax,w,
     /        pred,xtw,s)
         call mkxtwx(ntot,r,p,pcol,xcol,st,fin,patt,
     /        ni,nmax,pred,xtw,xtwx)
         call mkxtwey(ntot,r,p,st,fin,nmax,ni,patt,
     /        xtw,ey,xtwey)
 100  continue
      call chfce(p*r,p*r,xtwx,err)
      if(err.eq.1) goto 999
      call bkslv(p*r,p*r,xtwx)
      call mm(p*r,p*r,xtwx,xtwxinv)
      ia=0
      do 200 i=1,r
         do 150 j=1,p
            ia=ia+1
            sum=dble(0.)
            do 110 ja=1,ia
               sum=sum+xtwxinv(ja,ia)*xtwey(ja)
 110        continue
            do 120 ja=ia+1,p*r
               sum=sum+xtwxinv(ia,ja)*xtwey(ja)
 120        continue
            beta(j,i)=sum
 150     continue
 200  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkxtw(ntot,r,p,m,pcol,xcol,patt,ni,st,fin,nmax,w,
     /     pred,xtw,s)
C calculates t(I_r Otimes X)%*%W for subject s 
      implicit none
      integer ntot,r,p,m,pcol,xcol(p),patt(ntot),ni,st,fin,
     /     nmax,s,i,j,k,ia,ja,la,ll,l,indi,indl
      double precision pred(ntot,pcol),w(r*nmax,r*nmax,m),
     /     xtw(p*r,nmax*r),sum
      do 500 i=1,r
         do 450 j=1,r
            do 400 k=1,p
               ia=(i-1)*p+k           
               indi=0
               do 300 l=st,fin
                  if(patt(l).ne.0) then
                     indi=indi+1
                     ja=(j-1)*ni+indi
                     sum=dble(0.)
                     indl=0
                     do 250 ll=st,fin
                        if(patt(ll).ne.0) then
                           indl=indl+1
                           la=(i-1)*ni+indl
                           if(la.le.ja) then
                              sum=sum+w(la,ja,s)*pred(ll,xcol(k))
                           else
                              sum=sum+w(ja,la,s)*pred(ll,xcol(k))
                           endif
                        endif
 250                 continue
                  endif
                  xtw(ia,ja)=sum
 300           continue
 400        continue
 450     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkxtwx(ntot,r,p,pcol,xcol,st,fin,patt,
     /     ni,nmax,pred,xtw,xtwx)
C increments t(I_r Otimes X)%*%W%*%(I_r Otimes X) for subject s
      implicit none
      integer ntot,r,p,pcol,xcol(p),st,fin,patt(ntot),ni,
     /     nmax,i,j,k,l,ii,jj,ll,lla,indi
      double precision pred(ntot,pcol),xtw(p*r,nmax*r),
     /     sum,xtwx(p*r,p*r)
      do 500 i=1,r
         do 400 j=1,r
            do 200 k=1,p
               do 100 l=1,p
                  ii=(i-1)*p+k
                  jj=(j-1)*p+l
                  sum=dble(0.)
                  indi=0
                  do 50 ll=st,fin
                     if(patt(ll).ne.0) then
                        indi=indi+1
                        lla=(i-1)*ni+indi
                        sum=sum+xtw(jj,lla)*pred(ll,xcol(k))
                     endif
 50               continue
                  if(ii.le.jj) then
                     xtwx(ii,jj)=xtwx(ii,jj)+sum
                  endif
 100           continue
 200        continue
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkxtwey(ntot,r,p,st,fin,nmax,ni,
     /     patt,xtw,ey,xtwey)
C increments t(I_r Otimes X)%*%W%*%E( vec(y) | for subject s
      implicit none
      integer ntot,r,p,st,fin,nmax,ni,patt(ntot),
     /     i,j,k,ka,indi
      double precision xtw(p*r,nmax*r),ey(ntot,r),
     /     xtwey(p*r),sum
      do 100 i=1,p*r
         sum=dble(0.)
         do 50 j=1,r
            indi=0
            do 30 k=st,fin
               if(patt(k).ne.0) then
                  indi=indi+1
                  ka=(j-1)*ni+indi
                  sum=sum+xtw(i,ka)*ey(k,j)
               endif
 30         continue
 50      continue
         xtwey(i)=xtwey(i)+sum
 100  continue
      return
      end
C***********************************************************************
      subroutine mkuszteeyxyxt(m,r,q,nmax,ni,s,wkqnm,eyxyxt,wkrqrq2)
C calculates 
C U_i%*%(inv(sigma) Otimes t(Z_i)) %*% 
C E( vec(y_i - X_i beta)%*%t(vec(y_i - X_i beta)) | y_i(obs),theta )%*%
C (inv(sigma) Otimes Z_i)%*%U_i, FOR SUBJECT S 
C note that U_i%*%(inv(sigma) Otimes t(Z_i)) is already contained 
C in wkqnm.
C result is stored in wkqq2(the upper tri part)
      integer m,r,q,nmax,ni,s,i,j,k,l
      double precision wkqnm(r*q,r*nmax,m),eyxyxt(r*nmax,r*nmax),
     /     wkrqrq2(r*q,r*q),sum
      do 2 i=1,r*q
         do 1 j=i,r*q
            wkrqrq2(i,j)=dble(0.)
 1       continue
 2    continue
      do 100 i=1,r*q
         do 50 j=i,r*q
            sum=dble(0.)
            do 40 k=1,r*ni
               do 30 l=1,k
                  sum=sum+wkqnm(i,l,s)*eyxyxt(l,k)*wkqnm(j,k,s)
 30            continue
               do 20 l=k+1,r*ni
                  sum=sum+wkqnm(i,l,s)*eyxyxt(k,l)*wkqnm(j,k,s)
 20            continue
 40         continue
            wkrqrq2(i,j)=sum
 50      continue
 100  continue
      return
      end
C***********************************************************************
      subroutine fscov(ntot,nmax,npatt,m,r,pcol,q,zcol,
     /     ist,ifin,nstari,patt,rmat,lmc,loc,mc,mc1,oc,oc2,
     /     pred,u,ztz,psi,opsi,wkrr1,wkrr2,wkqnm,wkqnm1,wkeyxyxt,
     /     wkqqu,wkqq1,wkqq2,wkqq3,varb,wkrrpt,wkrrb21,
     /     wxbeta,ggs,wkg,wkgg,w,sigma,osigma,
     /     y,ey,eyyt,eyxyxt,trwex,msg,ldsi,ldps,lduu)
C Fisher scoring on the variance components given beta. 
C scoring is performed on eta^* rather than eta:
C eta^* is the transformed vector from eta, in which
C log of the diag elements of inv(psi) and inv(sigma) are used.
C     opsi = old value of psi
C      psi = new value
C     sigma = old value of sigma, sigma = new value of sigma
C Error messages returned through msg:
C     0 = no error
C     5 = non-positive definite psi
C     6 = non-positive definite sigma
C     7 = non-positive definite wkgg
      implicit none
      integer ntot,nmax,npatt,m,r,pcol,q,zcol(q),
     /     ist(m),ifin(m),nstari(m),patt(ntot),rmat(npatt,r),
     /     lmc,loc,mc(lmc),mc1(lmc),oc(loc),oc2(loc),
     /     ggs,msg,i,j,st,
     /     fin,g,gs,s,ii,jj,gi,gj,jjmin,err,ni,ind
      double precision pred(ntot,pcol),u(r*q,r*q,m),
     /     ztz(q,q,m),psi(r*q,r*q),opsi(r*q,r*q),wkrr1(r,r),
     /     wkrr2(r,r),wkqnm(r*q,r*nmax,m),wkqnm1(r*nmax,r*nmax),
     /     wkeyxyxt(r*nmax,r*nmax),
     /     wkqqu(r*q,r*q,m),wkqq1(r*q,r*q),wkqq2(r*q,r*q),
     /     wkqq3(r*q,r*q),varb(r*q,r*q,m),
     /     wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),wxbeta(ntot,r),
     /     wkg(ggs),w(r*nmax,r*nmax,m),
     /     wkgg(ggs,ggs),sigma(r,r),osigma(r,r),
     /     y(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),ldsi,ldps,lduu,trwex,sum,deflate
      double precision trahah,trahaj,trajaj,treyxyxti,
     /     treyxyxtkl,truztzh,truztzhk,trhshoztzu,truztzhuztzh,
     /     trhsjoztzu,trjsjoztzu,truztzjuztzj,truztzhuztzj,
     /     truiulztz,truiulkztz,truijuztzk,truijuztzlk
      g=(r*q)*((r*q)+1)/2
      gs=r*(r+1)/2
      msg=0
      trwex=dble(0.)
C*****
      do 2 i=1,r
         do 1 j=i+1,r
            osigma(j,i)=osigma(i,j)
 1       continue
 2    continue
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=1,ggs
         wkg(i)=dble(0.)
         do 7 j=i,ggs
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue      
C     *** main loop to accumulate wkg and wkgg *********
      do 390 s=1,m
C        *** put U_i%*%t(Z_i)%*%inv(V_i)%*%Z_i%*%U_i into wkqq2 ***
C        *** for some reason we don't need to calculate this!!!!!***
C        *** in the univariate case we do need it****************
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
         call mkyyt(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,y,eyyt)
         call mkeyyt(ntot,nmax,npatt,m,r,st,fin,ni,patt,
     /        rmat,s,pcol,q,zcol,lmc,loc,mc,mc1,oc,oc2,
     /        pred,varb,wkrrpt,wkrrb21,ey,eyyt)
         call mkeyxyxt(ntot,nmax,r,st,fin,ni,patt,wxbeta,ey,eyyt,eyxyxt)
         call lltrwex(nmax,m,r,ni,s,w,eyxyxt,trwex)
         do 25 i=1,r*q
            do 24 j=i+1,r*q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
         do 27 i=1,q
            do 26 j=i+1,q
               ztz(j,i,s)=ztz(i,j,s)
 26         continue
 27      continue
C*********************
C******************
         call mkuszteeyxyxt(m,r,q,nmax,ni,s,wkqnm,eyxyxt,wkqq2)
C        *** store (psi-U_i) in wkqq1
         do 45 i=1,r*q
            do 44 j=i,r*q
               wkqq1(i,j)=opsi(i,j)-u(i,j,s)
               if(i.ne.j) then
                  wkqq1(j,i)=wkqq1(i,j)
                  wkqq2(j,i)=wkqq2(i,j)
               endif
 44         continue
 45      continue
         do 47 i=1,r*q
            do 46 j=1,r*q
               wkqq3(i,j)=wkqq1(i,j)-wkqq2(i,j)
 46         continue
 47      continue
         call mkwkeyxyxt(ntot,nmax,m,r,q,pcol,zcol,st,fin,s,
     /     patt,ni,pred,wkqnm,wkqnm1,eyxyxt,wkeyxyxt)
C     *** now we're ready to accumulate wkg and wkgg
         gi=0
         ind=r*q
         do 200 i=1,r*q
            do 190 j=i,r*q
               gi=gi+1
               if(i.eq.j) then
                  wkg(gi)=wkg(gi)+0.5*wkqq3(i,i)
               else
                  wkg(gi)=wkg(gi)+wkqq3(i,j)
               endif
               gj=gi-1
               do 180 ii=i,r*q
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 170 jj=jjmin,r*q
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahah(ind,wkqq1,i,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(ind,wkqq1,i,ii,jj)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trahaj(ind,wkqq1,ii,i,j)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          trajaj(ind,wkqq1,i,j,ii,jj)
                        endif
                     endif
 170              continue
 180           continue
 190        continue
 200     continue
         gi=g
         do 300 i=1,r
            do 250 j=i,r
               gi=gi+1
               if(i.eq.j) then
                  wkg(gi)=wkg(gi)+0.5*ni*osigma(i,i)-
     /                 0.5*truztzh(s,r,q,m,i,ztz,u)-
     /                 0.5*treyxyxti(r,nmax,ni,i,wkeyxyxt)
               else
                  wkg(gi)=wkg(gi)+ni*osigma(i,j)-
     /                 0.5*truztzhk(s,r,q,m,i,j,ztz,u)-
     /                 0.5*treyxyxtkl(r,nmax,ni,i,j,wkeyxyxt)
               endif
               gj=gi-1
               do 240 ii=i,r 
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 230 jj=jjmin,r
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahah(r,osigma,i,ii)-
     /                          dble(2.)*trhshoztzu(s,r,q,m,i,ii,
     /                          osigma,ztz,u)+
     /                          truztzhuztzh(s,r,q,m,i,ii,ztz,u)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahaj(r,sigma,i,ii,jj)-
     /                          dble(2.)*trhsjoztzu(s,r,q,m,i,
     /                          ii,jj,osigma,ztz,u)+
     /                          truztzhuztzj(s,r,q,m,i,ii,jj,ztz,u)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahaj(r,osigma,ii,i,j)-
     /                          dble(2.)*trhsjoztzu(s,r,q,m,ii,i,j,
     /                          osigma,ztz,u)+
     /                          truztzhuztzj(s,r,q,m,ii,i,j,ztz,u)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trajaj(r,osigma,i,j,ii,jj)-
     /                          dble(2.)*trjsjoztzu(s,r,q,m,i,j,ii,
     /                          jj,osigma,ztz,u)+
     /                          truztzjuztzj(s,r,q,m,i,j,ii,jj,ztz,u)
                        endif
                     endif
 230              continue
 240           continue
 250        continue
 300     continue
C  ***** so far we filled the sub-diagonal matrices of wkgg, now
C  **** we need to fill in upper right sub-matrix (g by gs, cross derivatives)
         gi=0
         do 385 i=1,r*q
            do 380 j=i,r*q 
               gi=gi+1
               gj=g
               do 370 ii=1,r
                  do 360 jj=ii,r
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          truiulztz(s,r,q,m,ztz,u,i,ii)
                         else
                          wkgg(gi,gj)=wkgg(gi,gj)+
     /                          truiulkztz(s,r,q,m,ztz,u,i,ii,jj)
                       endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          truijuztzk(s,r,q,m,ztz,u,i,j,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          truijuztzlk(s,r,q,m,ztz,u,i,j,ii,jj)
                        endif
                     endif
 360              continue
 370           continue
 380        continue
 385     continue
 390  continue
C***put the inverse of opsi and osigma into wkqq1 and wkrr1,respectively:
      do 397 i=1,r*q
         do 396 j=i,r*q
            wkqq2(i,j)=opsi(i,j)
 396     continue
 397  continue
      do 399 i=1,r
         do 398 j=i,r
            wkrr2(i,j)=osigma(i,j)
 398     continue
 399  continue
      call chfce(r*q,r*q,wkqq2,err)
      if(err.eq.1) then
         msg=5
         goto 999
      endif
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) then
         msg=6
         goto 999
      endif
      call bkslv(r*q,r*q,wkqq2)
      call mm(r*q,r*q,wkqq2,wkqq1)
      call bkslv(r,r,wkrr2)
      call mm(r,r,wkrr2,wkrr1)
      gi=0
C********
      do 403 i=1,r*q
         do 402 j=i,r*q
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=wkg(gi)*wkqq1(i,j)
            else
               wkg(gi)=wkg(gi)
            endif
            gj=gi-1
            do 401 ii=i,r*q
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 400 jj=jjmin,r*q
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)*
     /                       wkqq1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(ii,ii)
                     endif
                  endif
                  wkgg(gi,gj)=wkgg(gi,gj)
 400           continue
 401        continue
 402     continue
 403  continue
      gi=g
      do 407 i=1,r
         do 406 j=i,r
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=wkg(gi)*wkrr1(i,j)
            else
               wkg(gi)=wkg(gi)
            endif
            gj=gi-1
            do 405 ii=i,r 
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 404 jj=jjmin,r
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(i,i)*
     /                       wkrr1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)
                     endif
                  endif
                  wkgg(gi,gj)=wkgg(gi,gj)
 404           continue
 405        continue
 406     continue
 407  continue
      gi=0
      do 411 i=1,r*q
         do 410 j=i,r*q 
            gi=gi+1
            gj=g
            do 409 ii=1,r
               do 408 jj=ii,r
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)*wkqq1(i,i)
                    else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkqq1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)
                     endif
                  endif
 408           continue
 409        continue
 410     continue
 411  continue
C now compute wkg
      do 600 i=1,ggs
         do 590 j=i+1,ggs
            wkgg(j,i)=wkgg(i,j)
 590     continue
 600  continue
      gi=0
      do 417 i=1,r*q
         do 416 j=i,r*q
            gi=gi+1
            gj=0
            sum=dble(0.)
            do 415 ii=1,r*q
               do 414 jj=ii,r*q
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkqq1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
                  endif
 414           continue
 415        continue
            do 413 ii=1,r
               do 412 jj=ii,r
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkrr1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkrr1(ii,jj)
                  endif
 412           continue
 413        continue
            wkg(gi)=wkg(gi)+sum
 416     continue
 417  continue
C****
      do 423 i=1,r
         do 422 j=i,r
            gi=gi+1
            gj=0
            sum=dble(0.)
            do 421 ii=1,r*q
               do 420 jj=ii,r*q
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkqq1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkqq1(ii,jj)
                  endif
 420           continue
 421        continue
            do 419 ii=1,r
               do 418 jj=ii,r
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkrr1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkrr1(ii,jj)
                  endif
 418           continue
 419        continue
            wkg(gi)=wkg(gi)+sum
 422     continue
 423  continue
C********
C     *** solve the linear system, storing the result in wkg ******
      call chfce(ggs,ggs,wkgg,err)
      if(err.eq.1) then
         msg=7
         goto 999
      endif
      call bkslv(ggs,ggs,wkgg)
      do 510 i=ggs,1,-1
         sum=dble(0.)
         do 519 j=1,i
            sum=sum+wkgg(j,i)*wkg(j)
 519     continue
         wkg(i)=sum
 510  continue
      do 530 i=1,ggs
         sum=dble(0.)
         do 529 j=i,ggs
            sum=sum+wkgg(i,j)*wkg(j)
 529     continue
         wkg(i)=sum
 530  continue
C     *** invert wkg, putting the result into psi **************
C     *** step-halving is used here if wkg is not pos.def. ****
      deflate=dble(1.)      
      do 540 i=1,r*q
         wkqq1(i,i)=dlog(wkqq1(i,i))
 540  continue
      do 542 i=1,r
         wkrr1(i,i)=dlog(wkrr1(i,i))
 542  continue
 425  continue
      gi=0
      do 430 i=1,r*q
         do 429 j=i,r*q
            gi=gi+1
            wkqq2(i,j)=wkqq1(i,j)+deflate*(wkg(gi)-wkqq1(i,j))
 429     continue
         wkqq2(i,i)=dexp(wkqq2(i,i))
 430  continue      
      do 428 i=1,r
         do 427 j=i,r
            gi=gi+1
            wkrr2(i,j)=wkrr1(i,j)+deflate*(wkg(gi)-wkrr1(i,j))
 427     continue
         wkrr2(i,i)=dexp(wkrr2(i,i))
 428  continue
      call chfce(r*q,r*q,wkqq2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 425
      endif
      call bkslv(r*q,r*q,wkqq2)
      call mm(r*q,r*q,wkqq2,psi)
      call bkslv(r,r,wkrr2)
      call mm(r,r,wkrr2,sigma)
      do 440 i=1,r*q
         do 439 j=i+1,r*q
            psi(j,i)=psi(i,j)
 439     continue
 440  continue
      do 450 i=1,r
         do 449 j=i+1,r
            sigma(j,i)=sigma(i,j)
 449     continue
 450  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine fscovbd(ntot,nmax,npatt,m,r,pcol,q,zcol,
     /     ist,ifin,nstari,patt,rmat,lmc,loc,mc,mc1,oc,oc2,
     /     pred,u,ztz,psi,opsi,wkrr1,wkrr2,wkqnm,wkqnm1,wkeyxyxt,
     /     wkrqrq1,wkrqrq2,wkqq1bd,wkqq2bd,wkqq3,varb,wkrrpt,wkrrb21,
     /     wxbeta,ggs,wkg,wkgg,w,sigma,osigma,
     /     y,ey,eyyt,eyxyxt,trwex,msg)
C Block-diagonal version of the Fisher scoring on the variance components
C Returned error messages are same. 
      implicit none
      integer ntot,nmax,npatt,m,r,pcol,q,zcol(q),
     /     ist(m),ifin(m),nstari(m),patt(ntot),rmat(npatt,r),
     /     lmc,loc,mc(lmc),mc1(lmc),oc(loc),oc2(loc),
     /     ggs,msg,i,j,st,l,ll,
     /     fin,g,gs,s,ii,jj,gi,gj,jjmin,err,ni,ia,ja,iia,jja
      double precision pred(ntot,pcol),u(r*q,r*q,m),
     /     ztz(q,q,m),psi(q,q,r),opsi(q,q,r),wkrr1(r,r),
     /     wkrr2(r,r),wkqnm(r*q,r*nmax,m),wkqnm1(r*nmax,r*nmax),
     /     wkeyxyxt(r*nmax,r*nmax),wkrqrq1(r*q,r*q),wkrqrq2(r*q,r*q),
     /     wkqq1bd(q,q,r),wkqq2bd(q,q),wkqq3(r*q,r*q),varb(r*q,r*q,m),
     /     wkrrpt(r,r,npatt),wkrrb21(r,r,npatt),wxbeta(ntot,r),
     /     wkg(ggs),w(r*nmax,r*nmax,m),
     /     wkgg(ggs,ggs),sigma(r,r),osigma(r,r),
     /     y(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),trwex,sum,deflate
      double precision trahah,trahaj,trajaj,trahahbd,trahajbd1,
     /     trahajbd2,trajajbd,treyxyxti,
     /     treyxyxtkl,truztzh,truztzhk,trhshoztzu,truztzhuztzh,
     /     trhsjoztzu,trjsjoztzu,truztzjuztzj,truztzhuztzj,
     /     truiulztzbd,truiulkztzbd,truijuztzkbd,truijuztzlkbd
      g=r*q*(q+1)/2
      gs=r*(r+1)/2
      msg=0
      trwex=dble(0.)
C*****
      do 2 i=1,r
         do 1 j=i+1,r
            osigma(j,i)=osigma(i,j)
 1       continue
 2    continue
      do 5 l=1,r
         do 4 i=1,q
            do 3 j=i+1,q
               opsi(j,i,l)=opsi(i,j,l)
 3          continue
 4       continue
 5    continue
C     *** initialize the workspaces wkg and wkgg *******
      do 8 i=1,ggs
         wkg(i)=dble(0.)
         do 7 j=i,ggs
            wkgg(i,j)=dble(0.)
 7       continue
 8    continue      
C     *** main loop to accumulate wkg and wkgg *********
      do 390 s=1,m
C        *** put U_i%*%t(Z_i)%*%inv(V_i)%*%Z_i%*%U_i into wkqq2bd ***
C        *** for some reason we don't need to calculate this!!!!!***
C        *** in the univariate case we do need it****************
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
         call mkyyt(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,y,eyyt)
         call mkeyyt(ntot,nmax,npatt,m,r,st,fin,ni,patt,
     /        rmat,s,pcol,q,zcol,lmc,loc,mc,mc1,oc,oc2,
     /        pred,varb,wkrrpt,wkrrb21,ey,eyyt)
         call mkeyxyxt(ntot,nmax,r,st,fin,ni,patt,wxbeta,ey,eyyt,eyxyxt)
         call lltrwex(nmax,m,r,ni,s,w,eyxyxt,trwex)
         do 25 i=1,r*q
            do 24 j=i+1,r*q
               u(j,i,s)=u(i,j,s)
 24         continue
 25      continue
         do 27 i=1,q
            do 26 j=i+1,q
               ztz(j,i,s)=ztz(i,j,s)
 26         continue
 27      continue
C*********************
C******************
         call mkuszteeyxyxt(m,r,q,nmax,ni,s,wkqnm,eyxyxt,wkrqrq2)
C        *** store (psi-U_i) in wkrqrq1
         do 36 l=1,r
            do 35 ll=l,r
               do 34 i=1,q
                  do 33 j=1,q
                     ia=(l-1)*q+i
                     ja=(ll-1)*q+j
                     if(l.eq.ll) then
                        wkrqrq1(ia,ja)=opsi(i,j,l)-u(ia,ja,s)
                     else
                        wkrqrq1(ia,ja)=-u(ia,ja,s)
                     endif
                     if(ia.ne.ja) then
                        wkrqrq1(ja,ia)=wkrqrq1(ia,ja)
                     endif
 33               continue
 34            continue
 35         continue
 36      continue
C note that we will only use the block diagonal (q by q) matrices in the 
C calculations for the first derivatives wrt elements of psi
CCC         do 40 l=1,r
         do 39 i=1,r*q
            do 38 j=i,r*q
               wkqq3(i,j)=wkrqrq1(i,j)-wkrqrq2(i,j)
               if(i.ne.j) then
                  wkqq3(j,i)=wkqq3(i,j)
               endif
 38         continue
 39      continue
         call mkwkeyxyxt(ntot,nmax,m,r,q,pcol,zcol,st,fin,s,
     /     patt,ni,pred,wkqnm,wkqnm1,eyxyxt,wkeyxyxt)
C     *** now we're ready to accumulate wkg and wkgg
C   NOTE THAT IN BLOCK DIAG. VERSION THE DIMENSIONS OF WKGG WILL BE:
C   r*(r+1)/2 + r*q*(q+1)/2 instead of r*(r+1)/2 + r*q*(r*q+1)/2 
         gi=0
         do 200 l=1,r
            do 190 i=1,q
               do 185 j=i,q
                  gi=gi+1
                  ia=(l-1)*q+i
                  ja=(l-1)*q+j
                  if(i.eq.j) then
                     wkg(gi)=wkg(gi)+0.5*wkqq3(ia,ia)
                  else
                     wkg(gi)=wkg(gi)+wkqq3(ia,ja)
                  endif
                  gj=gi-1
                  do 180 ll=l,r
                     if(l.eq.ll) then                        
                        do 170 ii=i,q
                           if(ii.eq.i) then
                              jjmin=j
                           else
                              jjmin=ii
                           endif
                           do 165 jj=jjmin,q
                              gj=gj+1
                              if(i.eq.j) then
                                 if(ii.eq.jj) then
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahahbd(r*q,wkrqrq1,q,
     /                                   l,ll,i,ii)
                                 else
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahajbd1(r*q,wkrqrq1,q,l,
     /                                   ll,i,ii,jj)
                                 endif
                              else
                                 if(ii.eq.jj) then
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahajbd2(r*q,wkrqrq1,q,l,ll,
     /                                   i,j,ii)
                                 else
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trajajbd(r*q,wkrqrq1,q,l,ll,i,
     /                                   j,ii,jj)
                                 endif
                              endif
 165                       continue
 170                    continue
                     else
                        do 175 ii=1,q
                           do 174 jj=ii,q
                              gj=gj+1
                              if(i.eq.j) then
                                 if(ii.eq.jj) then
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahahbd(r*q,wkrqrq1,q,
     /                                   l,ll,i,ii)
                                 else
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahajbd1(r*q,wkrqrq1,q,l,
     /                                   ll,i,ii,jj)
                                 endif
                              else
                                 if(ii.eq.jj) then
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trahajbd2(r*q,wkrqrq1,q,l,ll,
     /                                   i,j,ii)
                                 else
                                    wkgg(gi,gj)=wkgg(gi,gj)+
     /                                   trajajbd(r*q,wkrqrq1,q,l,ll,i,
     /                                   j,ii,jj)
                                 endif
                              endif
 174                       continue
 175                    continue
                     endif
 180              continue
 185           continue
 190        continue
 200     continue
         gi=g
         do 300 i=1,r
            do 250 j=i,r
               gi=gi+1
               if(i.eq.j) then
                  wkg(gi)=wkg(gi)+0.5*ni*osigma(i,i)-
     /                 0.5*truztzh(s,r,q,m,i,ztz,u)-
     /                 0.5*treyxyxti(r,nmax,ni,i,wkeyxyxt)
               else
                  wkg(gi)=wkg(gi)+ni*osigma(i,j)-
     /                 0.5*truztzhk(s,r,q,m,i,j,ztz,u)-
     /                 0.5*treyxyxtkl(r,nmax,ni,i,j,wkeyxyxt)
               endif
               gj=gi-1
               do 240 ii=i,r 
                  if(ii.eq.i) then
                     jjmin=j
                  else
                     jjmin=ii
                  endif
                  do 230 jj=jjmin,r
                     gj=gj+1
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahah(r,osigma,i,ii)-
     /                          dble(2.)*trhshoztzu(s,r,q,m,i,ii,
     /                          osigma,ztz,u)+
     /                          truztzhuztzh(s,r,q,m,i,ii,ztz,u)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahaj(r,sigma,i,ii,jj)-
     /                          dble(2.)*trhsjoztzu(s,r,q,m,i,
     /                          ii,jj,sigma,ztz,u)+
     /                          truztzhuztzj(s,r,q,m,i,ii,jj,ztz,u)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trahaj(r,osigma,ii,i,j)-
     /                          dble(2.)*trhsjoztzu(s,r,q,m,ii,i,j,
     /                          osigma,ztz,u)+
     /                          truztzhuztzj(s,r,q,m,ii,i,j,ztz,u)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)+
     /                          ni*trajaj(r,osigma,i,j,ii,jj)-
     /                          dble(2.)*trjsjoztzu(s,r,q,m,i,j,ii,
     /                          jj,osigma,ztz,u)+
     /                          truztzjuztzj(s,r,q,m,i,j,ii,jj,ztz,u)
                        endif
                     endif
 230              continue
 240           continue
 250        continue
 300     continue
C  ***** so far we filled the sub-diagonal matrices of wkgg, now
C  **** we need to fill in upper right sub-matrix (g by gs, cross derivatives)
         gi=0
         do 386 l=1,r
            do 385 i=1,q
               do 384 j=i,q 
                  gi=gi+1
                  gj=g
                  do 370 ii=1,r
                     do 360 jj=ii,r
                        gj=gj+1
                        if(i.eq.j) then
                           if(ii.eq.jj) then
                              wkgg(gi,gj)=wkgg(gi,gj)+
     /                             truiulztzbd(s,r,q,m,ztz,u,l,i,ii)
                           else
                              wkgg(gi,gj)=wkgg(gi,gj)+
     /                             truiulkztzbd(s,r,q,m,ztz,u,l,i,ii,jj)
                           endif
                        else
                           if(ii.eq.jj) then
                              wkgg(gi,gj)=wkgg(gi,gj)+
     /                             truijuztzkbd(s,r,q,m,ztz,u,l,i,j,ii)
                           else
                              wkgg(gi,gj)=wkgg(gi,gj)+
     /                             truijuztzlkbd(s,r,q,m,ztz,u,l,
     /                             i,j,ii,jj)
                           endif
                        endif
 360                 continue
 370              continue
 384           continue
 385        continue
 386     continue
 390  continue
C***  put the inverse of opsi and osigma into wkrqrq1(using the fact that psi
C is block diagonal) and wkrr1,respectively:
      do 397 i=1,r*q
         do 396 j=i,r*q
            wkrqrq1(i,j)=dble(0.)
 396     continue
 397  continue
      do 395 l=1,r
         do 394 i=1,q
            do 393 j=i,q
               wkqq1bd(i,j,l)=opsi(i,j,l)
 393        continue
 394     continue
         call chle(q,q,r,wkqq1bd,l,err)
         if(err.eq.1) then
            msg=4
            goto 999
         endif
         call bkslvl(q,q,r,wkqq1bd,l)
         call mmul(q,q,r,wkqq1bd,l,wkqq2bd)
         do 392 i=1,q
            do 391 j=i,q
               ii=(l-1)*q+i
               jj=(l-1)*q+j
               wkrqrq1(ii,jj)=wkqq2bd(i,j)
 391        continue
 392     continue
 395  continue
      do 399 i=1,r
         do 398 j=i,r
            wkrr2(i,j)=osigma(i,j)
 398     continue
 399  continue
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) then
         msg=5
         goto 999
      endif
      call bkslv(r,r,wkrr2)
      call mm(r,r,wkrr2,wkrr1)
      gi=0
C********
      do 407 l=1,r
         do 406 i=1,q
            do 405 j=i,q
               ia=(l-1)*q+i
               ja=(l-1)*q+j
               gi=gi+1
               if(i.eq.j) then
                  wkg(gi)=wkg(gi)*wkrqrq1(ia,ja)
               else
                  wkg(gi)=wkg(gi)
               endif
               gj=gi-1
               do 404 ll=l,r
                  if(l.eq.ll) then
                     do 401 ii=i,q
                        if(ii.eq.i) then
                           jjmin=j
                        else
                           jjmin=ii
                        endif
                        do 400 jj=jjmin,q
                           gj=gj+1
                           iia=(ll-1)*q+jj
                           wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                           if(i.eq.j) then
                              if(ii.eq.jj) then
                                 wkgg(gi,gj)=wkgg(gi,gj)*wkrqrq1(ia,ia)*
     /                                wkrqrq1(iia,iia)
                              else
                                 wkgg(gi,gj)=wkgg(gi,gj)*wkrqrq1(ia,ia)
                              endif
                           else
                              if(ii.eq.jj) then
                                 wkgg(gi,gj)=wkgg(gi,gj)*
     /                                wkrqrq1(iia,iia)
                              endif
                           endif
                           wkgg(gi,gj)=wkgg(gi,gj)
 400                    continue
 401                 continue
                  else
                     do 403 ii=1,q
                        do 402 jj=ii,q
                           gj=gj+1
                           iia=(ll-1)*q+jj
                           wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                           if(i.eq.j) then
                              if(ii.eq.jj) then
                                 wkgg(gi,gj)=wkgg(gi,gj)*wkrqrq1(ia,ia)*
     /                                wkrqrq1(iia,iia)
                              else
                                 wkgg(gi,gj)=wkgg(gi,gj)*wkrqrq1(ia,ia)
                              endif
                           else
                              if(ii.eq.jj) then
                                 wkgg(gi,gj)=wkgg(gi,gj)*
     /                                wkrqrq1(iia,iia)
                              endif
                           endif
                           wkgg(gi,gj)=wkgg(gi,gj)
 402                    continue
 403                 continue
                  endif
 404           continue
 405        continue
 406     continue
 407  continue
      gi=g
      do 700 i=1,r
         do 690 j=i,r
            gi=gi+1
            if(i.eq.j) then
               wkg(gi)=wkg(gi)*wkrr1(i,j)
            else
               wkg(gi)=wkg(gi)
            endif
            gj=gi-1
            do 680 ii=i,r 
               if(ii.eq.i) then
                  jjmin=j
               else
                  jjmin=ii
               endif
               do 670 jj=jjmin,r
                  gj=gj+1
                  wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                  if(i.eq.j) then
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(i,i)*
     /                       wkrr1(ii,ii)
                     else
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(i,i)
                     endif
                  else
                     if(ii.eq.jj) then
                        wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)
                     endif
                  endif
                  wkgg(gi,gj)=wkgg(gi,gj)
 670           continue
 680        continue
 690     continue
 700  continue
      gi=0
      do 413 l=1,r
         do 412 i=1,q
            do 411 j=i,q
               ia=(l-1)*q+i
               gi=gi+1
               gj=g
               do 410 ii=1,r
                  do 409 jj=ii,r
                     gj=gj+1
                     wkgg(gi,gj)=wkgg(gi,gj)/dble(2.)
                     if(i.eq.j) then
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)*
     /                          wkrqrq1(ia,ia)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)*wkrqrq1(ia,ia)
                        endif
                     else
                        if(ii.eq.jj) then
                           wkgg(gi,gj)=wkgg(gi,gj)*wkrr1(ii,ii)
                        else
                           wkgg(gi,gj)=wkgg(gi,gj)
                        endif
                     endif
 409              continue
 410           continue
 411        continue
 412     continue
 413  continue
CC now compute wkg
      do 600 i=1,ggs
         do 590 j=i+1,ggs
            wkgg(j,i)=wkgg(i,j)
 590     continue
 600  continue
      gi=0
      do 430 l=1,r
         do 429 i=1,q
            do 428 j=i,q
               ia=(l-1)*q+i
               ja=(l-1)*q+j
               gi=gi+1
               gj=0
               sum=dble(0.)
               do 427 ll=1,r
                  do 426 ii=1,q
                     do 425 jj=ii,q
                        iia=(ll-1)*q+ii
                        jja=(ll-1)*q+jj
                        gj=gj+1
                        if(ii.eq.jj) then
                           sum=sum+wkgg(gi,gj)*dlog(wkrqrq1(iia,iia))
                        else
                           sum=sum+wkgg(gi,gj)*wkrqrq1(iia,jja)
                        endif
 425                 continue
 426              continue
 427           continue
               do 424 ii=1,r
                  do 423 jj=ii,r
                     gj=gj+1
                     if(ii.eq.jj) then
                        sum=sum+wkgg(gi,gj)*dlog(wkrr1(ii,ii))
                     else
                        sum=sum+wkgg(gi,gj)*wkrr1(ii,jj)
                     endif
 423              continue
 424           continue
               wkg(gi)=wkg(gi)+sum
 428        continue
 429     continue
 430  continue
C****
      do 440 i=1,r
         do 439 j=i,r
            gi=gi+1
            gj=0
            sum=dble(0.)
            do 438 l=1,r
               do 437 ii=1,q
                  do 436 jj=ii,q
                     ia=(l-1)*q+ii
                     ja=(l-1)*q+jj
                     gj=gj+1
                     if(ii.eq.jj) then
                        sum=sum+wkgg(gi,gj)*dlog(wkrqrq1(ia,ia))
                     else
                        sum=sum+wkgg(gi,gj)*wkrqrq1(ia,ja)
                     endif
 436              continue
 437           continue
 438        continue
            do 435 ii=1,r
               do 434 jj=ii,r
                  gj=gj+1
                  if(ii.eq.jj) then
                     sum=sum+wkgg(gi,gj)*dlog(wkrr1(ii,ii))
                  else
                     sum=sum+wkgg(gi,gj)*wkrr1(ii,jj)
                  endif
 434           continue
 435        continue
            wkg(gi)=wkg(gi)+sum
 439     continue
 440  continue
C********
C     *** solve the linear system, storing the result in wkg ******
      call chfce(ggs,ggs,wkgg,err)
      if(err.eq.1) then
         msg=3
         goto 999
      endif
      call bkslv(ggs,ggs,wkgg)
      do 442 i=ggs,1,-1
         sum=dble(0.)
         do 441 j=1,i
            sum=sum+wkgg(j,i)*wkg(j)
 441     continue
         wkg(i)=sum
 442  continue
      do 445 i=1,ggs
         sum=dble(0.)
         do 444 j=i,ggs
            sum=sum+wkgg(i,j)*wkg(j)
 444     continue
         wkg(i)=sum
 445  continue
C     *** invert wkg, putting the result into psi **************
C     *** step-halving is used here if wkg is not pos.def. ****
      deflate=dble(1.)      
      do 455 i=1,r*q
         wkrqrq1(i,i)=dlog(wkrqrq1(i,i))
 455  continue
      do 460 i=1,r
         wkrr1(i,i)=dlog(wkrr1(i,i))
 460  continue
 475  continue
      gi=0
      do 470 l=1,r
         do 465 i=1,q
            do 464 j=i,q
               ia=(l-1)*q+i
               ja=(l-1)*q+j
               gi=gi+1
               wkqq1bd(i,j,l)=wkrqrq1(ia,ja)+deflate*(wkg(gi)-
     /              wkrqrq1(ia,ja))
 464        continue
            wkqq1bd(i,i,l)=dexp(wkqq1bd(i,i,l))
 465     continue
 470  continue      
      do 485 i=1,r
         do 480 j=i,r
            gi=gi+1
            wkrr2(i,j)=wkrr1(i,j)+deflate*(wkg(gi)-wkrr1(i,j))
 480     continue
         wkrr2(i,i)=dexp(wkrr2(i,i))
 485  continue
      do 490 l=1,r
         call chle(q,q,r,wkqq1bd,l,err)
         if(err.eq.1) then
            deflate=deflate/dfloat(2)
            goto 475
         endif
         call bkslvl(q,q,r,wkqq1bd,l)
         call mmul(q,q,r,wkqq1bd,l,wkqq2bd)
         do 484 i=1,q
            do 483 j=i,q
               psi(i,j,l)=wkqq2bd(i,j)
 483        continue
 484     continue
 490  continue
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) then
         deflate=deflate/dfloat(2)
         goto 475
      endif
      call bkslv(r,r,wkrr2)
      call mm(r,r,wkrr2,sigma)
      do 497 l=1,r
         do 496 i=1,q
            do 495 j=i+1,q
               psi(j,i,l)=psi(i,j,l)
 495        continue
 496     continue
 497  continue
      do 499 i=1,r
         do 498 j=i+1,r
            sigma(j,i)=sigma(i,j)
 498     continue
 499  continue
 999  continue
      return
      end
C***********************************************************************
      subroutine mkwkeyxyxt(ntot,nmax,m,r,q,pcol,zcol,st,fin,s,
     /     patt,ni,pred,wkqnm,wkqnm1,eyxyxt,wkeyxyxt)
C let 
C eyxyxt = E( vec(y_i - X_i beta)%*%t(vec(y_i - X_i beta)) | y_i(obs),theta ). 
C  This subroutine calculates 
C  eyxyxt - 2*(I_r Otimes Z_i) %*% U_i%*%(inv(sigma) Otimes t(Z_i)) %*% eyxyxt
C + (I_r Otimes Z_i) %*% U_i%*%(inv(sigma) Otimes t(Z_i))%*%eyxyxt%*%
C  (inv(sigma) Otimes Z_i)%*%U_i%*%(I_r Otimes t(Z_i)), for subject s,
C  storing in wkeyxyxt.
C  wkqnm1 is used to store 
C     (I_r Otimes Z_i) %*% U_i%*%(inv(sigma) Otimes t(Z_i))
C  for subject s
      implicit none
      integer ntot,nmax,m,r,q,pcol,zcol(q),st,fin,s,
     /     patt(ntot),ni,i,ja,indi,j,ii,
     /     indj,k,l,ia,iia
      double precision pred(ntot,pcol),wkqnm(r*q,r*nmax,m),
     /     wkqnm1(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),wkeyxyxt(r*nmax,r*nmax),sum
C first initialize wkqnm1 and wkeyxyxt:
      do 2 i=1,r*nmax
         do 1 j=1,r*nmax
            wkqnm1(i,j)=dble(0.)
            wkeyxyxt(i,j)=dble(0.)
 1       continue
 2    continue
      do 100 i=1,r
         indi=0
         do 96 k=st,fin
            if(patt(k).ne.0) then
               indi=indi+1
               ii=(i-1)*ni+indi
               do 90 j=1,r
                  indj=0
                  do 80 l=st,fin
                     if(patt(l).ne.0) then
                        indj=indj+1
                        ja=(j-1)*ni+indj
                        sum=dble(0.)
                        do 70 ia=1,q
                           iia=(i-1)*q+ia
                           sum=sum+pred(k,zcol(ia))*
     /                          wkqnm(iia,ja,s)
 70                     continue
                        wkqnm1(ii,ja)=sum
                     endif
 80               continue
 90            continue
            endif
 96      continue
 100  continue
      do 150 i=1,r*ni
         do 140 j=1,r*ni
            sum=dble(0.)
            do 130 k=1,j
               sum=sum+wkqnm1(i,k)*eyxyxt(k,j)
 130        continue
            do 135 k=j+1,r*ni
               sum=sum+wkqnm1(i,k)*eyxyxt(j,k)
 135        continue
            if(i.le.j) then
               wkeyxyxt(i,j)=eyxyxt(i,j)-2*sum
            else
               wkeyxyxt(i,j)=eyxyxt(j,i)-2*sum
            endif
 140     continue
 150  continue
      do 300 i=1,r*ni
         do 250 j=1,r*ni
            sum=dble(0.)
            do 240 k=1,r*ni
               do 230 l=1,k
                  sum=sum+wkqnm1(i,l)*eyxyxt(l,k)*wkqnm1(j,k)
 230           continue
               do 235 l=k+1,r*ni
                  sum=sum+wkqnm1(i,l)*eyxyxt(k,l)*wkqnm1(j,k) 
 235           continue
 240        continue
            wkeyxyxt(i,j)=wkeyxyxt(i,j)+sum
 250     continue
 300  continue
      return
      end
C***********************************************************************
      function treyxyxti(r,nmax,ni,j,wkeyxyxt)
C calculates the trace of wkeyxyxt %*% (F_j Otimes I_ni), which is
C needed for calculating the first derivative wrt sigma_ii
      integer r,nmax,ni,j,l,ia
      double precision wkeyxyxt(r*nmax,r*nmax),treyxyxti,sum
      sum=dble(0.)
      do 10 l=1,ni
         ia=(j-1)*ni+l
          sum=sum+wkeyxyxt(ia,ia)
 10    continue
       treyxyxti=sum
       return
       end
C***********************************************************************
      function treyxyxtkl(r,nmax,ni,k,l,wkeyxyxt)
C calculates the trace of wkeyxyxt %*% (F_kl Otimes I_ni), which is
C needed for calculating the first derivative wrt sigma_kl
      integer r,nmax,ni,k,l,i,ia,ja
      double precision wkeyxyxt(r*nmax,r*nmax),treyxyxtkl,sum
      sum=dble(0.)
      do 10 i=1,ni
         ia=(k-1)*ni+i
         ja=(l-1)*ni+i
          sum=sum+wkeyxyxt(ia,ja)+wkeyxyxt(ja,ia)
 10    continue
       treyxyxtkl=sum
       return
       end
C***********************************************************************
      function trajaj(b,a,j,k,l,i)
C Calculates trace of A%*%J_jk%*%A%*%J_li, where A is a symmetric matrix
C and J_jk is the matrix with ones in positions (j,k) and (k,j) 
C and zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,j,k,l,i
      double precision a(b,b),trajaj
      trajaj=dble(2.)*(a(j,l)*a(k,i)+a(k,l)*a(j,i))
      return
      end
C***********************************************************************
      function trahaj(b,a,i,j,k)
C Calculates trace of A%*%H_i%*%A%*%J_jk, where A is a symmetric matrix,
C H_i is the matrix with a one in position (i,i) and zeroes elsewhere, 
C and J_jk is the matrix with  ones in positions (j,k) and (k,j) and 
C zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,i,j,k
      double precision a(b,b),trahaj
      trahaj=dble(2.)*a(i,j)*a(i,k)
      return
      end
C***********************************************************************
      function trahah(b,a,i,j)
C Calculates trace of A%*%H_i%*%A%*%H_j, where A is a symmetric matrix,
C and H_i is the matrix with a one in position (i,i) and zeroes 
C elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,i,j
      double precision a(b,b),trahah
      trahah=a(i,j)*a(i,j)
      return
      end
C***********************************************************************
      function trajajbd(b,a,c,l,ll,i,j,ii,jj)
C Calculates trace of A%*%J_jk%*%A%*%J_li, where A is a symmetric matrix
C and J_jk is the matrix with ones in positions (j,k) and (k,j) 
C and zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,c,l,ll,i,j,ii,jj,ia,iia,ja,jja
      double precision a(b,b),trajajbd
      ia=(l-1)*c+i
      iia=(l-1)*c+j
      ja=(ll-1)*c+ii
      jja=(ll-1)*c+jj
      trajajbd=dble(2.)*(a(ia,jja)*a(iia,ja)+a(iia,jja)*a(ia,ja))
      return
      end
C***********************************************************************
      function trahajbd1(b,a,c,l,ll,i,ii,jj)
C Calculates trace of A%*%H_i%*%A%*%J_jk, where A is a symmetric matrix,
C H_i is the matrix with a one in position (i,i) and zeroes elsewhere, 
C and J_jk is the matrix with  ones in positions (j,k) and (k,j) and 
C zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,c,l,ll,i,ii,jj,ia,ja,jja
      double precision a(b,b),trahajbd1
      ia=(l-1)*c+i
      ja=(ll-1)*c+ii
      jja=(ll-1)*c+jj
      trahajbd1=dble(2.)*a(ia,ja)*a(ia,jja)
      return
      end
C***********************************************************************
      function trahajbd2(b,a,c,l,ll,i,j,ii)
C Calculates trace of A%*%H_i%*%A%*%J_jk, where A is a symmetric matrix,
C H_i is the matrix with a one in position (i,i) and zeroes elsewhere, 
C and J_jk is the matrix with  ones in positions (j,k) and (k,j) and 
C zeroes elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,c,l,ll,i,j,ii,ia,iia,ja
      double precision a(b,b),trahajbd2
      ia=(l-1)*c+i
      iia=(l-1)*c+j
      ja=(ll-1)*c+ii
      trahajbd2=dble(2.)*a(ia,ja)*a(iia,ja)
      return
      end
C***********************************************************************
      function trahahbd(b,a,c,l,ll,i,ii)
C Calculates trace of A%*%H_i%*%A%*%H_j, where A is a symmetric matrix,
C and H_i is the matrix with a one in position (i,i) and zeroes 
C elsewhere.
C Note that A must be filled in above and below the diagonal.
      integer b,c,l,ll,i,ii
      double precision a(b,b),trahahbd
      trahahbd=a((l-1)*c+i,(ll-1)*c+ii)*a((l-1)*c+i,(ll-1)*c+ii)
      return
      end
C***********************************************************************
      function treps2h(ntot,r,h,st,fin,eps2)
C calculates the trace of 
C vec(y_i-X_ibeta-Z_ib_i)%*%t(vec(y_i-X_ibeta-Z_ib_i))%*%(F_k otimes I_ni),
C where F_k is the matrix with a one in position (k,k) and zeroes elsewhere
      implicit none
      integer ntot,r,h,st,fin,ii
      double precision eps2(ntot,r),sum,treps2h
      sum=dble(0.)
      do 10 ii=st,fin
         sum=sum+eps2(ii,h)*eps2(ii,h)
 10   continue
      treps2h=sum
      return
      end
C***********************************************************************
      function treps2hj(ntot,r,h,j,st,fin,eps2)
C calculates the trace of 
C vec(y_i-X_ibeta-Z_ib_i)%*%t(vec(y_i-X_ibeta-Z_ib_i))%*%(F_hj Otimes I_ni),
C where F_k is the matrix with a one in position (h,j) and (j,h) and 
C zeroes elsewhere
      implicit none
      integer ntot,r,h,j,st,fin,ii
      double precision eps2(ntot,r),sum,treps2hj
      sum=dble(0.)
      do 10 ii=st,fin
         sum=sum+dble(2.)*eps2(ii,h)*eps2(ii,j)
 10   continue
      treps2hj=sum
      return
      end
C***********************************************************************
      function truztzh(s,r,q,m,h,ztz,u)
C calculates the trace of 
C U_i %*% (F_h Otimes Z^T%*%Z), where F_h is the matrix with a one in 
C position (h,h) and zeroes elsewhere.
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,h,ii,jj,ia,ja
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truztzh
      sum=dble(0.)
      do 10 ii=1,q
         do 5 jj=1,q
            ia=(h-1)*q+ii
            ja=(h-1)*q+jj
            sum=sum+u(ia,ja,s)*ztz(jj,ii,s)
 5       continue
 10   continue
      truztzh=sum
      return
      end
C***********************************************************************
      function truztzhk(s,r,q,m,h,k,ztz,u)
C calculates the trace of 
C U_i %*% (F_hk Otimes Z^T%*%Z), where F_hk is the matrix with a one in 
C position (h,k) and (k,h) and zeroes elsewhere.
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,h,k,ii,jj,ia,ja
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truztzhk
      sum=dble(0.)
      do 10 ii=1,q
         do 5 jj=1,q
            ia=(h-1)*q+ii
            ja=(k-1)*q+jj
            sum=sum+dble(2.)*u(ia,ja,s)*ztz(jj,ii,s)
 5       continue
 10   continue
      truztzhk=sum
      return
      end
C************************************************************************
      function trhshoztzu(s,r,q,m,i,j,sigma,ztz,u)
C calculates the trace of 
C (F_i%*%Sigma%*%F_j Otimes Z^T%*%Z)%*%U_i, where F_i is the matrix 
C with a one in position (i,i) 
C Note that sigma,ztz and u must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,ii,jj,ia,ja
      double precision sigma(r,r),ztz(q,q,m),u(r*q,r*q,m),
     /      trhshoztzu,sum
      sum=dble(0.)
      do 10 ii=1,q
         do 5 jj=1,q
            ia=(i-1)*q+ii
            ja=(j-1)*q+jj
            sum=sum+sigma(i,j)*ztz(ii,jj,s)*u(ja,ia,s)
 5       continue
 10   continue
      trhshoztzu=sum
      return
      end
C************************************************************************
      function trhsjoztzu(s,r,q,m,i,j,k,sigma,ztz,u)
C calculates the trace of 
C (F_i%*%Sigma%*%F_jk Otimes Z^T%*%Z)%*%U_i, where F_i is the matrix 
C with a one in position (i,i), and F_jk is the matrix with ones in positions 
C (j,k) and (k,j) and zeroes elsewhere.
C Note that sigma,ztz and u must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,k,ii,jj,ia,ja,ka
      double precision sigma(r,r),ztz(q,q,m),u(r*q,r*q,m),
     /      trhsjoztzu,sum
      sum=dble(0.)
      do 10 ii=1,q
         do 5 jj=1,q
            ia=(i-1)*q+ii
            ja=(j-1)*q+jj
            ka=(k-1)*q+jj
            sum=sum+sigma(k,i)*ztz(ii,jj,s)*u(ja,ia,s)+
     /           sigma(j,i)*ztz(ii,jj,s)*u(ka,ia,s)
 5       continue
 10   continue
      trhsjoztzu=sum
      return
      end
C************************************************************************
      function trjsjoztzu(s,r,q,m,j,k,l,mm,sigma,ztz,u)
C calculates the trace of 
C (F_jk%*%Sigma%*%F_lmm Otimes Z^T%*%Z)%*%U_i, where F_jk is the matrix 
C with ones in positions (j,k) and (k,j) and zeroes elsewhere.
C Note that sigma,ztz and u must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,j,k,l,mm,ii,jj,ia,ja,ka
      double precision sigma(r,r),ztz(q,q,m),u(r*q,r*q,m),
     /      trjsjoztzu,sum
      sum=dble(0.)
      do 10 ii=1,q
         do 5 jj=1,q
            ia=(l-1)*q+ii
            ja=(j-1)*q+jj
            ka=(k-1)*q+jj
            sum=sum+sigma(k,mm)*ztz(ii,jj,s)*u(ja,ia,s)+
     /           sigma(j,mm)*ztz(ii,jj,s)*u(ka,ia,s)
 5       continue
 10   continue
      do 20 ii=1,q
         do 15 jj=1,q
            ia=(mm-1)*q+ii
            ja=(j-1)*q+jj
            ka=(k-1)*q+jj
            sum=sum+sigma(k,l)*ztz(ii,jj,s)*u(ja,ia,s)+
     /           sigma(j,l)*ztz(ii,jj,s)*u(ka,ia,s)
 15      continue
 20   continue
      trjsjoztzu=sum
      return
      end
C***********************************************************************
      function truztzhuztzh(s,r,q,m,i,j,ztz,u)
C calculates the trace of
C U_i%*%(F_j Otimes Z^T%*%Z)%*%U_i%*%(F_i Otimes Z^T%*%Z), where F_j is the 
C matrix with a one in (j,j) and zeroes elsewhere.
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,ii,jj,iia,jja,iu,ju,iua,jua
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truztzhuztzh
      sum=dble(0.)
      do 50 ii=1,q
         do 40 jj=1,q
            iia=(i-1)*q+ii
            jja=(j-1)*q+jj
            do 30 iu=1,q
               do 20 ju=1,q
                  iua=(j-1)*q+iu
                  jua=(i-1)*q+ju
                  sum=sum+u(iia,iua,s)*ztz(iu,jj,s)*
     /                 u(jja,jua,s)*ztz(ju,ii,s)
 20            continue
 30         continue
 40      continue
 50   continue
      truztzhuztzh=sum
      return
      end
C***********************************************************************
      function truztzhuztzj(s,r,q,m,i,j,k,ztz,u)
C calculates the trace of
C U-i%*%(F_i Otimes Z^T%*%Z)%*%U_i%*%(F_jk Otimes Z^T%*%Z), where F_i is the 
C matrix with a one in (i,i) and zeroes elsewhere, and F_jk is the matrix with
C ones in positions (j,k) and (k,j) and zeroes elsewhere.
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,k,ii,jj,iia,jja,iu,ju,iua,jua
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truztzhuztzj
      sum=dble(0.)
      do 50 ii=1,q
         do 40 jj=1,q
            iia=(i-1)*q+ii
            jja=(j-1)*q+jj
            do 30 iu=1,q
               do 20 ju=1,q
                  iua=(k-1)*q+iu
                  jua=(i-1)*q+ju
                  sum=sum+dble(2.)*u(iia,iua,s)*ztz(iu,jj,s)*
     /                 u(jja,jua,s)*ztz(ju,ii,s)
 20            continue
 30         continue
 40      continue
 50   continue
      truztzhuztzj=sum
      return
      end
C************************************************************************
      function truztzjuztzj(s,r,q,m,j,k,l,mm,ztz,u)
C calculates the trace of
C U_i%*%(F_jk Otimes Z^T%*%Z)%*%U_i%*%(F_lmm Otimes Z^T%*%Z), where F_jk is 
C the matrix with ones in positions (j,k) and (k,j) and zeroes elsewhere.
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,j,k,l,mm,ii,jj,iia,jja,iu,ju,iua,jua
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truztzjuztzj
      sum=dble(0.)
      do 50 ii=1,q
         do 40 jj=1,q
            iia=(l-1)*q+ii
            jja=(k-1)*q+jj
            do 30 iu=1,q
               do 20 ju=1,q
                  iua=(j-1)*q+iu
                  jua=(mm-1)*q+ju
                  sum=sum+dble(2.)*u(iia,iua,s)*ztz(iu,jj,s)*
     /                 u(jja,jua,s)*ztz(ju,ii,s)
 20            continue
 30         continue
 40      continue
 50   continue
      do 100 ii=1,q
         do 90 jj=1,q
            iia=(l-1)*q+ii
            jja=(j-1)*q+jj
            do 80 iu=1,q
               do 70 ju=1,q
                  iua=(k-1)*q+iu
                  jua=(mm-1)*q+ju
                  sum=sum+dble(2.)*u(iia,iua,s)*ztz(iu,jj,s)*
     /                 u(jja,jua,s)*ztz(ju,ii,s)
 70            continue
 80         continue
 90      continue
 100  continue
      truztzjuztzj=sum
      return
      end
C************************************************************************
      function truiulztz(s,r,q,m,ztz,u,b,a)
C calculates the trace of
C U_i%*%(F_j Otimes Z^T%*%Z)%*%U_i%*%G_i
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,b,a,i,j,ai,aj
      double precision ztz(q,q,m),u(r*q,r*q,m),sum,truiulztz
      sum=dble(0.)
      do 10 i=1,q
         do 5 j=1,q
            aj=(a-1)*q+j
            ai=(a-1)*q+i
            sum=sum+u(aj,b,s)*u(b,ai,s)*ztz(i,j,s)
 5       continue
 10   continue
      truiulztz=sum
      return
      end
C************************************************************************
      function truiulztzbd(s,r,q,m,ztz,u,l,i,j)
C calculates the trace of
C U_i%*%(F_j Otimes Z^T%*%Z)%*%U_i%*%G_i,l
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,l,i,j,ii,a,b,ia,ja
      double precision ztz(q,q,m),u(r*q,r*q,m),sum,truiulztzbd
      sum=dble(0.)
      ii=(l-1)*q+i
      do 10 a=1,q
         do 5 b=1,q
            ia=(j-1)*q+b
            ja=(j-1)*q+a
            sum=sum+u(ia,ii,s)*u(ii,ja,s)*ztz(a,b,s)
 5       continue
 10   continue
      truiulztzbd=sum
      return
      end
C************************************************************************
      function truiulkztz(s,r,q,m,ztz,u,i,j,k)
C calculates the trace of
C U_i%*%(F_jk Otimes Z^T%*%Z)%*%U_i%*%G_i
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,k,kk,ll,kki,kkj,kkii,kkjj
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truiulkztz
      sum=dble(0.)
      do 10 kk=1,q
         do 5 ll=1,q
            kki=(k-1)*q+ll
            kkj=(j-1)*q+kk
            kkii=(j-1)*q+ll
            kkjj=(k-1)*q+kk
            sum=sum+u(kki,i,s)*u(kkj,i,s)*ztz(kk,ll,s)+
     /           u(kkii,i,s)*u(kkjj,i,s)*ztz(kk,ll,s)
 5       continue
 10   continue
      truiulkztz=sum
      return
      end
C************************************************************************
      function truiulkztzbd(s,r,q,m,ztz,u,l,i,ii,jj)
C calculates the trace of
C U_i%*%(F_jk Otimes Z^T%*%Z)%*%U_i%*%G_i
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,l,i,ii,jj,ja,jja,a,b,ia,iia,iaa,jaa
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truiulkztzbd
      sum=dble(0.)
      jja=(l-1)*q+i
      iaa=(ii-1)*q+i
      do 10 a=1,q
         do 5 b=1,q
            ia=(ii-1)*q+a
            ja=(jj-1)*q+b
            iia=(ii-1)*q+b
            jaa=(jj-1)*q+a
            sum=sum+ztz(a,b,s)*(u(ja,jja,s)*u(jja,ia,s)+
     /           u(iia,jja,s)*u(jja,jaa,s))
 5       continue
 10   continue
      truiulkztzbd=sum
      return
      end
C************************************************************************
      function truijuztzk(s,r,q,m,ztz,u,i,j,k)
C calculates the trace of
C U_i%*%(F_k Otimes Z^T%*%Z)%*%U_i%*%G_ij
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,k,kk,ll,kki,kkj
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truijuztzk
      sum=dble(0.)
      do 10 kk=1,q
         do 5 ll=1,q
            kki=(k-1)*q+kk
            kkj=(k-1)*q+ll
            sum=sum+ztz(ll,kk,s)*(u(kki,j,s)*u(i,kkj,s)+
     /           u(j,kkj,s)*u(kki,i,s))
 5       continue
 10   continue
      truijuztzk=sum
      return
      end
C************************************************************************
      function truijuztzkbd(s,r,q,m,ztz,u,l,i,j,ii)
C calculates the trace of
C U_i%*%(F_k Otimes Z^T%*%Z)%*%U_i%*%G_ij
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,l,i,j,ii,ja,jja,a,b,ia,iia
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truijuztzkbd
      sum=dble(0.)
      ja=(l-1)*q+j
      jja=(l-1)*q+i
      do 10 a=1,q
         do 5 b=1,q
            ia=(ii-1)*q+a
            iia=(ii-1)*q+b
            sum=sum+ztz(b,a,s)*(u(ia,ja,s)*u(jja,iia,s)+
     /           u(ia,jja,s)*u(ja,iia,s))
 5       continue
 10   continue
      truijuztzkbd=sum
      return
      end
C************************************************************************
      function truijuztzlk(s,r,q,m,ztz,u,i,j,k,l)
C calculates the trace of
C U_i%*%(F_kl Otimes Z^T%*%Z)%*%U_i%*%G_ij
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,i,j,k,l,kk,ll,kki,kkj
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truijuztzlk
      sum=dble(0.)
      do 10 kk=1,q
         do 5 ll=1,q
            kki=(l-1)*q+ll
            kkj=(k-1)*q+kk
            sum=sum+dble(2.)*ztz(kk,ll,s)*(u(kki,j,s)*u(i,kkj,s)+
     /           u(j,kkj,s)*u(kki,i,s))
 5       continue
 10   continue
      truijuztzlk=sum
      return
      end
C***********************************************************************
      function truijuztzlkbd(s,r,q,m,ztz,u,l,i,j,ii,jj)
C calculates the trace of
C U_i%*%(F_kl Otimes Z^T%*%Z)%*%U_i%*%G_ij
C Note that u and ztz must be filled in above and belove the diagonal.
      implicit none
      integer s,r,q,m,l,i,j,ii,jj,ja,ia,a,b,ki,kj,kki,kkj
      double precision u(r*q,r*q,m),ztz(q,q,m),sum,truijuztzlkbd
      ja=(l-1)*q+j
      ia=(l-1)*q+i
      sum=dble(0.)
      do 10 a=1,q
         do 5 b=1,q
            ki=(ii-1)*q+a
            kj=(jj-1)*q+b
            kki=(jj-1)*q+a
            kkj=(ii-1)*q+b
            sum=sum+ztz(b,a,s)*(u(ki,ja,s)*u(ia,kj,s)+
     /           u(ki,ia,s)*u(ja,kj,s)+u(kki,ja,s)*u(ia,kkj,s)+
     /           u(kki,ia,s)*u(ja,kkj,s))
 5       continue
 10   continue
      truijuztzlkbd=sum
      return
      end
C***********************************************************************
      subroutine lltrwex(nmax,m,r,ni,s,w,eyxyxt,trwex)
C calculates trace of 
C W_i %*% E( vec(y_i - X_i beta)%*%t(vec(y_i - X_i beta)) | y_i(obs),theta )
C for subject i
      implicit none
      integer nmax,m,r,ni,s,i,j
      double precision w(r*nmax,r*nmax,m),eyxyxt(r*nmax,r*nmax),
     /     trwex,sum
      sum=dble(0.)
      do 100 i=1,r*ni
         do 50 j=1,r*ni
            if(i.le.j) then
               sum=sum+w(i,j,s)*eyxyxt(j,i)
            else
               sum=sum+w(j,i,s)*eyxyxt(j,i)
            endif
 50      continue
 100  continue
      trwex=trwex+sum
      return
      end
C************************************************************************
      subroutine trdelwdel(nmax,r,ntot,st,fin,npatt,patt,rmat,
     /     p,xcol,pcol,pdwo,pred,beta,y,wxbeta1,vdel,wo1,trdel)
C calculates t(vec(y_i - X_i beta))%*%W_i%*%vec(y_i - X_i beta)
      implicit none
      integer nmax,r,ntot,st,fin,npatt,patt(ntot),rmat(npatt,r),
     /     p,xcol(p),pcol,pdwo,posn,i,j,k
      double precision pred(ntot,pcol),beta(p,r),y(ntot,r),
     /     wxbeta1(ntot,r),vdel(r*nmax),wo1(r*nmax,r*nmax),
     /     trdel,sum
      do 2 i=1,r*nmax
         vdel(i)=dble(0.)
 2    continue
      do 100 i=st,fin
         if(patt(i).ne.0) then
            do 90 j=1,r
               sum=dble(0.)
               do 80 k=1,p
                  sum=sum+pred(i,xcol(k))*beta(k,j)
 80            continue
               wxbeta1(i,j)=sum
 90         continue
         endif
 100  continue
      posn=0
      do 150 j=1,r
         do 140 i=st,fin
            if(patt(i).ne.0) then
               if(rmat(patt(i),j).eq.1) then
                  posn=posn+1
                  vdel(posn)=y(i,j)-wxbeta1(i,j)
               endif
            endif
 140     continue
 150  continue
      do 200 i=1,pdwo
         sum=dble(0.)
         do 195 j=1,i
            sum=sum+vdel(j)*wo1(j,i)
 195     continue
         do 190 j=i+1,pdwo
            sum=sum+vdel(j)*wo1(i,j)
 190     continue
         trdel=trdel+sum*vdel(i)
 200  continue
      return
      end
C***********************************************************************
      function var1(ntot,r,q,m,s,i,ii,mcj,mck,zcol,pcol,pred,varb)
C calculates
C Cov( ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j,mc(k)) , 
C      ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j',mc(k)) )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,s,i,ii,mcj,mck,zcol(q),pcol,l,ll,ia,ja
      double precision pred(ntot,pcol),varb(r*q,r*q,m),var1,sum
      sum=dble(0.)
      do 20 l=1,q
         do 10 ll=1,q
            ia=(mcj-1)*q+l
            ja=(mck-1)*q+ll
            sum=sum+pred(i,zcol(l))*pred(ii,zcol(ll))*varb(ia,ja,s)
 10      continue
 20   continue
      var1=sum
      return
      end
C***********************************************************************
      function var2(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mcj,mck,i,
     /     zcol,pcol,wkrrb21,varb,pred)                            
C calculates
C Cov( ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j,mc(k)) , 
C ( (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j(obs) )
C  _mc(k)) )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,pt,s,noc,loc,oc(loc),npatt,mcj,mck,i,
     /     zcol(q),pcol,a,l,ll,ia,ja
      double precision wkrrb21(r,r,npatt),varb(r*q,r*q,m),
     /     pred(ntot,pcol),var2,sum
      sum=dble(0.)
      do 30 a=1,noc
         do 20 l=1,q
            do 10 ll=1,q
               ia=(mcj-1)*q+l
               ja=(oc(a)-1)*q+ll
               sum=sum+wkrrb21(mck,oc(a),pt)*pred(i,zcol(l))*
     /              pred(i,zcol(ll))*varb(ia,ja,s)
C               sum=sum+wkrrb21(oc(a),mck,pt)*pred(i,zcol(l))*
C     /              pred(ii,zcol(ll))*varb(ia,ja,s)
 10         continue
 20      continue
 30   continue
      var2=sum
      return
      end
C***********************************************************************
      function var3(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mcj,mck,i,
     /     zcol,pcol,wkrrb21,varb,pred)                            
C calculates
C Cov( ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j,mc(k)) , 
C ( (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j(obs) )
C  _mc(k)) )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,pt,s,noc,loc,oc(loc),npatt,mcj,mck,i,
     /     zcol(q),pcol,a,l,ll,ia,ja
      double precision wkrrb21(r,r,npatt),varb(r*q,r*q,m),
     /     pred(ntot,pcol),var3,sum
      sum=dble(0.)
      do 30 a=1,noc
         do 20 l=1,q
            do 10 ll=1,q
               ja=(mck-1)*q+ll
               ia=(oc(a)-1)*q+l
               sum=sum+wkrrb21(mcj,oc(a),pt)*pred(i,zcol(l))*
     /              pred(i,zcol(ll))*varb(ia,ja,s)
C               sum=sum+wkrrb21(oc(a),mck,pt)*pred(i,zcol(l))*
C     /              pred(ii,zcol(ll))*varb(ia,ja,s)
 10         continue
 20      continue
 30   continue
      var3=sum
      return
      end
C*********************************************************************
      function var22(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mcj,mck,i,j,
     /     zcol,pcol,wkrrb21,varb,pred)                            
C calculates
C Cov( ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j,mc(k)) , 
C ( (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j(obs) )
C  _mc(k)) )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,pt,s,noc,loc,oc(loc),npatt,mcj,mck,i,j,
     /     zcol(q),pcol,a,l,ll,ia,ja
      double precision wkrrb21(r,r,npatt),varb(r*q,r*q,m),
     /     pred(ntot,pcol),var22,sum
      sum=dble(0.)
      do 30 a=1,noc
         do 20 l=1,q
            do 10 ll=1,q
               ia=(mcj-1)*q+l
               ja=(oc(a)-1)*q+ll
               sum=sum+wkrrb21(mck,oc(a),pt)*pred(i,zcol(l))*
     /              pred(j,zcol(ll))*varb(ia,ja,s)
C               sum=sum+wkrrb21(oc(a),mck,pt)*pred(i,zcol(l))*
C     /              pred(ii,zcol(ll))*varb(ia,ja,s)
 10         continue
 20      continue
 30   continue
      var22=sum
      return
      end
C***********************************************************************
      function var32(ntot,r,q,m,pt,s,noc,oc,loc,npatt,mcj,mck,i,j,
     /     zcol,pcol,wkrrb21,varb,pred)                            
C calculates
C Cov( ( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_(j,mc(k)) , 
C ( (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j(obs) )
C  _mc(k)) )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,pt,s,noc,loc,oc(loc),npatt,mcj,mck,i,j,
     /     zcol(q),pcol,a,l,ll,ia,ja
      double precision wkrrb21(r,r,npatt),varb(r*q,r*q,m),
     /     pred(ntot,pcol),var32,sum
      sum=dble(0.)
      do 30 a=1,noc
         do 20 l=1,q
            do 10 ll=1,q
               ja=(mck-1)*q+ll
               ia=(oc(a)-1)*q+l
               sum=sum+wkrrb21(mcj,oc(a),pt)*pred(i,zcol(l))*
     /              pred(j,zcol(ll))*varb(ia,ja,s)
C               sum=sum+wkrrb21(oc(a),mck,pt)*pred(i,zcol(l))*
C     /              pred(ii,zcol(ll))*varb(ia,ja,s)
 10         continue
 20      continue
 30   continue
      var32=sum
      return
      end
C***********************************************************************
      function var4(ntot,r,q,m,s,pt1,pt2,noc1,noc2,loc,oc1,oc2,
     /     npatt,mcj,mck,i,ii,zcol,pcol,wkrrb21,varb,pred)
C calculates
C Cov( (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j(obs) )
C  _mc(k),
C  (Sigma_21%*%Sigma_11)( Z_i( b_i- E(b_i | y_i(obs),theta ) ) )_j'(obs) )
C  _mc(k') )
C note that varb must be filled in above and below the diagonal
      integer ntot,r,q,m,s,pt1,pt2,noc1,noc2,
     /     loc,oc1(loc),oc2(loc),
     /     npatt,mcj,mck,i,ii,zcol(q),pcol,a,aa,l,ll,ia,ja
      double precision wkrrb21(r,r,npatt),varb(r*q,r*q,m),
     /     pred(ntot,pcol),var4,sum
      sum=dble(0.)
      do 40 a=1,noc1
         do 30 aa=1,noc2
            do 20 l=1,q
               do 10 ll=1,q
                  ia=(oc1(a)-1)*q+l
                  ja=(oc2(aa)-1)*q+ll
                  sum=sum+wkrrb21(mcj,oc1(a),pt1)*
     /                 wkrrb21(mck,oc2(aa),pt2)*
     /                 pred(i,zcol(l))*pred(ii,zcol(ll))*varb(ia,ja,s)
C                  sum=sum+wkrrb21(oc1(a),mcj,pt1)*
C     /                 wkrrb21(oc2(aa),mck,pt2)*
C     /                 pred(i,zcol(l))*pred(ii,zcol(ll))*varb(ia,ja,s)
 10            continue
 20         continue
 30      continue
 40   continue
      var4=sum
      return
      end
C***********************************************************************
      subroutine obsll(ntot,m,r,nmax,p,pcol,ist,ifin,xcol,nstari,ormat,
     /     nor,mrmat,nmr,npatt,rmat,patt,y,pred,beta,wxbeta,w,wkwmm1,
     /     wkwmm2,wom,wo,wo1,wm,vdel,trdet,trdel,msg)
C calculates observed log-likelihood for given theta. Note that 
C we use workspace wkwmm1&wkwmm2 to store the inv & det of the 
C Cov(y_i(obs)) 
      implicit none
      integer ntot,m,r,nmax,p,pcol,ist(m),ifin(m),xcol(p),nstari(m),
     /     ormat(r,nmax),nor(r),mrmat(r,nmax),nmr(r),npatt,
     /     rmat(npatt,r),patt(ntot),msg,st,fin,ni,pdwm,pdwo,s,
     /     i,j,k,l,err
      double precision y(ntot,r),pred(ntot,pcol),beta(p,r),
     /     wxbeta(ntot,r),w(r*nmax,r*nmax,m),
     /     wkwmm1(r*nmax,r*nmax),wkwmm2(r*nmax,r*nmax),
     /     wom(r*nmax,r*nmax),wo(r*nmax,r*nmax),
     /     wo1(r*nmax,r*nmax),wm(r*nmax,r*nmax),
     /     vdel(r*nmax),sum,trdel,trdet
      trdel=dble(0.)
      trdet=dble(0.)
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
         ni=nstari(s)
         do 2 i=1,r*nmax
            do 1 j=i,r*nmax
               wkwmm1(i,j)=dble(0.)
               wkwmm2(i,j)=dble(0.)
               wo1(i,j)=dble(0.)
 1          continue
 2       continue
         call getormat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nor,ormat)
         call getmrmat(ntot,r,nmax,npatt,patt,
     /     rmat,st,fin,nmr,mrmat)
         call mkwo(m,r,nmax,s,ormat,nor,ni,pdwo,w,wo)
         call mkwom(m,r,nmax,s,mrmat,ormat,nmr,nor,
     /        ni,pdwo,pdwm,w,wom)
         call mkwm(m,r,nmax,s,mrmat,nmr,ni,pdwm,w,wm)
         do 4 i=1,pdwm
            do 3 j=i,pdwm
               wkwmm1(i,j)=wm(i,j)
 3          continue
 4       continue
         call chfce(r*nmax,pdwm,wkwmm1,err)
         if(err.eq.1) then 
            msg=99
            goto 999
         endif
         call bkslv(r*nmax,pdwm,wkwmm1)
         call mm(r*nmax,pdwm,wkwmm1,wkwmm2)
         do 100 k=1,pdwo
            do 90 l=k,pdwo
              sum=dble(0.)
              do 80 j=1,pdwm
                 do 70 i=1,j
                    sum=sum+wom(k,i)*wkwmm2(i,j)*wom(l,j)
 70              continue
                 do 75 i=j+1,pdwm
                    sum=sum+wom(k,i)*wkwmm2(j,i)*wom(l,j)
 75              continue
 80           continue
              wo1(k,l)=wo(k,l)-sum
 90        continue
 100    continue
        call trdelwdel(nmax,r,ntot,st,fin,npatt,patt,rmat,
     /     p,xcol,pcol,pdwo,pred,beta,y,wxbeta,vdel,wo1,trdel)
        call chfce(r*nmax,pdwo,wo1,err)
        if(err.eq.1) then
           msg=90
           goto 999
        endif
        do 110 i=1,pdwo
           trdet=trdet+dlog(wo1(i,i))
 110    continue
 500  continue
 999  continue
      return
      end
C************************************************************************


         
         
