C***********************************************************************
      subroutine istfin2(ntot,subj,m,ist,ifin)
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
C*********************************************************************
       subroutine mlem(ntot,m,r,p,subj,ist,ifin,nmax,iposn,npatt,pstfin,
     /     patt,nstar,nstari,rmat,pcol,xcol,pred,wxbeta,wkrrpt,y,ey,
     /     eyyt,eyxyxt,iter,msg,sigma,beta,xtx,xtw,xtwx,xtwy,xtwxinv,
     /     wkrr1,wkrr2,wkeyxyxt,wkqnm1,cvgd,obeta,osigma,
     /     bigm,maxits,llovec,epsil,sflag,epsi,wkpr,wkpp,xtxinv)
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
      integer ntot,m,r,p,subj(ntot),ist(m),ifin(m),nmax,iposn(ntot),
     /     npatt,pstfin(npatt,2),patt(ntot),nstar,nstari(m),
     /     rmat(npatt,r),pcol,xcol(p),err,iter,msg,cvgd,maxits,
     /     c1,c2,i,j,oc(100),oc2(100),mc(100),mc1(100),
     /     sflag
      double precision pred(ntot,pcol),wxbeta(ntot,r),
     /     wkrrpt(r,r,npatt),y(ntot,r),ey(ntot,r),eyyt(r*nmax,r*nmax),
     /     eyxyxt(r*nmax,r*nmax),sigma(r,r),
     /     beta(p,r),xtx(p,p,m),xtw(p*r,nmax*r),xtwx(p*r,p*r),
     /     xtwy(p*r),xtwxinv(p*r,p*r),wkrr1(r,r),wkrr2(r,r),
     /     wkqnm1(r*nmax,r*nmax),wkeyxyxt(r*nmax,r*nmax),obeta(p,r),
     /     bigm(r,r),osigma(r,r),epsil,ll,epsi,
     /     llovec(maxits),eps(ntot,r),wkpr(p,r),wkpp(p,p),xtxinv(p,p)
      msg=0
      iter=0
      call prefem2(ntot,subj,m,ist,ifin,pcol,pred,p,
     /     xcol,patt,nstar,nstari,xtx,wkpp,xtxinv,err)
      if(err.eq.1) then
         msg=1
         goto 999
      endif
      if(sflag.ne.1) then
         call mimpy(ntot,r,y,patt,npatt,rmat)
         call mkxty(ntot,r,y,pcol,pred,p,xcol,patt,wkpr)
         call mkbeta(p,r,xtxinv,wkpr,beta)
         call mkeps1(ntot,r,y,pcol,pred,p,xcol,beta,eps,patt)
         call mksigma(ntot,r,eps,nstar,sigma,patt)
      endif
      cvgd=0
C****************START OF MAIN ITERATION ****************************
 1    continue
      iter=iter+1
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
C******* caculate the quantities required for E-step
      call mkxbeta(ntot,m,ist,ifin,p,r,pcol,xcol,patt,pred,beta,wxbeta)
      call mkey2(100,100,oc,mc,m,r,ntot,nstari,iposn,npatt,pstfin,
     /     rmat,patt,p,xcol,pcol,ist,ifin,pred,y,beta,sigma,wkrr1,
     /     wxbeta,ll,bigm,wkrrpt,ey,eps)
      llovec(iter)=ll
C      call mkxbeta(ntot,m,ist,ifin,p,r,pcol,xcol,patt,pred,beta,wxbeta)
      call sigmaem2(ntot,nmax,m,r,pcol,ist,ifin,nstari,100,100,mc,mc1,
     /     oc,oc2,nstar,npatt,patt,rmat,pred,wxbeta,y,ey,eyyt,
     /     wkrrpt,sigma,err)
      if(err.eq.1) then
         msg=3
         goto 999
      endif
      call gls2(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,nstari,patt,pred,
     /     sigma,wkrr1,wkrr2,ey,beta,xtx,xtw,xtwx,xtwy,xtwxinv,err)
      if(err.eq.1) then
         msg=4
         goto 999
      endif
C********  CHECK CONVERGENCE  *************************
      c1=0
      do 30 i=1,p
         do 27 j=1,r
            if(dabs(beta(i,j)-obeta(i,j)).gt.(epsil*
     /           dabs(obeta(i,j)))) c1=1
 27      continue
 30   continue
      c2=0
      do 50 i=1,r
         do 45 j=i,r
            if(dabs(sigma(i,j)-osigma(i,j)).gt.(epsil*
     /           dabs(osigma(i,j)))) c2=1
 45      continue
 50   continue
      if((c1.eq.0).and.(c2.eq.0)) cvgd=1
      if((cvgd.eq.0).and.(iter.lt.maxits)) goto 1
C********* end of main iteration *****************
 999  continue
      iter=iter
      msg=msg
      return
      end
C***********************************************************************
      subroutine prefem2(ntot,subj,m,ist,ifin,pcol,pred,p,
     /     xcol,patt,nstar,nstari,xtx,wkpp,xtxinv,err)
C Preliminary manipulations for pan. After execution, 
C     ist  = vector of starting positions for subject i, i=1,...,m
C     ifin = vector of finishing positions for subject i, i=1,...,m
C     xtx  = t(x_i)%*%x, i=1,...,m
C     nstar = total number of rows in y containing data
C     nstari = total number of rows in y_i containing data
C     wkpp = t(pred)%*%pred
C     xtxinv = inv(t(pred)%*%pred)
      implicit none
      integer ntot,subj(ntot),m,ist(m),ifin(m),pcol,
     /     p,xcol(p),patt(ntot),nstar,nstari(m),err,
     /     st,fin,s,i,j,k
      double precision pred(ntot,pcol),xtx(p,p,m),
     /     wkpp(p,p),xtxinv(p,p),sum
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
         do 90 i=1,p
            do 80 j=i,p
               sum=dble(0.)
               do 60 k=ist(s),ifin(s)
                  if(patt(k).ne.0) then
                     sum=sum+pred(k,xcol(i))*pred(k,xcol(j))
                  endif
 60            continue
               xtx(i,j,s)=sum
               if(i.ne.j) xtx(j,i,s)=sum
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
C*********************************************************************
      subroutine mimpy2(ntot,r,y,patt,npatt,rmat)
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
      subroutine mkbeta2(p,r,xtxinv,xty,beta)
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
      subroutine mkxty2(ntot,r,y,pcol,pred,p,xcol,patt,xty)
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
      subroutine mkeps12(ntot,r,y,pcol,pred,p,xcol,beta,epsi,patt)
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
      subroutine mksigma2(ntot,r,epsi,nstar,sigma,patt)
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
      subroutine mkyyt2(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,
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
      subroutine mkey2(loc,lmc,oc,mc,m,r,ntot,nstari,iposn,npatt,pstfin,
     /     rmat,patt,p,xcol,pcol,ist,ifin,pred,y,beta,sigma,wkrr1,
     /     wxbeta,ll,bigm,wkrrpt,ey,epsi)
      implicit none
      integer loc,lmc,oc(loc),mc(lmc),m,r,ntot,nstari(m),iposn(ntot),
     /     npatt,pstfin(npatt,2),rmat(npatt,r),patt(ntot),
     /     p,xcol(p),pcol,ist(m),ifin(m),
     /     i,j,pt,nmc,noc,k,l,s,st,fin,ni,gi
      double precision pred(ntot,pcol),y(ntot,r),beta(p,r),
     /     sigma(r,r),wkrr1(r,r),wkrrpt(r,r,npatt),wxbeta(ntot,r),ll,
     /     bigm(r,r),epsi(ntot,r),ey(ntot,r),sum,t,d
C below epsi is like ystar
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
                  if(rmat(patt(i),j).eq.1) then
                     epsi(i,j)=y(i,j)-sum
                  else
                     epsi(i,j)=y(i,j)
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
 85      continue
         do 100 i=pstfin(pt,1),pstfin(pt,2)
            do 95 k=1,nmc
               sum=dble(0.)
               do 90 l=1,noc
                  sum=sum+wkrr1(oc(l),mc(k))*epsi(iposn(i),oc(l))
 90            continue
               epsi(iposn(i),mc(k))=sum
 95         continue
 100     continue
 110  continue
      d=dble(0.)
      ll=dble(0.)
      do 200 pt=1,npatt
         do 190 j=1,r
            if((rmat(pt,j).eq.1).and.(sigma(j,j).gt.0))then
               d=d+dlog(sigma(j,j))
               call swp(r,sigma,j)
            elseif((rmat(pt,j).eq.0).and.(sigma(j,j).lt.0))then
               call rsw(r,sigma,j)
               d=d-dlog(sigma(j,j))
            endif
 190     continue
         do 195 i=1,r
            do 194 j=i,r
               bigm(i,j)=dble(0.)
 194        continue
 195     continue
         call getoc(r,npatt,rmat,pt,loc,oc,noc)
         do 185 i=pstfin(pt,1),pstfin(pt,2)
            do 180 j=1,noc
               do 170 k=j,noc
                  bigm(oc(j),oc(k))=bigm(oc(j),oc(k))+(y(iposn(i),
     /                 oc(j))-wxbeta(iposn(i),oc(j)))*(y(iposn(i),
     /                 oc(k))-wxbeta(iposn(i),oc(k)))
 170           continue
 180        continue
 185     continue
         t=dble(0.)
         do 160 j=1,noc
            do 150 k=1,noc
               t=t-sigma(oc(j),oc(k))*bigm(oc(j),oc(k))
C               ll=ll-((pstfin(pt,2)-pstfin(pt,1)+1)*d+t)/2
 150        continue
 160     continue
         ll=ll-((pstfin(pt,2)-pstfin(pt,1)+1)*d+t)/2
 200  continue
C now calculate E(y_i/y_i(obs),theta) 
      do 500 s=1,m
         st=ist(s)
         fin=ifin(s)
         do 450 i=st,fin
            gi=0
            if(patt(i).ne.0) then
               do 400 j=1,r
                  sum=dble(0.)
                  do 350 k=1,p
                     sum=sum+pred(i,xcol(k))*beta(k,j)
 350              continue
                  ey(i,j)=sum+epsi(i,j)
 400           continue
            endif
 450     continue
 500  continue
      return
      end
C*********************************************************************
      subroutine mkeyyt2(ntot,nmax,npatt,m,r,st,fin,ni,patt,
     /     rmat,s,pcol,lmc,loc,mc,mc1,oc,oc2,pred,wkrrpt,ey,eyyt)
C this subroutine fills only the missing portions of the 
C E(vec(y_i) %*%vec(y_i)^T  | y_i(obs),theta )
      implicit none
      integer ntot,nmax,npatt,m,r,st,fin,ni,patt(ntot),
     /     rmat(npatt,r),s,pcol,lmc,loc,
     /     mc(lmc),mc1(lmc),oc(loc),nmc,nmc1,noc,i,j,ri,rj,k,
     /     jj,ii,pt,pt1,pt2,gi,gj,noc2,oc2(loc)
      double precision pred(ntot,pcol),wkrrpt(r,r,npatt),ey(ntot,r),
     /     eyyt(r*nmax,r*nmax)
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
     /                 ey(i,mc(j))*ey(i,mc(k))
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
                           eyyt(ri,rj)=ey(i,mc1(ii))*ey(j,mc(k))
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
C*********************************************************************
      subroutine sigmaem2(ntot,nmax,m,r,pcol,ist,ifin,nstari,lmc,loc,mc,
     /     mc1,oc,oc2,nstar,npatt,patt,rmat,pred,wxbeta,
     /     y,ey,eyyt,wkrrpt,emsigma,err)
C replaces sigma with its EM estimate
      implicit none
      integer ntot,nmax,m,r,pcol,ist(m),ifin(m),nstari(m),
     /     nstar,npatt,patt(ntot),rmat(npatt,r),
     /     err,lmc,loc,mc(lmc),mc1(lmc),oc(loc),oc2(loc),
     /     st,fin,ni,i,j,s,ii,gi,jj,ia
      double precision pred(ntot,pcol),wxbeta(ntot,r),y(ntot,r),
     /     ey(ntot,r),eyyt(r*nmax,r*nmax),wkrrpt(r,r,npatt),
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
               call mkyyt(ntot,nmax,r,st,fin,ni,patt,npatt,rmat,y,eyyt)
               call mkeyyt2(ntot,nmax,npatt,m,r,st,fin,ni,patt,rmat,s,
     /              pcol,lmc,loc,mc,mc1,oc,oc2,pred,wkrrpt,ey,eyyt)
               gi=0
               do 390 ii=st,fin
                  if(patt(ii).ne.0) then
                     gi=gi+1
                     sum=sum+eyyt((i-1)*ni+gi,(j-1)*ni+gi)
                  endif
 390           continue
               do 365 jj=st,fin
                  if(patt(jj).ne.0) then
                     sum2=sum2+wxbeta(jj,j)*ey(jj,i)
                  endif
 365           continue
               do 355 jj=st,fin
                  if(patt(jj).ne.0) then
                     sum2=sum2+wxbeta(jj,i)*ey(jj,j)
                  endif
 355           continue
               do 300 ia=st,fin
                  if(patt(ia).ne.0) then
                     sum4=sum4+wxbeta(ia,i)*wxbeta(ia,j)
                  endif
 300           continue
 400        continue
            emsigma(i,j)=(sum-sum1-sum2+sum4)/dfloat(nstar)
            if(i.ne.j) emsigma(j,i)=emsigma(i,j)
 450     continue
 500  continue
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
      subroutine getmc2(r,npatt,rmat,pt,lmc,mc,nmc)
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
      subroutine getoc2(r,npatt,rmat,pt,loc,oc,noc)
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
      subroutine swpobs2(r,sigma,npatt,rmat,pt)
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
      subroutine swp2(n,mat,p)
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
      subroutine rsw2(n,mat,p)
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
C***********************************************************************
      subroutine chfce2(p,pw,s,err)
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
      subroutine bkslv2(p,pw,s)
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
      subroutine mm2(p,pw,wm,cm)
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
C*********************************************************************
      subroutine gls2(ntot,m,r,ist,ifin,nmax,pcol,p,xcol,
     /     nstari,patt,pred,sigma,wkrr1,wkrr2,ey,beta,xtx,xtw,xtwx,
     /     xtwey,xtwxinv,err)
C calculates gls estimate of beta, using weights in 
C w= inv(sigma Otimes I_(n_i))
C for the model without random effects
      implicit none
      integer ntot,m,r,ist(m),ifin(m),nmax,pcol,p,
     /     xcol(p),nstari(m),patt(ntot),ni,st,fin,s,
     /     i,j,err,ia,ja,msg
      double precision pred(ntot,pcol),sigma(r,r),wkrr1(r,r),wkrr2(r,r),
     /     ey(ntot,r),beta(p,r),xtx(p,p,m),xtw(p*r,nmax*r),
     /     xtwx(p*r,p*r),xtwey(p*r),xtwxinv(p*r,p*r),sum
      err=0
C initialize ( inv(sigma) Otimes X_i^t ) E(vec(y_i) |yobs, theta) and
C ( inv(sigma) Otimes X_i^t%*%X_i )
      do 10 i=1,r*p
         xtwey(i)=dble(0.)
         do 5 j=i,r*p
            xtwx(i,j)=dble(0.)
 5       continue
 10   continue
      do 20 i=1,r
         do 15 j=i,r
            wkrr2(i,j)=sigma(i,j)
 15      continue
 20   continue
      call chfce(r,r,wkrr2,err)
      if(err.eq.1) then
         msg=6
         goto 999
      endif
      call bkslv(r,r,wkrr2)
      call mm(r,r,wkrr2,wkrr1)
      do 100 s=1,m
         ni=nstari(s)
         st=ist(s)
         fin=ifin(s)
         call mkxtw2(ntot,r,p,m,pcol,xcol,patt,ni,st,fin,nmax,wkrr1,
     /        pred,xtw,s)
         call mkxtwx2(ntot,m,r,p,pcol,xcol,st,fin,patt,s,
     /        ni,nmax,pred,xtx,wkrr1,xtw,xtwx)
         call mkxtwey2(ntot,r,p,st,fin,nmax,ni,patt,
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
      subroutine mkxtw2(ntot,r,p,m,pcol,xcol,patt,ni,st,fin,nmax,wkrr1,
     /     pred,xtw,s)
C calculates t(inv(sigma) Otimes X_i) for subject s 
      implicit none
      integer ntot,r,p,m,pcol,xcol(p),patt(ntot),ni,st,fin,
     /     nmax,s,i,j,k,ia,ja,l,indi
      double precision pred(ntot,pcol),wkrr1(r,r),xtw(p*r,nmax*r)
      do 10 i=1,r
         do 5 j=i+1,r
            wkrr1(j,i)=wkrr1(i,j)
 5       continue
 10   continue
      do 500 i=1,r
         do 450 j=1,r
            do 400 k=1,p
               ia=(i-1)*p+k           
               indi=0
               do 300 l=st,fin
                  if(patt(l).ne.0) then
                     indi=indi+1
                     ja=(j-1)*ni+indi
                     xtw(ia,ja)=wkrr1(i,j)*pred(l,xcol(k))
                  endif
 300           continue
 400        continue
 450     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkxtwx2(ntot,m,r,p,pcol,xcol,st,fin,patt,
     /     s,ni,nmax,pred,xtx,wkrr1,xtw,xtwx)
C increments inv(sigma) Otimes t(X_i) %*% X_i  for subject s
      implicit none
      integer ntot,m,r,p,pcol,xcol(p),st,fin,patt(ntot),s,ni,
     /     nmax,i,j,k,l,ii,jj
      double precision pred(ntot,pcol),xtx(p,p,m),wkrr1(r,r),
     /     xtw(p*r,nmax*r),xtwx(p*r,p*r),tmp
      do 500 i=1,r
         do 400 j=i,r
            do 200 k=1,p
               do 100 l=1,p
                  ii=(i-1)*p+k
                  jj=(j-1)*p+l
                  tmp=wkrr1(i,j)*xtx(k,l,s)
                  xtwx(ii,jj)=xtwx(ii,jj)+tmp
 100           continue
 200        continue
 400     continue
 500  continue
      return
      end
C***********************************************************************
      subroutine mkxtwey2(ntot,r,p,st,fin,nmax,ni,
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
