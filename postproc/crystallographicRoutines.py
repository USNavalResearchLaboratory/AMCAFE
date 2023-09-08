#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 15:03:18 2018

@author: kteferra
"""

import numpy as np

def EulerAnglesToQuats(EulerAngles):
    '''
    Euler angles (varphi1, Phi,varphi2) must be in radians and in Bunge Euler
    description (zxz). They are converted to quats using
    the approach in appendix A of:
    Rowenhorst, David, et al. "Consistent representations of and conversions
    between 3D rotations.    " Modelling and Simulation in Materials Science 
    and Engineering 23.8 (2015): 083501.
    '''
    P=1.0
    if (int(EulerAngles.size/3)-np.finfo(float).eps >1):
        n = EulerAngles.shape[0]
        sigma = .5*(EulerAngles[:,0]+EulerAngles[:,2])
        delta = .5*(EulerAngles[:,0]-EulerAngles[:,2])
        c = np.cos(EulerAngles[:,1]/2.0)
        s = np.sin(EulerAngles[:,1]/2.0)
        quats = np.concatenate((c*np.cos(sigma),-P*s*np.cos(delta), \
                                -P*s*np.sin(delta),-P*c*np.sin(sigma))).\
                                reshape([4,n]).transpose()
        z=np.where(quats[:,0]<0)[0]
        quats[z,:] = -quats[z,:]                          
    else:
        sigma = .5*(EulerAngles[0]+EulerAngles[2])
        delta = .5*(EulerAngles[0]-EulerAngles[2])
        c = np.cos(EulerAngles[1]/2.0)
        s = np.sin(EulerAngles[1]/2.0)
        quats = np.array([c*np.cos(sigma),-P*s*np.cos(delta), \
                                -P*s*np.sin(delta),-P*c*np.sin(sigma)])
        if (quats[0]<0): quats = -quats                          
        
    return quats


def AxisAngleToQuat(axis,angle):
    '''
    angle is in radians
    axis is unit norm
    '''
    a=np.shape(axis)
    if (len(a)==1):
        quat = np.append(np.cos(angle/2.0),axis*np.sin(angle/2.0))
    else:
        N = a[0]
        quat = np.zeros([N,4])
        for j in range(0,N):
            quat[j,:] = np.append(np.cos(angle[j]/2.0),
                axis[j,:]*np.sin(angle[j]/2.0))            
    return quat

def QuatToAxisAngle(quat):
    a=np.shape(quat)
    if (len(a)==1):
        angle = 2*np.arccos(quat[0])
        s = 1/np.linalg.norm(quat[1:])
        axis = quat[1:]*s
    else:
        N=a[0]
        angle=np.zeros(N)
        axis = np.zeros([N,3])
        for j in range(0,N):
            angle[j] = 2*np.arccos(quat[j,0])
            s = 1.0/np.linalg.norm(quat[j,1:])
            axis[j,:] = quat[j,1:]*s
    return (axis,angle)

def QuatToEulerAngle(quat):
    
    '''
    converts quats to Euler angles in ZXZ convention. Note that
    the returned Euler Angles are active rotations (meaning if plug into
    rotation matrix will get active rotation)
    note that if passive rotation eu are a,b,c then active is -c,-b,-a 
    '''
    a=np.shape(quat)
    if (len(a)==1):
        q0=quat[0]
        q1=quat[1]
        q2=quat[2]
        q3=quat[3]
        q03=q0**2+q3**2
        q12=q1**2+q2**2
        chi=np.sqrt( q03*q12)
        eulerAngles=np.zeros(3)
        eulerAngles[0]=-np.arctan2((q0*q2+q1*q3)/chi,(q2*q3-q1*q0)/chi)
        eulerAngles[1]=-np.arctan2(2*chi,q03-q12)
        eulerAngles[2]=-np.arctan2((q1*q3-q0*q2)/chi,(-q0*q1-q2*q3)/chi)
    else:
        N=a[0]
        eulerAngles=np.zeros([N,3])
        q0=quat[:,0]
        q1=quat[:,1]
        q2=quat[:,2]
        q3=quat[:,3]
        q03=q0**2+q3**2
        q12=q1**2+q2**2
        chi=np.sqrt( q03*q12)
        eulerAngles[:,0]=-np.arctan2((q0*q2+q1*q3)/chi,(q2*q3-q1*q0)/chi)
        eulerAngles[:,1]=-np.arctan2(2*chi,q03-q12)
        eulerAngles[:,2]=-np.arctan2((q1*q3-q0*q2)/chi,(-q0*q1-q2*q3)/chi)
    return eulerAngles        

def AxisAngleToEulerAngle(axis,ang):
    q=AxisAngleToQuat(axis,ang)
    eu=QuatToEulerAngle(q)
    return eu


def QuatMultiply(quat1,quat2):
    ''' see Eq. (23) in 
    Rowenhorst, David, et al. "Consistent representations of and conversions
    between 3D rotations.    " Modelling and Simulation in Materials Science 
    and Engineering 23.8 (2015): 083501.
    '''
    P=1.0
    qprod=np.append(quat1[0]*quat2[0]-np.dot(quat1[1:],quat2[1:]), \
                          quat2[0]*quat1[1:]+quat1[0]*quat2[1:]+ \
                          P*np.cross(quat1[1:],quat2[1:]))
    return qprod

def QuatRotate(quat1,r):
    '''
    See Eq. (24) in 
    Rowenhorst, David, et al. "Consistent representations of and conversions
    between 3D rotations.    " Modelling and Simulation in Materials Science 
    and Engineering 23.8 (2015): 083501.

    this rotates vector r by quat1. Can do a passive rotation by using the 
    conjugate of quat1. Note that quat1 should be a unit quaternion
    '''
    qDim=quat1.shape
    rDim=r.shape
    P=1.0
    if (len(qDim)==1 and len(rDim)==1):
        r2 = (quat1[0]**2-np.linalg.norm(quat1[1:])**2)*r + \
            2*np.dot(quat1[1:],r)*quat1[1:] + \
            2*P*quat1[0]*np.cross(quat1[1:],r)
    if (len(qDim)>1 and len(rDim)==1):
        N = qDim[0]
        r2=np.zeros([N,3])
        b=np.linalg.norm(quat1[:,1:],axis=1)**2
        v1=np.matmul(np.reshape(quat1[:,0]**2-b,[N,1]),
            np.reshape(r,[1,3]))
        b=2*np.matmul(quat1[:,1:],np.reshape(r,[3,1]))
        v2 = b*quat1[:,1:]
        v3 = 2*P*np.reshape(quat1[:,0],[N,1])*np.cross(quat1[:,1:],r)
        r2 = v1+v2+v3
    if (len(qDim)==1 and len(rDim)>1):        
        N = rDim[0]
        r2=np.zeros([N,3])
        b=quat1[0]**2-np.linalg.norm(quat1[1:])**2
        v1=b*r
        b=np.reshape(np.tile(quat1[1:],N),[N,3])
        b2=2*np.reshape(np.matmul(r,quat1[1:]),[N,1])
        v2=b2*b
        v3=2*P*np.reshape(np.tile(quat1[0],N),[N,1])*np.cross(quat1[1:],r)        
        r2=v1+v2+v3
    if (len(qDim)>1 and len(rDim)>1):                
        Nq=qDim[0]
        Nr=rDim[0]
        Nt=Nq*Nr
        r2=np.zeros([Nt,3])
        qv=np.tile(quat1,(Nr,1))
        rv=np.repeat(r,Nq,axis=0)
        b=qv[:,0]**2 - np.linalg.norm(qv[:,1:],axis=1)**2
        v1=np.reshape(b,[Nt,1])*rv
        v2=2*np.reshape(np.sum(qv[:,1:]*rv,axis=1),[Nt,1])*qv[:,1:]
        v3=2*P*np.reshape(qv[:,0],[Nt,1])*np.cross(qv[:,1:],rv)
        r2=v1+v2+v3
    '''        
        for j in range(0,N):
            r2[j,:] = (quat1[j,0]**2-np.linalg.norm(quat1[j,1:])**2)*r + \
                2*np.dot(quat1[j,1:],r)*quat1[j,1:] + \
                2*P*quat1[j,0]*np.cross(quat1[j,1:],r)
    '''
    return r2



def QuatConjugate(quat1):
    a=quat1.shape
    if (len(a)==1):
        quatconj = np.append(quat1[0],-quat1[1:])
    else:
        N = a[0]
        quatconj=np.zeros([N,4])
        for j in range(0,N):
            quatconj[j,:] = np.append(quat1[j,0],-quat1[j,1:])
    return quatconj


def EulerAngleXZXToRmat(phi1,psi,phi2):
    ' constructs rotation matrix from Euler angles in XZX convention'
    rRot=np.zeros([3,3])
    rRot[0,0] = np.cos(psi)
    rRot[0,1] = -np.cos(phi2)*np.sin(psi)
    rRot[0,2] = np.sin(psi)*np.sin(phi2)
    rRot[1,0] = np.cos(phi1)*np.sin(psi)
    rRot[1,1] = np.cos(phi1)*np.cos(psi)*np.cos(phi2) \
        - np.sin(phi1)*np.sin(phi2)
    rRot[1,2] = -np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*np.cos(psi)*np.sin(phi2)
    rRot[2,0] = np.sin(phi1)*np.sin(psi)
    rRot[2,1] = np.cos(phi1)*np.sin(phi2)+np.cos(psi)*np.cos(phi2)*np.sin(phi1)
    rRot[2,2] = np.cos(phi1)*np.cos(phi2) - np.cos(psi)*np.sin(phi1)*np.sin(phi2)
    return rRot

def EulerAngleZYZToRmat(phi1,psi,phi2):
    rRot=np.zeros([3,3])
    rRot[0,0] = np.cos(phi1)*np.cos(psi)*np.cos(phi2) - \
        np.sin(phi1)*np.sin(phi2)
    rRot[0,1] = -np.cos(phi2)*np.sin(phi1) - np.cos(phi1)*np.cos(psi)*np.sin(phi2)
    rRot[0,2] = np.cos(phi1)*np.sin(psi)
    rRot[1,0] = np.cos(phi1)*np.sin(phi2)+np.cos(psi)*np.cos(phi2)*np.sin(phi1)
    rRot[1,1] = np.cos(phi1)*np.cos(phi2)-np.cos(psi)*np.sin(phi1)*np.sin(phi2)
    rRot[1,2] = np.sin(phi1)*np.sin(psi)
    rRot[2,0] = -np.cos(phi2)*np.sin(psi)
    rRot[2,1] = np.sin(psi)*np.sin(phi2)
    rRot[2,2] = np.cos(psi)
    return rRot

def EulerAngleZXZToRmat(phi1,psi,phi2):
    # active rotation
    rRot=np.zeros([3,3])
    c1=np.cos(phi1)
    c2=np.cos(psi)
    c3=np.cos(phi2)
    s1=np.sin(phi1)
    s2=np.sin(psi)
    s3=np.sin(phi2)
    rRot[0,0] = c1*c3-c2*s1*s3
    rRot[0,1] = -c1*s3 -c2*c3*s1
    rRot[0,2] = s1*s2
    rRot[1,0] = c3*s1+c1*c2*s3
    rRot[1,1] = c1*c2*c3 - s1*s3
    rRot[1,2] = -c1*s2
    rRot[2,0] = s2*s3
    rRot[2,1] = c3*s2
    rRot[2,2] = c2
    return rRot


def AxisAngleToRmat(ax,ang):
    rRot =np.zeros([3,3])
    rRot[0][0] = np.cos(ang) + (ax[0]**2.0)*(1-np.cos(ang))
    rRot[0][1] = ax[0]*ax[1]*(1-np.cos(ang)) - ax[2]*np.sin(ang)
    rRot[0][2] = ax[0]*ax[2]*(1-np.cos(ang)) + ax[1]*np.sin(ang)
    rRot[1][0] = ax[0]*ax[1]*(1-np.cos(ang)) + ax[2]*np.sin(ang)
    rRot[1][1] = np.cos(ang) + (ax[1]**2.0)*(1-np.cos(ang))
    rRot[1][2] = ax[1]*ax[2]*(1-np.cos(ang)) - ax[0]*np.sin(ang)
    rRot[2][0] = ax[2]*ax[0]*(1-np.cos(ang)) - ax[1]*np.sin(ang)
    rRot[2][1] = ax[2]*ax[1]*(1-np.cos(ang)) + ax[0]*np.sin(ang)
    rRot[2][2] = np.cos(ang) + (ax[2]**2.0)*(1-np.cos(ang))
    return rRot    

def QuatToRmat(quat):
    q0=quat[0]
    q1=quat[1]
    q2=quat[2]
    q3=quat[3]
    rRot =np.zeros([3,3])
    rRot[0,0] = 2 * (q0 * q0 + q1 * q1) - 1
    rRot[0,1] = 2 * (q1 * q2 - q0 * q3)
    rRot[0,2] = 2 * (q1 * q3 + q0 * q2)
    rRot[1,0] = 2 * (q1 * q2 + q0 * q3)
    rRot[1,1] = 2 * (q0 * q0 + q2 * q2) - 1
    rRot[1,2] = 2 * (q2 * q3 - q0 * q1)
    rRot[2,0] = 2 * (q1 * q3 - q0 * q2)
    rRot[2,1] = 2 * (q2 * q3 + q0 * q1)
    rRot[2,2] = 2 * (q0 * q0 + q3 * q3) - 1
    return rRot
    


def WignerSmallD(l,m,n,beta):
    '''
    This computes the wigner small "d_l^{mn}" based on 
    Hielscher, Ralf. "Kernel density estimation on the
    rotation group and its application to crystallographic
    texture analysis." Journal of Multivariate Analysis
    119 (2013): 119-143.
    '''
    from scipy.special import factorial
    from scipy.special import comb    
    epsm=np.finfo(float).eps
    x=np.cos(beta)
    if (len(x.shape)==0):
        x=np.reshape(x,1)
    
    i1 = np.where(np.abs(x-1.0)<2*epsm)[0]
    i2 = np.where(np.abs(x+1.0)<2*epsm)[0] 
    ind1=np.concatenate((i1,i2))
    if (len(ind1)>0):
        x[ind1]=0.0
    coeff=(-1)**(l-m)/(2**l)*np.sqrt(factorial(l+m)/ \
               (factorial(l-m)*factorial(l+n)*factorial(l-n)) )* \
               np.sqrt( (1-x)**(n-m)/(1+x)**(n+m) )

    nd = l - m
    coeff1 = np.zeros(nd+1)
    for j in range(0,(nd+1)):
        coeff1[j] = comb(nd,j)
    i1 = np.zeros(nd+1)
    i2 = np.zeros(nd+1)
    for j in range(0,(nd+1)):
        i1[j]=nd-j
        i2[j]=j
        
    coeff2=0.0
    for j in range(0,(nd+1)):
        y1=1       
        a=l-n        
        for j1 in range(0,int(i1[j])):
            y1 *= (-1)*(a-j1)
        y1*=(1-x)**(a-i1[j])
        y2=1
        a=l+n
        for j1 in range(0,int(i2[j])):
            y2*=(a-j1)
        y2*=(1+x)**(a-i2[j])
        coeff2+= coeff1[j]*y1*y2
    ld=coeff*coeff2
    if (len(ind1)>0):
        ld[ind1]=0.0
    return ld

def legendreA(l,m,beta):
    '''
    This computes the associated legendre polynomials based on 
    note that this associated legendre uses the (-1)**m coefficient.
    Denoting it as P_l^m, it is the same as
    P_l^m(beta) = d_l^{m0}(beta)
    Also, if you for spherical harmonic s_l^m(th,beta)
    s_l^m(th,beta) = P_l^m(beta)*sqrt((2*l+1)/4/np.pi)*exp(1j*m*th)

    and that s_l^m forms an orthonormal basis
    
    Note that s_l^m(th,beta) can be obtain from
    scipy.special.sph_harm(m,l,th,beta)    
        
    '''
    from scipy.special import factorial
    from scipy.special import gamma
    from scipy.special import comb
    x=np.cos(beta)
    coeff1=2**l*(1-x**2)**(m/2.0)
    coeff2=0
    for j in range(m,(l+1)):
        alph=(l+j-1)/2.0
        coeff2+= factorial(j)/factorial(j-m)*x**(j-m)*comb(l,j)* \
            gamma(alph+1)/gamma(l+1)/gamma(alph-l+1)
        
    Plm=(-1)**m*coeff1*coeff2*np.sqrt(factorial(l-m)/factorial(l+m))
    return Plm

def computeSPHCoeff(qIn,isw,lU):

    '''
    Computes the coefficients of the spherical harmonics
    P(h,y) = \sum_{l=0}^{L}\sum_{m=-l}^l \sum_{n=-1}^l
        c_l^{mn} \overline{Y_l^m(h)} Y_l^n(y)


    qIn ... N x 4 array of samples of quats for polycrstal
    isw ... = 0 for pole figure or 1 for texture plot
    lU .... interger that is maximum degree of polynomials
    
    '''
    from scipy.special import sph_harm

    if (isw==0):
        Nth=64
        Nphi=32
        Nsym=48
        Nq=qIn.shape[0]
        th=np.linspace(0,2*np.pi,Nth)
        dth=th[1]-th[0]
        phi=np.linspace(0,np.pi,Nphi)
        dphi=phi[1]-phi[0]
        thv=np.repeat(th,Nphi)
        phiv=np.tile(phi,Nth)
        Nt=Nth*Nphi
        rv0=np.zeros([Nt,3])
        rv0[:,0]=1.0*np.cos(thv)*np.sin(phiv)
        rv0[:,1]=1.0*np.sin(thv)*np.sin(phiv)
        rv0[:,2]=1.0*np.cos(phiv)
        rv=ComputeCubicSymmetries(rv0)
        rvR=QuatRotate(qIn,rv)
        a=np.arctan2(rvR[:,1],rvR[:,0])                
        thvh=a*(a>0) + (a+2*np.pi)*(a<0)
        phivh=np.arccos(rvR[:,2])
        nlm=np.sum( (np.arange(3,(2*lU+1)+2,2))**2)+1
        clmn=np.zeros(nlm,dtype=complex)
        cc=0
        for l in range(0,(lU+1)):
            for m in range(-l,(l+1)):
                z=sph_harm(m,l,thvh,phivh)
                sph1=np.conjugate(np.sum(np.reshape(z,[Nt,Nsym*Nq]),axis=1))/(Nq*Nsym)
                sph2=sph_harm(m,l,thv,phiv)
                dv=np.sin(phiv)*dth*dphi
                z=np.sum(np.conjugate(sph1)*sph2*dv)
                ni1=2*l+1
                clmn[cc:(cc+ni1)]=z*np.ones(ni1)
                cc+=ni1
    return clmn

def evaluateSPHPF(clmn,lU,isw,hhat,yhat):
    '''
    evaluates the pole figure (or texture plot)
    P(h,y) = \sum_{l=0}^{L}\sum_{m=-l}^l \sum_{n=-1}^l
        c_l^{mn} \overline{Y_l^m(h)} Y_l^n(y)


    clmn ... array of coefficients
    isw .... = 0 for pole figure or 1 for texture plot
    lU ..... interger that is maximum degree of polynomials
    hhat ... values of the crystallography orientations if pole figure
             Nx3. if texture plot 1x3
    yhat ... values of reference frame direction. if pole figure then
             1x3. if texture plot Nx3
    
    '''
    from scipy.special import sph_harm
    
    if (isw==0):
        a=np.arctan2(yhat[1],yhat[0])                
        thyhat=a*(a>0) + (a+2*np.pi)*(a<0)
        phiyhat=np.arccos(yhat[2])  
        Nsym=48
        Nt=hhat.shape[0]
        hs=ComputeCubicSymmetries(hhat)
        a=np.arctan2(hs[1],hs[0])                
        thvh=a*(a>0) + (a+2*np.pi)*(a<0)
        phivh=np.arccos(hs[:,2])
        '''
        a=np.arctan2(hhat[1],hhat[0])                
        thvh=a*(int(a>0)) + (a+2*np.pi)*(int(a<0))
        phivh=np.arccos(hhat[:,2])
        '''
        cc=0
        N=hhat.shape[0]
        pf=np.zeros(N,dtype=complex)
        for l in range(0,(lU+1)):
            for m in range(-l,(l+1)):
                z=sph_harm(m,l,thvh,phivh)
                sph1=np.conjugate(np.sum(np.reshape(z,[Nt,Nsym]),axis=1))/Nsym
                for n in range(-l,(l+1)):
                    sph2=sph_harm(n,l,thyhat,phiyhat)
                    pf+=clmn[cc]*sph1*sph2

    return pf
    
    

def ComputeMisorientationsCubicQuats(quat1,quat2):
    '''
    This approach follows section 1.3.3.4 in
    Sutton, Adrian P. "Interfaces in crystalline materials." Monographs
    on the Physice and Chemistry of Materials (1995): 414-423.    
    
    note that we use the format: quat = (q_0, \textbf{q}), where q_0 is the
    scalar and \textbf{q} is the vector. Therefore, i need to do a sort to
    go from the formula in this book back to the format we use
    
    '''
    
    q12 = QuatMultiply(quat1,QuatConjugate(quat2))
    q12 = np.abs(q12)
    q12 = np.sort(q12)
    c = np.zeros([3,4])
    c[0,:] = q12
    c[1,:] = np.array([q12[0]-q12[1],q12[0]+q12[1],\
                     q12[2]-q12[3],q12[2]+q12[3]])*1.0/np.sqrt(2.0)
    c[2,:] = np.array([q12[0]-q12[1]+q12[2]-q12[3], q12[0]+q12[1]-q12[2]-q12[3], \
                  -q12[0]+q12[1]+q12[2]-q12[3], q12[0]+q12[1]+q12[2]+q12[3]])*\
                  1.0/2.0
    i1 = np.argmax(c[:,3])
    q12 = c[i1,[3,2,1,0]]           
    return q12       
    
def ComputeMisorientationFCCTwin(q12):
    '''
    Computes misorientation angle of quaternion q12 (which is ought to be
    the misorientation between grains) with 60 degrees about <1,1,1>. Note that
    if q12 is computed using ComputeMisorientationCubicQuats, then 
    symmetry has already been accounted for and only 60<1 1 1> needs to be 
    considered 
    '''
    rho = np.array([1.0,1.0,1.0])
    theta = np.pi/3.0
    qtwin = np.array([np.cos(theta/2.0),rho[0]*np.sin(theta/2.0), \
                      rho[1]*np.sin(theta/2.0),rho[2]*np.sin(theta/2.0)])
        
def ComputeElasticTensorCubic(C11,C12,C44,eulerAngles):
    '''
    computes the 6x6 elasticity tensor for a cubic material in the global
    reference frame    
    '''
    cmap2 = np.array([[1,2,3,4,5,6],\
                      [2,7,8,9,10,11],\
                      [3,8,12,13,14,15],\
                      [4,9,13,16,17,18], \
                      [5,10,14,17,19,20],\
                      [6,11,15,18,20,21]],dtype=int)
    cmap = np.array([[1,2,3,1,2,3],[1,2,3,2,3,1]],dtype=int).transpose()
    cmap2 = cmap2-1 # note that above is fortran counting convention
    cmap = cmap - 1 # note that above is fortran counting convention    
    th1 = eulerAngles[0]
    phi = eulerAngles[1]
    th2 = eulerAngles[2]
    Q = np.zeros([3,3])
    Q[0,0] = np.cos(th2)*np.cos(th1)-np.sin(th2)*np.cos(phi)*np.sin(th1)
    Q[0,1] = np.cos(th2)*np.sin(th1)+np.sin(th2)*np.cos(phi)*np.cos(th1)
    Q[0,2] = np.sin(th2)*np.sin(phi)
    Q[1,0] = -np.sin(th2)*np.cos(th1) - np.cos(th2)*np.cos(phi)*np.sin(th1)
    Q[1,1] = -np.sin(th2)*np.sin(th1)+np.cos(th2)*np.cos(phi)*np.cos(th1)
    Q[1,2] = np.cos(th2)*np.sin(phi)
    Q[2,0] = np.sin(phi)*np.sin(th1)
    Q[2,1] = -np.sin(phi)*np.cos(th1)
    Q[2,2] = np.cos(phi)
    chet = np.zeros(21)
    c=0
    for j2 in range(0,6):
        for j3 in range(j2,6):
            chet[c] = C11*( \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],0]*Q[cmap[j3,0],0]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],1]*Q[cmap[j3,0],1]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],2]*Q[cmap[j3,0],2]*Q[cmap[j3,1],2] ) + C12*( \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],0]*Q[cmap[j3,0],1]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],0]*Q[cmap[j3,0],2]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],1]*Q[cmap[j3,0],2]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],1]*Q[cmap[j3,0],0]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],2]*Q[cmap[j3,0],0]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],2]*Q[cmap[j3,0],1]*Q[cmap[j3,1],1] ) + C44*( \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],2]*Q[cmap[j3,0],1]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],1]*Q[cmap[j3,0],1]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],2]*Q[cmap[j3,0],2]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],1]*Q[cmap[j3,0],2]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],0]*Q[cmap[j3,0],2]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],2]*Q[cmap[j3,0],2]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],2]*Q[cmap[j2,1],0]*Q[cmap[j3,0],0]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],2]*Q[cmap[j3,0],0]*Q[cmap[j3,1],2] + \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],1]*Q[cmap[j3,0],0]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],0]*Q[cmap[j3,0],0]*Q[cmap[j3,1],1] + \
             Q[cmap[j2,0],0]*Q[cmap[j2,1],1]*Q[cmap[j3,0],1]*Q[cmap[j3,1],0] + \
             Q[cmap[j2,0],1]*Q[cmap[j2,1],0]*Q[cmap[j3,0],1]*Q[cmap[j3,1],0] )
            c=c+1
    cvoigt=np.zeros([6,6])
    for j1 in range(0,6):
        for j2 in range(0,6):
            cvoigt[j1,j2] = chet[cmap2[j1,j2]]
    return cvoigt        
    
def StereographicProjectionCubic(vecIn):
    ''' 
        computes stereographic projection of vecIn and returns x,y values.
        It is assumed vecIn values are (1) unit length 
        (2) is matrix of N x 3 (even if N=1, must be matrix). Output is a
        N x 2 matrix of x and y position on pole
    '''            
    N = vecIn.shape[0]
    posXY = np.zeros([N,2])

    th = np.arccos(vecIn[:,2])
    phi = np.arctan2(vecIn[:,1],vecIn[:,0])
    posXY[:,0] = np.tan(th/2.0)*np.cos(phi)
    posXY[:,1] = np.tan(th/2.0)*np.sin(phi)
    return posXY

def ComputeCubicSymmetries(vecIn):
    '''
        Computes the 48 equivalent vectors for a cubic system.
        Cubic symmetry includes (1) 4 fold symmetry about axis
        connecting centers of opposite faces (3x3=9), (2) 3 fold symmetry
        of axis connecting a vertex and the center (4x2=8), (3) 2 fold
        symmetry of axis midpoint of edge with center (6x1=6), (4) the 
        identity (1), and all the rotations involved followed by an
        inversion
    '''
    if (len(vecIn.shape)==1):
        vecIn=np.array([vecIn])
    N=vecIn.shape[0]
    quatSymm = np.zeros([24,4])
    quatSymm[0,:] = np.append(np.cos(0/2),np.sin(0/2)*np.array([0,0,1]))
    ca = np.cos(np.pi/2/2); sa = np.sin(np.pi/2/2)
    quatSymm[1,:] = np.append(ca,sa*np.array([0,0,1]))
    quatSymm[4,:] = np.append(ca,sa*np.array([1,0,0]))
    quatSymm[7,:] = np.append(ca,sa*np.array([0,1,0]))        
    ca = np.cos(np.pi/2); sa = np.sin(np.pi/2)
    quatSymm[2,:] = np.append(ca,sa*np.array([0,0,1]))
    quatSymm[5,:] = np.append(ca,sa*np.array([1,0,0]))
    quatSymm[8,:] = np.append(ca,sa*np.array([0,1,0]))            
    ca = np.cos(3*np.pi/2/2); sa = np.sin(3*np.pi/2/2)
    quatSymm[3,:] = np.append(ca,sa*np.array([0,0,1]))
    quatSymm[6,:] = np.append(ca,sa*np.array([1,0,0]))
    quatSymm[9,:] = np.append(ca,sa*np.array([0,1,0]))        
    a=1.0/np.sqrt(3)    
    ca = np.cos(2*np.pi/3/2); sa = np.sin(2*np.pi/3/2)
    quatSymm[10,:] = np.append(ca,sa*np.array([a,a,a]))
    quatSymm[12,:] = np.append(ca,sa*np.array([a,a,-a]))
    quatSymm[14,:] = np.append(ca,sa*np.array([a,-a,a]))
    quatSymm[16,:] = np.append(ca,sa*np.array([a,-a,-a]))
    ca = np.cos(4*np.pi/3/2); sa = np.sin(4*np.pi/3/2)
    quatSymm[11,:] = np.append(ca,sa*np.array([a,a,a]))
    quatSymm[13,:] = np.append(ca,sa*np.array([a,a,-a]))
    quatSymm[15,:] = np.append(ca,sa*np.array([a,-a,a]))
    quatSymm[17,:] = np.append(ca,sa*np.array([a,-a,-a]))
    a=1.0/np.sqrt(2)
    ca = np.cos(np.pi/2); sa = np.sin(np.pi/2)
    quatSymm[18,:] = np.append(ca,sa*np.array([a,0,-a]))
    quatSymm[19,:] = np.append(ca,sa*np.array([a,0,a]))
    quatSymm[20,:] = np.append(ca,sa*np.array([0,a,-a]))
    quatSymm[21,:] = np.append(ca,sa*np.array([0,a,a]))
    quatSymm[22,:] = np.append(ca,sa*np.array([a,a,0]))
    quatSymm[23,:] = np.append(ca,sa*np.array([a,-a,0]))
    v2=np.repeat(vecIn,2,axis=0)
    i1=np.reshape(np.tile(np.array([1,-1]),N),[2*N,1])
    v2=i1*v2
    vecOut = QuatRotate(quatSymm,v2)
    '''
    for j in range(0,N):
        for j1 in range(0,24):
            vecOut[48*j+j1,:] = QuatRotate(quatSymm[j1,:],vecIn[j,:])
        for j1 in range(0,24):
            vecOut[48*j+24+j1,:] = QuatRotate(quatSymm[j1,:],-vecIn[j,:])
    '''
    return vecOut

def ComputeMRD(vecIn,weight,meshR):
    '''
    vecIn .... unit normal vectors on sphere
    weight ... weight (or volume) associated with each vector
    meshR .... resolution of meshgrid for square for computing histogram
          .... the number of grid cells is (2meshR)**2
    This program computes the multiples of random distribution of the 
    vector distribution on the sphere. The approach taken follows:
        
        Roşca, D. "New uniform grids on the sphere." Astronomy & 
        Astrophysics 520 (2010): A63.APA	
    
    First the vectors are mapped onto a grid (via Lambert projection), the 
    histogram is computed as well as the empirical pdf normalized with 
    respect to a completely random distribution of vector orientations.
    The centroids of the gridpoints are then inverted back to the points on
    the sphere (vecOut). Along with these points, the multiples of random
    distribution (i.e., intensities) are passed back (MRD).
    '''

    N = vecIn.shape[0]
    r = 1.0
    # FYI radius of Lambert disc = np.sqrt(2)*r (but not used below)
    L = np.sqrt(np.pi/2.0)*r
    DL = L/meshR
    Nc = 2*meshR
    hBinNorth = np.zeros([Nc,Nc])
    hBinSouth = np.zeros([Nc,Nc])    
    vecOut = np.zeros([2*Nc**2,3])
    MRD = np.zeros(2*Nc**2)
    wSouth = 0
    wNorth = 0
    for j in range(0,N):
        x = vecIn[j,0] 
        y = vecIn[j,1]
        z = vecIn[j,2] 
        (A,B) = (np.sqrt(2*r/(r+np.abs(z)))*x,np.sqrt(2*r/(r+np.abs(z)))*y)
        if (np.abs(B) < np.abs(A)):
            (xs,ys) = (np.sign(A)*np.sqrt(A**2+B**2)*np.sqrt(np.pi)/2, \
                np.sign(A)*np.sqrt(A**2+B**2)*2/np.sqrt(np.pi)*np.arctan(B/A))
        if (np.abs(A) <= np.abs(B)):
            if ( (A==0) and (B==0)):
                (xs,ys)=0,0
            else:
                (xs,ys) = (np.sign(B)*np.sqrt(A**2+B**2)*\
                 2/np.sqrt(np.pi)*np.arctan(A/B),np.sign(B)*np.sqrt(A**2+B**2)*\
                 np.sqrt(np.pi)/2)
        jx = np.min([int(np.floor( (1-2*np.finfo(float).eps)*(xs+L)/DL)),Nc-1])
        jy = np.min([int(np.floor( (1-2*np.finfo(float).eps)*(ys+L)/DL)),Nc-1])
        if (z<0):
            hBinSouth[jx,jy] = hBinSouth[jx,jy]+weight[j]
            wSouth = wSouth + weight[j]
        else:
            hBinNorth[jx,jy] = hBinNorth[jx,jy]+weight[j]    
            wNorth = wNorth + weight[j]            
    cc=0
    for jy in range(0,Nc):
        for jx in range(0,Nc):
            xs = (jx + 0.5)*DL - L
            ys = (jy + 0.5)*DL - L        
            if (np.abs(ys)<=np.abs(xs)):
                (A,B) = (2*xs/np.sqrt(np.pi)*np.cos(ys*np.pi/4/xs),
                         2*xs/np.sqrt(np.pi)*np.sin(ys*np.pi/4/xs))
            if (np.abs(xs)<=np.abs(ys)):
                (A,B) = (2*ys/np.sqrt(np.pi)*np.sin(xs*np.pi/4/ys),
                         2*ys/np.sqrt(np.pi)*np.cos(xs*np.pi/4/ys))
            x2 = np.sqrt(1-(A**2+B**2)/4/r**2)*A
            y2 = np.sqrt(1-(A**2+B**2)/4/r**2)*B
            z2 = -r + (A**2+B**2)/2/r
            vecOut[cc,:] = np.array([x2,y2,z2])
            MRD[cc] = hBinSouth[jx,jy]/wSouth*Nc**2
            z2 = r - (A**2+B**2)/2/r
            vecOut[Nc**2+cc,:] = np.array([x2,y2,z2]) 
            MRD[Nc**2+cc] = hBinNorth[jx,jy]/wNorth*Nc**2            
            cc=cc+1

    return (vecOut,MRD)

def MakeContourCubic(posXY,MRD,Np):
    '''
    posXY ....... stereographic xy position (Nx2)
    MRD ......... multiples of random intensities (N)
    Np .......... number of points to plot in contour
    '''
    import scipy.interpolate
    import matplotlib.pyplot as plt
#    xv = np.linspace(np.min(posXY[:,0]),np.max(posXY[:,0]),Np)
#    yv = np.linspace(np.min(posXY[:,1]),np.max(posXY[:,1]),Np)
    xv = np.linspace(-1,1,Np)
    yv = np.linspace(-1,1,Np)

    zv=scipy.interpolate.griddata((posXY[:,0],posXY[:,1]),MRD,(xv[None,:],yv[:,None]),method='linear')
    i1 = np.zeros([Np,Np])
    for j in range(0,Np):
        for j1 in range(0,Np):
            r = np.sqrt(xv[j]**2 + yv[j1]**2)
            if (r>1):
                i1[j,j1] = 1
    zvMask=np.ma.masked_array(zv,mask=i1)
    #make lines for great circles
    vecGCNodes = np.zeros([40,6])
    ap=np.zeros([17,3])
    ap[0,:]=np.array([-1,0,0])
    ap[1,:]=np.array([-1,1,1])/np.sqrt(3)
    ap[2,:]=np.array([0,1,1])/np.sqrt(2)
    ap[3,:]=np.array([1,1,1])/np.sqrt(3)
    ap[4,:]=np.array([1,0,0])
    ap[5,:]=np.array([-1,-1,1])/np.sqrt(3)
    ap[6,:]=np.array([0,-1,1])/np.sqrt(2)
    ap[7,:]=np.array([1,-1,1])/np.sqrt(3)
    ap[8,:]=np.array([-1,0,1])/np.sqrt(2)
    ap[9,:]=np.array([1,0,1])/np.sqrt(2)
    ap[10,:]=np.array([1,1,0])/np.sqrt(2)
    ap[11,:]=np.array([0,1,0])
    ap[12,:]=np.array([-1,1,0])/np.sqrt(2)
    ap[13,:]=np.array([-1,-1,0])/np.sqrt(2)
    ap[14,:]=np.array([0,-1,0])
    ap[15,:]=np.array([1,-1,0])/np.sqrt(2)
    ap[16,:]=np.array([0,0,1])    
    vecGCNodes[0,:] = np.append(ap[0,:],ap[1,:])
    vecGCNodes[1,:] = np.append(ap[1,:],ap[2,:])    
    vecGCNodes[2,:] = np.append(ap[2,:],ap[3,:])
    vecGCNodes[3,:] = np.append(ap[3,:],ap[4,:])    
    vecGCNodes[4,:] = np.append(ap[0,:],ap[5,:])        
    vecGCNodes[5,:] = np.append(ap[5,:],ap[6,:])        
    vecGCNodes[6,:] = np.append(ap[6,:],ap[7,:])        
    vecGCNodes[7,:] = np.append(ap[7,:],ap[4,:])            
    vecGCNodes[8,:] = np.append(ap[14,:],ap[7,:])        
    vecGCNodes[9,:] = np.append(ap[7,:],ap[9,:])        
    vecGCNodes[10,:] = np.append(ap[9,:],ap[3,:])
    vecGCNodes[11,:] = np.append(ap[3,:],ap[11,:])    
    vecGCNodes[12,:] = np.append(ap[14,:],ap[5,:])    
    vecGCNodes[13,:] = np.append(ap[5,:],ap[8,:])    
    vecGCNodes[14,:] = np.append(ap[8,:],ap[1,:])    
    vecGCNodes[15,:] = np.append(ap[1,:],ap[11,:])    
    vecGCNodes[16,:] = np.append(ap[12,:],ap[1,:])    
    vecGCNodes[17,:] = np.append(ap[1,:],ap[16,:])    
    vecGCNodes[18,:] = np.append(ap[16,:],ap[7,:])    
    vecGCNodes[19,:] = np.append(ap[7,:],ap[15,:])    
    vecGCNodes[20,:] = np.append(ap[13,:],ap[5,:])        
    vecGCNodes[21,:] = np.append(ap[5,:],ap[16,:])            
    vecGCNodes[22,:] = np.append(ap[16,:],ap[3,:])            
    vecGCNodes[23,:] = np.append(ap[3,:],ap[10,:])            
    vecGCNodes[24,:] = np.append(ap[0,:],ap[8,:])            
    vecGCNodes[25,:] = np.append(ap[8,:],ap[16,:])            
    vecGCNodes[26,:] = np.append(ap[16,:],ap[9,:])                
    vecGCNodes[27,:] = np.append(ap[9,:],ap[4,:])                
    vecGCNodes[28,:] = np.append(ap[14,:],ap[6,:])            
    vecGCNodes[29,:] = np.append(ap[6,:],ap[16,:])            
    vecGCNodes[30,:] = np.append(ap[16,:],ap[2,:])            
    vecGCNodes[31,:] = np.append(ap[2,:],ap[11,:])            
    vecGCNodes[32,:] = np.append(ap[14,:],ap[13,:])            
    vecGCNodes[33,:] = np.append(ap[13,:],ap[0,:])            
    vecGCNodes[34,:] = np.append(ap[0,:],ap[12,:])                
    vecGCNodes[35,:] = np.append(ap[12,:],ap[11,:])                
    vecGCNodes[36,:] = np.append(ap[11,:],ap[10,:])                
    vecGCNodes[37,:] = np.append(ap[10,:],ap[4,:])            
    vecGCNodes[38,:] = np.append(ap[4,:],ap[15,:])            
    vecGCNodes[39,:] = np.append(ap[15,:],ap[14,:])                
    Ninc = 10
    vecGC = np.zeros([40,(Ninc+1)*3])
    vecGC[:,0:3] = vecGCNodes[:,0:3]
    vecGC[:,(Ninc*3):((Ninc+1)*3)]=vecGCNodes[:,3:6]
    rRot = np.zeros([3,3])
    for j in range(0,40):
        ep1 = vecGCNodes[j,0:3]
        ep3 = np.cross(ep1,vecGCNodes[j,3:6])
        ep3 =ep3/np.linalg.norm(ep3)
        ep2 = np.cross(ep3,ep1)
        rRot = np.array([ep1,ep2,ep3])
        thMax = np.arccos(np.dot(vecGCNodes[j,0:3],vecGCNodes[j,3:6]))
        for j1 in range(1,Ninc):
            thj = thMax/Ninc*j1
            vp = np.array([np.cos(thj),np.sin(thj),0])
            vecGC[j,(3*j1):(3*(j1+1))] = np.matmul(np.transpose(rRot),vp)
    xyGC = np.zeros([40,(Ninc+1)*2])
    for j in range(0,(Ninc+1)):
            xyGC[:,(2*j):(2*(j+1))]= \
            StereographicProjectionCubic(vecGC[:,(3*j):(3*(j+1))])    
            #xyGC[:,(2*j):(2*(j+1))]= \
            #crystallographicRoutines.StereographicProjectionCubic(vecGC[:,(3*j):(3*(j+1))])    

    plt.figure(figsize=[10,8])
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)
    #plt.contourf(xv,yv,zv,cmap=plt.cm.jet,corner_mask=True)
    plt.contourf(xv,yv,zvMask,cmap=plt.cm.jet,corner_mask=True)    
    for j in range(0,40):
        for j1 in range(0,Ninc):
            plt.plot([xyGC[j,2*j1],xyGC[j,2*j1+2]],[xyGC[j,2*j1+1],xyGC[j,2*j1+3]], \
                     'w',linewidth=1)
    #ax1=fig.add_axes([0.0, 0.8,.80, .5])    
    plt.colorbar(orientation='horizontal')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    
def ComputeSchmidFactorsFCC(sigma,eulerAngles):

    th1 = eulerAngles[0]
    phi = eulerAngles[1]
    th2 = eulerAngles[2]
    Q = np.zeros([3,3])
    Q[0,0] = np.cos(th2)*np.cos(th1)-np.sin(th2)*np.cos(phi)*np.sin(th1)
    Q[0,1] = np.cos(th2)*np.sin(th1)+np.sin(th2)*np.cos(phi)*np.cos(th1)
    Q[0,2] = np.sin(th2)*np.sin(phi)
    Q[1,0] = -np.sin(th2)*np.cos(th1) - np.cos(th2)*np.cos(phi)*np.sin(th1)
    Q[1,1] = -np.sin(th2)*np.sin(th1)+np.cos(th2)*np.cos(phi)*np.cos(th1)
    Q[1,2] = np.cos(th2)*np.sin(phi)
    Q[2,0] = np.sin(phi)*np.sin(th1)
    Q[2,1] = -np.sin(phi)*np.cos(th1)
    Q[2,2] = np.cos(phi)
    slocal=np.zeros([12,3])
    slocal[0,:] = [-1.0,1.0,0.0]
    slocal[1,:] = [-1.0,0.0,1.0]
    slocal[2,:] = [0.0,-1.0,1.0]
    slocal[3,:] = [1.0,1.0,0.0]
    slocal[4,:] = [1.0,0.0,1.0]
    slocal[5,:] = [0.0,-1.0,1.0]
    slocal[6,:] = [-1.0,1.0,0.0]
    slocal[7,:] = [1.0,0.0,1.0]
    slocal[8,:] = [0.0,1.0,1.0]
    slocal[9,:] = [1.0,1.0,0.0]
    slocal[10,:] = [0.0,1.0,1.0]
    slocal[11,:] = [-1.0,0.0,1.0]
    slocal = slocal/np.sqrt(2.0)
    mlocal=np.zeros([12,3])
    mlocal[0,:] = [1.0,1.0,1.0]
    mlocal[1,:] = [1.0,1.0,1.0]
    mlocal[2,:] = [1.0,1.0,1.0]
    mlocal[3,:] = [-1.0,1.0,1.0]
    mlocal[4,:] = [-1.0,1.0,1.0]
    mlocal[5,:] = [-1.0,1.0,1.0]
    mlocal[6,:] = [1.0,1.0,-1.0]
    mlocal[7,:] = [1.0,1.0,-1.0]
    mlocal[8,:] = [1.0,1.0,-1.0]
    mlocal[9,:] = [1.0,-1.0,1.0]
    mlocal[10,:] = [1.0,-1.0,1.0]
    mlocal[11,:] = [1.0,-1.0,1.0]
    mlocal = mlocal/np.sqrt(3.0)
    schmidFactors = np.zeros(12)
    ml=np.zeros([12,3,3])
    for j in range(0,12):
        ml[j,:,:] = np.matmul(np.reshape(mlocal[j,:],[3,1]),np.reshape(slocal[j,:],[1,3]))
    mg=np.zeros([12,3,3])
    for j in range(0,12):
        mg[j,:,:] = np.matmul(Q,np.matmul(ml[j,:,:],np.transpose(Q)))
        schmidFactors[j] = float(np.tensordot(sigma,mg[j,:,:]))
    return schmidFactors


def getFSCoeffPlm(l,m,**kwargs):
    
    '''
    Computes the fourier series coefficients for the associated
    legendre polynomial P_l^m(x). It either computes as many as necessary
    such that the error is less than tol or computes a fixed number
    
    You must call this program either specifying tol or nS (if both, then
    nS will prevail)
    
    1) cR.getFSCoeffPlm(l,m,tol=xtol,nS=xNs)
    2) cR.getFSCoeffPlm(l,m,tol=xtol)
    3) cR.getFSCoeffPlm(l,m,nS=xNs)
    where xNs and xtol are variables you've defined
    
    '''
    for key,value in kwargs.items():
        if key=="nS":
            nS = value
            toltmp=0.0
        if key=="tol":
            tol=value
            nStmp=1000
    try:
        if (nS < 1000):
            tol=toltmp
    except: 
        nS=nStmp   
    N=8192
    beta=np.linspace(0,np.pi,N)
    db=beta[1]-beta[0]
    x=np.linspace(-1,1,N)
    la=legendreAX(l,m,x)
    a2=np.zeros(nS)
    a2[0]=np.sum(la)*db/np.pi
    las2 = a2[0]*np.cos(0*beta)    
    lal2=np.linalg.norm(la)    
    nF=1
    weOut=1
    if (np.mod(m,2)==0):
        while (nF<nS and weOut==1):
            cosb1=np.cos(nF*beta)
            a2[nF]=2/np.pi*np.sum(la*cosb1)*db
            las2 += a2[nF]*np.cos(nF*beta)    
            err=np.linalg.norm(la-las2)/lal2
            nF+=1
            if (err<tol):
                weOut=0
    else:
        while (nF<nS and weOut==1):
            sinb1=np.sin(nF*beta)
            a2[nF]=2/np.pi*np.sum(la*sinb1)*db
            las2 += a2[nF]*np.sin(nF*beta)    
            err=np.linalg.norm(la-las2)/lal2
            nF+=1
            if (err<tol):
                weOut=0
            
    a=a2[0:nF]
    return a

def computeRGB(vecIn):
    '''                                                                                                                                                                
    function takes normal vector(or list of them), vecIn, from the cubic unit 
    triangle(please apply all symetry before hand) and give a color for it in
    RGB.  
    '''
    if (len(vecIn.shape)==1):
        vecIn = np.reshape(vecIn,[1,3])
    N=vecIn.shape[0]

    xp=(2*vecIn[:,0])/(1+vecIn[:,2])
    yp=(2*vecIn[:,1])/(1+vecIn[:,2])
    triPts=np.transpose(np.array([[0,0],[2./np.sqrt(2.)/(1.+1./np.sqrt(2.)),0], \
                     [2./np.sqrt(3.)/(1.+1./np.sqrt(3)),2./np.sqrt(3)/(1. +1./np.sqrt(3.))]]))
    
    m=np.tan(1./2.*np.arctan2(triPts[1,2],triPts[0,2]))
    a=np.sqrt( (triPts[1,2]-triPts[1,1])**2 + (triPts[0,2]-triPts[0,1])**2)
    b=np.sqrt( triPts[0,1]**2 + triPts[1,1]**2 )
    c=np.sqrt( triPts[0,2]**2 + triPts[1,2]**2 )
    y0=1./2.*np.sqrt( (b+c-a)*(c+a-b)*(a+b-c)/(a+b+c) )
    x0=y0/m
    S=np.sqrt( (xp-x0)**2 + (yp-y0)**2)
    H=np.arctan2( (yp-y0),(xp-x0))*180./np.pi
    V=np.ones(N)
    i1=np.where(xp<x0)[0]
    H[i1]=H[i1]+180
    sMax=np.sqrt(x0**2+y0**2)
    S=S/sMax*0.8 + 0.2
    HSV=np.transpose(np.reshape(np.concatenate((H,S,V)),[3,N]))
    RGB=hsv2rgb(HSV)
    return RGB

def computeIPF(quats,r):
    if (len(quats.shape)==1):
        quats=np.reshape(quats,[1,4])
    #N=quats.shape[0]
    vec1= np.abs(QuatRotate(quats,r))
    vec = np.sort(vec1,axis=1)
    vec[:,[0,1]] = vec[:,[1,0]]
    RGB=computeRGB(vec)
    return RGB
    
def hsv2rgb(hsv):
    if (len(hsv.shape)==1):
        hsv=np.reshape(hsv,[1,3])
    N=hsv.shape[0]
    rgb=np.zeros([N,3])
    hh=hsv[:,0]
    i1=np.where(hh>=360)[0]
    hh[i1]=0.0
    hh/=60.0
    i=int(hh)
    ff=hh-i
    p=hsv[:,2]*(1.-hsv[:,1])
    q=hsv[:,2]*(1.-hsv[:,1]*ff)
    t=hsv[:,2]*(1.-(hsv[:,1]*(1.-ff)))
    if (i==0):
        rgb[:,0]=hsv[:,2]
        rgb[:,1]=t
        rgb[:,2]=p
    if (i==1):
        rgb[:,0]=q
        rgb[:,1]=hsv[:,2]
        rgb[:,2]=p
    if (i==2):
        rgb[:,0]=p
        rgb[:,1]=hsv[:,2]
        rgb[:,2]=t
    if (i==3):
        rgb[:,0]=p
        rgb[:,1]=q
        rgb[:,2]=hsv[:,2]
    if (i==4):
        rgb[:,0]=t
        rgb[:,1]=p
        rgb[:,2]=hsv[:,2]
    if (i==5):
        rgb[:,0]=hsv[:,2]
        rgb[:,1]=p
        rgb[:,2]=q
        
    i1=np.where(hsv[:,1]<=0.0)[0]
    rgb[i1,0]= hsv[i1,2]
    rgb[i1,1]= hsv[i1,2]
    rgb[i1,2]= hsv[i1,2]    
    return rgb

