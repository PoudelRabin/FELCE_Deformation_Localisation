//********************NOTE BY RABIN POUDEL******************************
/*Note for the fourth order operation
dyad operation of the tens4d are diffrent to the one done in ...
matlab (ie Mikhail Itskov Book). The equivalent are

TensorProd_AoB(A,B) = dyad1(A,B) 
TensorProd_AxB(A,B) = dyad2(A,B.transpose())
TensorProd_AxBt(A,B) = dyad3(A, B.transpose())

AoB(i,j,k,l) = A(i,j)*B(k,l)----
AxB(i,j,k,l) = A(i,k)*B(l,j)    |--Defination in Mikhail Itskov Book 
AxBt(i,j,k,l) = A(i,l)*B(k,j)---

A(dyad1)B(i,j,k,l) = A(i,j)*B(k,l)----
A(dyad2)B(i,j,k,l) = A(i,k)*B(j,l)    |--Defination in tens4d.hpp
A(dyad3)B(i,j,k,l) = A(i,l)*B(j,k)----


The 9*9 matrix obtaine by tensor3333_99 from 4th order is not same
extract 9*9 used in tens4d. The order of indices for voigt notation
are not in same order (see tens4d.h and tens4d.hpp). Also 9*9 
matrix function is to be first transposed for our use(since tens4d.h
first fills the single column). we use the voigt order as in tens4d.h*/
//**********************************************************************




#include "stdafx.h"
#include "FELCEMaterial444.h"
#include<FECore/FESolidDomain.h>
#include "FECore/log.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FELCEMaterial444, FEElasticMaterial)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER(0.0), "mu");
	ADD_PARAMETER(m_kfrank, FE_RANGE_GREATER(0.0), "K");
    ADD_PARAMETER(m_a, FE_RANGE_GREATER(0.0), "a");
    ADD_PARAMETER(m_bulk, FE_RANGE_GREATER(0.0), "B");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FELCEMaterial444::FELCEMaterial444(FEModel* pfem) : FEElasticMaterial(pfem) {}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
mat3ds FELCEMaterial444::Stress(FEMaterialPoint& mp)
{
    mat3ds s;
    s.zero();
    return s;
}
//-----------------------------------------------------------------------------
tens4ds FELCEMaterial444::Tangent(FEMaterialPoint& mp)
{
    tens4ds T;
    T.zero();
    return T;
}

//-----------------------------------------------------------------------------
mat3d FELCEMaterial444::StressLCE(FESolidElement el, int ngauss, double theta,  
        double theta0, mat3d F, matrix Gradt, double pressmix, matrix &PK1_vgt, 
        matrix &Pi0, matrix &gt)
            
{
    //Declare some temporary double, matrix and tensors 
    double temp1,temp2,temp3,temp4,temp5; 
    mat3d tempmat1, tempmat2, tempmat3, tempmat4, tempmat5;
    tens4d tempA, tempB, tempC, tempD, tempE;

    
    
    //-----------------------------------------------------
    //Get material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

    //Get Current time
    double cTime = GetFEModel()->GetTime().currentTime;

    //GEt material parameter value
    //todo is writing (0) good? Change it for material point, here assumed parameter is same everywhere
    double mu = m_mu(mp);
    double mu_neo = 1.0*mu;

    //mu = ((1.0e-10 - 2.0e-2)/(1.20-1.0))*(cTime+1.0-1.20) + 1.0e-10;
    //if (cTime>0.20) {mu = 1.0e-10;}
    double kfrank = m_kfrank(mp);

    //Get deformation gradient and its determinant
    double detF = F.det();
    mat3d Finv = F.inverse();
    mat3d Fbar = F* pow(detF, -1.0 / 3.0);


    //double a = (1.0+2.0*Qorder)/(1.0-Qorder);
    double a = m_a(mp);

    double a0 = a;
    double bulkmod = m_bulk(mp);
    //--------------------------
    
    
    
    //-------------
    mat3d eye;
    eye.zero();
    eye[0][0] = 1.0;
    eye[1][1] = 1.0;
    eye[2][2] = 1.0;
    
    //-------------
    matrix n0(3, 1); //intial director vector at current gauss point
    n0(0,0) = cos(theta0);
    n0(1,0) = sin(theta0);
    n0(2,0) = 0.0;
    
    matrix n0Tn0_temp(3, 3);
    mat3d n0Tn0;
    n0Tn0_temp = n0*n0.transpose(); //dyadic product of n0 and n0 
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            n0Tn0(i, j) = n0Tn0_temp(i, j);
        }
    }
    
    temp1 = pow(a0, 1.0/3.0);
    temp2 = pow(a0, -1.0/6.0);
    mat3d G0 = n0Tn0 * (temp1-temp2)  + eye*temp2; //Initial Spontanious deformation tensor
   //mat3d G0 = eye;
    mat3d G02 = G0*G0;
    //-------------------
    
    matrix n(3, 1); //director vector at current gauss point
    n(0,0) = cos(theta);
    n(1,0) = sin(theta);
    n(2,0) = 0.0;
    
    matrix nTn_temp(3, 3);
    mat3d nTn;
    nTn_temp = n*n.transpose(); //dyadic product of n and n 
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            nTn(i, j) = nTn_temp(i, j);
        }
    }
      
    mat3d dnTndt; //derivative of nTn wrt theta
    dnTndt.zero();
    dnTndt[0][0] = -sin(theta*2.0 );
    dnTndt[0][1] = cos(theta * 2.0);
    dnTndt[1][0] = cos(theta * 2.0);
    dnTndt[1][1] = sin(theta * 2.0);
   
    temp1 = pow(a, 1.0/3.0);
    temp2 = pow(a, -1.0/6.0);
    mat3d G = nTn * (temp1-temp2)  + eye*temp2; //Spontanious deformation tensor
    //mat3d G = (eye + Qtensor*2.0) * pow(Qtensor.det(), -1.0 / 3.0);
    
    //temp1 = pow(m_a, -1.0/3.0);
    //temp2 = pow(m_a, 1.0/6.0);
    //matrix (3,3) Ginv = nTn * (temp1-temp2)  + eye*temp2; //inverse of G
    mat3d Ginv = G.inverse();
    
    //temp1 = pow(m_a, -2.0/3.0);
    //temp2 = pow(m_a, 2.0/6.0);
    //matrix (3,3) Ginv2 = nTn * (temp1-temp2)  + eye*temp2; //Ginv*Ginv
    
    mat3d Ginv2 = Ginv*Ginv;
    
    temp1 = pow(a, -2.0/3.0);
    temp2 = pow(a, 2.0/6.0);
    mat3d dGinv2dt = dnTndt * (temp1-temp2);
    //mat3d dGinv2dt = dnTndt * 0.0;

    //Derivative of J=detF
    mat3d dJdF = Finv.transpose() * detF;
    
    //--------First order derivatives of free energies-----------------
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F) 
    tempmat1 =  (Finv.transpose()*temp1)*(-1.0 / 3.0);
    tempmat2 = (Ginv2*F)*G02;
    mat3d delas2dF = (tempmat1 + tempmat2) * mu * pow(detF,(-2.0/3.0));

    //Semi-Soft elasticity, Pure neohookean
    temp1 = (F ).dotdot( F); //doublecontraction(F,F) 
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = F;
    mat3d delas1dF = (tempmat1 + tempmat2) * mu_neo * pow(detF, (-2.0 / 3.0));

    mat3d delasdF = delas1dF + delas2dF;

 
    //------
    matrix dfrankdGradt(3, 1);
    dfrankdGradt = Gradt * kfrank;
    
    //-----
    matrix delasdt(1, 1);
    tempmat1 = (G02*Fbar.transpose()) ;
    tempmat2 = (dGinv2dt * Fbar);
    tempmat3 = tempmat1 * tempmat2;
    temp1 = tempmat3.trace();
    delasdt(0, 0) = 0.5 * mu * temp1;

     
    
    //--------Calculation of Stresses---------------
    mat3d PK1 = delasdF + dJdF*pressmix;
    Pi0 = dfrankdGradt;

    ////--------Later Addition to check------///
    //TODO Testing only as gt = skewsymmetric part of cauchy stress
    // mat3d Cauchy = PK1*F.transpose();
    // matrix dndt(3, 1); //director vector at current gauss point
    // dndt(0,0) = -sin(theta);
    // dndt(1,0) = cos(theta);
    // dndt(2,0) = 0.0;
    // matrix nTdndt_temp(3,3);
    // nTdndt_temp = n*dndt.transpose();
    // mat3d nTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         nTdndt(i, j) = nTdndt_temp(i, j);
    //     }
    // }
    // temp1 = (Cauchy.transpose() - Cauchy).dotdot(nTdndt);
    // delasdt(0,0) = -temp1;
    ///--------Later addition to be end here------////

    gt = delasdt*(-1.0);



    
    //Get PK1 in voigt notation
    tens33_9(PK1,PK1_vgt);
    
    return PK1;
}



//==============================================================================================
//==============================================================================================
void FELCEMaterial444::TangentLCE(FESolidElement el, int ngauss, double theta, 
          double theta0, mat3d F, matrix Gradt, double pressmix, matrix &dPdF_vgt, 
          matrix &dPdt_vgt, matrix &dgdF_vgt, mat3d &dPidGradt, matrix &dgdt,
          matrix &dJdF_vgt)          
{
    //Declare some temporary double, matrix and tensors 
    double temp1,temp2,temp3,temp4,temp5; 
    mat3d tempmat1, tempmat2, tempmat3, tempmat4, tempmat5;
    tens4d tempA, tempB, tempC, tempD, tempE;
    
    

    //-----------------------------------------------------
    //Get material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(ngauss);
    FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

    //Get Current time
    double cTime = GetFEModel()->GetTime().currentTime;


    //GEt material parameter value
    //todo is writing (0) good? Change it for material point, here assumed parameter is same everywhere
    double mu = m_mu(mp);
    double mu_neo = 1.0*mu;
    //mu = ((1.0e-10 - 2.0e-2)/(1.20-1.0))*(cTime+1.0-1.20) + 1.0e-10;
    //if (cTime>0.20) {mu = 1.0e-10;}
    double kfrank = m_kfrank(mp);

    //Get deformation gradient and its determinant
    double detF = F.det();
    mat3d Finv = F.inverse();
    mat3d Fbar = F* pow(detF, -1.0 / 3.0);

    //double a = (1.0+2.0*Qorder)/(1.0-Qorder);
    double a = m_a(mp);

    double a0 = a;
    double bulkmod = m_bulk(mp);

    //--------------------------
    
    
    
    //-------------
    mat3d eye;
    eye.zero();
    eye[0][0] = 1.0;
    eye[1][1] = 1.0;
    eye[2][2] = 1.0;
    
    //-------------
    matrix n0(3, 1); //intial director vector at current gauss point
    n0(0,0) = cos(theta0);
    n0(1,0) = sin(theta0);
    n0(2,0) = 0.0;
    
    matrix n0Tn0_temp(3, 3);
    mat3d n0Tn0;
    n0Tn0_temp = n0*n0.transpose(); //dyadic product of n0 and n0 
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            n0Tn0(i, j) = n0Tn0_temp(i, j);
        }
    }
    
    temp1 = pow(a0, 1.0/3.0);
    temp2 = pow(a0, -1.0/6.0);
    mat3d G0 = n0Tn0 * (temp1-temp2)  + eye*temp2; //Initial Spontanious deformation tensor
    //mat3d G0 = eye;
    mat3d G02 = G0*G0;
    
    //-------------------
    
    matrix n(3, 1); //director vector at current gauss point
    n(0,0) = cos(theta);
    n(1,0) = sin(theta);
    n(2,0) = 0.0;
    
    matrix nTn_temp(3, 3);
    mat3d nTn;
    nTn_temp = n*n.transpose(); //dyadic product of n and n 
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            nTn(i, j) = nTn_temp(i, j);
        }
    }
      
    mat3d dnTndt; //derivative of nTn wrt theta
    dnTndt.zero();
    dnTndt[0][0] = -sin(theta*2.0 );
    dnTndt[0][1] = cos(theta * 2.0);
    dnTndt[1][0] = cos(theta * 2.0);
    dnTndt[1][1] = sin(theta * 2.0);
    
    mat3d ddnTnddt; //second derivative of nTn wrt theta   
    ddnTnddt.zero();
    ddnTnddt[0][0] = -2.0 * cos(theta * 2.0);
    ddnTnddt[0][1] = -2.0 * sin(theta * 2.0);
    ddnTnddt[1][0] = -2.0 * sin(theta * 2.0);
    ddnTnddt[1][1] = 2.0*cos(theta * 2.0);
   
    temp1 = pow(a, 1.0/3.0);
    temp2 = pow(a, -1.0/6.0);
    mat3d G = nTn * (temp1-temp2)  + eye*temp2; //Spontanious deformation tensor
    //mat3d G = (eye + Qtensor*2.0) * pow(Qtensor.det(), -1.0 / 3.0);

    //temp1 = pow(m_a, -1.0/3.0);
    //temp2 = pow(m_a, 1.0/6.0);
    //matrix (3,3) Ginv = nTn * (temp1-temp2)  + eye*temp2; //inverse of G
    mat3d Ginv = G.inverse();
    
    //temp1 = pow(m_a, -2.0/3.0);
    //temp2 = pow(m_a, 2.0/6.0);
    //matrix (3,3) Ginv2 = nTn * (temp1-temp2)  + eye*temp2; //Ginv*Ginv
    
    mat3d Ginv2 = Ginv*Ginv;
    
    temp1 = pow(a, -2.0/3.0);
    temp2 = pow(a, 2.0/6.0);
    mat3d dGinv2dt = dnTndt * (temp1-temp2);
    mat3d ddGinv2ddt = ddnTnddt * (temp1-temp2);

    
    //Derivative of J=detF
    mat3d dJdF = Finv.transpose() * detF;
    tempA = dyad1(Finv.transpose(),dJdF);
    tempB = dyad3(-Finv.transpose(),Finv);
    tens4d ddJddF = tempA + tempB*detF;
       
    //---------Second order derivate of free energies----------------
    //First, we need delasdF to calculate ddelasddf 
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F) 
    tempmat1 =  (Finv.transpose()*temp1)*(-1.0 / 3.0);
    tempmat2 = (Ginv2*F)*G02;
    mat3d delas2dF = (tempmat1 + tempmat2) * mu * pow(detF,(-2.0/3.0)); 

    //Semi-Soft elasticity, Pure neohookean
    temp1 = (F).dotdot(F); //doublecontraction(F,F) 
    tempmat1 = (Finv.transpose() * temp1) * (-1.0 / 3.0);
    tempmat2 = F;
    mat3d delas1dF = (tempmat1 + tempmat2) * mu_neo * pow(detF, (-2.0 / 3.0));


    mat3d delasdF = delas1dF + delas2dF;
    
    
    tempA = dyad2(Ginv2,G02.transpose());
    tempB= dyad1(Finv.transpose(),(Ginv2*(F*G02)))*2.0;
    temp1 = (F*G02).dotdot(Ginv2*F); //doublecontraction(F*G02,Ginv2*F);
    tempC = dyad3(-Finv.transpose(),Finv) * temp1;
    tempD = dyad1(delas2dF,Finv.transpose())*(2.0/3.0);
    tens4d ddelas2ddF = ((tempB + tempC)*(-1.0/3.0) + tempA)* 
                             mu*pow(detF,(-2.0/3.0)) - tempD;

    //Semi-Soft elasticity, Pure neohookean
    tempA = dyad2(eye, eye.transpose());
    tempB = dyad1(Finv.transpose(), F) * 2.0;
    temp1 = (F ).dotdot( F); //doublecontraction(F,F);
    tempC = dyad3(-Finv.transpose(), Finv) * temp1;
    tempD = dyad1(delas1dF, Finv.transpose()) * (2.0 / 3.0);
    tens4d ddelas1ddF = ((tempB + tempC) * (-1.0 / 3.0) + tempA) *
        mu_neo * pow(detF, (-2.0 / 3.0)) - tempD;
    
    tens4d ddelasddF = ddelas1ddF + ddelas2ddF;


                             
    //----
    temp1 = (F*G02).dotdot(dGinv2dt*F); //doublecontraction(F*G02,dGinv2dt*F);
    tempmat1 = Finv.transpose()*(1.0/3.0)*temp1;
    mat3d ddelasdFdt = ((dGinv2dt*(F*G02))- tempmat1)*mu*pow(detF,(-2.0/3.0));



    
    //----
    mat3d ddelasdtdF = ddelasdFdt;
    
    //-----
    mat3d ddfrankddGradt = eye*kfrank;
    
    //-----
    matrix ddelasdtdt(1, 1);
    tempmat1 = (G02*Fbar.transpose()) * (ddGinv2ddt*Fbar);
    temp1 = tempmat1.trace();
    ddelasdtdt(0, 0) = 0.5 * mu * temp1;

    
    //--------Calculation of Stress Tangents---------------
    tens4d dPdF = ddJddF*pressmix + ddelasddF;
    mat3d dPdt = ddelasdFdt;
    mat3d dgdF = -ddelasdtdF;
    dPidGradt = ddfrankddGradt;

    ////--------Later Addition to check------///
    //TODO Testing only as gt = skewsymmetric part of cauchy stress
    // mat3d PK1 = delasdF + dJdF*pressmix;
    // mat3d Cauchy = PK1*F.transpose();
    // mat3d dCauchydt = dPdt * F.transpose();
    // matrix dndt(3, 1); //director vector at current gauss point
    // dndt(0,0) = -sin(theta);
    // dndt(1,0) = cos(theta);
    // dndt(2,0) = 0.0;
    // matrix nTdndt_temp(3,3);
    // nTdndt_temp = n*dndt.transpose();
    // mat3d nTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         nTdndt(i, j) = nTdndt_temp(i, j);
    //     }
    // }
    // matrix dndtTdndt_temp(3,3);
    // dndtTdndt_temp = dndt*dndt.transpose();
    // mat3d dndtTdndt;
    //  for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < 3; ++j) {
    //         dndtTdndt(i, j) = dndtTdndt_temp(i, j);
    //     }
    // }
    // temp1 = (dCauchydt.transpose() - dCauchydt).dotdot(nTdndt);
    // temp2 = (Cauchy.transpose() - Cauchy).dotdot(dndtTdndt);
    // ddelasdtdt(0,0) = -temp1 - temp2;
    ///--------Later addition to be end here------////

    dgdt = ddelasdtdt*(-1.0);
    
    //------Transforming Tensor to voigt notation
    //matrix(9,9) dPdF_vgt = tens3333_99(dPdF);
    double dPdF_vgt_temp[9][9];
    dPdF.extract(dPdF_vgt_temp);
    matrix dPdF_vgt_temp2(9, 9);
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {

            dPdF_vgt_temp2(i, j) = dPdF_vgt_temp[i][j];
        }
    }
    dPdF_vgt = dPdF_vgt_temp2.transpose();
    
    
    tens33_9(dPdt,dPdt_vgt);
    tens33_9(dgdF,dgdF_vgt);
    tens33_9(dJdF,dJdF_vgt);
    
}

//----------------------------Function definations----------------------
//----------------------------------------------------------------------
void FELCEMaterial444::tens33_9(mat3d A, matrix& B)
{
    //Notice the voigt order is different than I used in matlab tens33_9

    B(0, 0) = A(0, 0);
    B(1, 0) = A(1, 1);
    B(2, 0) = A(2, 2);
    B(3, 0) = A(0, 1);
    B(4, 0) = A(1, 2);
    B(5, 0) = A(2, 0);
    B(6, 0) = A(1, 0);
    B(7, 0) = A(2, 1);
    B(8, 0) = A(0, 2);
}
   
