#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEModelParam.h>
#include <FECore/FENode.h>



//-----------------------------------------------------------------------------
//! LCE material with Neo-classical + frank energy

class FEBIOMECH_API FELCEMaterial444 : public FEElasticMaterial
{


public:
	FELCEMaterial444(FEModel* pfem);

public:
	FEParamDouble		m_mu;	//!< Shear modulus
	FEParamDouble		m_kfrank;	//!< frank constant
    FEParamDouble		m_a;        //!< stretch parameter
    FEParamDouble		m_bulk;     //!< bulk modulus
    
public:
    //! calculate stress at material point
    virtual mat3ds Stress(FEMaterialPoint& pt) override;

    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! Calculate LCE stresses PK1_vgt, Pi0, gt and returns mat3d PK1
    mat3d StressLCE(FESolidElement el, int ngauss, double theta, double theta0, 
        mat3d F, matrix Gradt, double pressmix, matrix &PK1_vgt, matrix &Pi0,
        matrix &gt);
            
    //! Calculate all the tangents needed for the LCE stiffness calculation
    void TangentLCE(FESolidElement el, int ngauss, double theta, double theta0,
          mat3d F, matrix Gradt, double pressmix, matrix &dPdF_vgt, 
          matrix &dPdt_vgt, matrix &dgdF_vgt, mat3d &dPidGradt, matrix &dgdt,
          matrix &dJdF_vgt);

    //Change 3by3 tensor to 9by1 voigt matrix
    void tens33_9(mat3d A, matrix &B);
               
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};    