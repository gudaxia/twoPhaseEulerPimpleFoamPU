/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phasea_;
        const volVectorField& Ua_;
        const volVectorField& Ub_;
        const volScalarField& alpha_;
        const surfaceScalarField& phia_;

        const dragModel& draga_;

        const dimensionedScalar& rhoa_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nua_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        //- use generation == dissipation
        Switch equilibrium_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        autoPtr<frictionalStressModel> frictionalStressModel_;

        //- coefficient of restitution
        const dimensionedScalar e_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- min value for which the frictional stresses are zero
        const dimensionedScalar alphaMinFriction_;

        //- material constant for frictional normal stress
        const dimensionedScalar Fr_;

        //- material constant for frictional normal stress
        const dimensionedScalar eta_;

        //- material constant for frictional normal stress
        const dimensionedScalar p_;

        //- angle of internal friction
        const dimensionedScalar phi_;

        //- The granular energy/temperature
        volScalarField Theta_;

        //- The granular viscosity
        volScalarField mua_;

        //- The granular bulk viscosity
        volScalarField lambda_;

        //- The granular pressure
        volScalarField pa_;
	
	//- The granular pressure derivative
        volScalarField ppMagf_;

        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;
	
	//- The derivative of the radial distribution function
        volScalarField gs0Prime_;
        
	//_AO_09/02/2014 - The shear stress ratio
        volScalarField upsilon_;
	
	//_AO_09/02/2014 - The shear stress 
	volSymmTensorField tau_;

        //_AO_09/02/2014 - effective restitution coefficient
        mutable dimensionedScalar eEff_;
	
        //_AO_09/02/2014 - friction coefficient
        mutable scalar muFric_;
		
        //_AO_09/02/2014 - max vol.frac for inertial regime
        //mutable dimensionedScalar alphaf_;
        mutable scalar alphaf_;

        //_AO_09/02/2014 - critical volume fraction
        mutable scalar alphac_;

         // YG 12/27/2014 - alphad_
        mutable scalar alphad_;

        //_AO_09/02/2014 - yield stress ratio
        mutable scalar upsilons_;
	
        //_AO_09/30/2014 - use modified kinetic theory by Chialvo-Sundaresan
        Switch mofidiedKineticTheoryPU_;
	
        //-AO, YG - Decompose particle pressure, Sundar's idea
        Switch decomposePp_;
	
	
	//-AO, YG - Berzi Model using SKT
	Switch Berzi_;

        //-AO, YG - Decompose particle pressure, Sundar's idea
        volScalarField paStar_;			

    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phasea,
            const volVectorField& Ub,
            const volScalarField& alpha,
            const dragModel& draga
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions	
	void solve(const volTensorField& gradUat);

        bool on() const
        {
            return kineticTheory_;
        }

        const volScalarField& mua() const
        {
            return mua_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }
	
 	//_AO_09/30/2014 modified kinetic theory		
        bool PU() const
        {
            return mofidiedKineticTheoryPU_;
        }

        //-AO, YG - Decompose particle pressure, Sundar's idea
        bool decomposePp() const
        {
            return decomposePp_;
        }
	
	
	bool Berzi() const
        {
            return Berzi_;
        }
	
	//_AO_09/02/2014_
        const volSymmTensorField& tau() const
        //const volTensorField& tau() const
	{
            return tau_; // Convert symmTensor to Tensor, used in UEqns.H
	}

	//_AO_09/02/2014 - 
	tmp<volScalarField> ppMagfCalc(const volScalarField& alphaUpdate);
	
	// There is a problem here !!! Herein, really stupid solution
	const volScalarField& ppMagf(const volScalarField& alphaUpdate)
	{
	  // Keep ppMagf previous value
	  ppMagf_.storePrevIter();
	  
	  const scalar constSMALL = 1.e-06;

	  if(mofidiedKineticTheoryPU_)
	  {
	  // Read from dictionary
	   alphaf_ = readScalar(kineticTheoryProperties_.lookup("alphaDiluteInertialUpperLimit"));
	   alphac_ = readScalar(kineticTheoryProperties_.lookup("alphaCritical"));
            alphad_ = readScalar(kineticTheoryProperties_.lookup("alphaDelta"));
	  }else
	  {
	   alphaf_ = alphaMax_.value();
	   alphac_ = alphaMax_.value();
           alphad_ = 0.01;
	  }
	  volScalarField alpha = alphaUpdate;

	  //gs0_ = radialModel_->g0(min(alpha, alphaMinFriction_), alphaMax_); 
	  //gs0Prime_ = radialModel_->g0prime(min(alpha, alphaMinFriction_), alphaMax_);

	  // 
	  gs0_ = radialModel_->g0jamming
	  (
	      Ua_.mesh(),
	     //max(alpha, scalar(constSMALL)),
              min(max(alpha_, scalar(constSMALL)), alphaMax_ - 0.01),
	      alphaMax_,
	      alphad_,  //changed by YG
	      alphac_    
	  );

     	  Info << "kinTheory: min(g0) = " << min(gs0_).value()
               << ", max(g0) = "          << max(gs0_).value() << endl;

	  gs0Prime_ = radialModel_->g0jammingPrime
	  (
	      Ua_.mesh(),
	     //max(alpha, scalar(constSMALL)),
              min(max(alpha_, scalar(constSMALL)), alphaMax_ - 0.01),
	      alphaMax_,
	      alphad_,  //changed by YG
	      alphac_
	  );	

     	  Info << "kinTheory: min(g0Prime) = " << min(gs0Prime_).value()
               << ", max(g0Prime) = "          << max(gs0Prime_).value() << endl;

	  volScalarField ppMagfFriction = frictionalStressModel_->frictionalPressurePrime
	  (
	      alpha, 
	      alphaMinFriction_, 
	      alphaMax_,
	      Fr_,
	      eta_,
	      p_
	  );

	  // Computing ppMagf
	  ppMagf_ = Theta_*granularPressureModel_->granularPressureCoeffPrime
	  (
	      alpha, 
	      gs0_, 
	      gs0Prime_, 
	      rhoa_, 
	      e_
	  );  

	  // NOTE: this might not be appropriate if J&J model is used (verify)
	  


	  forAll(alpha, cellI)
	  {
	      if(alpha[cellI] >= alphaMinFriction_.value())
	      {
		  ppMagf_[cellI] = ppMagf_[cellI] + ppMagfFriction[cellI];
	      }
	  }


	 // ppMagf_.relax();
	  
	  ppMagf_.correctBoundaryConditions();	  

	  Info<< "kinTheory: min(ppMagf) = " << min(ppMagf_).value()
	      << ", max(ppMagf) = "          << max(ppMagf_).value() << endl;         

	if(decomposePp_)
	{
		 // particle pressure star 
		 volScalarField PsCoeffStar	// -> rho_p * H 
		 (
        	     granularPressureModel_->granularPressureCoeff
        	     (
        		 alpha_,
        		 gs0_,
        		 rhoa_,
        		 e_
        	     )
		 );    

		 PsCoeffStar.max(1.0e-15);

		 // Particle pressure star
        	 paStar_ = PsCoeffStar * Theta_;

		 //paStar_ /= (fvc::average(alpha_) + scalar(0.001));
        	 paStar_ /= (alpha_ + scalar(0.001));
		 //paStar_.correctBoundaryConditions();
		 paStar_.max(1.0e-15);	  

		 Info<< "kinTheory: min(ppStar) = " << min(paStar_).value()
		     << ", max(ppStar) = "          << max(paStar_).value() << endl;  
	}

	  return ppMagf_;
	}
		
        const volScalarField& paStar() const
        {
            return paStar_;
        }	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
