{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = scalar(1) - alphaf;

    //volScalarField rUaA = 1.0/UaEqn().A();
    //volScalarField rUbA = 1.0/UbEqn().A();
 
    volScalarField rUaA = rAUa;
    volScalarField rUbA = rAUb; 

    phia == (fvc::interpolate(Ua) & mesh.Sf());
    phib == (fvc::interpolate(Ub) & mesh.Sf());

    rUaAf = rAUaf; //fvc::interpolate(rUaA);
    surfaceScalarField rUbAf = rAUbf; //fvc::interpolate(rUbA);

    //Ua = rUaA*UaEqn().H();
    volVectorField HbyA("HbyA", Ua);    
    HbyA = rUaA*UaEqn().H();
    
    volVectorField HbyB("HbyB", Ub);    
    //Ub = rUbA*UbEqn().H();
    HbyB = rUbA*UbEqn().H();
    
/*     // Clearing temporary momentum predictor objects
    if (nCorr <= 1)
    {
	UaEqn.clear();
	UbEqn.clear();
    } */

/*     surfaceScalarField phiDraga =
        fvc::interpolate(beta/rhoa*K*rUaA)*phib + rUaAf*(g & mesh.Sf());

    //_AO_09/08/2014 - if (g0.value() > 0.0 || kineticTheory.on())
    if (g0.value() > 0.0 || kineticTheory.on() || kineticTheory.PU())
    {
        phiDraga -= ppMagf*fvc::snGrad(alpha)*mesh.magSf();
 	//phiDraga -= rUaAf*fvc::snGrad(kineticTheory.pa()/rhoa)*mesh.magSf();
    } */

//     if (kineticTheory.on())
//     {
//         phiDraga -= rUaAf*fvc::snGrad(kineticTheory.pa()/rhoa)*mesh.magSf();
//     }

/*     surfaceScalarField phiDragb =
        fvc::interpolate(alpha/rhob*K*rUbA)*phia + rUbAf*(g & mesh.Sf()); */

    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            //phiDraga.boundaryField()[patchi] = 0.0;
            //phiDragb.boundaryField()[patchi] = 0.0;
            phiaSources.boundaryField()[patchi] = 0.0;
            phibSources.boundaryField()[patchi] = 0.0;        
	}
    }

    //phia = (fvc::interpolate(Ua) & mesh.Sf()) + fvc::ddtPhiCorr(rUaA, Ua, phia)
         //+ phiDraga;
    phia = (fvc::interpolate(HbyA) & mesh.Sf()) + fvc::ddtPhiCorr(rUaA, Ua, phia)
         + phiaSources;	
	  
    //phib = (fvc::interpolate(Ub) & mesh.Sf()) + fvc::ddtPhiCorr(rUbA, Ub, phib)
         //+ phiDragb;
    phib = (fvc::interpolate(HbyB) & mesh.Sf()) + fvc::ddtPhiCorr(rUbA, Ub, phib)
	 + phibSources;
	 
    phi = alphaf*phia + betaf*phib;

    surfaceScalarField Dp("(rho*(1|A(U)))", alphaf*rUaAf/rhoa + betaf*rUbAf/rhob);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
	
        pEqn.solve
	(
	    mesh.solver
	    (
		// Selecting linear solver
		p.select
		(
		     finalIter		 	// Final outer corrector
		  && corr == nCorr-1		// Final inner corrector
		  && nonOrth == nNonOrthCorr	// Final non-orhogonal corrector
		)
	    )
	);

        if (nonOrth == nNonOrthCorr)
        {
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phia -= rUaAf*SfGradp/rhoa;
            phib -= rUbAf*SfGradp/rhob;
            phi = alphaf*phia + betaf*phib;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

   	    //Ua += fvc::reconstruct(phiDraga - rUaAf*SfGradp/rhoa);
            Ua = HbyA + fvc::reconstruct(phiaSources - rUaAf*SfGradp/rhoa);
	    Ua.correctBoundaryConditions();

            //Ub += fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob);
            Ub = HbyB + fvc::reconstruct(phibSources - rUbAf*SfGradp/rhob);
	    Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

#include "continuityErrs.H"
