     Info<< "kinTheory: min(g0) = " << min(gs0_).value()
         << ", max(g0) = "          << max(gs0_).value() << endl;

     Info<< "kinTheory: min(PsCoeff) = " << min(PsCoeff).value()
         << ", max(PsCoeff) = "          << max(PsCoeff).value() << endl;
	 	 
     Info<< "kinTheory: min(mua) = " << min(mua_).value()
         << ", max(mua) = "          << max(mua_).value() << endl; 
         
     Info<< "kinTheory: min(lamdba) = " << min(lambda_).value()
         << ", max(lamdba) = "          << max(lambda_).value() << endl;

     Info<< "kinTheory: min(J) = " << min(J_).value()
         << ", max(J) = "          << max(J_).value() << endl;

     Info<< "kinTheory: min(K) = " << min(K_).value()
         << ", max(K) = "          << max(K_).value() << endl;

     Info<< "kinTheory: min(Kmod) = " << min(Kmod_).value()
         << ", max(Kmod) = "          << max(Kmod_).value() << endl;

     Info<< "kinTheory: min(M) = " << min(M_).value()
         << ", max(M) = "          << max(M_).value() << endl;
     
     Info<< "kinTheory: min(D) = " << min(D).value()
         << ", max(D) = "          << max(D).value() << endl;

     Info<< "kinTheory: min(gammadot) = " << min(gammaDot).value()
         << ", max(gammadot) = "          << max(gammaDot).value() << endl;

     Info<< "kinTheory: min(ThetaDil) = "   << min(ThetaDil_).value()
         << ", max(ThetaDil) = "            << max(ThetaDil_).value() << endl;
	 
     Info<< "kinTheory: min(ThetaDense) = " << min(ThetaDense_).value()
      << ", max(ThetaDense) = "          << max(ThetaDense_).value() << endl;		 

     Info<< "kinTheory: psi = " << psi.value() << endl;

     Info<< "kinTheory: min(InertiaNumber) = " << min(inertiaNumber).value()
       << ", max(InertiaNumber) = "          << max(inertiaNumber).value() << endl;     
     
     Info<< "kinTheory: min(modInertiaNumber) = " << min(modInertiaNumber).value()
         << ", max(modInertiaNumber) = "          << max(modInertiaNumber).value() << endl;

     Info<< "kinTheory: min(chi) = " << min(chi).value()
         << ", max(chi) = "          << max(chi).value() << endl;

     Info<< "kinTheory: min(beta) = " << min(beta).value()
         << ", max(beta) = "          << max(beta).value() << endl;
       
     Info<< "kinTheory: min(sigmaTau) = " << min(sigmaTau).value()
         << ", max(sigmaTau) = "          << max(sigmaTau).value() << endl;
	
     Info<< "kinTheory: min(sigmaGamma) = " << min(sigmaGamma).value()
     	 << ", max(sigmaGamma) = "          << max(sigmaGamma).value() << endl;

     Info<< "kinTheory: min(gammaCoeff) = " << min(gammaCoeff).value()
         << ", max(gammaCoeff) = "          << max(gammaCoeff).value() << endl;

     Info<< "kinTheory: min(func_B) = " << min(func_B).value()
         << ", max(func_B) = "          << max(func_B).value() << endl;
	  	     
     Info<< "kinTheory: min(S) = " << min(S).value()
         << ", max(S) = "          << max(S).value() << endl;

     Info<< "kinTheory: min(hatS) = " << min(hatS).value()
         << ", max(hatS) = "          << max(hatS).value() << endl;
    
     volScalarField ktn(mua_/rhoa_);

     Info<< "kinTheory: min(nu1) = " << min(ktn).value()
         << ", max(nu1) = "          << max(ktn).value() << endl;

     Info<< "kinTheory: min(upsilon) = " << min(upsilon_).value()
         << ", max(upsilon) = "          << max(upsilon_).value() << endl;	

     Info<< "kinTheory: min(Theta) = "   << min(Theta_).value()
         << ", max(Theta) = "           << max(Theta_).value() << endl;

     Info<< "kinTheory: min(pa) = " << min(alpha_*pa_).value()
         << ", max(pa) = "          << max(alpha_*pa_).value() << endl; 
      
     Info<< "kinTheory: min(tau) = " << min(tau_).value()
         << ", max(tau) = "          << max(tau_).value() << endl;

     Info << " " << endl;
