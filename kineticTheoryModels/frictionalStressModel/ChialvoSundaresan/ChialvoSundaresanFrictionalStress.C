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

\*---------------------------------------------------------------------------*/

#include "ChialvoSundaresanFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChialvoSundaresanFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        ChialvoSundaresanFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanFrictionalStress::ChialvoSundaresanFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanFrictionalStress::~ChialvoSundaresanFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return
        dimensionedScalar("1e24", dimensionSet(1, -1, -2, 0, 0), 0.083e5)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), 2/3);
}


Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return
        dimensionedScalar("1e25", dimensionSet(1, -1, -2, 0, 0), 0.055333e5)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), -1/3);
}


Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanFrictionalStress::muf
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const dimensionedScalar& phi
) const
{
    const scalar I2Dsmall = 1.0e-15;

    // Creating muf assuming it should be 0 on the boundary which may not be
    // true
    tmp<volScalarField> tmuf
    (
        new volScalarField
        (
            IOobject
            (
                "muf",
                alpha.mesh().time().timeName(),
                alpha.mesh()
            ),
            alpha.mesh(),
            dimensionedScalar("muf", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );

    volScalarField& muff = tmuf();

//    forAll (D, celli)
 //   {
   //     if (alpha[celli] > alphaMax.value() )
     //   {
       //     muff[celli] = pf * upsilon_ / max( gammaDot, gammaDotSmall );
//        }
  //  }

    return tmuf;
}


// ************************************************************************* //
