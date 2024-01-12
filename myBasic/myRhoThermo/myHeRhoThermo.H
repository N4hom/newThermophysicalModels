/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::myHeRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    myHeRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef myHeRhoThermo_H
#define myHeRhoThermo_H

#include "myRhoThermo.H"
#include "myHeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class myHeRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class myHeRhoThermo
:
    public myHeThermo<BasicPsiThermo, MixtureType>
{
    // Private data

    bool initialize_;
    // Private Member Functions

        //- Calculate the thermo variables
        // void calculate
        // (
        //     const volScalarField& p,
        //     volScalarField& T,
        //     volScalarField& he,
        //     volScalarField& psi,
        //     volScalarField& rho,
        //     volScalarField& mu,
        //     volScalarField& alpha,
        //     const bool doOldTimes
        // );

        //- Calculate the thermo variables
        void calculateFromRhoE();

        //- Construct as copy (not implemented)
        myHeRhoThermo(const myHeRhoThermo<BasicPsiThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("myHeRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        myHeRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        myHeRhoThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~myHeRhoThermo();


    // Member functions

        //- Update properties
        //virtual void correct();


        virtual void correctFromRhoE();
        //- Update properties from T
        //virtual void correctT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myHeRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
