/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::mySpecie

Description
    Base class of the thermophysical property types.

SourceFiles
    mySpecieI.H
    mySpecie.C

\*---------------------------------------------------------------------------*/

#ifndef mySpecie_H
#define mySpecie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mySpecie;

inline mySpecie operator+(const mySpecie&, const mySpecie&);
inline mySpecie operator*(const scalar, const mySpecie&);
inline mySpecie operator==(const mySpecie&, const mySpecie&);

Ostream& operator<<(Ostream&, const mySpecie&);


/*---------------------------------------------------------------------------*\
                           Class mySpecie Declaration
\*---------------------------------------------------------------------------*/

class mySpecie
{
    // Private Data

        //- Name of mySpecie
        word name_;

        //- Number of moles of this component in the mixture
        scalar Y_;

        //- Molecular weight of mySpecie [kg/kmol]
        scalar molWeight_;


public:

    //- Declare type-name (with debug switch)
    ClassName("mySpecie");


    // Generated Methods

        //- Copy construct
        mySpecie(const mySpecie&) = default;


    // Constructors

        //- Construct from components without name
        inline mySpecie(const scalar Y, const scalar molWeight);

        //- Construct from components with name
        inline mySpecie
        (
            const word& name,
            const scalar Y,
            const scalar molWeight
        );

        //- Construct as named copy
        inline mySpecie(const word& name, const mySpecie&);

        //- Construct from dictionary
        explicit mySpecie(const dictionary& dict);


    // Member Functions

        // Access

            //- Name
            inline const word& name() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- No of moles of this mySpecies in mixture
            inline scalar Y() const;

            //- Gas constant [J/(kg K)]
            inline scalar R() const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment, preserve original name
        inline void operator=(const mySpecie&);

        inline void operator+=(const mySpecie&);
        inline void operator*=(const scalar);


    // Friend Operators

        inline friend mySpecie operator+(const mySpecie&, const mySpecie&);
        inline friend mySpecie operator*(const scalar, const mySpecie&);
        inline friend mySpecie operator==(const mySpecie&, const mySpecie&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const mySpecie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mySpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
