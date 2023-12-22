/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::myAtomicWeightTable

Description
    A table of myAtomic weights for all the elements

SourceFiles
    myAtomicWeights.C

\*---------------------------------------------------------------------------*/

#ifndef myAtomicWeights_H
#define myAtomicWeights_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class myAtomicWeights Declaration
\*---------------------------------------------------------------------------*/

class myAtomicWeightTable
:
    public HashTable<scalar>
{

public:

    // Public types and data

        //- Structure to hold the element name and myAtomic weight pair
        struct myAtomicWeight
        {
            char name[3];
            scalar weight;
        };

        static const int nElements = 104;

        //- Static table of the weights of all known elements
        static const myAtomicWeight myAtomicWeights[nElements];


    // Constructors

        //- Construct from myAtomicWeights_
        myAtomicWeightTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// myAtomic weights table for every element in the periodic table
extern myAtomicWeightTable myAtomicWeights;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
