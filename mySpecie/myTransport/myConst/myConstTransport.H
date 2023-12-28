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
    Foam::myConstTransport

Group
    grpSpecieTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    myConstTransportI.H
    myConstTransport.C

\*---------------------------------------------------------------------------*/

#ifndef myConstTransport_H
#define myConstTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class myConstTransport;

template<class Thermo>
inline myConstTransport<Thermo> operator+
(
    const myConstTransport<Thermo>&,
    const myConstTransport<Thermo>&
);

template<class Thermo>
inline myConstTransport<Thermo> operator*
(
    const scalar,
    const myConstTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const myConstTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class myConstTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class myConstTransport
:
    public Thermo
{
    // Private Data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number [-]
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline myConstTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct as named copy
        inline myConstTransport(const word&, const myConstTransport&);

        //- Construct from dictionary
        explicit myConstTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<myConstTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<myConstTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "const<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const myConstTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend myConstTransport operator+ <Thermo>
        (
            const myConstTransport&,
            const myConstTransport&
        );

        friend myConstTransport operator* <Thermo>
        (
            const scalar,
            const myConstTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const myConstTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myConstTransportI.H"

#ifdef NoRepository
    #include "myConstTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
