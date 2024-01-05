/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::myPerfectGas

Group
    grpSpecieEquationOfState

Description
    Perfect gas equation of state.

SourceFiles
    myPerfectGasI.H
    myPerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef myPerfectGas_H
#define myPerfectGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class myPerfectGas;

template<class Specie>
inline myPerfectGas<Specie> operator+
(
    const myPerfectGas<Specie>&,
    const myPerfectGas<Specie>&
);

template<class Specie>
inline myPerfectGas<Specie> operator*
(
    const scalar,
    const myPerfectGas<Specie>&
);

template<class Specie>
inline myPerfectGas<Specie> operator==
(
    const myPerfectGas<Specie>&,
    const myPerfectGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const myPerfectGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class myPerfectGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class myPerfectGas
:
    public Specie
{

public:

    // Constructors

        //- Construct from components
        inline myPerfectGas(const Specie& sp);

        //- Construct from dictionary
        myPerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline myPerfectGas(const word& name, const myPerfectGas&);

        //- Construct and return a clone
        inline autoPtr<myPerfectGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<myPerfectGas> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "myPerfectGas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3] - Not needed anymore as the density is now the independent variable
            inline scalar rho(scalar p, scalar T) const;

            //- Return pressure []
            inline scalar p(scalar rho, scalar e, scalar T) const;

            //- Return speed of sound []
            inline scalar cSqr(scalar rho, scalar e, scalar T, scalar cv) const;

            inline scalar gamma(const scalar rho, const scalar e, const scalar T, const scalar cv) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const myPerfectGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend myPerfectGas operator+ <Specie>
        (
            const myPerfectGas&,
            const myPerfectGas&
        );

        friend myPerfectGas operator* <Specie>
        (
            const scalar s,
            const myPerfectGas&
        );

        friend myPerfectGas operator== <Specie>
        (
            const myPerfectGas&,
            const myPerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const myPerfectGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myPerfectGasI.H"

#ifdef NoRepository
    #include "myPerfectGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
