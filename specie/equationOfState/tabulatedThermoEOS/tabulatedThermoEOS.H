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
    Foam::tabulatedThermoEOS

Group
    grpSpecieEquationOfState

Description
    Perfect gas equation of state.

SourceFiles
    tabulatedThermoEOSI.H
    tabulatedThermoEOS.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedThermoEOS_H
#define tabulatedThermoEOS_H

#include "autoPtr.H"
#include "lookupTables2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class tabulatedThermoEOS;

template<class Specie>
inline tabulatedThermoEOS<Specie> operator+
(
    const tabulatedThermoEOS<Specie>&,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
inline tabulatedThermoEOS<Specie> operator*
(
    const scalar,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
inline tabulatedThermoEOS<Specie> operator==
(
    const tabulatedThermoEOS<Specie>&,
    const tabulatedThermoEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const tabulatedThermoEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class tabulatedThermoEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabulatedThermoEOS
:
    public Specie
{

    // Density lookup table
    scalarLookupTable2D rhoTable_;

    // Internal energy lookup table
    scalarLookupTable2D eTable_;

protected:
    scalarLookupTable2D cvTable_;

    scalarLookupTable2D cpTable_;

private:
    scalarLookupTable2D cSqrTable_;

public:

    // Constructors

        //- Construct from components
        inline tabulatedThermoEOS(const Specie& sp);

        //- Construct from dictionary
        tabulatedThermoEOS(const dictionary& dict);

        //- Construct as named copy
        inline tabulatedThermoEOS(const word& name, const tabulatedThermoEOS&);

        //- Construct and return a clone
        inline autoPtr<tabulatedThermoEOS> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedThermoEOS> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulatedThermoEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

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


            /////////////////////////////////////
            //- Return speed of sound
            inline scalar cSqr
            (
                const scalar p,
                const scalar T
            ) const;

            //- Return pressure
            inline scalar p
            (
                const scalar rho,
                const scalar T
            ) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const tabulatedThermoEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedThermoEOS operator+ <Specie>
        (
            const tabulatedThermoEOS&,
            const tabulatedThermoEOS&
        );

        friend tabulatedThermoEOS operator* <Specie>
        (
            const scalar s,
            const tabulatedThermoEOS&
        );

        friend tabulatedThermoEOS operator== <Specie>
        (
            const tabulatedThermoEOS&,
            const tabulatedThermoEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const tabulatedThermoEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedThermoEOSI.H"

#ifdef NoRepository
    #include "tabulatedThermoEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
