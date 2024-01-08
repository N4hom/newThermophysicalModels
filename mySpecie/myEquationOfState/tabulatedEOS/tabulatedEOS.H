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
    Foam::tabulatedEOS

Group
    grpSpecieEquationOfState

Description
    Perfect gas equation of state.

SourceFiles
    tabulatedEOSI.H
    tabulatedEOS.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedEOS_H
#define tabulatedEOS_H

#include "autoPtr.H"
#include "lookupTables2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class tabulatedEOS;

template<class Specie>
inline tabulatedEOS<Specie> operator+
(
    const tabulatedEOS<Specie>&,
    const tabulatedEOS<Specie>&
);

template<class Specie>
inline tabulatedEOS<Specie> operator*
(
    const scalar,
    const tabulatedEOS<Specie>&
);

template<class Specie>
inline tabulatedEOS<Specie> operator==
(
    const tabulatedEOS<Specie>&,
    const tabulatedEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const tabulatedEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class tabulatedEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabulatedEOS
:
    public Specie
{
protected:
    // Density lookup table
    scalarLookupTable2D pTable_;

    // Internal energy lookup table
    scalarLookupTable2D TTable_;

/*protected:
    scalarLookupTable2D cvTable_;

    scalarLookupTable2D cpTable_;

private:
    scalarLookupTable2D cSqrTable_;
*/
public:

    // Constructors

        //- Construct from components
        inline tabulatedEOS(const Specie& sp);

        //- Construct from dictionary
        tabulatedEOS(const dictionary& dict);

        //- Construct as named copy
        inline tabulatedEOS(const word& name, const tabulatedEOS&);

        //- Construct and return a clone
        inline autoPtr<tabulatedEOS> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedEOS> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulatedEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;


            //- Return pressure from lookup table
            inline scalar p(const scalar rho, const scalar e, const scalar T) const;

            //- Return speed of sound from lookup table
            inline scalar cSqr(const scalar rho, const scalar e, const scalar T, const scalar cv) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(const scalar rho, const scalar e, const scalar T) const;

            inline scalar Cv(const scalar rho, const scalar e, const scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar rho, const scalar e, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(const scalar rho, const scalar e, scalar T) const;

            inline scalar gamma(const scalar rho, const scalar e, const scalar T, const scalar cv) const;

            //- Return derivative of pressure w.r.t. Temperature
            inline scalar dpdT( const scalar rho, const scalar e, const scalar T) const;

            //- Return derivative of pressure w.r.t. internal energy
            inline scalar dpde( const scalar rho, const scalar e, const scalar T) const;

            //- Return derivative of pressure w.r.t. specific volume
            inline scalar dpdv
            (
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;
            
            


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const tabulatedEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend tabulatedEOS operator+ <Specie>
        (
            const tabulatedEOS&,
            const tabulatedEOS&
        );

        friend tabulatedEOS operator* <Specie>
        (
            const scalar s,
            const tabulatedEOS&
        );

        friend tabulatedEOS operator== <Specie>
        (
            const tabulatedEOS&,
            const tabulatedEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const tabulatedEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedEOSI.H"

#ifdef NoRepository
    #include "tabulatedEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
