/*---------------------------------------------------------------------------* \
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
    Foam::eTabulatedThermo

Group
    grpSpecieThermo

Description
    Constant properties thermodynamics package templated on an equation of
    state

SourceFiles
    eTabulatedThermoI.H
    eTabulatedThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eTabulatedThermo_H
#define eTabulatedThermo_H

#include "thermo.H"
#include "lookupTables2D.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class eTabulatedThermo;

/*template<class EquationOfState>
inline eTabulatedThermo<EquationOfState> operator+
(
    const eTabulatedThermo<EquationOfState>&,
    const eTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline eTabulatedThermo<EquationOfState> operator*
(
    const scalar,
    const eTabulatedThermo<EquationOfState>&
);

template<class EquationOfState>
inline eTabulatedThermo<EquationOfState> operator==
(
    const eTabulatedThermo<EquationOfState>&,
    const eTabulatedThermo<EquationOfState>&
);*/

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const eTabulatedThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class eTabulatedThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eTabulatedThermo
:
    public EquationOfState
{
    // Private data

        /*//- Heat capacity at constant volume
        //  Note: input in [J/(kg K)], but internally uses [J/(kmol K)]
        scalar Cv_;

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;

        //- Reference temperature around which to linearise [K]
        scalar Tref_;

        //- Reference sensible enthalpy around which to linearise [J/kg]
        scalar Esref_;*/

        /////////////////
        //- cvTable_ is inherited from <EquationOfState>
       //scalarLookupTable2D cvTable_;
protected:

        scalarLookupTable2D eTable_;

        //- Temperature bounds
        scalar Tlow_;
        scalar Thigh_;

        //- Heat of formation
        //scalar Hf_;
       //scalar banana;

private:
    // Private Member Functions

        //- Construct from components
        inline eTabulatedThermo
        (
            const EquationOfState& st,
            const scalarLookupTable2D& cvTable,
            const scalarLookupTable2D& eTable,
            const scalar Hf
        );


public:

    // Constructors

        //- Construct from dictionary
        eTabulatedThermo(const dictionary& dict);

        /*//- Construct as named copy
        inline eTabulatedThermo(const word&, const eTabulatedThermo&);*/

        //- Construct and return a clone
        inline autoPtr<eTabulatedThermo> clone() const;

        // Selector from dictionary
        inline static autoPtr<eTabulatedThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "eTabulated<" + EquationOfState::typeName() + '>';
        }

        //- Limit temperature to be within the range
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant volume [J/(kg K)]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return table
            inline scalarLookupTable2D eTable() const
            {
                return eTable_;
            }

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            //inline scalar Gstd(const scalar T) const;

            #include "EtoHthermo.H"


        // Derivative term used for Jacobian


            //- Temperature derivative of heat capacity at constant pressure
            //inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

      //  inline void operator+=(const eTabulatedThermo&);


    // Friend operators

       /* friend eTabulatedThermo operator+ <EquationOfState>
        (
            const eTabulatedThermo&,
            const eTabulatedThermo&
        );

        friend eTabulatedThermo operator* <EquationOfState>
        (
            const scalar,
            const eTabulatedThermo&
        );

        friend eTabulatedThermo operator== <EquationOfState>
        (
            const eTabulatedThermo&,
            const eTabulatedThermo&
        );*/


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const eTabulatedThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eTabulatedThermoI.H"

#ifdef NoRepository
    #include "eTabulatedThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
