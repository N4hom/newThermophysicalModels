/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::tabulatedInternalEnergy

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the tabulated internal energy
    functions.

\*---------------------------------------------------------------------------*/

#ifndef tabulatedInternalEnergy_H
#define tabulatedInternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class tabulatedInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedInternalEnergy
{

public:

    // Constructors

        //- Construct
        tabulatedInternalEnergy()
        {
            Info << "Constructing " << typeName() <<  endl;
        }


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulatedInternalEnergy";
        }

        // Fundamental properties

            static word energyName()
            {
                return "ea";
            }

            // Heat capacity at constant volume [J/(kg K)]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar cv = thermo.Cv(p, T);
                return cv;
                #else
                return thermo.Cv(p, T);
                #endif
            }

            //- Cp/Cv []
            scalar CpByCpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar gamma = thermo.gamma(p, T);
                return gamma;
                #else
                return thermo.gamma(p, T);
                #endif
            }

            // Absolute internal energy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar ea = thermo.Ea(p, T);
                return ea;
                #else
                return thermo.Ea(p, T);
                #endif
            }

            //- Temperature from tabulated internal energy
            //  given an initial temperature T0
            scalar THE
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                Info << "Here " << endl;

                const scalar& energy(e);
                Info << "Internal energy " << energy << endl;
                Info <<" Temperature = " <<  thermo.eTable().reverseLookupY(energy,p) << endl;
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar tea = thermo.TEa(e, p, T0);
                return tea;
                #else

                return thermo.eTable().reverseLookupY(energy,p);
                #endif
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
