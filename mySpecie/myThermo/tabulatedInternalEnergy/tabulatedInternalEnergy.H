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

            scalar Cpv
            (
                const Thermo& thermo,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                //Info << "Cpv " << thermo.Cv(rho, e, T) << endl;
                return thermo.Cv(rho, e, T);
            }


            


            //- Ratio of specific heats Cp/Cv []
            // scalar CpByCpv
            // (
            //     const Thermo& thermo,
            //     const scalar p,
            //     const scalar T
            // ) const
            // {
            //     #ifdef __clang__
            //     // Using volatile to prevent compiler optimisations leading to
            //     // a sigfpe
            //     volatile const scalar gamma = thermo.gamma(p, T);
            //     return gamma;
            //     #else
            //     return thermo.gamma(p, T);
            //     #endif
            // }

            scalar CpByCpv
            (
                const Thermo& thermo,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                
                return thermo.gamma(rho, e, T);
            }

            //- Sensible internal energy [J/kg]
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
                volatile const scalar es = thermo.Es(p, T);
                return es;
                #else
                return thermo.Es(p, T);
                #endif
            }



             //- Sensible internal energy (rho,T) [J/kg]
            scalar HErhoT
            (
                const Thermo& thermo,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                //Info << "In HErhoT " << endl;
                //Info << "rho = " << rho << endl;
                return thermo.EsRhoT(rho, e, T);
            }

            //- Temperature from sensible internal energy given an initial
            //- temperature T0 [K]
            scalar THE
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar tes = thermo.TEs(e, p, T0);
                return tes;
                #else
                return thermo.TEs(e, p, T0);
                #endif
            }

            //- Temperature from sensible internal energy given an initial
            //- temperature T0 [K]
            scalar TRhoE
            (
                const Thermo& thermo,
                const scalar T,
                const scalar rho,
                const scalar e
            ) const
            {
                
                return thermo.tabulatedTRhoEs(rho, e, T);
                
            }



            scalar gamma
            (
                const Thermo& thermo,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const
            {
                
                return CpByCpv(thermo, rho, e, T);
                
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
