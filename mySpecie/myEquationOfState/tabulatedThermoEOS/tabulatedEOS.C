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

\*---------------------------------------------------------------------------*/

#include "tabulatedEOS.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::tabulatedEOS<Specie>::tabulatedEOS(const dictionary& dict)
:
    Specie(dict),
    rhoTable_(dict.subDict("equationOfState"), "p", "T", "rho"),
    eTable_(dict.subDict("thermodynamics"), "p", "T", "e"),
    cvTable_(dict.subDict("thermodynamics"), "p", "T", "Cv"),
    cpTable_(dict.subDict("thermodynamics"), "p", "T", "Cp"),
    cSqrTable_(dict.subDict("thermodynamics"), "p", "T", "cSqr")
{
    Info << "Constructing equationOfState tabulatedEOS " << endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::tabulatedEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const tabulatedEOS<Specie>& pg)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
