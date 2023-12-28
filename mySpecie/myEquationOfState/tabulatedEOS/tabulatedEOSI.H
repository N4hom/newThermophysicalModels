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
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::tabulatedEOS<Specie>::tabulatedEOS(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::tabulatedEOS<Specie>::tabulatedEOS
(
    const word& name,
    const tabulatedEOS<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::tabulatedEOS<Specie>>
Foam::tabulatedEOS<Specie>::clone() const
{
    return autoPtr<tabulatedEOS<Specie>>::New(*this);
}


template<class Specie>
inline Foam::autoPtr<Foam::tabulatedEOS<Specie>> Foam::tabulatedEOS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedEOS<Specie>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::rho(scalar p, scalar T) const
{
    Info << "Density from lookup " << rhoTable_.lookup(p,T) << endl;
    return rhoTable_.lookup(p,T);
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::H(scalar p, scalar T) const
{
    
    return eTable_.lookup(p,T) + p/rhoTable_.lookup(p,T);
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::Cp(scalar p, scalar T) const
{
    return cpTable_.lookup(p,T);
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::E(scalar p, scalar T) const
{
    return eTable_.lookup(p,T);
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::Cv(scalar p, scalar T) const
{
    return cvTable_.lookup(p,T);
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::S(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::psi(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::Z(scalar p, scalar T) const
{
    return p/rho(p,T)/this->R()/T;
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::CpMCv(scalar p, scalar T) const
{
    return Cp(p,T) - Cv(p,T);
}



//////////////////////////////

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::cSqr
(
    const scalar p,
    const scalar T
) const
{
    return cSqrTable_.lookup(p,T);
}

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::cSqr
(
    const scalar rho,
    const scalar e
) const
{
    return pTable_.dFdX(rho, e) + Gamma(rho, e, T)*p/max(rho, 1e-10);
}

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::p
(
    const scalar rho,
    const scalar T
) const
{
    Info << "reverse lookup in rhoTable_ " << endl;
    return rhoTable_.reverseLookupX(rho, T);
}

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::p
(
    const scalar rho,
    const scalar e
) const
{
    Info << "Lookup in pTable_ " << endl;
    return pTable_.lookup(rho, e);
}




// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::tabulatedEOS<Specie>::operator+=(const tabulatedEOS<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::tabulatedEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::tabulatedEOS<Specie> Foam::operator+
(
    const tabulatedEOS<Specie>& pg1,
    const tabulatedEOS<Specie>& pg2
)
{
    return tabulatedEOS<Specie>
    (
        static_cast<const Specie&>(pg1) + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::tabulatedEOS<Specie> Foam::operator*
(
    const scalar s,
    const tabulatedEOS<Specie>& pg
)
{
    return tabulatedEOS<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::tabulatedEOS<Specie> Foam::operator==
(
    const tabulatedEOS<Specie>& pg1,
    const tabulatedEOS<Specie>& pg2
)
{
    return tabulatedEOS<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
    );
}


// ************************************************************************* //
