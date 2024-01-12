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
Foam::scalar Foam::tabulatedEOS<Specie>::p(const scalar rho, const scalar e, const scalar T) const
{
    return pTable_.lookup(rho, e);
}

template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::rho(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::H(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::Cp(scalar rho, scalar e, scalar T) const
{
    NotImplemented;
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::Cv(scalar rho, scalar e, scalar T) const
{
    NotImplemented;
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::E(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}




template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::S(scalar rho, scalar e, scalar T) const
{
    return 1.0/max(TTable_.dFdY(rho, e), SMALL);
   
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
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::tabulatedEOS<Specie>::CpMCv(scalar rho, scalar e, scalar T) const
{
    return dpdT(rho, e, T)/max(rho, 1e-10);
}



//////////////////////////////


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpde
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return pTable_.dFdY(rho, e);
}

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::gamma
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return max(dpdT(rho, e, T)/(cv*max(rho, 1e-10)), SMALL) + 1.0;
}

template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::cSqr
(
    const scalar rho,
    const scalar e,
    const scalar T,
    const scalar cv
) const
{
    return
        (sqr(dpdT(rho, e, T))*T/cv - dpdv(rho, e, T))
       /sqr(max(rho, 1e-10));
}



template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpdv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return -sqr(rho)*pTable_.dFdX(rho, e);
}


template<class Specie>
Foam::scalar Foam::tabulatedEOS<Specie>::dpdT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    scalar dTdRho(TTable_.dFdX(rho, e));
    scalar dTde(TTable_.dFdY(rho, e));
    scalar dpdt = 0;

    if (mag(dTdRho) > SMALL)
    {
        dpdt += pTable_.dFdX(rho, e)/dTdRho;
    }
    if (mag(dTde) > SMALL)
    {
        dpdt += pTable_.dFdY(rho, e)/dTde;
    }
    return dpdt;
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
