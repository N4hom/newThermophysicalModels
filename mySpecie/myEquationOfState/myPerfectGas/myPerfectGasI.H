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

#include "myPerfectGas.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::myPerfectGas<Specie>::myPerfectGas(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::myPerfectGas<Specie>::myPerfectGas
(
    const word& name,
    const myPerfectGas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::myPerfectGas<Specie>>
Foam::myPerfectGas<Specie>::clone() const
{
    return autoPtr<myPerfectGas<Specie>>::New(*this);
}


template<class Specie>
inline Foam::autoPtr<Foam::myPerfectGas<Specie>> Foam::myPerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<myPerfectGas<Specie>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::rho(scalar p, scalar T) const
{
   
    return p/(this->R()*T);
}

template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::p(scalar rho , scalar e, scalar T) const
{
    
    return rho * this->R() * T;
}

template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::cSqr(scalar rho , scalar e, scalar T, scalar cv) const
{
    // hardcoded for now

    // scalar gamma(gamma(rho, e, T, cv));
    // Info << "gamma = " << gamma;

    // Info << "gamma(rho, e, T, cv) \n" << gamma(rho, e, T, cv) << endl;
    // Info << "R \n" << this->R() << endl;
    // Info << "T \n" << T << endl;

    // Info << "here " << endl;
    // Info << "gamma(rho, e, T, cv) " << gamma(rho, e, T, cv)<< endl;

    return gamma(rho, e, T, cv) * this->R() * T;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::H(scalar p, scalar T) const
{
    
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
Foam::scalar Foam::myPerfectGas<Specie>::gamma(const scalar rho, const scalar e, const scalar T, const scalar cv) const
{

    return this->R()/cv + 1.0;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::Z(scalar p, scalar T) const
{
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::myPerfectGas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::myPerfectGas<Specie>::operator+=(const myPerfectGas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::myPerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::myPerfectGas<Specie> Foam::operator+
(
    const myPerfectGas<Specie>& pg1,
    const myPerfectGas<Specie>& pg2
)
{
    return myPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1) + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::myPerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const myPerfectGas<Specie>& pg
)
{
    return myPerfectGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::myPerfectGas<Specie> Foam::operator==
(
    const myPerfectGas<Specie>& pg1,
    const myPerfectGas<Specie>& pg2
)
{
    return myPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
    );
}


// ************************************************************************* //
