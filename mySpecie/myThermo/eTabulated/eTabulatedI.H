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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eTabulated<EquationOfState>::eTabulated
(
    const EquationOfState& st
)
:
    EquationOfState(st)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*template<class EquationOfState>
inline Foam::eTabulated<EquationOfState>::eTabulated
(
    const word& name,
    const eTabulated<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Esref_(ct.Esref_)
{}*/


template<class EquationOfState>
inline Foam::autoPtr<Foam::eTabulated<EquationOfState>>
Foam::eTabulated<EquationOfState>::clone() const
{
    return autoPtr<eTabulated<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eTabulated<EquationOfState>>
Foam::eTabulated<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eTabulated<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::limit
(
    const scalar T
) const
{
    // ok for now
    return T;
}


template<class EquationOfState>
Foam::scalar Foam::eTabulated<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return 1.0/max(this->TTable_.dFdY(rho, e), SMALL);
}


template<class EquationOfState>
Foam::scalar Foam::eTabulated<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    const scalar cv = Cv(rho, e, T);
    return CpMCv(rho, e, T) + cv;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    //Info << "returned Es " << eTable_.lookup(rho,T) << endl;

    //return this->TTable_.reverseLookupY(T, rho);

    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::EsRhoT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    //Info << "returned Es " << eTable_.lookup(rho,T) << endl;
    // Info << "reverseLookupY(T, rho),  rho = " << rho << endl;
    // Info << "T = " << T << endl; 
    return this->TTable_.reverseLookupY(T, rho);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::Ea
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->TTable_.reverseLookupY(T, rho);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::tabulatedTRhoEs
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Info << "rho " << rho << endl;
    // Info << "e " << e << endl;
    // Info << "T " << this->TTable_.lookup(rho, e) << endl;
    return this->TTable_.lookup(rho, e);
}

template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::gamma
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    //Info << "gamma " << EquationOfState::gamma(rho, e, T, Cv(rho, e, T));
    return EquationOfState::gamma(rho, e, T, Cv(rho, e, T)); 
}

template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::CpMCv(scalar rho, scalar e, scalar T) const
{
    return EquationOfState::dpdT(rho, e, T)/max(rho, 1e-10);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::Hc() const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::S
(
    const scalar rho,
    const scalar e
) const
{
    NotImplemented;
    return 0;
}


/*template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return
        Cv_*(T - Tref_) + Esref_ + Hc() + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}*/


/*template<class EquationOfState>
inline Foam::scalar Foam::eTabulated<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}*/

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*template<class EquationOfState>
inline void Foam::eTabulated<EquationOfState>::operator+=
(
    const eTabulated<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        Esref_ = Y1*Esref_ + Y2*ct.Esref_;
    }
}*/


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

/*template<class EquationOfState>
inline Foam::eTabulated<EquationOfState> Foam::operator+
(
    const eTabulated<EquationOfState>& ct1,
    const eTabulated<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return eTabulated<EquationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.Hf_,
            ct1.Tref_,
            ct1.Esref_
        );
    }
    else
    {
        return eTabulated<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Tref_,
            ct1.Y()/eofs.Y()*ct1.Esref_
          + ct2.Y()/eofs.Y()*ct2.Esref_
        );
    }
}
*/

/*template<class EquationOfState>
inline Foam::eTabulated<EquationOfState> Foam::operator*
(
    const scalar s,
    const eTabulated<EquationOfState>& ct
)
{
    return eTabulated<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_,
        ct.Tref_,
        ct.Esref_
    );
}
*/

/*template<class EquationOfState>
inline Foam::eTabulated<EquationOfState> Foam::operator==
(
    const eTabulated<EquationOfState>& ct1,
    const eTabulated<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return eTabulated<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.Esref_
      - ct1.Y()/eofs.Y()*ct1.Esref_
    );
}
*/

// ************************************************************************* //
