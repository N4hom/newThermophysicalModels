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
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const EquationOfState& st,
    const scalar cv,
    const scalar hf,
    const scalar Tref,
    const scalar Esref
)
:
    EquationOfState(st),
    Cv_(cv),
    Hf_(hf),
    Tref_(Tref),
    Esref_(Esref)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const word& name,
    const eConstThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Esref_(ct.Esref_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState>>
Foam::eConstThermo<EquationOfState>::clone() const
{
    return autoPtr<eConstThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState>>
Foam::eConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eConstThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


/*template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + EquationOfState::Cv(p, T);
}*/

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Cv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Cv_;
}

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Cp
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    // Just for testing
    //Info << "Cp in eConstThermo " << endl;

    
    
    return Cv(rho, e, T) * EquationOfState::gamma(rho, e, T, Cv(rho, e, T)); 
}

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::gamma
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
   // Info << "gamma " << EquationOfState::gamma(rho, e, T, Cv(rho, e, T));
    return EquationOfState::gamma(rho, e, T, Cv(rho, e, T)); 
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    return Cv_*(T - Tref_) + Esref_ + EquationOfState::E(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::EsRhoT
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    //Info << "EsRhoT " << endl;    
    return Cv_*(T - Tref_) + Esref_ ;
}




template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Ea
(
    const scalar p,
    const scalar T
) const
{
    return Es(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T)*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return
        Cv_*(T - Tref_) + Esref_ + Hc() + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eConstThermo<EquationOfState>::operator+=
(
    const eConstThermo<EquationOfState>& ct
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
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator+
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return eConstThermo<EquationOfState>
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
        return eConstThermo<EquationOfState>
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


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eConstThermo<EquationOfState>& ct
)
{
    return eConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_,
        ct.Tref_,
        ct.Esref_
    );
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator==
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return eConstThermo<EquationOfState>
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


// ************************************************************************* //
