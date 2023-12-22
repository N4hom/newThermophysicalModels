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
inline Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const EquationOfState& st,
    const scalarLookupTable2D& cvTable,
    const scalarLookupTable2D& eTable,
    const scalar Hf
)
:
    EquationOfState(st),
    eTable_(eTable),
    Tlow_(min(this->cvTable_.x())),
    Thigh_(max(this->cvTable_.x()))
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*template<class EquationOfState>
inline Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const word& name,
    const eTabulatedThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Esref_(ct.Esref_)
{}*/


template<class EquationOfState>
inline Foam::autoPtr<Foam::eTabulatedThermo<EquationOfState>>
Foam::eTabulatedThermo<EquationOfState>::clone() const
{
    return autoPtr<eTabulatedThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eTabulatedThermo<EquationOfState>>
Foam::eTabulatedThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eTabulatedThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    // ok for now
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return this->cvTable_.lookup(p,T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Es
(
    const scalar p,
    const scalar T
) const
{
    Info << "returned Es " << eTable_.lookup(p,T) << endl;
    return eTable_.lookup(p,T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Ea
(
    const scalar p,
    const scalar T
) const
{
    Info << "Es(p, T) @ " << "p="  << p << " T=" << T << " " << Es(p, T) << endl;
    return Es(p, T) ;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Hc() const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


/*template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return
        Cv_*(T - Tref_) + Esref_ + Hc() + Pstd/EquationOfState::rho(Pstd, T)
      - S(Pstd, T)*T;
}*/


/*template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dCpdT
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
inline void Foam::eTabulatedThermo<EquationOfState>::operator+=
(
    const eTabulatedThermo<EquationOfState>& ct
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
inline Foam::eTabulatedThermo<EquationOfState> Foam::operator+
(
    const eTabulatedThermo<EquationOfState>& ct1,
    const eTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return eTabulatedThermo<EquationOfState>
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
        return eTabulatedThermo<EquationOfState>
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
inline Foam::eTabulatedThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eTabulatedThermo<EquationOfState>& ct
)
{
    return eTabulatedThermo<EquationOfState>
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
inline Foam::eTabulatedThermo<EquationOfState> Foam::operator==
(
    const eTabulatedThermo<EquationOfState>& ct1,
    const eTabulatedThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return eTabulatedThermo<EquationOfState>
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
