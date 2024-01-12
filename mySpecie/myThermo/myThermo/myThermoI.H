/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include "myThermo.H"


// Type parameter defines wether the state function is absolute or relative to a standard state.

// Thermo parameter defines how to convert from the state function to temperature. 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::mySpecies::thermo<Thermo, Type>::thermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{
        Info << "Constructing thermo as named copy " << endl;

}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::T
// (
//     scalar f,
//     scalar p,
//     scalar T0,
//     scalar (thermo<Thermo, Type>::*F)(const scalar, const scalar) const,
//     scalar (thermo<Thermo, Type>::*dFdT)(const scalar, const scalar)
//         const,
//     scalar (thermo<Thermo, Type>::*limit)(const scalar) const
// ) const
// {
//     if (T0 < 0)
//     {
//         FatalErrorInFunction
//             << "Negative initial temperature T0: " << T0
//             << abort(FatalError);
//     }

//     scalar Test = T0;
//     scalar Tnew = T0;
//     scalar Ttol = T0*tol_;
//     int    iter = 0;

//     do
//     {
//         Test = Tnew;
//         Tnew =
//             (this->*limit)
//             (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

//         if (iter++ > maxIter_)
//         {
//             FatalErrorInFunction
//                 << "Maximum number of iterations exceeded: " << maxIter_
//                 << " when starting from T0:" << T0
//                 << " old T:" << Test << " new T:" << Tnew
//                 << " f:" << f
//                 << " p:" << p
//                 << " tol:" << Ttol
//                 << abort(FatalError);
//         }

//     } while (mag(Tnew - Test) > Ttol);

//     return Tnew;
// }

// *F, *dFdT and *limit are function pointer and all require three arguments 
// instead of two or one
template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::TFromRhoE
(
    const scalar f,
    const scalar e,
    const scalar rho,
    const scalar T0,
    scalar (thermo<Thermo, Type>::*F)(const scalar rho, const scalar e, const scalar T0) const,
    scalar (thermo<Thermo, Type>::*dFdT)(const scalar rho, const scalar e, const scalar T0) 
        const,
    scalar (thermo<Thermo, Type>::*limit)(const scalar T0) const,
    const bool diagnostics
) const
{
    Info << "TFromRhoE " << endl;


    // I need to find a way to do a lookup without this whole thing below.
    // It should be something like


    // if(Thermo::typeName() == "eTabulated")
    // {   
    //     return Thermo::TFromRhoE(rho, e, T);
    // }

    if (rho < SMALL)
    {
         Info << "Here " << endl;
        return T0;
    }

    scalar Test = (this->*limit)(T0);
    scalar Tnew = Test;

    //Info << "Test = " << Test << endl; 


    int    iter = 0;
    scalar relError = GREAT;
    scalar absError = GREAT;

    if (diagnostics)
    {
        const unsigned int width = IOstream::defaultPrecision() + 8;

        InfoInFunction
            << "Energy -> temperature conversion failed to converge:" << endl;
        // Pout<< setw(width) << "iter"
        //     << setw(width) << "rho"
        //     << setw(width) << "e"
        //     << setw(width) << "Test"
        //     << setw(width) << "e/h Actual"
        //     << setw(width) << "e/h Test"
        //     << setw(width) << "Cv/Cp"
        //     << setw(width) << "Tnew"
        //     << endl;
    }

    do
    {
        Test = Tnew;
        //Info << "Test = " << Test << endl; 

        // limit has only one parameter
        Tnew =
            (this->*limit)
            (
                Test
              - ((this->*F)(rho, f, Test) - f)/(this->*dFdT)(rho, f, Test)
                //,Test     
            );

        //Info << "Test = " << Test << endl; 

        absError = mag((Tnew - Test));
        relError = absError/max(Tnew, SMALL);

        if (diagnostics)
        {
            // It looks like a debugging thing so I commenting out and leave it for the future
            // const unsigned int width = IOstream::defaultPrecision() + 8;

            // Pout<< setw(width) << iter
            //     << setw(width) << rho
            //     << setw(width) << e
            //     << setw(width) << Test
            //     << setw(width) << f
            //     << setw(width) << ((this->*F)(rho, f, Test))
            //     << setw(width) << ((this->*dFdT)(rho, f, Test))
            //     << setw(width) << Tnew
            //     << endl;
        }

        if (iter++ > maxIter_)
        {
            if (!diagnostics)
            {
                TFromRhoE(f, e, rho, T0, F, dFdT, limit, true);
            }

            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << abort(FatalError);
        }
    } while (relError > relTol_ && absError > absTol_);
    
    //Hardcoding just for testing
    
    //scalar Cv(705.0);

    //scalar Tnew = e/Cv;
    //Info << "Returning Tnew " << endl;
    return Tnew;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::mySpecies::thermo<Thermo, Type>::thermo
(
    const word& name,
    const thermo& st
)
:
    Thermo(name, st)
{
    Info << "Constructing thermo " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::mySpecies::thermo<Thermo, Type>::heName()
{
    return Type<thermo<Thermo, Type>>::energyName();
}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar
// Foam::mySpecies::thermo<Thermo, Type>::Cp(const scalar rho, const scalar e, const scalar T) const
// {
//     Info << "Cp in myThermoI " << endl;
//     return Type<thermo<Thermo, Type>>::Cp(*this, rho, e, T);
// }

// template<class Thermo, template<class> class Type>
// inline Foam::scalar
// Foam::mySpecies::thermo<Thermo, Type>::Cv(const scalar rho, const scalar e, const scalar T) const
// {
//     return Type<thermo<Thermo, Type>>::Cv(*this, rho, e, T);
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar
// Foam::mySpecies::thermo<Thermo, Type>::Cpv(const scalar p, const scalar T) const
// {
//     return Type<thermo<Thermo, Type>>::Cpv(*this, p, T);
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::Cpv(const scalar rho, const scalar e, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::Cpv(*this, rho, e, T);
}

// Modified to be a function of (rho, e)
// template<class Thermo, template<class> class Type>
// inline Foam::scalar
// Foam::mySpecies::thermo<Thermo, Type>::gamma(const scalar p, const scalar T) const
// {
//     #ifdef __clang__
//     volatile const scalar Cp = this->Cp(p, T);
//     #else
//     const scalar Cp = this->Cp(p, T);
//     #endif

//     return Cp/(Cp - this->CpMCv(p, T));
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::gamma(const scalar rho, const scalar e, const scalar T) const
{
    
    //const scalar Cp = this->Cp(rho, e, T);
    //return Cp/(Cp - this->CpMCv(p, T));

    // Hardcoding just for testing purposes
    //const scalar gamma(5/3);
    return Type<thermo<Thermo, Type>>::gamma(*this, rho, e, T);
}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar
// Foam::mySpecies::thermo<Thermo, Type>::CpByCpv
// (
//     const scalar p,
//     const scalar T
// ) const
// {
//     return Type<thermo<Thermo, Type>>::CpByCpv(*this, p, T);
// }

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::CpByCpv
(
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return Type<thermo<Thermo, Type>>::CpByCpv(*this, rho, e, T);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::HE(const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::HE(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::HErhoT(const scalar rho, const scalar e, const scalar T) const
{
    // Implementation of Type<thermo<Thermo, Type>>::HErhoT in sensibleInternalEnergy.H
    Info << "In myThermo \n" << endl; 
    Info << "rho = " << rho << endl;
    return Type<thermo<Thermo, Type>>::HErhoT(*this, rho, e, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T*this->S(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T*this->S(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T)*this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*Pstd/p, nm);
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    //Info << "Returning temperature from internal energy " << endl;
    return Type<thermo<Thermo, Type>>::THE(*this, he, p, T0);
}

template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::TRhoE
(
    const scalar T,
    const scalar rho,
    const scalar e
) const
{
    //Info << "Returning temperature from internal energy " << endl;
    return Type<thermo<Thermo, Type>>::TRhoE(*this, T, rho, e);
}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::THs
// (
//     const scalar hs,
//     const scalar p,
//     const scalar T0
// ) const
// {
//     return T
//     (
//         hs,
//         p,
//         T0,
//         &thermo<Thermo, Type>::Hs,
//         &thermo<Thermo, Type>::Cp,
//         &thermo<Thermo, Type>::limit
//     );
// }


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::THa
// (
//     const scalar ha,
//     const scalar p,
//     const scalar T0
// ) const
// {
//     return T
//     (
//         ha,
//         p,
//         T0,
//         &thermo<Thermo, Type>::Ha,
//         &thermo<Thermo, Type>::Cp,
//         &thermo<Thermo, Type>::limit
//     );
// }


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    NotImplemented;
    return 0;
}



template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::TRhoEs
(
    const scalar T,
    const scalar rho,
    const scalar es
) const
{
    // thermo<Thermo, Type>::EsRhoT can be found in eConstThermo.H et simila
    //Info << "TRhoE in myThermoI.H " << endl;
    // Info << "TFromRhoE = " << TFromRhoE(
    //     es,
    //     es,
    //     rho,
    //     T,
    //     &thermo<Thermo, Type>::EsRhoT,
    //     &thermo<Thermo, Type>::Cv,
    //     &thermo<Thermo, Type>::limit,
    //     false
    // )  << endl;
    
    return TFromRhoE
    (
        es,
        es,
        rho,
        T,
        &thermo<Thermo, Type>::EsRhoT,
        &thermo<Thermo, Type>::Cv,
        &thermo<Thermo, Type>::limit,
        false
    ); 
}


// template<class Thermo, template<class> class Type>
// inline Foam::scalar Foam::mySpecies::thermo<Thermo, Type>::TEa
// (
//     const scalar ea,
//     const scalar p,
//     const scalar T0
// ) const
// {
//     return T
//     (
//         ea,
//         p,
//         T0,
//         &thermo<Thermo, Type>::Ea,
//         &thermo<Thermo, Type>::Cv,
//         &thermo<Thermo, Type>::limit
//     );
// }


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::dKcdTbyKc
(
    const scalar p,
    const scalar T
) const
{
    const scalar dKcdTbyKc =
        (this->S(Pstd, T) + this->Gstd(T)/T)*this->Y()/(RR*T);

    const scalar nm = this->Y()/this->W();
    if (equal(nm, SMALL))
    {
        return dKcdTbyKc;
    }
    else
    {
        return dKcdTbyKc - nm/T;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::mySpecies::thermo<Thermo, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::mySpecies::thermo<Thermo, Type>::operator+=
(
    const thermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::mySpecies::thermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::mySpecies::thermo<Thermo, Type> Foam::mySpecies::operator+
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::mySpecies::thermo<Thermo, Type> Foam::mySpecies::operator*
(
    const scalar s,
    const thermo<Thermo, Type>& st
)
{
    return thermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::mySpecies::thermo<Thermo, Type> Foam::mySpecies::operator==
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) == static_cast<const Thermo&>(st2)
    );
}


// ************************************************************************* //
