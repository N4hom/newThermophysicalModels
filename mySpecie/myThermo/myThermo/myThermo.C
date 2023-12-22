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

#include "myThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class Thermo, template<class> class Type>
const Foam::scalar Foam::species::myThermo<Thermo, Type>::tol_ = 1.0e-4;

template<class Thermo, template<class> class Type>
const int Foam::species::myThermo<Thermo, Type>::maxIter_ = 100;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
Foam::species::myThermo<Thermo, Type>::myThermo(const dictionary& dict)
:
    Thermo(dict)
{
    Info << "Constructing myThermo derived from " << 
    Thermo::typeName() <<  " and " << Type<myThermo<Thermo, Type>>::typeName() << endl;
    Info << "myThermo typeName " << typeName() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
void Foam::species::myThermo<Thermo, Type>::write(Ostream& os) const
{
    Thermo::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
Foam::Ostream& Foam::species::operator<<
(
    Ostream& os, const myThermo<Thermo, Type>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
