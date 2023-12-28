/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "myBasicThermo.H"
#include "stringOps.H"
#include "wordIOList.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myBasicThermo, 0);
    defineRunTimeSelectionTable(myBasicThermo, fvMesh);
    defineRunTimeSelectionTable(myBasicThermo, fvMeshDictPhase);
}

const Foam::word Foam::myBasicThermo::dictName("thermophysicalProperties");

const Foam::wordList Foam::myBasicThermo::componentHeader4
({
    "type",
    "mixture",
    "properties",
    "energy"
});

const Foam::wordList Foam::myBasicThermo::componentHeader7
({
    "type",
    "mixture",
    "transport",
    "thermo",
    "equationOfState",
    "specie",
    "energy"
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::Ostream& Foam::myBasicThermo::printThermoNames
(
    Ostream& os,
    const wordList& cmptNames,
    const wordList& thermoNames
)
{
    const int nCmpt = cmptNames.size();

    // Build a table of constituent parts by split name into constituent parts
    // - remove incompatible entries from the list
    // - note: row-0 contains the names of constituent parts (ie, the header)

    DynamicList<wordList> outputTbl;
    outputTbl.resize(thermoNames.size()+1);

    label rowi = 0;

    // Header
    outputTbl[rowi] = cmptNames;
    if (!outputTbl[rowi].empty())
    {
        ++rowi;
    }

    for (const word& thermoName : thermoNames)
    {
        outputTbl[rowi] = myBasicThermo::splitThermoName(thermoName, nCmpt);
        if (!outputTbl[rowi].empty())
        {
            ++rowi;
        }
    }

    if (rowi > 1)
    {
        outputTbl.resize(rowi);
        Foam::printTable(outputTbl, os);
    }

    return os;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::myBasicThermo::makeThermoName
(
    const dictionary& thermoTypeDict,
    const wordList*& cmptHeaderPtr
)
{
    if (thermoTypeDict.found("properties"))
    {
        if (cmptHeaderPtr)
        {
            cmptHeaderPtr = &(componentHeader4);
        }

        return word
        (
            thermoTypeDict.get<word>("type") + '<'
          + thermoTypeDict.get<word>("mixture") + '<'
          + thermoTypeDict.get<word>("properties") + ','
          + thermoTypeDict.get<word>("energy") + ">>"
        );
    }
    else
    {
        if (cmptHeaderPtr)
        {
            cmptHeaderPtr = &(componentHeader7);
        }

        return word
        (
            thermoTypeDict.get<word>("type") + '<'
          + thermoTypeDict.get<word>("mixture") + '<'
          + thermoTypeDict.get<word>("transport") + '<'
          + thermoTypeDict.get<word>("thermo") + '<'
          + thermoTypeDict.get<word>("equationOfState") + '<'
          + thermoTypeDict.get<word>("specie") + ">>,"
          + thermoTypeDict.get<word>("energy") + ">>>"
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::myBasicThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf = this->T_.boundaryField();

    wordList hbt(tbf.size());

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const auto& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const auto& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList Foam::myBasicThermo::heBoundaryTypes()
{
    const volScalarField::Boundary& tbf = this->T_.boundaryField();

    wordList hbt(tbf.types());

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::volScalarField& Foam::myBasicThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& fieldName,
    bool& isOwner
)
{
    auto* ptr = mesh.objectRegistry::getObjectPtr<volScalarField>(fieldName);

    isOwner = !ptr;


    if (!ptr)
    {
        
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        // Transfer ownership of this object to the objectRegistry
        ptr->store();
    }

    Info << " Constructed field in myBasicThermo constructor : " << fieldName << endl;
    Info << *ptr << endl;
    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myBasicThermo::myBasicThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phaseName_(phaseName),

    pOwner_(false),
    TOwner_(false),
    dpdt_(getOrDefault<bool>("dpdt", true)),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),
    T_(lookupOrConstruct(mesh, phasePropertyName("T"), TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),
    e_
    (
        IOobject
        (
            myBasicThermo::phasePropertyName
            (
                "e", phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0.0),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            IOobject::groupName("Cv", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    Info << "Reading dictionary named " << this->name() << endl;
    Info << "Constructing myBasicThermo " << endl;

    this->readIfPresent("updateT", TOwner_);  // Manual override

}


Foam::myBasicThermo::myBasicThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),

    phaseName_(phaseName),

    pOwner_(false),
    TOwner_(false),
    dpdt_(getOrDefault<bool>("dpdt", true)),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),
    T_(lookupOrConstruct(mesh, phasePropertyName("T"), TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),
    e_
    (
        IOobject
        (
            myBasicThermo::phasePropertyName
            (
                "e", phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0.0),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            IOobject::groupName("Cv", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    this->readIfPresent("updateT", TOwner_);  // Manual override
}


Foam::myBasicThermo::myBasicThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    IOdictionary
    (
        IOobject
        (
            dictionaryName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phaseName_(phaseName),

    pOwner_(false),
    TOwner_(false),
    dpdt_(getOrDefault<bool>("dpdt", true)),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),
    T_(lookupOrConstruct(mesh, "T", TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),
    e_
    (
        IOobject
        (
            myBasicThermo::phasePropertyName
            (
                "e", phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0.0),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            IOobject::groupName("Cv", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    this->readIfPresent("updateT", TOwner_);  // Manual override

    if (debug)
    {
        Pout<< "Constructed shared thermo : mesh:" << mesh.name()
            << " phase:" << phaseName
            << " dictionary:" << dictionaryName
            << " T:" << T_.name()
            << " updateT:" << TOwner_
            << " alphaName:" << alpha_.name()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::myBasicThermo> Foam::myBasicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<myBasicThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myBasicThermo::~myBasicThermo()
{
    if (pOwner_)
    {
        db().checkOut(p_.name());
    }

    if (TOwner_)
    {
        db().checkOut(T_.name());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::myBasicThermo& Foam::myBasicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    const myBasicThermo* thermo = pf.db().findObject<myBasicThermo>(dictName);

    if (thermo)
    {
        return *thermo;
    }

    HashTable<const myBasicThermo*> thermos =
        pf.db().lookupClass<myBasicThermo>();

    forAllConstIters(thermos, iter)
    {
        thermo = iter.val();
        if
        (
            &(thermo->he().internalField())
         == &(pf.internalField())
        )
        {
            return *thermo;
        }
    }

    return pf.db().lookupObject<myBasicThermo>(dictName);
}


void Foam::myBasicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::myBasicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types: " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::myBasicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types: " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::myBasicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
         || he().name() == phasePropertyName(d)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types: " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::wordList Foam::myBasicThermo::splitThermoName
(
    const std::string& thermoName,
    const int nExpectedCmpts
)
{
    // Split on ",<>" but include space for good measure.
    // Splits things like
    // "hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,enthalpy>>>"

    const auto parsed = stringOps::splitAny<std::string>(thermoName, " ,<>");
    const int nParsed(parsed.size());

    wordList cmpts;

    if (!nExpectedCmpts || nParsed == nExpectedCmpts)
    {
        cmpts.resize(nParsed);

        auto iter = cmpts.begin();
        for (const auto& sub : parsed)
        {
            *iter = word(sub.str());
            ++iter;
        }
    }

    return cmpts;
}


Foam::volScalarField& Foam::myBasicThermo::p()
{
    Info << "returning p from myBasicThermo" << p_ << endl;

    return p_;
}


const Foam::volScalarField& Foam::myBasicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::myBasicThermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::myBasicThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::myBasicThermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::myBasicThermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


bool Foam::myBasicThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
