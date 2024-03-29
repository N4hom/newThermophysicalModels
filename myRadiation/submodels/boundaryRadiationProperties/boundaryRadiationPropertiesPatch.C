/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "boundaryRadiationPropertiesPatch.H"
#include "mappedPatchBase.H"
#include "myRadiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //
namespace Foam
{
    namespace myRadiation
    {
        defineTypeNameAndDebug(boundaryRadiationPropertiesPatch, 0);
        defineRunTimeSelectionTable
        (
            boundaryRadiationPropertiesPatch,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::myRadiation::boundaryRadiationPropertiesPatch>
Foam::myRadiation::boundaryRadiationPropertiesPatch::New
(
    const dictionary& dict,
    const polyPatch& pp
)
{
    const word modelType(dict.getCompat<word>("type", {{"mode", 1812}}));

    Info<< "Selecting boundary radiation Model: "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "radiationModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, pp);
}


// * * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * //

Foam::label
Foam::myRadiation::boundaryRadiationPropertiesPatch::nbrPatchIndex() const
{
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_);

    return (mpp.samplePolyPatch().index());
}


const Foam::fvMesh&
Foam::myRadiation::boundaryRadiationPropertiesPatch::nbrRegion() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_);

     return (refCast<const fvMesh>(mpp.sampleMesh()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myRadiation::boundaryRadiationPropertiesPatch::
boundaryRadiationPropertiesPatch
(
    const dictionary& dict,
    const polyPatch& p
)
:
    dict_(dict),
    patch_(p),
    absorptionEmission_(nullptr),
    transmissivity_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::myRadiation::wallAbsorptionEmissionModel&
Foam::myRadiation::boundaryRadiationPropertiesPatch::absorptionEmission() const
{
    return *absorptionEmission_;
}


const Foam::myRadiation::wallTransmissivityModel&
Foam::myRadiation::boundaryRadiationPropertiesPatch::transmissiveModel() const
{
    return *transmissivity_;
}


void Foam::myRadiation::boundaryRadiationPropertiesPatch::write(Ostream& os) const
{
    NotImplemented;
}


// ************************************************************************* //
