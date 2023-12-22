/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "myHeRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& speedOfSound,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    Info << "Calculating thermo quantities " << endl;
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            speedOfSound.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& speedOfSoundCells = speedOfSound.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT() && !initialize_)
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        speedOfSoundCells[celli] = mixture_.cSqr(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
        Info << "rhoCells[celli] " << rhoCells[celli] << endl;

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);

       
    }


    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& speedOfSoundBf = speedOfSound.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pspeedOfSound = speedOfSoundBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                Info << "pT.fixesValue() is true " << endl; 
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                pspeedOfSound[facei] = mixture_.cSqr(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                Info << "rho boundaryField " << prho[facei] << endl;
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT() && !initialize_)
                {
                    if (!initialize_)
                    {
                        
                        Info << "!initialize_ " << !initialize_ << endl;
                    }
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                pspeedOfSound[facei] = mixture_.cSqr(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::myHeRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    myHeThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    initialize_(true)
{
    Info << "Constructing myHeRhoThermo " << endl;
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->speedOfSound_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
    Info << "myHeRhoThermo constructed " << endl;

    initialize_ = false;
}


template<class BasicPsiThermo, class MixtureType>
Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::myHeRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    myHeThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->speedOfSound_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::~myHeRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->speedOfSound_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );


    Info << "p to be updated " << this->p_ << endl;
    updateP(this->rho_, this->T_);

    Info << "updated p " << this->p_ << endl;

    DebugInFunction << "Finished" << endl;
}


//////////////////////////////////

template<class BasicPsiThermo, class MixtureType>
void Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::updateP(volScalarField& rho, volScalarField& T)
{
    scalarField& pCells = this->p_.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& TCells = T.primitiveFieldRef();

    forAll(pCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);   
    
        pCells[celli] = mixture_.p(rhoCells[celli], TCells[celli]);
        Info << "celli " << celli << endl;
        Info << "rhoCells " << rhoCells[celli] << endl;
        Info << "TCells " << TCells[celli] << endl;
    }


    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    
    forAll(pBf, patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        
        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pp[facei] = mixture_.p(prho[facei], pT[facei]);
                
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);
                
                pp[facei] = mixture_.p(prho[facei], pT[facei]);
            }
        }
    }
}


// ************************************************************************* //
