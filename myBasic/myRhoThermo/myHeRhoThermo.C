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
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
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

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::calculateFromRhoE()
{
    Info << "calculateFromRhoE " << endl;
    scalarField& eCells = this->heRef().primitiveFieldRef();
    scalarField& TCells = this->TRef().primitiveFieldRef();
    scalarField& pCells = this->pRef().primitiveFieldRef();
    scalarField& CpCells = this->CpRef().primitiveFieldRef();
    scalarField& CvCells = this->CvRef().primitiveFieldRef();
    scalarField& speedOfSoundCells = this->speedOfSoundRef().primitiveFieldRef();

    const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(0);

    Info << "mixture_.Cp(this->rho_[0], eCells[0], TCells[0]) " << mixture_.Cp(this->rho_[0], eCells[0], TCells[0])  << endl;


    forAll(this->rho_, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        const scalar& rhoi(this->rho_[celli]);
        //Info << "scalar& rhoi(this->rho_[celli]) " << rhoi << endl;
        scalar& ei(eCells[celli]);
        scalar& Ti(TCells[celli]);

        TCells[celli] = mixture_.TRhoE(Ti, rhoi, ei);

       

        scalar pi = mixture_.p(rhoi, ei, Ti);
        //Info << "calculated pi  " << pi << endl; 
        scalar Cpi = mixture_.Cp(rhoi, ei, Ti);
        scalar Cvi = mixture_.Cv(rhoi, ei, Ti);

        pCells[celli] = pi;
        CpCells[celli] = Cpi;
        CvCells[celli] = Cvi;

        speedOfSoundCells[celli] = sqrt(max(mixture_.cSqr(rhoi, ei, Ti, Cvi), SMALL));
    }

     this->TRef().correctBoundaryConditions();
     this->pRef().correctBoundaryConditions();
     this->eRef().correctBoundaryConditions();
     // this->eRef().correctBoundaryConditions(); gives a runtime error due to an error during lookup



    volScalarField::Boundary& be  = this->eRef().boundaryFieldRef();   
    volScalarField::Boundary& bCp = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& bCv = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& bmu = this->muRef().boundaryFieldRef();
    volScalarField::Boundary& balpha = this->alphaRef().boundaryFieldRef();
    volScalarField::Boundary& bspeedOfSound =
        this->speedOfSoundRef().boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {

        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->TRef().boundaryField()[patchi];
        const fvPatchScalarField& phe = this->heRef().boundaryField()[patchi];
        const fvPatchScalarField& pp = this->pRef().boundaryField()[patchi];

        fvPatchScalarField& pCp = bCp[patchi];
        fvPatchScalarField& pCv = bCv[patchi];
        fvPatchScalarField& pmu = bmu[patchi];
        fvPatchScalarField& palpha = balpha[patchi];
        fvPatchScalarField& pspeedOfSound = bspeedOfSound[patchi];

        forAll(prho, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                        this->patchFaceMixture(patchi, facei);

            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);

            const scalar Cpi = mixture_.Cp(rhoi, ei, Ti);
            const scalar Cvi = mixture_.Cv(rhoi, ei, Ti);
            pCp[facei] = Cpi;
            pCv[facei] = Cvi;
            //pmu[facei] = mixture_.mu(rhoi, ei, Ti);
            //palpha[facei] = mixture_.kappa(rhoi, ei, Ti)/Cpi;
            pspeedOfSound[facei] =
                sqrt(max(mixture_.cSqr(rhoi, ei, Ti, Cvi), SMALL));
        }
    }

    Info << "End of calculateFromRhoE " << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::myHeRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    myHeThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    Info << "Constructing myHeRhoThermo " << endl;

    Info << this->name() << endl;
    
    calculateFromRhoE();

    // calculate
    // (
    //     this->p_,
    //     this->T_,
    //     this->he_,
    //     this->psi_,
    //     this->rho_,
    //     this->mu_,
    //     this->alpha_,
    //     true                    // Create old time fields
    // );

    Info << "myHeRhoThermo constructed" << endl;
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
    calculateFromRhoE();

    // calculate
    // (
    //     this->p_,
    //     this->T_,
    //     this->he_,
    //     this->psi_,
    //     this->rho_,
    //     this->mu_,
    //     this->alpha_,
    //     true                    // Create old time fields
    // );
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
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::myHeRhoThermo<BasicPsiThermo, MixtureType>::correctFromRhoE()
{

    calculateFromRhoE();

}
// ************************************************************************* //
