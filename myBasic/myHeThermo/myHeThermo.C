/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "myHeThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::myHeThermo<BasicThermo, MixtureType>::
heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hBf[patchi]).gradient()
                = hBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hBf[patchi]).refGrad()
                = hBf[patchi].fvPatchField::snGrad();
        }
    }
}


// template<class BasicThermo, class MixtureType>
// void Foam::myHeThermo<BasicThermo, MixtureType>::init
// (
//     const volScalarField& p,
//     const volScalarField& T,
//     volScalarField& he
// )
// {
//     scalarField& heCells = he.primitiveFieldRef();
//     const scalarField& pCells = p.primitiveField();
//     const scalarField& TCells = T.primitiveField();

//     forAll(heCells, celli)
//     {
//         heCells[celli] =
//             this->cellMixture(celli).HE(pCells[celli], TCells[celli]);
//         Info << "heCells[celli] " << heCells[celli] << endl;
//     }

//     volScalarField::Boundary& heBf = he.boundaryFieldRef();

//     forAll(heBf, patchi)
//     {
//         heBf[patchi] == this->he
//         (
//             p.boundaryField()[patchi],
//             T.boundaryField()[patchi],
//             patchi
//         );

//         heBf[patchi].useImplicit(T.boundaryField()[patchi].useImplicit());
//     }

//     this->heBoundaryCorrection(he);

//     // Note: T does not have oldTime
//     if (p.nOldTimes() > 0)
//     {
//         init(p.oldTime(), T.oldTime(), he.oldTime());
//     }
// }


template<class BasicThermo, class MixtureType>
void Foam::myHeThermo<BasicThermo, MixtureType>::initFromRhoT
(
    const volScalarField& rho,
    const volScalarField& p,    
    const volScalarField& T,
    volScalarField& he
)
{
    // Fill internal field values for the internal energy

    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& rhoCells = rho.primitiveField();
    const scalarField& TCells = T.primitiveField();

    forAll(heCells, celli)
    {
        //Info << "rhoCells[celli] " << rhoCells[celli] << endl; 
        //Info << "heCells[celli] " << this->cellMixture(celli).HErhoT(rhoCells[celli], heCells[celli], TCells[celli]) << endl;
        heCells[celli] =
            this->cellMixture(celli).HErhoT(rhoCells[celli], heCells[celli], TCells[celli]);
    }

    // Fill all the boundary values for every patch

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        heBf[patchi] == this->heRhoT
        (
            rho.boundaryField()[patchi],
            T.boundaryField()[patchi],
            patchi
        );

        heBf[patchi].useImplicit(T.boundaryField()[patchi].useImplicit());
    }

    this->heBoundaryCorrection(he);

    // Note: T does not have oldTime
    // if (p.nOldTimes() > 0)
    // {
    //     init(p.oldTime(), T.oldTime(), he.oldTime());
    // }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::myHeThermo<BasicThermo, MixtureType>::myHeThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BasicThermo(mesh, phaseName),
    MixtureType(*this, mesh, phaseName),
    he_(this->e_)
{
    Info << "Constructing myHeThermo \n" << endl;
    
    Info << "Initializing internal energy " << endl;
    initFromRhoT(this->rho_, this->p_, this->T_, he_);
}


template<class BasicThermo, class MixtureType>
Foam::myHeThermo<BasicThermo, MixtureType>::myHeThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh, phaseName),
    he_(this->e_)
{
    initFromRhoT(this->rho_, this->p_, this->T_, he_);
}


template<class BasicThermo, class MixtureType>
Foam::myHeThermo<BasicThermo, MixtureType>::myHeThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    BasicThermo(mesh, phaseName, dictionaryName),
    MixtureType(*this, mesh, phaseName),
    he_(this->e_)
{
    initFromRhoT(this->rho_, this->p_, this->T_, he_);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::myHeThermo<BasicThermo, MixtureType>::~myHeThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::myHeThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> the
    (
        new volScalarField
        (
            IOobject
            (
                "he",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& he = the.ref();
    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TCells = T;

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        scalarField& hep = heBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(hep, facei)
        {
            hep[facei] =
                this->patchFaceMixture(patchi, facei).HE(pp[facei], Tp[facei]);
        }
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, celli)
    {
        he[celli] = this->cellMixture(cells[celli]).HE(p[celli], T[celli]);
    }

    return the;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::he
// (
//     const scalarField& rho,
//     const scalarField& T,
//     const labelList& cells
// ) const
// {
//     tmp<scalarField> the(new scalarField(T.size()));
//     scalarField& he = the.ref();

//     forAll(T, celli)
//     {
//         he[celli] = this->cellMixture(cells[celli]).HE(p[celli], T[celli]);
//     }

//     return the;
// }


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, facei)
    {
        he[facei] =
            this->patchFaceMixture(patchi, facei).HE(p[facei], T[facei]);
    }

    return the;
}

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::heRhoT
(
    const scalarField& rho,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    forAll(T, facei)
    {
        he[facei] =
            this->patchFaceMixture(patchi, facei).HErhoT(rho[facei], he[facei] ,T[facei]);
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& hcf = thc.ref();
    scalarField& hcCells = hcf.primitiveFieldRef();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    volScalarField::Boundary& hcfBf = hcf.boundaryFieldRef();

    forAll(hcfBf, patchi)
    {
        scalarField& hcp = hcfBf[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}

// Commenting out because the method uses Cp(p,T) which doesn't exist anymore
// To avoid compilation issues, the method Cp(rho,e,T) must be renamed.

// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::Cp
// (
//     const scalarField& p,
//     const scalarField& T,
//     const label patchi
// ) const
// {
//     tmp<scalarField> tCp(new scalarField(T.size()));
//     scalarField& cp = tCp.ref();

//     forAll(T, facei)
//     {
//         cp[facei] =
//             this->patchFaceMixture(patchi, facei).Cp(p[facei], T[facei]);
//     }

//     return tCp;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp.ref();

    forAll(T, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp(rho[facei], e[facei], T[facei]);
    }

    return tCp;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::Cp
// (
//     const scalarField& p,
//     const scalarField& T,
//     const labelList& cells
// ) const
// {
//     auto tCp = tmp<scalarField>::New(T.size());
//     auto& Cp = tCp.ref();

//     forAll(cells, i)
//     {
//         const label celli = cells[i];
//         Cp[i] = this->cellMixture(celli).Cp(p[i], T[i]);
//     }

//     return tCp;
// }


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& cells
) const
{
    auto tCp = tmp<scalarField>::New(T.size());
    auto& Cp = tCp.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        Cp[i] = this->cellMixture(celli).Cp(rho[i], e[i], T[i]);
    }

    return tCp;
}

// Modified to be a function of (rho, e)
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp.ref();

    forAll(this->T_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp(this->rho_[celli], this->e_[celli], this->T_[celli]);
    }

    volScalarField::Boundary& cpBf = cp.boundaryFieldRef();

    forAll(cpBf, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cpBf[patchi];

        forAll(prho, facei)
        {
            pCp[facei] =
                this->patchFaceMixture(patchi, facei).Cp(prho[facei], pe[facei] , pT[facei]);
        }
    }

    return tCp;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::Cv
// (
//     const scalarField& p,
//     const scalarField& T,
//     const label patchi
// ) const
// {
//     tmp<scalarField> tCv(new scalarField(T.size()));
//     scalarField& cv = tCv.ref();

//     forAll(T, facei)
//     {
//         cv[facei] =
//             this->patchFaceMixture(patchi, facei).Cv(p[facei], T[facei]);
//     }

//     return tCv;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv.ref();

    forAll(T, facei)
    {
        cv[facei] =
            this->patchFaceMixture(patchi, facei).Cv(rho[facei], e[facei], T[facei]);
    }

    return tCv;
}

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::rhoEoS
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    auto tRho = tmp<scalarField>::New(T.size());
    auto& rho = tRho.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        rho[i] = this->cellMixture(celli).rho(p[i], T[i]);
    }

    return tRho;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv.ref();

    forAll(this->T_, celli)
    {
        cv[celli] =
            this->cellMixture(celli).Cv(this->rho_[celli], this->e_[celli], this->T_[celli]);
    }

    volScalarField::Boundary& cvBf = cv.boundaryFieldRef();

    forAll(cvBf, patchi)
    {
        // cvBf[patchi] = Cv
        // (
        //     this->rho_.boundaryField()[patchi],
        //     this->e_.boundaryField()[patchi],
        //     this->T_.boundaryField()[patchi],
        //     patchi
        // );

        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCv = cvBf[patchi];

        forAll(prho, facei)
        {
            pCv[facei] =
                this->patchFaceMixture(patchi, facei).Cv(prho[facei], pe[facei] ,pT[facei]);
        }
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tgamma(new scalarField(T.size()));
    scalarField& gamma = tgamma.ref();

    forAll(T, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma(rho[facei], e[facei] , T[facei]);
    }

    return tgamma;
}

// Modified to be function of (rho,e)
// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::volScalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::gamma() const
// {
//     const fvMesh& mesh = this->T_.mesh();

//     tmp<volScalarField> tgamma
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "gamma",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh,
//             dimless
//         )
//     );

//     volScalarField& gamma = tgamma.ref();

//     forAll(this->T_, celli)
//     {
//         gamma[celli] =
//             this->cellMixture(celli).gamma(this->p_[celli], this->T_[celli]);
//     }

//     volScalarField::Boundary& gammaBf = gamma.boundaryFieldRef();

//     forAll(gammaBf, patchi)
//     {
//         const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//         const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//         fvPatchScalarField& pgamma = gammaBf[patchi];

//         forAll(pT, facei)
//         {
//             pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
//             (
//                 pp[facei],
//                 pT[facei]
//             );
//         }
//     }

//     return tgamma;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tgamma
    (
        new volScalarField
        (
            IOobject
            (
                "gamma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless
        )
    );

    volScalarField& gamma = tgamma.ref();

    forAll(this->T_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma(this->rho_[celli], this->e_[celli], this->T_[celli]);
    }

    volScalarField::Boundary& gammaBf = gamma.boundaryFieldRef();

    forAll(gammaBf, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gammaBf[patchi];

        forAll(prho, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                prho[facei],
                pe[facei],
                pT[facei]
            );
        }
    }

    return tgamma;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::Cpv
// (
//     const scalarField& p,
//     const scalarField& T,
//     const label patchi
// ) const
// {
//     tmp<scalarField> tCpv(new scalarField(T.size()));
//     scalarField& Cpv = tCpv.ref();

//     forAll(T, facei)
//     {
//         Cpv[facei] =
//             this->patchFaceMixture(patchi, facei).Cpv(p[facei], T[facei]);
//     }

//     return tCpv;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::Cpv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpv(new scalarField(T.size()));
    scalarField& Cpv = tCpv.ref();

    forAll(T, facei)
    {
        Cpv[facei] =
            this->patchFaceMixture(patchi, facei).Cpv(rho[facei], e[facei], T[facei]);
    }

    return tCpv;
}

// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::volScalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::Cpv() const
// {
//     const fvMesh& mesh = this->T_.mesh();

//     tmp<volScalarField> tCpv
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "Cpv",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh,
//             dimEnergy/dimMass/dimTemperature
//         )
//     );

//     volScalarField& Cpv = tCpv.ref();

//     forAll(this->T_, celli)
//     {
//         Cpv[celli] =
//             this->cellMixture(celli).Cpv(this->p_[celli], this->T_[celli]);
//     }

//     volScalarField::Boundary& CpvBf = Cpv.boundaryFieldRef();

//     forAll(CpvBf, patchi)
//     {
//         const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//         const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//         fvPatchScalarField& pCpv = CpvBf[patchi];

//         forAll(pT, facei)
//         {
//             pCpv[facei] =
//                 this->patchFaceMixture(patchi, facei).Cpv(pp[facei], pT[facei]);
//         }
//     }

//     return tCpv;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::Cpv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCpv
    (
        new volScalarField
        (
            IOobject
            (
                "Cpv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& Cpv = tCpv.ref();

    forAll(this->T_, celli)
    {
        Cpv[celli] =
            this->cellMixture(celli).Cpv(this->rho_[celli], this->e_[celli], this->T_[celli]);
    }

    volScalarField::Boundary& CpvBf = Cpv.boundaryFieldRef();

    forAll(CpvBf, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpv = CpvBf[patchi];

        forAll(prho, facei)
        {
            pCpv[facei] =
                this->patchFaceMixture(patchi, facei).Cpv(prho[facei], pe[facei], pT[facei]);
        }
    }

    return tCpv;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::CpByCpv
// (
//     const scalarField& p,
//     const scalarField& T,
//     const label patchi
// ) const
// {
//     tmp<scalarField> tCpByCpv(new scalarField(T.size()));
//     scalarField& CpByCpv = tCpByCpv.ref();

//     forAll(T, facei)
//     {
//         CpByCpv[facei] =
//             this->patchFaceMixture(patchi, facei).CpByCpv(p[facei], T[facei]);
//     }

//     return tCpByCpv;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::CpByCpv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpByCpv(new scalarField(T.size()));
    scalarField& CpByCpv = tCpByCpv.ref();

    forAll(T, facei)
    {
        CpByCpv[facei] =
            this->patchFaceMixture(patchi, facei).CpByCpv(rho[facei], e[facei], T[facei]);
    }

    return tCpByCpv;
}

// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::volScalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::CpByCpv() const
// {
//     const fvMesh& mesh = this->T_.mesh();

//     tmp<volScalarField> tCpByCpv
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "CpByCpv",
//                 mesh.time().timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh,
//             dimless
//         )
//     );

//     volScalarField& CpByCpv = tCpByCpv.ref();

//     forAll(this->T_, celli)
//     {
//         CpByCpv[celli] = this->cellMixture(celli).CpByCpv
//         (
//             this->p_[celli],
//             this->T_[celli]
//         );
//     }

//     volScalarField::Boundary& CpByCpvBf =
//         CpByCpv.boundaryFieldRef();

//     forAll(CpByCpvBf, patchi)
//     {
//         const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//         const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//         fvPatchScalarField& pCpByCpv = CpByCpvBf[patchi];

//         forAll(pT, facei)
//         {
//             pCpByCpv[facei] = this->patchFaceMixture(patchi, facei).CpByCpv
//             (
//                 pp[facei],
//                 pT[facei]
//             );
//         }
//     }

//     return tCpByCpv;
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::CpByCpv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCpByCpv
    (
        new volScalarField
        (
            IOobject
            (
                "CpByCpv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless
        )
    );

    volScalarField& CpByCpv = tCpByCpv.ref();

    forAll(this->T_, celli)
    {
        CpByCpv[celli] = this->cellMixture(celli).CpByCpv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }

    volScalarField::Boundary& CpByCpvBf =
        CpByCpv.boundaryFieldRef();

    forAll(CpByCpvBf, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpByCpv = CpByCpvBf[patchi];

        forAll(pT, facei)
        {
            pCpByCpv[facei] = this->patchFaceMixture(patchi, facei).CpByCpv
            (
                prho[facei],
                pe[facei],
                pT[facei]
            );
        }
    }

    return tCpByCpv;
}

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();

    forAll(h, celli)
    {
        T[celli] =
            this->cellMixture(cells[celli]).THE(h[celli], p[celli], T0[celli]);
    }

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{

    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();
    forAll(h, facei)
    {
        T[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).THE(h[facei], p[facei], T0[facei]);
    }

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::myHeThermo<BasicThermo, MixtureType>::W
(
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tW
    (
        new volScalarField
        (
            IOobject
            (
                "W",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimMass/dimMoles
        )
    );

    volScalarField& W = tW.ref();
    scalarField& WCells = W.primitiveFieldRef();

    forAll(WCells, celli)
    {
        WCells[celli] = this->cellMixture(celli).W();
    }

    volScalarField::Boundary& WBf = W.boundaryFieldRef();

    forAll(WBf, patchi)
    {
        scalarField& Wp = WBf[patchi];
        forAll(Wp, facei)
        {
            Wp[facei] = this->patchFaceMixture(patchi, facei).W();
        }
    }

    return tW;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::kappa() const
{
    tmp<Foam::volScalarField> kappa(Cp()*this->alpha_);
    kappa.ref().rename("kappa");
    return kappa;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::kappa
// (
//     const label patchi
// ) const
// {
//     return
//         Cp
//         (
//             this->p_.boundaryField()[patchi],
//             this->T_.boundaryField()[patchi],
//             patchi
//         )*this->alpha_.boundaryField()[patchi];
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::myHeThermo<BasicThermo, MixtureType>::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::alphahe() const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpv()*this->alpha_);
    alphaEff.ref().rename("alphahe");
    return alphaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::alphahe(const label patchi) const
{
    return
    this->CpByCpv
    (
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi],
        patchi
    )
   *this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(Cp()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::kappaEff
// (
//     const scalarField& alphat,
//     const label patchi
// ) const
// {
//     return
//         Cp
//         (
//             this->p_.boundaryField()[patchi],
//             this->T_.boundaryField()[patchi],
//             patchi
//         )
//        *(
//            this->alpha_.boundaryField()[patchi]
//          + alphat
//         );
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        Cp
        (
            this->rho_.boundaryField()[patchi],
            this->e_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(
           this->alpha_.boundaryField()[patchi]
         + alphat
        );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpv()*(this->alpha_ + alphat));
    alphaEff.ref().rename("alphaEff");
    return alphaEff;
}


// template<class BasicThermo, class MixtureType>
// Foam::tmp<Foam::scalarField>
// Foam::myHeThermo<BasicThermo, MixtureType>::alphaEff
// (
//     const scalarField& alphat,
//     const label patchi
// ) const
// {
//     return
//     this->CpByCpv
//     (
//         this->p_.boundaryField()[patchi],
//         this->T_.boundaryField()[patchi],
//         patchi
//     )
//    *(
//         this->alpha_.boundaryField()[patchi]
//       + alphat
//     );
// }

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::myHeThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
    this->CpByCpv
    (
        this->rho_.boundaryField()[patchi],
        this->e_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi],
        patchi
    )
   *(
        this->alpha_.boundaryField()[patchi]
      + alphat
    );
}

template<class BasicThermo, class MixtureType>
bool Foam::myHeThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }

    return false;
}


// ************************************************************************* //
