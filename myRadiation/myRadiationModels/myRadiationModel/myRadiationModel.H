/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Namespace
    Foam::myRadiationModel

Description
    Namespace for radiation modelling

Class
    Foam::myRadiation::myRadiationModel

Description
    Top level model for radiation modelling

SourceFiles
    myRadiation.C
    myRadiationNew.C

\*---------------------------------------------------------------------------*/

#ifndef myRadiation_H
#define myRadiation_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myBasicThermo;
class fvMesh;

namespace myRadiation
{

// Forward declaration of classes
//class absorptionEmissionModel;
class scatterModel;
class sootModel;

/*---------------------------------------------------------------------------*\
                       Class myRadiation Declaration
\*---------------------------------------------------------------------------*/

class myRadiationModel
:
    public IOdictionary
{
public:

    // Static data

        //- Static name external radiative fluxes
        static const word externalRadHeatFieldName_;

        //- Static name for primary solar fluxes
        static const word primaryFluxName_;

        //- Static name for reflected solar fluxes
        static const word relfectedFluxName_;


protected:



    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Radiation model on/off flag
        Switch radiation_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<absorptionEmissionModel> absorptionEmission_;

            //- Scatter model
            autoPtr<scatterModel> scatter_;

            //- Soot model
            //autoPtr<sootModel> soot_;

            //- Transmissivity model
            //autoPtr<transmissivityModel> transmissivity_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- No copy construct
        myRadiationModel(const myRadiationModel&) = delete;

        //- No copy assignment
        void operator=(const myRadiationModel&) = delete;


public:

    //- Runtime type information
    TypeName("myRadiation");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            myRadiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            myRadiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        myRadiationModel(const volScalarField& T);

        //- Construct from components
        myRadiationModel(const word& type, const volScalarField& T);

        //- Construct from components
        myRadiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<myRadiationModel> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<myRadiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~myRadiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Radiation model on/off flag
            const Switch radiation() const
            {
                return radiation_;
            }

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const myBasicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                tmp<volScalarField> rhoCp,
                volScalarField& T
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                volScalarField& T
            ) const;

            virtual label nBands() const = 0;

            //- Access to absorptionEmission model
            const absorptionEmissionModel& absorptionEmission() const;

            //- Access to soot Model
            const sootModel& soot() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        myRadiationModel,                                                        \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        myRadiationModel,                                                        \
        model,                                                                 \
        T                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myRadiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
