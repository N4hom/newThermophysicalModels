/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::myBasicThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

    An important concept is that the thermo can share an existing T
    (similar to p) in which case it will not try to update it. This gets
    triggered purely on construction order - the first one to register
    T is responsible for updating it. Note that the mechanism also means
    that if multiple thermos are running on the same mesh, only the first one
    will update the temperature. The behaviour can be overridden using the
    'updateT' dictionary entry.

SourceFiles
    myBasicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef myBasicThermo_H
#define myBasicThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class myBasicThermo Declaration
\*---------------------------------------------------------------------------*/

class myBasicThermo
:
    public IOdictionary
{
    // Private Data

        //- Components names/order
        static const wordList componentHeader4;

        //- Components names/order
        static const wordList componentHeader7;


    // Private Member Functions

        //- Construct name of thermo package from dictionary components
        static word makeThermoName
        (
            const dictionary& dict,
            const wordList*& cmptHeaderPtr
        );

        //- Look up field from registry or construct and store
        static volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const word& fieldName,
            bool& isOwner  //!< Stored to registry by this instance
        );

        //- Generic lookup for thermodynamics package thermoTypeName
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoTypeDict,
            ThermoConstructTable& thermoTable,
            const word& thermoTypeName,
            const wordList& cmptNames
        );

        //- Generic lookup for each of the related thermodynamics packages
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoDict,
            ThermoConstructTable& thermoTable
        );


protected:

    // Protected Data

        //- Phase-name
        const word& phaseName_;

        //- Pressure created and stored by this instance
        bool pOwner_;

        //- Temperature created and stored by this instance
        bool TOwner_;

        //- Include dpdt term in the enthalpy equation?
        bool dpdt_;


    // Fields

        //- Pressure [Pa]
        volScalarField& p_;

        //- Temperature [K]
        volScalarField& T_;

        //- Laminar thermal diffusivity [kg/m/s]
        volScalarField alpha_;
        
        //- Internal energy
        volScalarField e_;

        //- Density field [kg/m^3]
        volScalarField rho_;

        //- Heat capacity at constant pressure [J/kg/K]
        volScalarField Cp_;

        //- Heat capacity at constant volume [J/kg/K]
        volScalarField Cv_;



    // Protected Member Functions

        //- Return the enthalpy/internal energy field boundary types
        //- by interrogating the temperature field boundary types
        wordList heBoundaryTypes();

        //- Return the enthalpy/internal energy field boundary base types
        //- by interrogating the temperature field boundary types
        wordList heBoundaryBaseTypes();

         //- Temperature [K]
            //  Non-const access allowed for transport equations
            volScalarField& TRef()
            {
                return this->T_;
            }


            volScalarField& pRef()
            {
                return this->p_;
            }


            //- Heat capacity at constant pressure [J/kg/K]
            volScalarField& CpRef()
            {
                return Cp_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            volScalarField& CvRef()
            {
                return Cv_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            volScalarField& alphaRef()
            {
                return alpha_;
            }


            
        //- No copy construct
        myBasicThermo(const myBasicThermo&) = delete;


public:

    //- Runtime type information
    TypeName("myBasicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        myBasicThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        myBasicThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        myBasicThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        myBasicThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh, phase name and explicit naming of the
        //- dictionary (so it can be shared amongst phases).
        myBasicThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    // Selectors

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName=word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );

        //- Specialisation of the Generic New for myBasicThermo
        static autoPtr<myBasicThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~myBasicThermo();


    // Member Functions
        

        static const word dictName;

        static word phasePropertyName
        (
            const word& name,
            const word& phaseName
        )
        {
            return IOobject::groupName(name, phaseName);
        }

        word phasePropertyName(const word& name) const
        {
            return IOobject::groupName(name, phaseName_);
        }

        static const myBasicThermo& lookupThermo(const fvPatchScalarField& pf);

        //- Print (filtered) table of thermo names, splits on \c " ,<>"
        static Ostream& printThermoNames
        (
            Ostream& os,
            const wordList& cmptNames,
            const wordList& thermoNames
        );

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&,
            const word&
        ) const;

        //- Split thermo package name into a list of components names
        //  Splits on \c " ,<>"
        //  \return empty list if the split name does not have the
        //  expected number of components (non-zero).
        static wordList splitThermoName
        (
            const std::string& thermoName,
            const int nExpectedCmpts
        );


        //- Update properties
        virtual void correct() = 0;

        virtual void correctFromRhoE() = 0;

        //- Return the name of the thermo physics
        virtual word thermoName() const = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- True if dpdt term should be included in enthalpy equation
        bool dpdt() const noexcept
        {
            return dpdt_;
        }

        //- True if T should be updated
        bool updateT() const noexcept
        {
            return TOwner_;
        }


        // Access to thermodynamic state variables

            //- Return enthalpy/Internal energy [J/kg]
            
            volScalarField& eRef()
            {
                return this->e_;
            }

            const volScalarField& e()
            {
                return e_;
            }

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;


            //- Density [kg/m^3]
            //- Non-const access 
            virtual volScalarField& getRho();

            //const volScalarField& rhoRef();

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he() = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const = 0;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& T();

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            // virtual tmp<scalarField> Cp
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const label patchi
            // ) const = 0;

            virtual tmp<scalarField> Cp
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity using pressure and temperature [J/kg/K]
            // virtual tmp<scalarField> Cp
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const labelList& cells
            // ) const = 0;

            virtual tmp<scalarField> Cp
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            // virtual tmp<scalarField> Cv
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const label patchi
            // ) const = 0;

            virtual tmp<scalarField> Cv
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Density from pressure and temperature from EoS
            virtual tmp<scalarField> rhoEoS
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const = 0;

            //- Gamma = Cp/Cv for patch []
            // virtual tmp<scalarField> gamma
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const label patchi
            // ) const = 0;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const label patchi
            ) const = 0;



            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            // virtual tmp<scalarField> Cpv
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const label patchi
            // ) const = 0;

            virtual tmp<scalarField> Cpv
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const label patchi
            ) const = 0;

            

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const = 0;

            //- Heat capacity ratio for patch []
            // virtual tmp<scalarField> CpByCpv
            // (
            //     const scalarField& p,
            //     const scalarField& T,
            //     const label patchi
            // ) const = 0;

            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;


        // Access to transport state variables

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual const scalarField& alpha
            (
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const = 0;

             //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const = 0;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const = 0;

            //- Effective thermal diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;




        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myBasicThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
