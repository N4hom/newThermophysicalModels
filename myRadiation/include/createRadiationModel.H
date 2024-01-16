Info << "Creating myRadiationModel \n" << endl;
    autoPtr<myRadiation::myRadiationModel> radiation
    (
        myRadiation::myRadiationModel::New(thermo.T())
    );
