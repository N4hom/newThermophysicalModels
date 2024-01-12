Info << "Creating radiation model \n" << endl;
    autoPtr<myRadiation::myRadiationModel> radiation
    (
        myRadiation::myRadiationModel::New(thermo.T())
    );
