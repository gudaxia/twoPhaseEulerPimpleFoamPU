{
    DDtUa =
        fvc::ddt(Ua)
      + fvc::div(phia, Ua)
      - fvc::div(phia)*Ua;
    // mrfZones.addCoriolis(Ua, DDtUa);

    DDtUb =
        fvc::ddt(Ub)
      + fvc::div(phib, Ub)
      - fvc::div(phib)*Ub;
    // mrfZones.addCoriolis(Ub, DDtUb);
}
