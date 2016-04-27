void sgp::CreateOrbitSGPExperiment()
{
    // Load a previously saved SGP
    if (!sgp::LoadSGP(gRun.loadSceneFromFile)) {
        //ncc__error("Could not load SGP from file; experiment aborted.\a");
        ncc__warning("Could not load SGP from file; experiment aborted.");
        CreateRubbleGrain(PxVec3(1,0,0),eSPHERE_GRAIN,0.5,*gPhysX.mDefaultMaterial,100);
//        CreateRubbleGrain(PxVec3(-1,0,0),eSPHERE_GRAIN,0.5,*gPhysX.mDefaultMaterial,100);
        RecenterScene(); // put center-of-mass at origin
    }
    DeadStop(); // stomp any residual velocities
    PxVec3 d = sgp::FindSGPCenterOfMass();
    RelocateScene(-d);
    if (sgp::orbit.sgpMass > 0) sgp::ReMassSGP(sgp::orbit.sgpMass);

    if (sgp::orbit.bPregenOrbit)
    {
        // Try to load from file run_base_name.orb (or make our own, some day)
        sgp::orbit.orbFile = gRun.workingDirectory + "/" + gRun.baseName + ".orb";
        nr3::MatDoub raw;
        if (ncc::load(sgp::orbit.orbFile.c_str(), &raw, 5))
        {
            // Distribute to vectors (kind of pointless really)
            sgp::orbit.tvec.resize(raw.nrows());
            sgp::orbit.xvec.resize(raw.nrows());
            sgp::orbit.yvec.resize(raw.nrows());
            for (int k=0; k<raw.nrows(); k++) {
                sgp::orbit.tvec[k] = raw[k][0];
                sgp::orbit.xvec[k] = raw[k][1];
                sgp::orbit.yvec[k] = raw[k][2];
            }

            // Estimate simulation length
            double t0 = sgp::orbit.tvec[0];
            double tf = sgp::orbit.tvec[raw.nrows() - 1];
            cout << "Requested orbit integration time " << (tf-t0) << " (cu) in " << (int)((tf-t0)/gSim.timeStep) << " time steps." << endl;
            if ((tf-t0)/gSim.timeStep > 2e4)
                ncc__warning("Long orbit requested -- this might take a while\a");
            sgp::orbit.tStart = t0;
            sgp::orbit.tEnd = tf;

            // Save the initial location and orbit parameters
            sgp::orbit.X0.x = sgp::orbit.xvec[0];
            sgp::orbit.X0.y = sgp::orbit.yvec[0];
            double q2 = 1e200;
            for (int k=0; k<sgp::orbit.xvec.size(); k++)
            {
                double x = sgp::orbit.xvec[k];
                double y = sgp::orbit.yvec[k];
                double r2 = x*x + y*y;
                if (r2 < q2) q2 = r2;
            }
            sgp::orbit.periapse = sqrt(q2);

            sgp::orbit.type = sgp::orbit.eASCII;
        } 
        else
        {
            sgp::GenerateOrbit(); // implement some day
        }
    } 
    else
    {
        // Determine initial conditions, based on orbit type
        if (sgp::orbit.type == sgp::orbit.eBOUND)
        {
            // Determine orbital parameters
            double q = sgp::orbit.periapse;
            double e = sgp::orbit.eccentricity;
            double a = q/(1 - e); // semi-major axis
            double Q = a*(1 + e); // apoapse
            double bigG = sgp::cunits.bigG;
            double bigM = sgp::orbit.bigM;
            double P = 2*PxPi*PxSqrt(a*a*a/bigM/bigG);

            // Estimate simulation length
            double t0 = 0;
            double tf = sgp::orbit.nbOrbits*P;
            cout << "Requested orbit integration time " << (tf-t0) << " (cu) in " << (int)((tf-t0)/gSim.timeStep) << " time steps." << endl;
            if ((tf-t0)/gSim.timeStep < 1e4)
                ncc__warning("Forward Euler integrator: consider taking smalller time steps.\a");
            if ((tf-t0)/gSim.timeStep > 2e4)
                ncc__warning("Long orbit requested -- this might take a while\a");
            sgp::orbit.tStart = t0;
            sgp::orbit.tEnd = tf;

            // Set initial position and velocity (bound, counterclockwise, orbit; starts at apoapse on positive x-axis)
            double vap = PxSqrt(bigG*bigM/a)*PxSqrt((1 - e)/(1 + e)); // v at apoapse
            sgp::orbit.X0 = PxVec3(Q,0,0);
            sgp::orbit.V0 = PxVec3(0,vap,0);

        }
        else if (sgp::orbit.type == sgp::orbit.eHYPERBOLIC)
        {
            // Determine orbital parameters
            double bigG = sgp::cunits.bigG;
            double bigM = sgp::orbit.bigM;
            double q = sgp::orbit.periapse;
            double vinf = sgp::orbit.v_inf;
            double a = bigG*bigM/vinf/vinf;            // semi-major axis
            double E = 0.5*vinf*vinf;                  // orbital energy per-unit-mass
            double el = PxSqrt(2*E + 2*bigG*bigM/q)*q; // orbital angular momentum per-unit-mass
            double e = 1 + (q*vinf*vinf)/(bigG*bigM);  // orbital eccentricity

            // Estimate Roche limit and requested start/end distance
            double rhoBulk = sgp::SGPBulkDensity();
            double roche = 1.51*PxPow(bigM/rhoBulk,1.0/3.0);
            double dStart = roche*sgp::orbit.rocheFactorInitial;
            double dEnd = roche*sgp::orbit.rocheFactorFinal;
            dStart = dEnd = 2*q; //DEBUG DEBUG DEBUG
            if (dStart < q || dEnd < q)
            {
                ncc__error("Requested orbit inside periapsis; experiment aborted.\a");
            }
            if (roche < q)
            {
                ncc__warning("Requested orbit completely outside roche limit; boring!");
            }

            // Convert to start/end eccentric anomaly
            double coshFStart = (1 + dStart/a)/e;
            double coshFEnd   = (1 + dEnd/a)/e;
            double FStart     = -PxLog(coshFStart + PxSqrt(coshFStart*coshFStart - 1)); // pre-periapse negative angle
            double FEnd       =  PxLog(coshFEnd   + PxSqrt(coshFEnd*coshFEnd - 1)); // post-periapse positive angle

            // And to start/end time
            double t0 = PxSqrt(a*a*a/bigG/bigM)*(e*sinh(FStart) - FStart);
            double tf = PxSqrt(a*a*a/bigG/bigM)*(e*sinh(FEnd) - FEnd);
            cout << "Requested orbit integration time " << (tf-t0) << " (cu) in " << (int)((tf-t0)/gSim.timeStep) << " time steps." << endl;
            if ((tf-t0)/gSim.timeStep > 1e4)
                ncc__warning("Long orbit requested -- this might take a while\a");
            sgp::orbit.tStart = t0;
            sgp::orbit.tEnd = tf;

            // And finally to initial coordinates and velocity
            double costeta = (el*el/(bigG*bigM*dStart) - 1)/e;
            double sinteta = PxSqrt(1 - costeta*costeta);
            sgp::orbit.X0 = PxVec3(dStart*costeta, dStart*sinteta, 0);
            double rdot = PxSqrt(2*E - el*el/dStart/dStart + 2*bigG*bigM/dStart);
            double rtetadot = el/dStart;
            double vx = -costeta*rdot + sinteta*rtetadot;
            double vy = -sinteta*rdot - costeta*rtetadot;
            sgp::orbit.V0 = PxVec3(vx,vy,0);
        }

        // Launch the sgp!
        KickActors(sgp::orbit.V0);
    }

    // Create the gravitator (actor representing the primary) located at -X0
    PxRigidDynamic* center = CreateRubbleGrain(-sgp::orbit.X0,eSPHERE_GRAIN,1,*gPhysX.mDefaultMaterial);
    center->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);
    center->setMass(sgp::orbit.bigM); // we don't care about the inertia for this kinematic actor
    ColorActor(center, ncc::rgb::oOrange);
    sgp::VIPs.gravitator = center;

    // Move the camera to a good location
    SpyOnSGP();
    sgp::orbit.bTrackingCamera = true;

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: ORBIT_SGP (" << sgp::eExperimentType << ")" << endl;
        header << "# Orbit type: " << sgp::orbit.type << endl;
        header << "# Central mass = " << sgp::orbit.bigM << " (cu)" << endl;
        if (sgp::orbit.type == sgp::orbit.eBOUND)
            header << "# Periapse = " << sgp::orbit.periapse << " (cu); eccentricity = " << sgp::orbit.eccentricity << " (cu)" << endl;
        if (sgp::orbit.type == sgp::orbit.eHYPERBOLIC)
            header << "# Periapse = " << sgp::orbit.periapse << " (cu); V_inf = " << sgp::orbit.v_inf << " (cu)" << endl;
        header << "# Rubble elements = " << gExp.rubbleCount << ";  Bulk density ~ " << SGPBulkDensity() << " (cu)" << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [time]    [CoM X]    [CoM Y]    [CoM Z]    [SGP bulk density]" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gSim.isRunning=true;
    gSim.bPause=false;
    gSim.codeTime = 0.0;

}
