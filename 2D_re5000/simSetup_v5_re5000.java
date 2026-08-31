package macro;

import java.util.*;
import java.io.File;
import star.common.*;
import star.base.neo.*;
import star.base.report.*;
import star.material.*;
import star.metrics.*;
import star.coupledflow.*;
import star.turbulence.*;
import star.flow.*;
import star.kwturb.*;
import star.energy.*;
import star.meshing.*;
import star.vis.*;

public class simSetup_v5_re5000 extends StarMacro {

    private Simulation sim;
    private Map<String, FieldFunction> ffMap = new HashMap<>();

    public void execute() {
        sim = getActiveSimulation();

        sim.println("Step 1: Setting up Physics Continuum...");
        setupPhysics();

        sim.println("Step 2: Adding Global Parameters...");
        createParameters();

        sim.println("Step 3: Creating Point Probe & Freestream Reports...");
        createPointProbe();
        createFreestreamReports();

        sim.println("Step 4: Generating Universal Vector Formulas (Version B)...");
        createFormulas();

        sim.println("Step 5: Applying Parameters to Conditions (Manual 1D Flow)...");
        applyParametersToPhysics();

        sim.println("Step 6: Creating Derived Parts (Thresholds & Planes)...");
        createDerivedParts();

        sim.println("Step 7: Generating Grouped Directional Reports (2D Layout)...");
        createRemainingReports();

        sim.println("Step 8: Generating Report & Field Mean Monitors...");
        createMonitors();
        
        sim.println("Macro setup execution complete.");
    }

    private void setupPhysics() {
        ImportManager importManager = sim.getImportManager();
		// change to the mesh location directory
        importManager.importMeshFiles(new StringVector(new String[] {resolvePath("D:\\Cranfield Files\\thesis\\2D vortex shedding cylinder\\mesh\\mesh_cylinder_v3_std.ccm")}), NeoProperty.fromString("{\'FileOptions\': [{\'Sequence\': 42}]}"));
        Region region = sim.getRegionManager().getRegion("Fluid");

        PhysicsContinuum physics = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        physics.enable(SingleComponentGasModel.class);
        physics.enable(CoupledFlowModel.class);
        physics.enable(ConstantDensityModel.class);
        physics.enable(CoupledEnergyModel.class);
        physics.enable(TurbulentModel.class);
        physics.enable(RansTurbulenceModel.class);
        physics.enable(KOmegaTurbulence.class);
        physics.enable(ImplicitUnsteadyModel.class);
        physics.enable(SstKwTurbModel.class);
        physics.enable(KwAllYplusWallTreatment.class);
        physics.enable(TwoDimensionalModel.class);
        
        SingleComponentGasModel gasModel = physics.getModelManager().getModel(SingleComponentGasModel.class);
        Gas gas = ((Gas) gasModel.getMaterial());
        gas.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).setMethod(SutherlandLaw.class);
        Units units_kelvin = ((Units) sim.getUnitsManager().getObject("K"));
        SutherlandLaw suthSetting = ((SutherlandLaw) gas.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).getMethod());
        suthSetting.getSutherlandConstant().setValueAndUnits(110.56, units_kelvin);

        Boundary outlet = region.getBoundaryManager().getBoundary("outlet");
        PressureBoundary pBound = ((PressureBoundary) sim.get(ConditionTypeManager.class).get(PressureBoundary.class));
        outlet.setBoundaryType(pBound);
    }

    private void createParameters() {
        Units m_s = sim.getUnitsManager().getObject("m/s");
        Units dimlessUnits = sim.getUnitsManager().getPreferredUnits(Dimensions.Builder().build());
        Units kelvin = sim.getUnitsManager().getObject("K");
        Units units_blank = ((Units) sim.getUnitsManager().getObject(""));
		Units units_s = ((Units) sim.getUnitsManager().getObject("s"));

        ScalarGlobalParameter tinfParam = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        tinfParam.setPresentationName("Tinf");
        tinfParam.getQuantity().setValueAndUnits(300.0, units_blank);
        tinfParam.setDimensions(Dimensions.Builder().temperature(1).build());

        VectorGlobalParameter vinfParam = (VectorGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(VectorGlobalParameter.class, "Vector");
        vinfParam.setPresentationName("Vinf");
        vinfParam.getQuantity().setComponentsAndUnits(7.844, 0.0, 0.0, units_blank);
        vinfParam.setDimensions(Dimensions.Builder().length(1).time(-1).build());
		
		ScalarGlobalParameter dt = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        dt.setPresentationName("timeStep");
        dt.getQuantity().setValueAndUnits(5.0635e-05, units_s);
        dt.setDimensions(Dimensions.Builder().time(1).build());
		
		ScalarGlobalParameter max_time = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        max_time.setPresentationName("maxTime");
        max_time.getQuantity().setValueAndUnits(1.0, units_s);
        max_time.setDimensions(Dimensions.Builder().time(1).build());
    }

    private void createPointProbe() {
        Region region = sim.getRegionManager().getRegion("Fluid");
        Units mUnits = sim.getUnitsManager().getPreferredUnits(Dimensions.Builder().length(1).build());

        PointPart point = sim.getPartManager().createPointPart(new ArrayList<>(Collections.<NamedObject>emptyList()), new DoubleVector(new double[] {0.0, 0.0, 0.0}), null);
        point.setPresentationName("Freestream_Probe");
        point.getPointCoordinate().setCoordinate(mUnits, mUnits, mUnits, new DoubleVector(new double[] {-0.3, 0.3, 0.0}));
        point.getInputParts().setQuery(null);
        point.getInputParts().setObjects(region);
    }

    private void createFreestreamReports() {
        Part point = sim.getPartManager().getPart("Freestream_Probe");

        MaxReport pRep = (MaxReport) sim.getReportManager().create("star.base.report.MaxReport");
        pRep.setPresentationName("FreestreamAbsolutePressure");
        pRep.setFieldFunction(((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("AbsolutePressure")));
        pRep.getParts().setObjects(point);

        MaxReport dRep = (MaxReport) sim.getReportManager().create("star.base.report.MaxReport");
        dRep.setPresentationName("FreestreamDensity");
        dRep.setFieldFunction(((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Density")));
        dRep.getParts().setObjects(point);

        MaxReport tRep = (MaxReport) sim.getReportManager().create("star.base.report.MaxReport");
        tRep.setPresentationName("FreestreamTemperature");
        tRep.setFieldFunction(((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Temperature")));
        tRep.getParts().setObjects(point);

        MaxReport muRep = (MaxReport) sim.getReportManager().create("star.base.report.MaxReport");
        muRep.setPresentationName("FreestreamDynamicViscosity");
        muRep.setFieldFunction(((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("DynamicViscosity")));
        muRep.getParts().setObjects(point);
    }

    private void createFormulas() {
        // --- Thermodynamic & Reference Formulations ---
        buildFunc("hinf", 0, Dimensions.Builder().length(2).time(-2).build(), "${SpecificHeat} * (${Tinf} - 298.15)");
        buildFunc("einf", 0, Dimensions.Builder().length(2).time(-2).build(), "${hinf} - ${FreestreamAbsolutePressure}/${FreestreamDensity}");
        buildFunc("SpecificEnthalpy", 0, Dimensions.Builder().length(2).time(-2).build(), "${TotalEnthalpy} - 0.5 * mag2($${Velocity})");
        buildFunc("SpecificEnergy", 0, Dimensions.Builder().length(2).time(-2).build(), "${SpecificEnthalpy} - ${AbsolutePressure} / ${Density}");
        buildFunc("deltaEnergy", 0, Dimensions.Builder().length(2).time(-2).build(), "${SpecificHeat} * (${Temperature} - ${FreestreamTemperature})");
		buildFunc("deltaEntropy", 0, Dimensions.Builder().length(2).time(-2).temperature(-1).build(), "(${SpecificHeat} * log(${Temperature} / ${FreestreamTemperature}))");
        buildFunc("xs", 0, Dimensions.Builder().length(2).time(-2).build(), "${deltaEnergy} - ${Tinf} * (${deltaEntropy}) + ${FreestreamAbsolutePressure} * (1 / ${Density} - 1 / ${FreestreamDensity})");
        buildFunc("MassSpec", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "${Density} * ${xs}");
        buildFunc("h0_defect_rho", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "(${Density} * ${SpecificHeat} * (${Temperature} - ${FreestreamTemperature})) + (${AbsolutePressure} - ${FreestreamAbsolutePressure}) + (0.5 * ${Density} * mag2($${PerturbationVelocity}))");
		
		// --- Energy quantities ---
		buildFunc("ek_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2))");
		buildFunc("theta_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "0.0");
        buildFunc("phi", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "(${DynamicViscosity} + ${TurbulentViscosity}) * pow(${StrainRate}, 2)");
		buildFunc("phi_alt_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "${Density} * ${deltaEnergy}");
		buildFunc("phi_alt_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${Density} * ${deltaEnergy} * $${Velocity}, $${Normal})");
		buildFunc("aphi_alt_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "(${Tinf} / ${TotalTemperature}) * ${Density} * ${deltaEnergy}");
		buildFunc("aphi_alt_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((${Tinf} / ${TotalTemperature}) * ${Density} * ${deltaEnergy} * $${Velocity}, $${Normal})");
		
		// --- Enthalpy formulation (NEW! from Drew Sanders)
		buildFunc("enthalpy_storage_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "(${Density} * ${SpecificHeat} * (${Temperature} - ${FreestreamTemperature})) + (0.5 * ${Density} * mag2($${PerturbationVelocity}))");
		buildFunc("enth_flux_enthalpy", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${h0_defect_rho} * $${Velocity}, $${Normal})");
		buildFunc("enth_flux_pressure", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((-1.0) * (${AbsolutePressure} - ${FreestreamAbsolutePressure}) * $${Vinf}, $${Normal})");
		buildFunc("enth_flux_viscous", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((-1.0) * dotVector($$${Tshear}, $${PerturbationVelocity}), $${Normal})");
		
		// --- Exergy Quantities ---
		buildFunc("aphi", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "${Tinf} / ${TotalTemperature} * ${phi}");
        buildFunc("ANablaT", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "${Tinf} * ${ThermalConductivity} * (1 / ${Temperature} * (pow(grad(${Temperature})[0], 2) + pow(grad(${Temperature})[1], 2)))"); 
        
        // --- Tensors & Base Vectors ---
		buildFunc("PerturbationVelocity", 1, Dimensions.Builder().length(1).time(-1).build(), "$${Velocity} - $${Vinf}");
        buildFunc("Tshear", 2, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "[ (2.0 * ${EffectiveViscosity} * $${U_VelocityGrad}[0]) - ((2.0 / 3.0) * ${Density} * ${TurbulentKineticEnergy}) ; (2.0 * ${EffectiveViscosity} * 0.5 * ($${U_VelocityGrad}[1] + $${V_VelocityGrad}[0])), (2.0 * ${EffectiveViscosity} * $${V_VelocityGrad}[1]) - ((2.0 / 3.0) * ${Density} * ${TurbulentKineticEnergy}) ; 0.0, 0.0, 0.0 ]");
		
		// --- Mechanical energy components ---
        buildFunc("ep_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${StaticPressure} * $${PerturbationVelocity}, $${Normal})");
        buildFunc("ea_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(0.5 * ${Density} * pow($${PerturbationVelocity}[0], 2) * $${Velocity}, $${Normal})");
        buildFunc("ev_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(0.5 * ${Density} * pow($${PerturbationVelocity}[1], 2) * $${Velocity}, $${Normal})");
        buildFunc("xth_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${MassSpec} * $${Velocity}, $${Normal})");
        buildFunc("etau_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(dotVector($$${Tshear}, $${PerturbationVelocity}), $${Normal})");
		buildFunc("etp_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2)) * $${Velocity}) + (${StaticPressure} * $${PerturbationVelocity}) - [($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1]), ($$${Tshear}[0,1] * $${PerturbationVelocity}[0] + $$${Tshear}[1,1] * $${PerturbationVelocity}[1]), 0.0], $${Normal})");
		
		// --- X-component Mechanical energy components evaluated as volume (for mean fields) ---
		buildFunc("ep_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "${StaticPressure} * $${PerturbationVelocity}[0]");
        buildFunc("ea_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "0.5 * ${Density} * pow($${PerturbationVelocity}[0], 2) * $${Velocity}[0]");
        buildFunc("ev_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "0.5 * ${Density} * pow($${PerturbationVelocity}[1], 2) * $${Velocity}[0]");
		buildFunc("xth_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "${MassSpec} * $${Velocity}[0]");
        buildFunc("etau_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1])");
		buildFunc("etp_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "(0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2)) * $${Velocity}[0]) + (${StaticPressure} * $${PerturbationVelocity}[0]) - ($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1])");
		
        // --- Box Filters ---
        buildFunc("box_cv1", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.02 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv2", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.05 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv3", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.1 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv4", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.15 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv5", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.2 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv6", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.3 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_cv7", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.5 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");

        buildFunc("box_x_cv1", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_x_cv2", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.05) ? 1.0 : 0.0");
        buildFunc("box_x_cv3", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_x_cv4", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.15) ? 1.0 : 0.0");
        buildFunc("box_x_cv5", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.2) ? 1.0 : 0.0");
        buildFunc("box_x_cv6", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.3) ? 1.0 : 0.0");
        buildFunc("box_x_cv7", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.5) ? 1.0 : 0.0");

        buildFunc("box_y_cv1", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv2", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv3", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv4", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv5", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv6", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_y_cv7", 0, Dimensions.Builder().build(), "($$Position[1] >= -0.1 && $$Position[1] <= 0.1) ? 1.0 : 0.0");
    }

    private void applyParametersToPhysics() {
        PhysicsContinuum physics = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        Region region = sim.getRegionManager().getRegion("Fluid");
        Boundary inlet = region.getBoundaryManager().getBoundary("inlet");

        VelocityProfile initVel = physics.getInitialConditions().get(VelocityProfile.class);
        initVel.getMethod(ConstantVectorProfileMethod.class).getQuantity().setDefinition("$${Vinf}");

        VelocityMagnitudeProfile inletVel = inlet.getValues().get(VelocityMagnitudeProfile.class);
        inletVel.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("$$Vinf[0]");

        StaticTemperatureProfile initTemp = physics.getInitialConditions().get(StaticTemperatureProfile.class);
        initTemp.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${Tinf}");
        
        StaticTemperatureProfile inletTemp = inlet.getValues().get(StaticTemperatureProfile.class);
        inletTemp.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${Tinf}");

        ImplicitUnsteadySolver solver = ((ImplicitUnsteadySolver) sim.getSolverManager().getSolver(ImplicitUnsteadySolver.class));
        solver.getTimeStep().setDefinition("${timeStep}");
        solver.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);
        
        CoupledImplicitSolver cSolver = ((CoupledImplicitSolver) sim.getSolverManager().getSolver(CoupledImplicitSolver.class));
        cSolver.setLeaveTemporaryStorage(true);

        PhysicalTimeStoppingCriterion ptStop = ((PhysicalTimeStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));
        ptStop.getMaximumTime().setDefinition("${maxTime}");
        ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps")).setIsUsed(false);

        AbortFileStoppingCriterion abortFileStop = ((AbortFileStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Stop File"));
        abortFileStop.setIsUsed(false);
    }

    private void createDerivedParts() {
        Region region = null;
        try {
            region = sim.getRegionManager().getRegion("Fluid");
        } catch (Exception e) {
            sim.println("Error: 'Fluid' region not found!");
            return;
        }

        Units mUnits = sim.getUnitsManager().getPreferredUnits(Dimensions.Builder().length(1).build());
        Units dimlessUnits = sim.getUnitsManager().getPreferredUnits(Dimensions.Builder().build());
        LabCoordinateSystem cs = sim.getCoordinateSystemManager().getLabCoordinateSystem();
        
        double[] xBounds = {0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5};

        for (int i = 0; i < 7; i++) {
            int cvNum = i + 1;
            String cvName = "cv" + cvNum;

            ThresholdPart threshVol = sim.getPartManager().createThresholdPart(
                new ArrayList<>(Arrays.<NamedObject>asList(region)), 
                new DoubleVector(new double[] {1.0, 1.0}), 
                dimlessUnits, ffMap.get("box_" + cvName), 0, null);
            threshVol.setPresentationName(cvName);

            ThresholdPart threshX = sim.getPartManager().createThresholdPart(
                new ArrayList<>(Arrays.<NamedObject>asList(region)), 
                new DoubleVector(new double[] {1.0, 1.0}), 
                dimlessUnits, ffMap.get("box_x_" + cvName), 0, null);
            threshX.setPresentationName("thresh_x_" + cvName);

            ThresholdPart threshY = sim.getPartManager().createThresholdPart(
                new ArrayList<>(Arrays.<NamedObject>asList(region)), 
                new DoubleVector(new double[] {1.0, 1.0}), 
                dimlessUnits, ffMap.get("box_y_" + cvName), 0, null);
            threshY.setPresentationName("thresh_y_" + cvName);

            buildPlaneSection(threshX, cs, mUnits, "bot_" + cvName, new double[]{0.0, -0.1, 0.0}, new double[]{0.0, -1.0, 0.0});
            buildPlaneSection(threshX, cs, mUnits, "top_" + cvName, new double[]{0.0, 0.1, 0.0}, new double[]{0.0, 1.0, 0.0});
            buildPlaneSection(threshY, cs, mUnits, "right_" + cvName, new double[]{xBounds[i], 0.0, 0.0}, new double[]{1.0, 0.0, 0.0});
        }
        
        ThresholdPart threshY_cv7 = (ThresholdPart) sim.getPartManager().getObject("thresh_y_cv7");
        buildPlaneSection(threshY_cv7, cs, mUnits, "left", new double[]{-0.1, 0.0, 0.0}, new double[]{-1.0, 0.0, 0.0});
    }

    private void buildPlaneSection(Part inputPart, LabCoordinateSystem cs, Units units, String name, double[] origin, double[] normal) {
        PlaneSection plane = (PlaneSection) sim.getPartManager().createImplicitPart(
            new ArrayList<>(Collections.<NamedObject>emptyList()), 
            new DoubleVector(new double[] {0.0, 0.0, 1.0}), 
            new DoubleVector(new double[] {0.0, 0.0, 0.0}), 
            0, 1, new DoubleVector(new double[] {0.0}), null);
        
        plane.setPresentationName(name);
        plane.setCoordinateSystem(cs);
        plane.getInputParts().setQuery(null);
        plane.getInputParts().setObjects(inputPart);
        plane.getOriginCoordinate().setValue(new DoubleVector(origin));
        plane.getOriginCoordinate().setCoordinate(units, units, units, new DoubleVector(origin));
        plane.getOriginCoordinate().setCoordinateSystem(cs);
        plane.getOrientationCoordinate().setValue(new DoubleVector(normal));
        plane.getOrientationCoordinate().setCoordinate(units, units, units, new DoubleVector(normal));
        plane.getOrientationCoordinate().setCoordinateSystem(cs);
    }

    private void createRemainingReports() {
        
        // --- 1. 2D "VOLUME" INTEGRALS ---
        String[] volParts = {"cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7"};
        String[][] volMetrics = {
            {"massspec_vol", "MassSpec"},
			{"anabla_vol", "ANablaT"},
            {"ek_vol", "ek_vol"},
			{"theta_vol", "theta_vol"},
            {"phi_vol", "phi"},
			{"aphi_vol", "aphi"},
			{"enthalpy_storage_vol", "enthalpy_storage_vol"},
			{"phi_alt_vol", "phi_alt_vol"},
			{"aphi_alt_vol", "aphi_alt_vol"}
        };

        for (String[] metric : volMetrics) {
            createVolumeReport(metric[0], metric[1], volParts); 
        }

        // --- 2. 1D "SURFACE" INTEGRALS (Standard Fluxes) ---
        String[][] surfaceMetrics = {
            {"ep", "ep_flux"},
			{"ea", "ea_flux"},
			{"ev", "ev_flux"},
			{"etp", "etp_flux"},
            {"xth", "xth_flux"},
			{"Etau", "etau_flux"},
			{"enth_flux_enthalpy", "enth_flux_enthalpy"},
			{"enth_flux_pressure", "enth_flux_pressure"},
			{"enth_flux_viscous", "enth_flux_viscous"},
			{"phi_alt_flux", "phi_alt_flux"},
			{"aphi_alt_flux", "aphi_alt_flux"}
        };

        String[] botParts = {"bot_cv1", "bot_cv2", "bot_cv3", "bot_cv4", "bot_cv5", "bot_cv6", "bot_cv7"};
        String[] topParts = {"top_cv1", "top_cv2", "top_cv3", "top_cv4", "top_cv5", "top_cv6", "top_cv7"};
        String[] rightParts = {"right_cv1", "right_cv2", "right_cv3", "right_cv4", "right_cv5", "right_cv6", "right_cv7"};
        String[] leftParts = {"left"};

        for (String[] metric : surfaceMetrics) {
            String repPrefix = metric[0];
            String funcKey = metric[1];

            createLineReport(repPrefix + "_bot", funcKey, botParts); 
            createLineReport(repPrefix + "_top", funcKey, topParts);
            createLineReport(repPrefix + "_right", funcKey, rightParts);
            createLineReport(repPrefix + "_left", funcKey, leftParts);
        }

        // --- 3. DRAG & LIFT FORCE REPORTS ---
        try {
            Boundary boundary = sim.getRegionManager().getRegion("Fluid").getBoundaryManager().getBoundary("wall_cylinder");
            Units uArea = sim.getUnitsManager().getObject("m^2");
            Units uVel = sim.getUnitsManager().getObject("m/s");
            Units uDens = sim.getUnitsManager().getObject("kg/m^3");

            ForceReport drag_total = sim.getReportManager().createReport(ForceReport.class);
            drag_total.setPresentationName("drag");
            drag_total.getDirection().setComponents(1.0, 0.0, 0.0);
            drag_total.getParts().setObjects(boundary);

            ForceReport lift = sim.getReportManager().createReport(ForceReport.class);
            lift.setPresentationName("lift");
            lift.getDirection().setComponents(0.0, 1.0, 0.0);
            lift.getParts().setObjects(boundary);

        } catch (Exception e) {
			sim.println("Warning: Could not create Drag / Lift reports");
		}
    }

    private void createMonitors() {
        MonitorManager monitorManager = sim.getMonitorManager();
        Region fluidRegion = sim.getRegionManager().getRegion("Fluid");
        
        String[] meanVars = {
            "ea_vol", "ev_vol", "ep_vol", "etau_vol", "xth_vol", "aphi", 
            "ANablaT", "theta_vol", "phi", "etp_vol", "MassSpec"
        };

        for (String var : meanVars) {
            try {
                // 1. Correctly look up the parameters using getObject() and get their numeric values
				ScalarGlobalParameter maxTimeParam = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("maxTime");
				ScalarGlobalParameter timeStepParam = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("timeStep");

				double maxTimeVal = maxTimeParam.getQuantity().getValue();
				double timeStepVal = timeStepParam.getQuantity().getValue();

				// 2. Compute the exact integer time steps using Java math
				int startTimeStep = (int) Math.round((maxTimeVal - 0.6) / timeStepVal);
				int stopTimeStep = (int) Math.round(maxTimeVal / timeStepVal);

				// 3. Create the monitor and apply the integer values directly using .setValue()
				FieldMeanMonitor meanMon = (FieldMeanMonitor) monitorManager.create("star.base.report.FieldMeanMonitor");
				meanMon.setPresentationName("Mean_" + var);
				meanMon.getParts().setObjects(fluidRegion);
				meanMon.setFieldFunction((UserFieldFunction) sim.getFieldFunctionManager().getFunction(var));

				StarUpdate starUpdate = meanMon.getStarUpdate();
				TimeStepUpdateFrequency tsUpdateFreq = starUpdate.getTimeStepUpdateFrequency();

				// Set the start step
				IntegerValue startInt = tsUpdateFreq.getStartTimeStepQuantity();
				startInt.getQuantity().setValue(startTimeStep); 

				// Set the stop step
				IntegerValue stopInt = tsUpdateFreq.getStopTimeStepQuantity();
				stopInt.getQuantity().setValue(stopTimeStep);
                
            } catch (Exception e) {
                sim.println("Warning: Could not create Mean monitor for " + var);
            }
        }

        // --- 2. REPORT MONITORS ---
        int count = 0;
        for (Report report : sim.getReportManager().getObjects()) {
            try {
                ReportMonitor monitor = report.createMonitor();
                if (report instanceof SurfaceIntegralReport || report instanceof LineIntegralReport) {
                    try { monitor.getValueType().setSelected(ReportMonitorValueType.Type.PartValue); } catch (Exception ex) {}
                }
                monitor.setPlotLimit(200000);
                count++;
            } catch (Exception e) {}
        }
        
        // --- 3. RESIDUAL MONITORS ---
        try {
            ((ResidualMonitor) monitorManager.getMonitor("Continuity")).setPlotLimit(200000);
            ((ResidualMonitor) monitorManager.getMonitor("Energy")).setPlotLimit(200000);
            ((ResidualMonitor) monitorManager.getMonitor("X-momentum")).setPlotLimit(200000);
            ((ResidualMonitor) monitorManager.getMonitor("Y-momentum")).setPlotLimit(200000);
            ((IterationMonitor) monitorManager.getMonitor("Iteration")).setPlotLimit(200000);
            ((PhysicalTimeMonitor) monitorManager.getMonitor("Physical Time")).setPlotLimit(200000);
        } catch (Exception e) {}
        
        sim.println("Successfully generated " + count + " report monitors and all Field Mean monitors.");
    }

    // --- Utility Methods ---
    private NamedObject getPartByName(String name) {
        for (Part p : sim.getPartManager().getObjects()) if (p.getPresentationName().equals(name)) return p;
        return null;
    }
    
    private FieldFunction getFieldFunction(String key) {
        if (ffMap.containsKey(key)) return ffMap.get(key);
        try {
            return sim.getFieldFunctionManager().getFunction(key);
        } catch (Exception e) {
            return null;
        }
    }

    private void createVolumeReport(String reportName, String funcKey, String[] partNames) {
        try {
            SurfaceIntegralReport report = sim.getReportManager().createReport(SurfaceIntegralReport.class);
            report.setPresentationName(reportName); 
            FieldFunction ff = getFieldFunction(funcKey);
            if (ff != null) report.setFieldFunction(ff);
            setReportParts(report, partNames);
        } catch (Exception e) {
            sim.println("Error generating 2D Volume (Surface) report (" + reportName + "): " + e.getMessage());
        }
    }

    private void createLineReport(String reportName, String funcKey, String[] partNames) {
        try {
            LineIntegralReport report = sim.getReportManager().createReport(LineIntegralReport.class);
            report.setPresentationName(reportName); 
            FieldFunction ff = getFieldFunction(funcKey);
            if (ff != null) report.setFieldFunction(ff);
            setReportParts(report, partNames);
        } catch (Exception e) {
            sim.println("Error generating 1D Surface (Line) report (" + reportName + "): " + e.getMessage());
        }
    }

    private void setReportParts(Report report, String[] partNames) {
        List<NamedObject> partList = new ArrayList<>();
        for (String pName : partNames) {
            NamedObject p = getPartByName(pName);
            if (p != null) partList.add(p);
        }
        if (!partList.isEmpty()) {
            NamedObject[] partArray = partList.toArray(new NamedObject[0]);
            if (report instanceof SurfaceIntegralReport) {
                ((SurfaceIntegralReport) report).getParts().setObjects(partArray);
            } else if (report instanceof LineIntegralReport) {
                ((LineIntegralReport) report).getParts().setObjects(partArray);
            }
        }
    }

    private void buildFunc(String name, int type, Dimensions dim, String def) {
        try {
            UserFieldFunction ff = sim.getFieldFunctionManager().createFieldFunction();
            if (type == 0) ff.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);
            else if (type == 1) ff.getTypeOption().setSelected(FieldFunctionTypeOption.Type.VECTOR);
            else if (type == 2) ff.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);
            
            ff.setPresentationName(name);
            ff.setFunctionName(name);
            ff.setDimensions(dim);
            ff.setDefinition(def);
            ffMap.put(name, ff);
        } catch (Exception e) {}
    }
}