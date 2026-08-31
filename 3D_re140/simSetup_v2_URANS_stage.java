// Simcenter STAR-CCM+ macro: simSetup_COUPLED_IDDES_v3.java
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

public class simSetup_v2_URANS_stage extends StarMacro {

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

        sim.println("Step 4: Generating Universal Vector Formulas (3D Version)...");
        createFormulas();

        sim.println("Step 5: Applying Parameters to Conditions...");
        applyParametersToPhysics();

        sim.println("Step 6: Creating Derived Parts (Thresholds & Planes)...");
        createDerivedParts();

        sim.println("Step 7: Generating Grouped Directional Reports (3D Layout)...");
        createRemainingReports();

        sim.println("Step 8: Generating Report & Field Mean Monitors...");
        createMonitors();
        
        sim.println("URANS stage setup is complete. Run stage_setup_laminar.java to setup the multi-stage config.");
    }

    private void setupPhysics() {
        ImportManager importManager = sim.getImportManager();
        importManager.importMeshFiles(new StringVector(new String[] {resolvePath("D:\\Cranfield Files\\thesis\\3D vortex shedding cylinder\\mesh\\mesh.ccm")}), NeoProperty.fromString("{'FileOptions': [{'Sequence': 42}]}"));
        
        FvRepresentation fvRep = ((FvRepresentation) sim.getRepresentationManager().getObject("Volume Mesh"));
        Region region = sim.getRegionManager().getRegion("Fluid");
        
        fvRep.generateCompactMeshReport(new ArrayList<>(Arrays.<Region>asList(region)));
        
        // --- Rename Boundaries ---
        Boundary inletBound = region.getBoundaryManager().getBoundary("Inlet");
        inletBound.setPresentationName("inlet");
        
        Boundary outletBound = region.getBoundaryManager().getBoundary("Outlet");
        outletBound.setPresentationName("outlet");
        
        Boundary wallBound = region.getBoundaryManager().getBoundary("Wall_Cyl");
        wallBound.setPresentationName("wall_cylinder");
        
        // --- Setup Periodic BC ---
        Boundary periodic1 = region.getBoundaryManager().getBoundary("Periodic1");
        Boundary periodic2 = region.getBoundaryManager().getBoundary("Periodic2");
        
        BoundaryInterface periodicInterface = sim.getInterfaceManager().createBoundaryInterface(periodic2, periodic1, "Interface");
        periodicInterface.getTopology().setSelected(InterfaceConfigurationOption.Type.PERIODIC);
        
        InterfacePeriodicTransformSpecification periodicTransform = periodicInterface.getPeriodicTransform();
        periodicTransform.getPeriodicityOption().setSelected(PeriodicityOption.Type.TRANSLATION);

        // --- Physics Continuum Setup ---
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
        physics.enable(ThreeDimensionalModel.class);

        SingleComponentGasModel gasModel = physics.getModelManager().getModel(SingleComponentGasModel.class);
        Gas gas = ((Gas) gasModel.getMaterial());
        gas.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).setMethod(SutherlandLaw.class);
        
        Units units_kelvin = ((Units) sim.getUnitsManager().getObject("K"));
        SutherlandLaw suthSetting = ((SutherlandLaw) gas.getMaterialProperties().getMaterialProperty(DynamicViscosityProperty.class).getMethod());
        suthSetting.getSutherlandConstant().setValueAndUnits(110.56, units_kelvin);

        PressureBoundary pBound = ((PressureBoundary) sim.get(ConditionTypeManager.class).get(PressureBoundary.class));
        outletBound.setBoundaryType(pBound);
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
        vinfParam.getQuantity().setComponentsAndUnits(0.2196, 0.0, 0.0, units_blank);
        vinfParam.setDimensions(Dimensions.Builder().length(1).time(-1).build());
		
        ScalarGlobalParameter dt_DES = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        dt_DES.setPresentationName("timeStep_DES");
        dt_DES.getQuantity().setValueAndUnits(1.25E-4, units_s); // IDDES timestep
        dt_DES.setDimensions(Dimensions.Builder().time(1).build());
		
        ScalarGlobalParameter max_time_DES = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        max_time_DES.setPresentationName("maxTime_DES");
        max_time_DES.getQuantity().setValueAndUnits(7.0, units_s); // IDDES max time
        max_time_DES.setDimensions(Dimensions.Builder().time(1).build());
		
		ScalarGlobalParameter dt_URANS = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        dt_URANS.setPresentationName("timeStep_URANS");
        dt_URANS.getQuantity().setValueAndUnits(1.25E-3, units_s); // URANS timestep
        dt_URANS.setDimensions(Dimensions.Builder().time(1).build());
		
        ScalarGlobalParameter max_time_URANS = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "Scalar");
        max_time_URANS.setPresentationName("maxTime_URANS");
        max_time_URANS.getQuantity().setValueAndUnits(5.0, units_s); // URANS max time
        max_time_URANS.setDimensions(Dimensions.Builder().time(1).build());
    }

    private void createPointProbe() {
        Region region = sim.getRegionManager().getRegion("Fluid");
        Units mUnits = sim.getUnitsManager().getPreferredUnits(Dimensions.Builder().length(1).build());

        PointPart point = sim.getPartManager().createPointPart(new ArrayList<>(Collections.<NamedObject>emptyList()), new DoubleVector(new double[] {0.0, 0.0, 0.0}), null);
        point.setPresentationName("Freestream_Probe");
        point.getPointCoordinate().setCoordinate(mUnits, mUnits, mUnits, new DoubleVector(new double[] {-0.15, 0.15, 0.01})); // 3D center
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
		
        // --- 3D Energy & Entropy Quantities ---
        buildFunc("ek_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2) + pow($${PerturbationVelocity}[2], 2))");
        buildFunc("theta_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "0.0");
        buildFunc("phi", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "(${DynamicViscosity} + ${TurbulentViscosity}) * pow(${StrainRate}, 2)");
        
		// --- Enthalpy formulation (NEW! from Drew Sanders)
		buildFunc("enthalpy_storage_vol", 0, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "(${Density} * ${SpecificHeat} * (${Temperature} - ${FreestreamTemperature})) + (0.5 * ${Density} * mag2($${PerturbationVelocity}))");
		buildFunc("enth_flux_enthalpy", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${h0_defect_rho} * $${Velocity}, $${Normal})");
		buildFunc("enth_flux_pressure", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((-1.0) * (${AbsolutePressure} - ${FreestreamAbsolutePressure}) * $${Vinf}, $${Normal})");
		buildFunc("enth_flux_viscous", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((-1.0) * dotVector($$${Tshear}, $${PerturbationVelocity}), $${Normal})");
		
		buildFunc("aphi", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "${Tinf} / ${TotalTemperature} * ${phi}");
        buildFunc("ANablaT", 0, Dimensions.Builder().mass(1).length(-1).time(-3).build(), "${Tinf} * ${ThermalConductivity} * (1 / ${Temperature} * (pow(grad(${Temperature})[0], 2) + pow(grad(${Temperature})[1], 2) + pow(grad(${Temperature})[2], 2)))"); 
        
        // --- Tensors & Base Vectors (3D Extension) ---
        buildFunc("PerturbationVelocity", 1, Dimensions.Builder().length(1).time(-1).build(), "$${Velocity} - $${Vinf}");
        buildFunc("Tshear", 2, Dimensions.Builder().mass(1).length(-1).time(-2).build(), "[ (2.0 * ${EffectiveViscosity} * $${U_VelocityGrad}[0]) - (0.6667 * ${Density} * ${TurbulentKineticEnergy}) ; (2.0 * ${EffectiveViscosity} * 0.5 * ($${U_VelocityGrad}[1] + $${V_VelocityGrad}[0])), (2.0 * ${EffectiveViscosity} * $${V_VelocityGrad}[1]) - (0.6667 * ${Density} * ${TurbulentKineticEnergy}) ; (2.0 * ${EffectiveViscosity} * 0.5 * ($${U_VelocityGrad}[2] + $${W_VelocityGrad}[0])), (2.0 * ${EffectiveViscosity} * 0.5 * ($${V_VelocityGrad}[2] + $${W_VelocityGrad}[1])), (2.0 * ${EffectiveViscosity} * $${W_VelocityGrad}[2]) - (0.6667 * ${Density} * ${TurbulentKineticEnergy}) ]");
        
        // --- 3D Mechanical Energy Fluxes ---
        buildFunc("ep_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${StaticPressure} * $${PerturbationVelocity}, $${Normal})");
        buildFunc("ea_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(0.5 * ${Density} * pow($${PerturbationVelocity}[0], 2) * $${Velocity}, $${Normal})");
        buildFunc("ev_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(0.5 * ${Density} * pow($${PerturbationVelocity}[1], 2) * $${Velocity}, $${Normal})");
        buildFunc("ew_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(0.5 * ${Density} * pow($${PerturbationVelocity}[2], 2) * $${Velocity}, $${Normal})");
        buildFunc("xth_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(${MassSpec} * $${Velocity}, $${Normal})");
        buildFunc("etau_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot(dotVector($$${Tshear}, $${PerturbationVelocity}), $${Normal})");
        buildFunc("etp_flux", 0, Dimensions.Builder().mass(1).time(-3).build(), "dot((0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2) + pow($${PerturbationVelocity}[2], 2)) * $${Velocity}) + (${StaticPressure} * $${PerturbationVelocity}) - [($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1] + $$${Tshear}[0,2] * $${PerturbationVelocity}[2]), ($$${Tshear}[0,1] * $${PerturbationVelocity}[0] + $$${Tshear}[1,1] * $${PerturbationVelocity}[1] + $$${Tshear}[1,2] * $${PerturbationVelocity}[2]), ($$${Tshear}[0,2] * $${PerturbationVelocity}[0] + $$${Tshear}[1,2] * $${PerturbationVelocity}[1] + $$${Tshear}[2,2] * $${PerturbationVelocity}[2])], $${Normal})");
        
        // --- 3D X-Component Volumetric Extractions (For Field Means) ---
        buildFunc("ep_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "${StaticPressure} * $${PerturbationVelocity}[0]");
        buildFunc("ea_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "0.5 * ${Density} * pow($${PerturbationVelocity}[0], 2) * $${Velocity}[0]");
        buildFunc("ev_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "0.5 * ${Density} * pow($${PerturbationVelocity}[1], 2) * $${Velocity}[0]");
        buildFunc("ew_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "0.5 * ${Density} * pow($${PerturbationVelocity}[2], 2) * $${Velocity}[0]");
        buildFunc("xth_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "${MassSpec} * $${Velocity}[0]");
        buildFunc("etau_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1] + $$${Tshear}[0,2] * $${PerturbationVelocity}[2])");
        buildFunc("etp_vol", 0, Dimensions.Builder().mass(1).time(-3).build(), "(0.5 * ${Density} * (pow($${PerturbationVelocity}[0], 2) + pow($${PerturbationVelocity}[1], 2) + pow($${PerturbationVelocity}[2], 2)) * $${Velocity}[0]) + (${StaticPressure} * $${PerturbationVelocity}[0]) - ($$${Tshear}[0,0] * $${PerturbationVelocity}[0] + $$${Tshear}[0,1] * $${PerturbationVelocity}[1] + $$${Tshear}[0,2] * $${PerturbationVelocity}[2])");
        
        // --- Box Filters (Revised for Contiguous CVs) ---
        buildFunc("box_cv1", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.02 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv2", 0, Dimensions.Builder().build(), "($$Position[0] > 0.02 && $$Position[0] <= 0.05 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv3", 0, Dimensions.Builder().build(), "($$Position[0] > 0.05 && $$Position[0] <= 0.1 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv4", 0, Dimensions.Builder().build(), "($$Position[0] > 0.1 && $$Position[0] <= 0.15 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv5", 0, Dimensions.Builder().build(), "($$Position[0] > 0.15 && $$Position[0] <= 0.2 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv6", 0, Dimensions.Builder().build(), "($$Position[0] > 0.2 && $$Position[0] <= 0.3 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_cv7", 0, Dimensions.Builder().build(), "($$Position[0] > 0.3 && $$Position[0] <= 0.45 && $$Position[1] >= -0.1 && $$Position[1] <= 0.1 && $$Position[2] >= 0.0 && $$Position[2] <= 0.02) ? 1.0 : 0.0");

        buildFunc("box_x_cv1", 0, Dimensions.Builder().build(), "($$Position[0] >= -0.1 && $$Position[0] <= 0.02) ? 1.0 : 0.0");
        buildFunc("box_x_cv2", 0, Dimensions.Builder().build(), "($$Position[0] > 0.02 && $$Position[0] <= 0.05) ? 1.0 : 0.0");
        buildFunc("box_x_cv3", 0, Dimensions.Builder().build(), "($$Position[0] > 0.05 && $$Position[0] <= 0.1) ? 1.0 : 0.0");
        buildFunc("box_x_cv4", 0, Dimensions.Builder().build(), "($$Position[0] > 0.1 && $$Position[0] <= 0.15) ? 1.0 : 0.0");
        buildFunc("box_x_cv5", 0, Dimensions.Builder().build(), "($$Position[0] > 0.15 && $$Position[0] <= 0.2) ? 1.0 : 0.0");
        buildFunc("box_x_cv6", 0, Dimensions.Builder().build(), "($$Position[0] > 0.2 && $$Position[0] <= 0.3) ? 1.0 : 0.0");
        buildFunc("box_x_cv7", 0, Dimensions.Builder().build(), "($$Position[0] > 0.3 && $$Position[0] <= 0.45) ? 1.0 : 0.0");

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
        solver.getTimeStep().setDefinition("${timeStep_URANS}");
        solver.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);
        
        CoupledImplicitSolver cSolver = ((CoupledImplicitSolver) sim.getSolverManager().getSolver(CoupledImplicitSolver.class));
        cSolver.setLeaveTemporaryStorage(true);

        PhysicalTimeStoppingCriterion ptStop = ((PhysicalTimeStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));
        ptStop.getMaximumTime().setDefinition("${maxTime_DES}");
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
        
        double[] xBoundsEnd   = {0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45};

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

            // Build Top, Bottom, and Right planes for each CV
            buildPlaneSection(threshX, cs, mUnits, "bot_" + cvName, new double[]{0.0, -0.1, 0.0}, new double[]{0.0, -1.0, 0.0});
            buildPlaneSection(threshX, cs, mUnits, "top_" + cvName, new double[]{0.0, 0.1, 0.0}, new double[]{0.0, 1.0, 0.0});
            buildPlaneSection(threshY, cs, mUnits, "right_" + cvName, new double[]{xBoundsEnd[i], 0.0, 0.0}, new double[]{1.0, 0.0, 0.0});
        }
        
        // Single global left boundary plane based on CV1 threshold
        ThresholdPart threshY_cv1 = (ThresholdPart) sim.getPartManager().getObject("thresh_y_cv1");
        buildPlaneSection(threshY_cv1, cs, mUnits, "left", new double[]{-0.1, 0.0, 0.0}, new double[]{-1.0, 0.0, 0.0});
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
        
        // --- 1. 3D "VOLUME" INTEGRALS ---
        String[] volParts = {"cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7"};
        String[][] volMetrics = {
            {"massspec_vol", "MassSpec"},
            {"anabla_vol", "ANablaT"},
            {"ek_vol", "ek_vol"},
            {"theta_vol", "theta_vol"},
            {"phi_vol", "phi"},
            {"aphi_vol", "aphi"},
			{"enthalpy_storage_vol", "enthalpy_storage_vol"}
        };

        for (String[] metric : volMetrics) {
            createVolumeReport(metric[0], metric[1], volParts); 
        }

        // --- 2. 2D "SURFACE" INTEGRALS ---
        String[][] surfaceMetrics = {
            {"ep", "ep_flux"},
            {"ea", "ea_flux"},
            {"ev", "ev_flux"},
            {"ew", "ew_flux"},
            {"etp", "etp_flux"},
            {"xth", "xth_flux"},
            {"Etau", "etau_flux"},
			{"enth_flux_enthalpy", "enth_flux_enthalpy"},
			{"enth_flux_pressure", "enth_flux_pressure"},
			{"enth_flux_viscous", "enth_flux_viscous"}
        };

        String[] botParts = {"bot_cv1", "bot_cv2", "bot_cv3", "bot_cv4", "bot_cv5", "bot_cv6", "bot_cv7"};
        String[] topParts = {"top_cv1", "top_cv2", "top_cv3", "top_cv4", "top_cv5", "top_cv6", "top_cv7"};
        String[] rightParts = {"right_cv1", "right_cv2", "right_cv3", "right_cv4", "right_cv5", "right_cv6", "right_cv7"};
        String[] leftParts = {"left"}; // Reverted to evaluate the single left boundary

        for (String[] metric : surfaceMetrics) {
            String repPrefix = metric[0];
            String funcKey = metric[1];

            createSurfaceReport(repPrefix + "_bot", funcKey, botParts); 
            createSurfaceReport(repPrefix + "_top", funcKey, topParts);
            createSurfaceReport(repPrefix + "_right", funcKey, rightParts);
            createSurfaceReport(repPrefix + "_left", funcKey, leftParts); 
        }

        // --- 3. DRAG & LIFT FORCE REPORTS ---
        try {
            Boundary boundary = sim.getRegionManager().getRegion("Fluid").getBoundaryManager().getBoundary("wall_cylinder");
            
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
            "ea_vol", "ev_vol", "ep_vol", "ew_vol", "etau_vol", "xth_vol", "aphi", 
            "ANablaT", "theta_vol", "phi", "etp_vol", "MassSpec"
        };

        for (String var : meanVars) {
            try {
				ScalarGlobalParameter maxTimeParam_DES = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("maxTime_DES");
				ScalarGlobalParameter timeStepParam_DES = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("timeStep_DES");
				ScalarGlobalParameter maxTimeParam_URANS = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("maxTime_URANS");
				ScalarGlobalParameter timeStepParam_URANS = (ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("timeStep_URANS");

				double maxTimeVal_DES = maxTimeParam_DES.getQuantity().getValue();
				double timeStepVal_DES = timeStepParam_DES.getQuantity().getValue();
				double maxTimeVal_URANS = maxTimeParam_URANS.getQuantity().getValue();
				double timeStepVal_URANS = timeStepParam_URANS.getQuantity().getValue();

				int startTimeStep = (int) Math.round((maxTimeVal_URANS / timeStepVal_URANS) + (((maxTimeVal_DES - maxTimeVal_URANS) - 2.0) / timeStepVal_DES));
				int stopTimeStep = (int) Math.round((maxTimeVal_URANS / timeStepVal_URANS) + ((maxTimeVal_DES - maxTimeVal_URANS) / timeStepVal_DES));

				FieldMeanMonitor meanMon = (FieldMeanMonitor) monitorManager.create("star.base.report.FieldMeanMonitor");
				meanMon.setPresentationName("Mean_" + var);
				meanMon.getParts().setObjects(fluidRegion);
				meanMon.setFieldFunction((UserFieldFunction) sim.getFieldFunctionManager().getFunction(var));

				StarUpdate starUpdate = meanMon.getStarUpdate();
				TimeStepUpdateFrequency tsUpdateFreq = starUpdate.getTimeStepUpdateFrequency();

				IntegerValue startInt = tsUpdateFreq.getStartTimeStepQuantity();
				startInt.getQuantity().setValue(startTimeStep); 

				IntegerValue stopInt = tsUpdateFreq.getStopTimeStepQuantity();
				stopInt.getQuantity().setValue(stopTimeStep);
                
            } catch (Exception e) {
                sim.println("Warning: Could not create Mean monitor for " + var);
            }
        }
		
        int count = 0;
        for (Report report : sim.getReportManager().getObjects()) {
            try {
                ReportMonitor monitor = report.createMonitor();
                if (report instanceof VolumeIntegralReport || report instanceof SurfaceIntegralReport) {
                    try { monitor.getValueType().setSelected(ReportMonitorValueType.Type.PartValue); } catch (Exception ex) {}
                }
                monitor.setPlotLimit(200000);
                count++;
            } catch (Exception e) {}
        }
    }

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
            VolumeIntegralReport report = sim.getReportManager().createReport(VolumeIntegralReport.class);
            report.setPresentationName(reportName); 
            FieldFunction ff = getFieldFunction(funcKey);
            if (ff != null) report.setFieldFunction(ff);
            setReportParts(report, partNames);
        } catch (Exception e) {
            sim.println("Error generating 3D Volume report (" + reportName + "): " + e.getMessage());
        }
    }

    private void createSurfaceReport(String reportName, String funcKey, String[] partNames) {
        try {
            SurfaceIntegralReport report = sim.getReportManager().createReport(SurfaceIntegralReport.class);
            report.setPresentationName(reportName); 
            FieldFunction ff = getFieldFunction(funcKey);
            if (ff != null) report.setFieldFunction(ff);
            setReportParts(report, partNames);
        } catch (Exception e) {
            sim.println("Error generating 2D Surface report (" + reportName + "): " + e.getMessage());
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
            if (report instanceof VolumeIntegralReport) {
                ((VolumeIntegralReport) report).getParts().setObjects(partArray);
            } else if (report instanceof SurfaceIntegralReport) {
                ((SurfaceIntegralReport) report).getParts().setObjects(partArray);
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