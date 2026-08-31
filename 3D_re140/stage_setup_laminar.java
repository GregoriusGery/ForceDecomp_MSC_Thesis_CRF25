// Simcenter STAR-CCM+ macro: test.java
// Written by Simcenter STAR-CCM+ 20.06.007
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.common.objectview.*;
import star.turbulence.*;
import star.kwturb.*;
import star.walldistance.*;
import star.automation.*;

public class stage_setup_laminar extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = getActiveSimulation();
	
    Stage stage_0 = simulation_0.get(StageManager.class).createStage("Stage");
    stage_0.setPresentationName("URANS");
    Stage stage_1 = simulation_0.get(StageManager.class).createStage("Stage");
    StagedObjectTreeView stagedObjectTreeView_0 = simulation_0.get(StageManager.class).getStageTreeViewManager().getTreeView().get();
    stagedObjectTreeView_0.setViewMode(ObjectTreeViewMode.Flat);
    stage_1.setPresentationName("DES");

    simulation_0.get(StageManager.class).setActiveStage(stage_0);

    PhysicsContinuum physicsContinuum_0 = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));

    simulation_0.get(StageManager.class).stage(true, new ArrayList<>(Arrays.<ClientServerObject>asList(physicsContinuum_0.getModelManager())));

    simulation_0.get(StageManager.class).setActiveStage(stage_1);

    KwAllYplusWallTreatment kwAllYplusWallTreatment_0 = 
      physicsContinuum_0.getModelManager().getModel(KwAllYplusWallTreatment.class);

    physicsContinuum_0.disableModel(kwAllYplusWallTreatment_0);

    WallDistanceModel wallDistanceModel_0 = physicsContinuum_0.getModelManager().getModel(WallDistanceModel.class);

    physicsContinuum_0.disableModel(wallDistanceModel_0);

    SstKwTurbModel sstKwTurbModel_0 = physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    physicsContinuum_0.disableModel(sstKwTurbModel_0);

    ImplicitUnsteadyModel implicitUnsteadyModel_0 = physicsContinuum_0.getModelManager().getModel(ImplicitUnsteadyModel.class);

    physicsContinuum_0.disableModel(implicitUnsteadyModel_0);

    KOmegaTurbulence kOmegaTurbulence_0 = 
      physicsContinuum_0.getModelManager().getModel(KOmegaTurbulence.class);

    physicsContinuum_0.disableModel(kOmegaTurbulence_0);

    RansTurbulenceModel ransTurbulenceModel_0 = 
      physicsContinuum_0.getModelManager().getModel(RansTurbulenceModel.class);

    physicsContinuum_0.disableModel(ransTurbulenceModel_0);
    physicsContinuum_0.enable(ImplicitUnsteadyModel.class);
    physicsContinuum_0.enable(DesTurbulenceModel.class);
    physicsContinuum_0.enable(SstKwTurbDesModel.class);
    physicsContinuum_0.enable(KwAllYplusWallTreatment.class);

    ImplicitUnsteadySolver implicitUnsteadySolver_0 = 
      ((ImplicitUnsteadySolver) simulation_0.getSolverManager().getSolver(ImplicitUnsteadySolver.class));

    simulation_0.get(StageManager.class).stage(true, new ArrayList<>(Arrays.<ClientServerObject>asList(implicitUnsteadySolver_0)));

    simulation_0.get(StageManager.class).setActiveStage(stage_0);

    ImplicitUnsteadySolver implicitUnsteadySolver_1 = 
      ((ImplicitUnsteadySolver) simulation_0.getSolverManager().getSolver(ImplicitUnsteadySolver.class));

    PhysicalTimeStoppingCriterion physicalTimeStoppingCriterion_0 = 
      implicitUnsteadySolver_1.getSolverStoppingCriterionManager().create("star.common.PhysicalTimeStoppingCriterion");

    Units units_1 = 
      simulation_0.getUnitsManager().getPreferredUnits(Dimensions.Builder().time(1).build());

    physicalTimeStoppingCriterion_0.getMaximumTime().setDefinition("${maxTime_URANS}");

    InnerIterationStoppingCriterion innerIterationStoppingCriterion_0 = 
      implicitUnsteadySolver_1.getSolverStoppingCriterionManager().create("star.common.InnerIterationStoppingCriterion");

    simulation_0.get(StageManager.class).setActiveStage(stage_1);

    ImplicitUnsteadySolver implicitUnsteadySolver_2 = 
      ((ImplicitUnsteadySolver) simulation_0.getSolverManager().getSolver(ImplicitUnsteadySolver.class));

    implicitUnsteadySolver_2.getTimeStep().setDefinition("${timeStep_DES}");

    PhysicalTimeStoppingCriterion physicalTimeStoppingCriterion_1 = 
      implicitUnsteadySolver_2.getSolverStoppingCriterionManager().create("star.common.PhysicalTimeStoppingCriterion");

    physicalTimeStoppingCriterion_1.getMaximumTime().setDefinition("${maxTime_DES}");

    MinimumInnerIterationStoppingCriterion minimumInnerIterationStoppingCriterion_0 = 
      implicitUnsteadySolver_2.getSolverStoppingCriterionManager().create("star.common.MinimumInnerIterationStoppingCriterion");

    IntegerValue integerValue_0 = minimumInnerIterationStoppingCriterion_0.getMinIterations();

    integerValue_0.getQuantity().setValue(5.0);

    simulation_0.get(StageManager.class).setActiveStage(stage_0);

    SimDriverWorkflow simDriverWorkflow_0 = simulation_0.get(SimDriverWorkflowManager.class).createSimDriverWorkflow("Simulation Operations");

    ClearSolutionAutomationBlock clearSolutionAutomationBlock_0 = 
      (ClearSolutionAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.ClearSolutionAutomationBlock", "Clear Solution");

    InitializeSolutionAutomationBlock initializeSolutionAutomationBlock_0 = 
      (InitializeSolutionAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.InitializeSolutionAutomationBlock", "Initialize Solution");

    SetActiveStageAutomationBlock setActiveStageAutomationBlock_0 = 
      (SetActiveStageAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.SetActiveStageAutomationBlock", "Set Active Stage");

    setActiveStageAutomationBlock_0.setStage(stage_0);

    SolvePhysics solvePhysics_0 = (SolvePhysics) simDriverWorkflow_0.getBlocks().createBlock("star.common.SolvePhysics", "Solve Physics");
    solvePhysics_0.getSimulationObjects().setQuery(null);
    solvePhysics_0.getSimulationObjects().setObjects(physicsContinuum_0);

    SaveSimulationAutomationBlock saveSimulationAutomationBlock_0 = 
      (SaveSimulationAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.common.SaveSimulationAutomationBlock", "Save Simulation");

    saveSimulationAutomationBlock_0.setBaseFilename("data");
    saveSimulationAutomationBlock_0.setPostfixMode(FilenamePostfixMode.TIMESTEP);
    SetActiveStageAutomationBlock setActiveStageAutomationBlock_1 = (SetActiveStageAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.automation.SetActiveStageAutomationBlock", "Set Active Stage");

    setActiveStageAutomationBlock_1.setStage(stage_1);

    SolvePhysics solvePhysics_1 = (SolvePhysics) simDriverWorkflow_0.getBlocks().createBlock("star.common.SolvePhysics", "Solve Physics");
    solvePhysics_1.getSimulationObjects().setQuery(null);
    solvePhysics_1.getSimulationObjects().setObjects(physicsContinuum_0);

    SaveSimulationAutomationBlock saveSimulationAutomationBlock_1 = (SaveSimulationAutomationBlock) simDriverWorkflow_0.getBlocks().createBlock("star.common.SaveSimulationAutomationBlock", "Save Simulation");
    saveSimulationAutomationBlock_1.setBaseFilename("data");
    saveSimulationAutomationBlock_1.setPostfixMode(FilenamePostfixMode.TIMESTEP);
  }
}
