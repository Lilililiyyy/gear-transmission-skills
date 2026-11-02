"""锥齿轮宏观几何、啮合力与常幅疲劳强度评估。"""
from .api import FatigueCurve, SteelMaterial, LoadCase, RatingSettings, analyze_pair, analyze_design, sweep_parameter, analyze_agma_pair, analyze_agma_design
from .geometry import GeometryInput, GeometryResult, calculate_geometry
from .forces import MeshForces, calculate_mesh_forces
from .io import read_config, export_result, write_json
from .assembly import build_assembly_scene, build_tooth_assembly_scene, export_assembly_view
from .manufacturing import ManufacturingInput, CutterOverrides, calculate_manufacturing, export_manufacturing_result
from .agma_rating import AGMARatingSettings, evaluate_agma_rating_mm, evaluate_agma_rating_si
from .rating_b2 import evaluate_b2_rating_mm
from .tooth_surface import build_spherical_tooth_surface, export_surface_json, export_surface_obj
from .kinematics import ideal_pair_pose
from .generated_cutter import CutterMachineSettings, evaluate_generated_point
from .generated_surface import GeneratedGearInput, build_generated_surface_grid, build_generated_accuracy_grid
from .generated_io import generated_grid_mesh, export_generated_grid
from .cad import GeneratedCadShell, build_generated_shell, export_step, export_brep, triangulate, measure_surface_accuracy
from .rigid_contact import RigidPose, RigidTransform, query_clearance
from .cad_import import ImportedCadModel, load_step_model, triangulate_imported_model, export_model_step, export_model_brep
from .clearance_phase import ClearancePhaseSettings, ClearancePhaseResult, solve_clearance_phase
from .gear_pair import StraightBevelPair, SpiralBevelPair
from .straight_cad import StraightCadSettings, StraightCadModel, build_straight_gear, build_straight_tooth_sector, build_straight_pair, sample_straight_contact_points
from .straight_tca import solve_straight_tca, straight_pair_frames
from .elasticity import ElasticMaterial, TetraMesh, ElasticSystem, ContactLoadChannel, contact_elasticity_factor, mesh_step_tetra, condense_contact_flexibility
from .contact_equilibrium import ContactEquilibriumResult, solve_contact_equilibrium, solve_normal_contact_equilibrium
from .straight_ltca import PreparedStraightLtca, prepare_straight_ltca, solve_straight_fe_ltca
from .iso2014_nominal import NominalMaterial, evaluate_iso2014_nominal
from .generated_solid import GeneratedBlankSettings, GeneratedBlank, GeneratedSolidModel, build_generated_blank, build_generated_gear
from .generated_tca import PreparedGeneratedTca, prepare_generated_tca, generated_tca_frames, solve_generated_tca
from .convergence import assess_straight_ltca_convergence
from .generated_ltca import PreparedGeneratedLtca, prepare_generated_ltca, solve_generated_fe_ltca
from .manufacturing_binding import BindingTolerances, measure_manufacturing_geometry, derive_nominal_geometry, bind_manufacturing_geometry
from .strength_assessment import DesignStressLimits, assess_nominal_strength
from .workflow import load_bevel_case, analyze_bevel_case, compare_workflow_reference
from .contact_cycle import validate_mesh_cycle_angles, sampled_pressure_maximum, build_mesh_cycle_pressure_data
from .contact_cycle_workflow import plot_mesh_cycle_workflows, audit_mesh_cycle_pressure_workflow

__version__ = "0.10.0"
__all__ = ["FatigueCurve", "SteelMaterial", "LoadCase", "RatingSettings", "analyze_pair", "analyze_design",
           "sweep_parameter", "GeometryInput", "GeometryResult", "calculate_geometry", "MeshForces",
           "calculate_mesh_forces", "read_config", "export_result", "write_json",
           "build_assembly_scene", "build_tooth_assembly_scene", "export_assembly_view",
           "ManufacturingInput","CutterOverrides",
           "calculate_manufacturing","export_manufacturing_result","AGMARatingSettings","analyze_agma_pair",
           "analyze_agma_design","evaluate_agma_rating_mm","evaluate_agma_rating_si","evaluate_b2_rating_mm",
           "build_spherical_tooth_surface",
           "export_surface_json","export_surface_obj","ideal_pair_pose",
           "CutterMachineSettings","evaluate_generated_point","GeneratedGearInput","build_generated_surface_grid",
           "build_generated_accuracy_grid","generated_grid_mesh","export_generated_grid","GeneratedCadShell",
           "build_generated_shell","export_step","export_brep","triangulate","measure_surface_accuracy",
           "RigidPose","RigidTransform","query_clearance","ImportedCadModel","load_step_model",
           "triangulate_imported_model","export_model_step","export_model_brep",
           "ClearancePhaseSettings","ClearancePhaseResult","solve_clearance_phase",
           "StraightBevelPair","SpiralBevelPair","StraightCadSettings","StraightCadModel",
           "build_straight_gear","build_straight_tooth_sector","build_straight_pair","sample_straight_contact_points",
           "solve_straight_tca","straight_pair_frames","ElasticMaterial","TetraMesh","ElasticSystem",
           "ContactLoadChannel","contact_elasticity_factor","mesh_step_tetra","condense_contact_flexibility",
           "ContactEquilibriumResult","solve_contact_equilibrium","solve_normal_contact_equilibrium",
           "PreparedStraightLtca","prepare_straight_ltca","solve_straight_fe_ltca",
           "NominalMaterial","evaluate_iso2014_nominal",
           "GeneratedBlankSettings","GeneratedBlank","GeneratedSolidModel","build_generated_blank","build_generated_gear",
           "PreparedGeneratedTca","prepare_generated_tca","generated_tca_frames","solve_generated_tca",
           "assess_straight_ltca_convergence",
           "PreparedGeneratedLtca","prepare_generated_ltca","solve_generated_fe_ltca",
           "BindingTolerances","measure_manufacturing_geometry","derive_nominal_geometry","bind_manufacturing_geometry",
           "DesignStressLimits","assess_nominal_strength","load_bevel_case","analyze_bevel_case","compare_workflow_reference",
           "validate_mesh_cycle_angles","sampled_pressure_maximum","build_mesh_cycle_pressure_data",
           "plot_mesh_cycle_workflows","audit_mesh_cycle_pressure_workflow"]
