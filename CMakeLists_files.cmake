# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.
#
# ATTIC_FILES           Unmaintained files. This for the projects developers
#                       only. Don't expect these files to build.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/core/flowdiagnostics/AnisotropicEikonal.cpp # test_anisotropicekonal.cpp, compute_eikonal_from_files
	opm/core/flowdiagnostics/DGBasis.cpp # test_dgbasis.cpp, compute_tof_from_files.cpp, compute_tof.cpp
	opm/core/flowdiagnostics/FlowDiagnostics.cpp # test_flowdiagnostics.cpp
	opm/core/flowdiagnostics/TofDiscGalReorder.cpp # compute_tof.cpp, compute_tof_from_files.cpp
	opm/core/flowdiagnostics/TofReorder.cpp # compute_tof.cpp, compute_tof_from_files.cpp
	opm/core/linalg/LinearSolverFactory.cpp # test_parallel_linearsolver.cpp, test_linearsolver.cpp, opm-simulators
        opm/core/linalg/LinearSolverInterface.cpp # ditto
        opm/core/linalg/LinearSolverIstl.cpp # ditto
        opm/core/linalg/LinearSolverPetsc.cpp # ditto
        opm/core/linalg/LinearSolverUmfpack.cpp # ditto
				opm/core/linalg/call_umfpack.c # ditto
				opm/core/linalg/sparse_sys.c # opm-simulators, opm-upscaling/sim_blackoil_impes.cpp, opm-upscaling/sim_co_impes.cpp
				opm/core/pressure/CompressibleTpfa.cpp # opm-simulators
				opm/core/pressure/FlowBCManager.cpp # opm-simulators
				opm/core/pressure/IncompTpfa.cpp # opm-simulators
				opm/core/pressure/IncompTpfaSinglePhase.cpp # unused
				opm/core/pressure/cfsh.c # unused
				opm/core/pressure/flow_bc.c # opm-simulators, opm-upscaling/sim_blackoil_impes.cpp, opm-upscaling/sim_co2_impes.cpp
				opm/core/pressure/fsh.c # unused
				opm/core/pressure/fsh_common_impl.c # unused
				opm/core/pressure/ifsh.c # unused
				opm/core/pressure/legacy_well.c # unused
				opm/core/pressure/mimetic/hybsys.c # unused
				opm/core/pressure/mimetic/hybsys_global.c # unused
				opm/core/pressure/mimetic/mimetic.c # opm-simulators
				opm/core/pressure/msmfem/coarse_conn.c # unused
				opm/core/pressure/msmfem/coarse_sys.c # unused
				opm/core/pressure/msmfem/dfs.c # test_equil.cpp
				opm/core/pressure/msmfem/hash_set.c # unused
				opm/core/pressure/msmfem/ifsh_ms.c # unused
				opm/core/pressure/msmfem/partition.c # test_equil.cpp
				opm/core/pressure/tpfa/TransTpfa.cpp # unused
				opm/core/pressure/tpfa/cfs_tpfa.c # opm-upscaling
				opm/core/pressure/tpfa/cfs_tpfa_residual.c # opm-simulators
				opm/core/pressure/tpfa/compr_bc.c # unused
				opm/core/pressure/tpfa/compr_quant.c # opm-upscaling
				opm/core/pressure/tpfa/compr_quant_general.c # unused
				opm/core/pressure/tpfa/compr_source.c # unused
				opm/core/pressure/tpfa/ifs_tpfa.c # opm-simulators
				opm/core/pressure/tpfa/trans_tpfa.c # test_pinchprocessor.cpp, opm-upscaling, opm-simulators
				opm/core/props/BlackoilPropertiesBasic.cpp # test_equil.cpp, opm-simulators
				opm/core/props/BlackoilPropertiesFromDeck.cpp # test_equil.cpp, test_satfunc.cpp, ewoms, opm-simulators
				opm/core/props/IncompPropertiesBasic.cpp # tutorial2-4.cpp, test_shadow.cpp, opm-simulators
				opm/core/props/IncompPropertiesFromDeck.cpp # opm-simulators
				opm/core/props/IncompPropertiesSinglePhase.cpp # compute_tof.cpp
				opm/core/props/pvt/PvtPropertiesBasic.cpp # test_equil.cpp, test_shadow.cpp, tutorial2-4.cpp, opm-simulators
				opm/core/props/pvt/PvtPropertiesIncompFromDeck.cpp # opm-simulators
				opm/core/props/rock/RockBasic.cpp # test_equil.cpp, test_shadow.cpp, tutorial2-4.cpp, opm-simulators
				opm/core/props/rock/RockCompressibility.cpp # compute_tof.cpp, tutorial2-4.cpp, opm-simulators, opm-upscaling
				opm/core/props/rock/RockFromDeck.cpp # test_equil.cpp, test_pinchprocessor.cpp, test_satfunc.cpp, test_wellsprocessor.cpp, test_stoppedwells.cpp ewoms, opm-simulators
				opm/core/props/satfunc/RelpermDiagnostics.cpp # diagnose_relperm.cpp, test_relpermdiagnostics.cpp, opm-simulators
				opm/core/props/satfunc/SaturationPropsBasic.cpp # test_equil.cpp, test_shadow.cpp, tutorial2-4.cpp, opm-simulators
				opm/core/props/satfunc/SaturationPropsFromDeck.cpp # test_equil.cpp, test_satfunc.pp, ewoms, opm-simulators
				opm/core/simulator/BlackoilState.cpp # test_equil.cpp, test_blackoilstate.cpp, test_stoppedwells.cpp, ewoms, opm-simulators
				opm/core/simulator/TwophaseState.cpp # tutorial2-4.cpp, opm-simulators
				opm/core/simulator/SimulatorReport.cpp # opm-simulators
				opm/core/transport/TransportSolverTwophaseInterface.cpp # tutorial3-4.cpp, opm-simulators
				opm/core/transport/implicit/TransportSolverTwophaseImplicit.cpp # opm-simulators
				opm/core/transport/implicit/transport_source.c  # opm-simulators
				opm/core/transport/minimal/spu_explicit.c # unused
				opm/core/transport/minimal/spu_implicit.c # unused
				opm/core/transport/reorder/ReorderSolverInterface.cpp # tutorial3-4.cpp, compute_tof.cpp, compute_tof_from_files.cpp, opm-simulators
				opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.cpp # opm-simulators
				opm/core/transport/reorder/TransportSolverTwophaseReorder.cpp # tutorial3-4.cpp, opm-simulators
				opm/core/transport/reorder/reordersequence.cpp # tutorial3-4.cpp, compute_tof.cpp, compute_tof_from_files.cpp, opm-simulators
				opm/core/transport/reorder/tarjan.c # opm-simulators
				opm/core/utility/Event.cpp # test_event.cpp, opm-simulators
				opm/core/utility/MonotCubicInterpolator.cpp # test_cubic.cpp, opm-simulators, opm-upscaling
				opm/core/utility/NullStream.cpp # opm-simulators
				opm/core/utility/VelocityInterpolation.cpp # test_velocityinterpolation.cpp
				opm/core/utility/WachspressCoord.cpp # test_velocityinterpolation.cpp, test_wachspresscoord.cpp
				opm/core/utility/compressedToCartesian.cpp # compute_initial_state.cpp, compute_tof.cpp, test_equil.cpp, test_relpermdiagnostics.cpp, test_wellsmanager.cpp, ewoms, opm-simulators
				 opm/core/utility/extractPvtTableIndex.cpp # compute_initial_state.cpp, test_equil.cpp, test_satfunc.cpp, ewoms, opm-simulators
				 opm/core/utility/miscUtilities.cpp # tutorial3-4.cpp, opm-simulators, opm-upscaling
				 opm/core/utility/miscUtilitiesBlackoil.cpp # opm-simulators
        opm/core/utility/parameters/Parameter.cpp # rubbel og bit
        opm/core/utility/parameters/ParameterGroup.cpp # ditto
        opm/core/utility/parameters/ParameterTools.cpp # ditto
				opm/core/wells/InjectionSpecification.cpp # test_stoppedwells.cpp, test_wellcollection.cpp, test_wellsgroup.cpp, test_wellsmanager.cpp, tutorial4.cpp, opm-simulators
				opm/core/wells/ProductionSpecification.cpp # compute_tof.cpp, test_stoppedwells.cpp, test_wellcollection.cpp, test_wellsgroup.cpp, tutorial4.cpp, opm-simulators
				opm/core/wells/WellCollection.cpp # compute_tof.cpp, tutorial4.cpp, test_wellsmanager.cpp, test_stoppedwells.cpp, test_wellcollection.cpp, opm-simulators
				opm/core/wells/WellsGroup.cpp # compute_tof.cpp, tutorial4.cpp, test_stoppedwells.cpp, test_wellcollection.cpp, test_wellsgroup.cpp, test_wellsmanager.cpp, opm-simulators
				opm/core/wells/WellsManager.cpp # compute_tof.cpp, test_stoppedwells.cpp, test_wellsmanager.cpp, opm-simulators
				opm/core/wells/well_controls.c # compute_tof.cpp, tutorial2-4.cpp, test_wells.cpp, test_wellcollection.cpp, test_wellgroups.cpp test_wellcontrols.cpp, opm-simulators, opm-upscaling
				opm/core/wells/wells.c # compute_tof.cpp, tutorial2-4.cpp, test_stoppedwells.cpp, test_wellcollection.cpp, test_wellsgroup.cpp, test_wellsmanager.cpp, opm-simulators
	)

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	tests/test_compressedpropertyaccess.cpp
		tests/test_dgbasis.cpp
		tests/test_cubic.cpp
	tests/test_event.cpp
		tests/test_flowdiagnostics.cpp
	tests/test_nonuniformtablelinear.cpp
	tests/test_parallelistlinformation.cpp
	tests/test_sparsevector.cpp
	tests/test_velocityinterpolation.cpp
	tests/test_uniformtablelinear.cpp
	tests/test_wells.cpp
	tests/test_wachspresscoord.cpp
	tests/test_linearsolver.cpp
		tests/test_parallel_linearsolver.cpp
	tests/test_param.cpp
	tests/test_satfunc.cpp
	tests/test_shadow.cpp
	 tests/test_equil.cpp
	tests/test_regionmapping.cpp
	tests/test_blackoilstate.cpp
	tests/test_wellsmanager.cpp
	tests/test_wellcontrols.cpp
	tests/test_wellsgroup.cpp
	tests/test_wellcollection.cpp
	tests/test_pinchprocessor.cpp
	tests/test_anisotropiceikonal.cpp
	tests/test_stoppedwells.cpp
	tests/test_relpermdiagnostics.cpp
        tests/test_norne_pvt.cpp
  )

# originally generated with the command:
# find tests -name '*.param' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	tests/testdata.param
	tests/liveoil.DATA
	tests/capillary.DATA
	tests/capillary_overlap.DATA
        tests/capillarySwatinit.DATA
	tests/compressed_gridproperty.data
	tests/deadfluids.DATA
	tests/equil_livegas.DATA
	tests/equil_liveoil.DATA
	tests/equil_rsvd_and_rvvd.DATA
	tests/wetgas.DATA
	tests/satfuncStandard.DATA
	tests/satfuncEPSBase.DATA
	tests/satfuncEPS_A.DATA
	tests/satfuncEPS_B.DATA
	tests/satfuncEPS_C.DATA
	tests/satfuncEPS_D.DATA
	tests/testBlackoilState1.DATA
	tests/testBlackoilState2.DATA
  tests/testPinch1.DATA
	tests/wells_manager_data.data
	tests/wells_manager_data_expanded.data
	tests/wells_manager_data_wellSTOP.data
        tests/wells_group.data
        tests/wells_stopped.data
				tests/relpermDiagnostics.DATA
        tests/norne_pvt.data
        )

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/compute_eikonal_from_files.cpp
	examples/compute_initial_state.cpp
	examples/compute_tof.cpp
	examples/compute_tof_from_files.cpp
	examples/diagnose_relperm.cpp
	tutorials/tutorial1.cpp
	tutorials/tutorial2.cpp
	tutorials/tutorial3.cpp
	tutorials/tutorial4.cpp
	)

# originally generated with the command:
# find attic -name '*.c*' -printf '\t%p\n' | sort
list (APPEND ATTIC_FILES
	attic/test_cfs_tpfa.c
	attic/test_jacsys.cpp
	attic/test_lapack.cpp
	attic/test_read_grid.c
	attic/test_read_vag.cpp
	attic/test_writeVtkData.cpp
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
        opm/core/doxygen_main.hpp
        opm/core/flowdiagnostics/AnisotropicEikonal.hpp
        opm/core/flowdiagnostics/DGBasis.hpp
        opm/core/flowdiagnostics/FlowDiagnostics.hpp
        opm/core/flowdiagnostics/TofDiscGalReorder.hpp
        opm/core/flowdiagnostics/TofReorder.hpp
        opm/core/linalg/LinearSolverFactory.hpp
        opm/core/linalg/LinearSolverInterface.hpp
        opm/core/linalg/LinearSolverIstl.hpp
        opm/core/linalg/LinearSolverPetsc.hpp
        opm/core/linalg/LinearSolverUmfpack.hpp
        opm/core/linalg/ParallelIstlInformation.hpp
        opm/core/linalg/blas_lapack.h
        opm/core/linalg/call_umfpack.h
        opm/core/linalg/sparse_sys.h
        opm/core/pressure/CompressibleTpfa.hpp
        opm/core/pressure/FlowBCManager.hpp
        opm/core/pressure/IncompTpfa.hpp
        opm/core/pressure/IncompTpfaSinglePhase.hpp
        opm/core/pressure/flow_bc.h
        opm/core/pressure/fsh.h
        opm/core/pressure/fsh_common_impl.h
        opm/core/pressure/legacy_well.h
        opm/core/pressure/mimetic/hybsys.h
        opm/core/pressure/mimetic/hybsys_global.h
        opm/core/pressure/mimetic/mimetic.h
        opm/core/pressure/msmfem/coarse_conn.h
        opm/core/pressure/msmfem/coarse_sys.h
        opm/core/pressure/msmfem/dfs.h
        opm/core/pressure/msmfem/hash_set.h
        opm/core/pressure/msmfem/ifsh_ms.h
        opm/core/pressure/msmfem/partition.h
        opm/core/pressure/tpfa/TransTpfa.hpp
        opm/core/pressure/tpfa/TransTpfa_impl.hpp
        opm/core/pressure/tpfa/cfs_tpfa.h
        opm/core/pressure/tpfa/cfs_tpfa_residual.h
        opm/core/pressure/tpfa/compr_bc.h
        opm/core/pressure/tpfa/compr_quant.h
        opm/core/pressure/tpfa/compr_quant_general.h
        opm/core/pressure/tpfa/compr_source.h
        opm/core/pressure/tpfa/ifs_tpfa.h
        opm/core/pressure/tpfa/trans_tpfa.h
        opm/core/props/BlackoilPhases.hpp
        opm/core/props/BlackoilPropertiesBasic.hpp
        opm/core/props/BlackoilPropertiesFromDeck.hpp
        opm/core/props/BlackoilPropertiesInterface.hpp
        opm/core/props/IncompPropertiesBasic.hpp
        opm/core/props/IncompPropertiesFromDeck.hpp
        opm/core/props/IncompPropertiesInterface.hpp
        opm/core/props/IncompPropertiesShadow.hpp
        opm/core/props/IncompPropertiesShadow_impl.hpp
        opm/core/props/IncompPropertiesSinglePhase.hpp
        opm/core/props/phaseUsageFromDeck.hpp
        opm/core/props/pvt/PvtPropertiesBasic.hpp
        opm/core/props/pvt/PvtPropertiesIncompFromDeck.hpp
        opm/core/props/pvt/ThermalGasPvtWrapper.hpp
        opm/core/props/pvt/ThermalOilPvtWrapper.hpp
        opm/core/props/pvt/ThermalWaterPvtWrapper.hpp
        opm/core/props/rock/RockBasic.hpp
        opm/core/props/rock/RockCompressibility.hpp
        opm/core/props/rock/RockFromDeck.hpp
        opm/core/props/satfunc/RelpermDiagnostics.hpp
        opm/core/props/satfunc/SaturationPropsBasic.hpp
        opm/core/props/satfunc/SaturationPropsFromDeck.hpp
        opm/core/props/satfunc/SaturationPropsInterface.hpp
	opm/core/props/satfunc/RelpermDiagnostics_impl.hpp
        opm/core/simulator/BlackoilState.hpp
        opm/core/simulator/BlackoilStateToFluidState.hpp
        opm/core/simulator/EquilibrationHelpers.hpp
        opm/core/simulator/ExplicitArraysFluidState.hpp
        opm/core/simulator/ExplicitArraysSatDerivativesFluidState.hpp
        opm/core/simulator/SimulatorReport.hpp
        opm/core/simulator/TwophaseState.hpp
        opm/core/simulator/WellState.hpp
        opm/core/simulator/initState.hpp
        opm/core/simulator/initStateEquil.hpp
        opm/core/simulator/initStateEquil_impl.hpp
        opm/core/simulator/initState_impl.hpp
        opm/core/transport/TransportSolverTwophaseInterface.hpp
        opm/core/transport/implicit/CSRMatrixBlockAssembler.hpp
        opm/core/transport/implicit/CSRMatrixUmfpackSolver.hpp
        opm/core/transport/implicit/ImplicitAssembly.hpp
        opm/core/transport/implicit/ImplicitTransport.hpp
        opm/core/transport/implicit/JacobianSystem.hpp
        opm/core/transport/implicit/NormSupport.hpp
        opm/core/transport/implicit/SimpleFluid2pWrappingProps.hpp
        opm/core/transport/implicit/SimpleFluid2pWrappingProps_impl.hpp
        opm/core/transport/implicit/SinglePointUpwindTwoPhase.hpp
        opm/core/transport/implicit/TransportSolverTwophaseImplicit.hpp
        opm/core/transport/implicit/transport_source.h
        opm/core/transport/minimal/spu_explicit.h
        opm/core/transport/minimal/spu_implicit.h
        opm/core/transport/reorder/ReorderSolverInterface.hpp
        opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp
        opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp
        opm/core/transport/reorder/reordersequence.h
        opm/core/transport/reorder/tarjan.h
        opm/core/utility/Average.hpp
        opm/core/utility/CompressedPropertyAccess.hpp
        opm/core/utility/DataMap.hpp
        opm/core/utility/Event.hpp
        opm/core/utility/Event_impl.hpp
        opm/core/utility/Factory.hpp
        opm/core/utility/initHydroCarbonState.hpp
        opm/core/utility/MonotCubicInterpolator.hpp
        opm/core/utility/NonuniformTableLinear.hpp
        opm/core/utility/NullStream.hpp
        opm/core/utility/RegionMapping.hpp
        opm/core/utility/RootFinders.hpp
        opm/core/utility/SparseVector.hpp
        opm/core/utility/UniformTableLinear.hpp
        opm/core/utility/VelocityInterpolation.hpp
        opm/core/utility/WachspressCoord.hpp
        opm/core/utility/buildUniformMonotoneTable.hpp
        opm/core/utility/compressedToCartesian.hpp
        opm/core/utility/extractPvtTableIndex.hpp
        opm/core/utility/have_boost_redef.hpp
        opm/core/utility/linearInterpolation.hpp
        opm/core/utility/miscUtilities.hpp
        opm/core/utility/miscUtilitiesBlackoil.hpp
        opm/core/utility/miscUtilities_impl.hpp
        opm/core/utility/parameters/Parameter.hpp
        opm/core/utility/parameters/ParameterGroup.hpp
        opm/core/utility/parameters/ParameterGroup_impl.hpp
        opm/core/utility/parameters/ParameterMapItem.hpp
        opm/core/utility/parameters/ParameterRequirement.hpp
        opm/core/utility/parameters/ParameterStrings.hpp
        opm/core/utility/parameters/ParameterTools.hpp
        opm/core/utility/share_obj.hpp
        opm/core/well_controls.h
        opm/core/wells.h
        opm/core/wells/InjectionSpecification.hpp
        opm/core/wells/ProductionSpecification.hpp
        opm/core/wells/WellCollection.hpp
        opm/core/wells/WellsGroup.hpp
        opm/core/wells/WellsManager.hpp
        opm/core/wells/DynamicListEconLimited.hpp
        opm/core/wells/WellsManager_impl.hpp
	)
