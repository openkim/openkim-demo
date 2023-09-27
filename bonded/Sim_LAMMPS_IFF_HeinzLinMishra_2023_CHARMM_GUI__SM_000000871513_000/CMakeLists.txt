cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API-ITEMS 2.2 REQUIRED CONFIG)

kim_api_items_setup_before_project(ITEM_TYPE "simulatorModel")
project(Sim_LAMMPS_IFF_HeinzLinMishra_2023_CHARMM_GUI__SM_000000871513_000)
kim_api_items_setup_after_project(ITEM_TYPE "simulatorModel")

add_kim_api_model_library(
  NAME            ${PROJECT_NAME}
  SM_SPEC_FILE    "smspec.edn"
  PARAMETER_FILES "iff.params"
                  "iff.mass"
  )