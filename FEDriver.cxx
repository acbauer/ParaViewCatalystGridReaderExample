// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#include "CatalystAdaptor.h"
#include <iostream>
#include <mpi.h>

// Example of a C++ adaptor for a simulation code
// where the simulation code reads in data from a file.

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int numRanks(1), myRank(0);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  //std::string dataFile = "/home/acbauer/Data/Helios/boundary_contour_volumetric/volumetricgrid_000050.vtm"; // has bad blanking information
  //std::string dataFile = "/home/acbauer/Data/Helios/SimpleGeometryVolumetricData/ugrid.vtm";
  //std::string dataFile = "/home/acbauer/Data/Helios/SimpleGeometryVolumetricData/amr.vtm";
  std::string dataFile = "/home/acbauer/Data/Helios/SimpleGeometryVolumetricData/volumetricgrid_000200.vtm";
  for (int i=0;i<argc-1;i++)
  {
    std::string arg = argv[i];
    if (arg == "-i")
    {
      dataFile = argv[i+1];
      break;
    }
  }
  std::cerr << "Input coming from: " << dataFile << std::endl;

  // The first argument is the program name
  CatalystAdaptor::Initialize(argc, argv, dataFile);
  unsigned int numberOfTimeSteps = 10;
  for (unsigned int timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
  {
    // use a time step length of 0.1
    double time = timeStep * 0.1;
    CatalystAdaptor::Execute(timeStep, time);
  }
  CatalystAdaptor::Finalize();

  MPI_Finalize();
  return 0;
}
