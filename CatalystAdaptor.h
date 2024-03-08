// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.hpp>

#include "vtkCellData.h"
#include <vtkCompositeDataIterator.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataReader.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include <vtkPoints.h>
#include <vtkSetGet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <iostream>
#include <math.h>
#include <mpi.h>
#include <string>
#include <vector>

//#include <ctime>
#include <chrono>

#define vtkTemplateMacroConduit(call)                                                                     \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call);                                                  \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);                                                    \
  vtkTemplateMacroCase(VTK_INT, int, call);                          \
  vtkTemplateMacroCase(VTK_UNSIGNED_INT, unsigned int, call);        \
  vtkTemplateMacroCase(VTK_SHORT, short, call);                                                    \
  vtkTemplateMacroCase(VTK_UNSIGNED_SHORT, unsigned short, call);                                  \
  vtkTemplateMacroCase(VTK_SIGNED_CHAR, signed char, call);                                        \
  vtkTemplateMacroCase(VTK_LONG, long, call);                           \
  vtkTemplateMacroCase(VTK_UNSIGNED_LONG, unsigned long, call);       \
  vtkTemplateMacroCase(VTK_UNSIGNED_CHAR, unsigned char, call)

// These are the data types that Conduit doesn't handle with zero-copy:
//   vtkTemplateMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, call);                          \
//  vtkTemplateMacroCase(VTK_LONG_LONG, long long, call);       \
//  vtkTemplateMacroCase(VTK_CHAR, char, call); \
// vtkTemplateMacroCase(VTK_ID_TYPE, vtkIdType, call);  \

namespace CatalystAdaptor
{
  vtkSmartPointer<vtkXMLMultiBlockDataReader> Reader;

  // used to print out what types of grids are being processed just once
  bool AMR_Processed = false;
  bool UnstructuredGrid_Processed = false;

// struct AMR
// {
//   int NumberOfAMRLevels;
//   std::array<int,6> GetLevelIndices(int i)
//     {
//       std::array<int,6> a = {i, i, i, i, i, i};
//       return a;
//     }
//   std::array<double,3> GetLevelOrigin(int i)
//     {
//       std::array<double,3> a = {double(i), double(i), double(i)};
//       return a;
//     }
//   std::vector<int> BlockId;

// };

void ProcessFieldData(vtkFieldData* field_data, std::string association, conduit_cpp::Node& fields)
{
  for (int i=0;i<field_data->GetNumberOfArrays();i++)
  {
    vtkDataArray* array = field_data->GetArray(i);
    std::string name = array->GetName();
    conduit_cpp::Node conduit_field = fields[name];
    conduit_field["association"] = association;
    conduit_field["topology"] = "mesh"; // ACB may be topo instead of mesh
    //conduit_field["topology"] = "topo"; // ACB may be topo instead of mesh
    conduit_field["volume_dependent"].set("false");
    if (array->GetNumberOfComponents() == 1)
    {
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples()));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,0);
        }
        conduit_field["values"] = data;
      }
    }
    else if (array->GetNumberOfComponents() == 3)
    {
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/x"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/0, /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,0);
        }
        conduit_field["values/x"] = data;
      }
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/y"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/sizeof(VTK_TT), /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,1);
        }
        conduit_field["values/y"] = data;
      }
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/z"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/2 * sizeof(VTK_TT), /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,2);
        }
        conduit_field["values/z"] = data;
      }
    }
    else
    {
      std::cerr << "ERROR: can't add " << name << " array because of having " << array->GetNumberOfComponents() << " number of components\n";
    }
  }
}
/**
 * In this example, we show how we can use Catalysts's C++
 * wrapper around conduit's C API to create Conduit nodes.
 * This is not required. A C++ adaptor can just as
 * conveniently use the Conduit C API to setup the
 * `conduit_node`. However, this example shows that one can
 * indeed use Catalyst's C++ API, if the developer so chooses.
 */
void Initialize(int argc, char* argv[], std::string& data_file)
{
  conduit_cpp::Node node;
  for (int cc = 1; cc < argc; ++cc)
  {
    if (strcmp(argv[cc], "-i") == 0)
    { // skip data file arguments
      cc += 1;
    }
    else
    {
      node["catalyst/scripts/script" + std::to_string(cc - 1)].set_string(argv[cc]);
      std::cerr << "Adding ParaView Catalyst script " << argv[cc] << std::endl;
    }
  }
  node["catalyst_load/implementation"] = "paraview";
  node["catalyst_load/search_paths/paraview"] = PARAVIEW_IMPL_DIR;
  std::cerr << "default search location for ParaView Catalyst is " << PARAVIEW_IMPL_DIR << std::endl;
  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
  }
  Reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
  Reader->SetFileName(data_file.c_str());
  std::cerr << "reading " << data_file << std::endl;
  Reader->Update();
}

void Execute(unsigned int cycle, double time)
{
  vtkMultiBlockDataSet* data_set = vtkMultiBlockDataSet::SafeDownCast(Reader->GetOutput());
  // vtkCompositeDataIterator* iter = data_set->NewIterator();
  // iter->SkipEmptyNodesOn();
  // iter->InitTraversal();
  // for(iter->GoToFirstItem();iter->IsDoneWithTraversal()==false;iter->GoToNextItem())
  // {
  //   vtkDataObject* data_object = iter->GetCurrentDataObject();
  //   unsigned int flat_index = iter->GetCurrentFlatIndex();
  //   std::cerr << "AT " << flat_index << " which is a " << data_object->GetClassName() << std::endl;
  // }

  //auto start = std::time(nullptr);
  auto start = std::chrono::high_resolution_clock::now();

  int numRanks(1), myRank(0);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  conduit_cpp::Node exec_params;

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(cycle);
  state["time"].set(time);
  state["multiblock"].set(1);

  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Since this example is using Conduit Mesh Blueprint to define the mesh,
  // we set the channel's type to "mesh".
  channel["type"].set("multimesh");

  // now create the mesh.
  conduit_cpp::Node mesh = channel["data"];

  vtkCompositeDataIterator* iter = data_set->NewIterator();
  iter->SkipEmptyNodesOn();
  iter->InitTraversal();
  for(iter->GoToFirstItem();iter->IsDoneWithTraversal()==false;iter->GoToNextItem())
  {
    unsigned int flat_index = iter->GetCurrentFlatIndex();
    std::string patch_name = "domain_" + std::to_string(flat_index);
    conduit_cpp::Node patch = mesh[patch_name];
    // add basic state info
    patch["state/domain_id"] = flat_index;
    patch["state/cycle"] = cycle;
    patch["state/time"] = time;
    //patch["state/level"] = level;

    vtkDataSet* data_set = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    if (vtkImageData* image_data = vtkImageData::SafeDownCast(iter->GetCurrentDataObject()))
    {
      if (!AMR_Processed)
      {
        std::cerr << "Processing AMR\n";
        AMR_Processed = true;
      }
      patch["coordsets/coords/type"] = "uniform";
      int extent[6];
      image_data->GetExtent(extent);
      patch["coordsets/coords/dims/i"] = extent[1] - extent[0] + 1;
      patch["coordsets/coords/dims/j"] = extent[3] - extent[2] + 1;
      patch["coordsets/coords/dims/k"] = extent[5] - extent[4] + 1;

      double spacing[3];
      image_data->GetSpacing(spacing);
      patch["coordsets/coords/spacing/dx"] = spacing[0];
      patch["coordsets/coords/spacing/dy"] = spacing[1];
      patch["coordsets/coords/spacing/dz"] = spacing[2];

      double origin[3];
      image_data->GetOrigin(origin);
      patch["coordsets/coords/origin/x"] = origin[0];
      patch["coordsets/coords/origin/y"] = origin[1];
      patch["coordsets/coords/origin/z"] = origin[2];

      // create a rectilinear topology that refs our coordset
      patch["topologies/mesh/type"] = "uniform";
      patch["topologies/mesh/coordset"] = "coords";

      // add logical elements origin
      patch["topologies/mesh/elements/origin/i0"] = extent[0];
      patch["topologies/mesh/elements/origin/j0"] = extent[2];
      patch["topologies/mesh/elements/origin/k0"] = extent[4];
    }
    else if (vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject()))
    {
      if (!UnstructuredGrid_Processed)
      {
        std::cerr << "Processing UnstructuredGrid\n";
        UnstructuredGrid_Processed = true;
      }
      patch["coordsets/coords/type"] = "explicit";

      conduit_cpp::Node conduit_field = patch["coordsets/coords"];
      vtkDataArray* array = ugrid->GetPoints()->GetData();
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/x"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/0, /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,0);
        }
        conduit_field["values/x"] = data;
      }
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/y"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/sizeof(VTK_TT), /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,1);
        }
        conduit_field["values/y"] = data;
      }
      switch(array->GetDataType())
      {
        vtkTemplateMacroConduit(conduit_field["values/z"].set_external(static_cast<VTK_TT*>(array->GetVoidPointer(0)), array->GetNumberOfTuples(), /*offset=*/2 * sizeof(VTK_TT), /*stride=*/3 * sizeof(VTK_TT)));
      default:
        std::cerr << "Doing deep-copy for array type " << array->GetDataType() << " for array " << array->GetName() << std::endl;
        std::vector<double> data(array->GetNumberOfTuples());
        for (vtkIdType i=0;i<array->GetNumberOfTuples();i++)
        {
          data[i] = array->GetComponent(i,2);
        }
        conduit_field["values/z"] = data;
      }
      // now put in cell connectivity
      patch["topologies/mesh/type"] = "unstructured";
      patch["topologies/mesh/coordset"] = "coords";

      patch["topologies/mesh/elements/shape"] = "mixed";
      patch["topologies/mesh/elements/shape_map/hex"] = VTK_HEXAHEDRON;
      patch["topologies/mesh/elements/shape_map/tet"] = VTK_TETRA;
      patch["topologies/mesh/elements/shape_map/wedge"] = VTK_WEDGE;
      patch["topologies/mesh/elements/shape_map/pyramid"] = VTK_PYRAMID;

      std::vector<int> shapes;
      std::vector<int> sizes;
      std::vector<int> offsets(1,0);
      std::vector<int> connectivity;
      vtkNew<vtkIdList> cell_ids;
      for (vtkIdType c=0;c<ugrid->GetNumberOfCells();c++)
      {
        ugrid->GetCellPoints(c, cell_ids);
        shapes.push_back(ugrid->GetCellType(c));
        sizes.push_back(cell_ids->GetNumberOfIds());
        offsets.push_back(offsets.size()-1+cell_ids->GetNumberOfIds());
        for (vtkIdType i=0;i<cell_ids->GetNumberOfIds();i++)
        {
          connectivity.push_back(cell_ids->GetId(i));
        }
      }
      offsets.pop_back(); // remove last value since it's not needed and will cause issues during verification
      patch["topologies/mesh/elements/shapes"] = shapes;
      patch["topologies/mesh/elements/sizes"] = sizes;
      patch["topologies/mesh/elements/offsets"] = offsets;
      patch["topologies/mesh/elements/connectivity"] = connectivity;

    } // if imagedata or unstructured grid

    // if (level > 0)
    // {
    //   conduit_cpp::Node nest_set;
    //   nest_set["association"] = "element";
    //   nest_set["topology"] = "topo";
    //   int parent_id = amr.BlockId[level-1];
    //   std::string parent_name = "windows/window_" + std::to_string(parent_id);
    //   conduit_cpp::Node parent = nest_set[parent_name];
    //   parent["domain_id"] = parent_id;
    //   parent["domain_type"] = "parent";
    //   std::array<int,6> parentLevelIndices = amr.GetLevelIndices(level-1);
    //   parent["origin/i"] = levelIndices[0]/2; //parentLevelIndices[0]+std::pow(2,level-1);
    //   parent["origin/j"] = parentLevelIndices[2];
    //   parent["origin/k"] = parentLevelIndices[4];
    //   parent["dims/i"] = parentLevelIndices[1] - levelIndices[0]/2 + 1;
    //   parent["dims/j"] = parentLevelIndices[3] - parentLevelIndices[2] + 1;;
    //   parent["dims/k"] = parentLevelIndices[5] - parentLevelIndices[4] + 1;;
    //   parent["ratio/i"] = 2;
    //   parent["ratio/j"] = 2;
    //   parent["ratio/k"] = 2;

    //   int child_id = amr.BlockId[level];
    //   std::string child_name = "windows/window_" + std::to_string(child_id);
    //   conduit_cpp::Node child = nest_set[child_name];
    //   child["domain_id"] = child_id;
    //   child["domain_type"] = "child";

    //   child["origin/i"] = levelIndices[0];
    //   child["origin/j"] = levelIndices[2];
    //   child["origin/k"] = levelIndices[4];

    //   child["dims/i"] = levelIndices[1] - levelIndices[0] + 1;
    //   child["dims/j"] = levelIndices[3] - levelIndices[2] + 1;
    //   child["dims/k"] = levelIndices[5] - levelIndices[4] + 1;

    //   child["ratio/i"] = 2;
    //   child["ratio/j"] = 2;
    //   child["ratio/k"] = 2;

    //   patch["nestsets/nest"].set(nest_set);
    // }

    // add fields
    conduit_cpp::Node fields = patch["fields"];
    ProcessFieldData(data_set->GetPointData(), "vertex", fields);
    ProcessFieldData(data_set->GetCellData(), "element", fields);
  }
  //auto end_construct = std::time(nullptr);
  auto end_construct = std::chrono::high_resolution_clock::now();

  //exec_params.print(); // for viewing the Conduit node information

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok)
  {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }
  iter->Delete();
  //auto end_processing = std::time(nullptr);
  auto end_processing = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> createdata = end_construct - start;
  std::chrono::duration<double> catalystuse = end_processing - end_construct;

  std::cerr << "For cycle " << cycle << " the time to generate the Conduit node is " << createdata.count() << " and the time to Catalyst::execute is " << catalystuse.count() << std::endl;
}

void Finalize()
{
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "Failed to finalize Catalyst: " << err << std::endl;
  }
  Reader = nullptr;
}
}

#endif
