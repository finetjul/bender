/*=========================================================================

  Program: Bender

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

// Bender includes
#include "DecimateModelCLP.h"
#include "benderIOUtils.h"

// VTK includes
#include <vtkDecimatePro.h>
#include <vtkFillHolesFilter.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>

int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  vtkSmartPointer<vtkPolyData> model;
  model.TakeReference(bender::IOUtils::ReadPolyData(InputModel.c_str()));

  if (model->GetNumberOfPolys() < 1 || Triangulate)
    {
    if (model->GetNumberOfPolys() < 1)
      {
      std::cout << "No triangle in model" << std::endl;
      }
    std::cout << "Triangulate..." << std::endl;
    vtkNew<vtkTriangleFilter> triangulate;
    triangulate->SetInput(model);
    triangulate->Update();
    model = triangulate->GetOutput();
    }

  if (Debug)
    {
    std::cout << "Input model:\n"
              << "  Points: " << model->GetNumberOfPoints() << "\n"
              << "  Cells: " << model->GetNumberOfCells() << "\n"
              << std::endl;
    }

  std::cout << "Decimate..." << std::endl;
  vtkNew<vtkDecimatePro> decimator;
  decimator->SetInput(model);
  decimator->SetPreSplitMesh(PreSplitMesh);
  decimator->SetTargetReduction(TargetReduction);
  decimator->SetPreserveTopology(PreserveTopology);
  decimator->SetFeatureAngle(FeatureAngle);
  decimator->SetSplitting( Splitting || PreSplitMesh );
  decimator->SetSplitAngle(SplitAngle);
  if (AbsoluteError)
    {
    decimator->SetAbsoluteError(MaximumError == 100000000000000 ? VTK_DOUBLE_MAX : MaximumError);
    }
  else
    {
    decimator->SetMaximumError(
      MaximumError == 100000000000000 ? VTK_DOUBLE_MAX : MaximumError);
    }
  decimator->SetAccumulateError(AccumulateError);
  decimator->SetBoundaryVertexDeletion(BoundaryVertexDeletion);
  decimator->SetDegree(Degree);

  if (Debug)
    {
    decimator->Print(std::cout);
    decimator->DebugOn();
    }
  decimator->Update();
  vtkSmartPointer<vtkPolyData> decimatedModel = decimator->GetOutput();

  if (Debug)
    {
    std::cout << "Decimated model:\n"
              << "  Points: " << decimatedModel->GetNumberOfPoints() << "\n"
              << "  Cells: " << decimatedModel->GetNumberOfCells() << "\n"
              << std::endl;
    }

  vtkSmartPointer<vtkPolyData> finalModel = decimatedModel;
  if (FillHoles)
    {
    std::cout << "Fill Holes..." << std::endl;
    vtkNew<vtkFillHolesFilter> fillHoles;
    fillHoles->SetInput(decimatedModel);
    fillHoles->Update();
    finalModel = fillHoles->GetOutput();
    }

  bender::IOUtils::WritePolyData(finalModel, DecimatedModel.c_str());

  return EXIT_SUCCESS;
}
