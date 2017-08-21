/*=========================================================================

  Program:   ParaView
  Module:    vtkCleanUnstructuredGrid.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class   vtkCleanUnstructuredGrid
 * @brief   merge duplicate points
 *
 *
 * vtkCleanUnstructuredGrid is a filter that takes unstructured grid data as
 * input and generates unstructured grid data as output. vtkCleanUnstructuredGrid can
 * merge duplicate points (with coincident coordinates) using the vtkMergePoints object
 * to merge points.
 *
 * @sa
 * vtkCleanPolyData
*/

#ifndef vtkCleanUnstructuredGrid_h
#define vtkCleanUnstructuredGrid_h

#include "vtkUnstructuredGridAlgorithm.h"

#if (__cplusplus >= 201103L) || ( defined(_MSC_VER) && _MSC_VER >= 1800 )
# define VTK_DELETE_FUNCTION =delete
#else
# define VTK_DELETE_FUNCTION
#endif

#ifndef VTK_OVERRIDE
#define VTK_OVERRIDE override
#endif

class vtkPointLocator;

class vtkCleanUnstructuredGrid
  : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkCleanUnstructuredGrid* New();

  vtkTypeMacro(vtkCleanUnstructuredGrid, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

protected:
  vtkCleanUnstructuredGrid();
  ~vtkCleanUnstructuredGrid();

  vtkPointLocator* Locator;

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;
  virtual int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

private:
  vtkCleanUnstructuredGrid(const vtkCleanUnstructuredGrid&) VTK_DELETE_FUNCTION;
  void operator=(const vtkCleanUnstructuredGrid&) VTK_DELETE_FUNCTION;
};
#endif
