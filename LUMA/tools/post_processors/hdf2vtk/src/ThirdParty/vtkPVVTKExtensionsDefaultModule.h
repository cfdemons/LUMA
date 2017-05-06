
#ifndef VTKPVVTKEXTENSIONSDEFAULT_EXPORT_H
#define VTKPVVTKEXTENSIONSDEFAULT_EXPORT_H

#ifdef VTKPVVTKEXTENSIONSDEFAULT_STATIC_DEFINE
#  define VTKPVVTKEXTENSIONSDEFAULT_EXPORT
#  define VTKPVVTKEXTENSIONSDEFAULT_NO_EXPORT
#else
#  ifndef VTKPVVTKEXTENSIONSDEFAULT_EXPORT
#    ifdef vtkPVVTKExtensionsDefault_EXPORTS
        /* We are building this library */
#      define VTKPVVTKEXTENSIONSDEFAULT_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define VTKPVVTKEXTENSIONSDEFAULT_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef VTKPVVTKEXTENSIONSDEFAULT_NO_EXPORT
#    define VTKPVVTKEXTENSIONSDEFAULT_NO_EXPORT 
#  endif
#endif

#ifndef VTKPVVTKEXTENSIONSDEFAULT_DEPRECATED
#  define VTKPVVTKEXTENSIONSDEFAULT_DEPRECATED __declspec(deprecated)
#  define VTKPVVTKEXTENSIONSDEFAULT_DEPRECATED_EXPORT VTKPVVTKEXTENSIONSDEFAULT_EXPORT __declspec(deprecated)
#  define VTKPVVTKEXTENSIONSDEFAULT_DEPRECATED_NO_EXPORT VTKPVVTKEXTENSIONSDEFAULT_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define VTKPVVTKEXTENSIONSDEFAULT_NO_DEPRECATED
#endif

/* AutoInit dependencies.  */
//#include "vtkFiltersParallelStatisticsModule.h"
//#include "vtkIOImportModule.h"
//#include "vtkIOParallelModule.h"
//#include "vtkIOParallelExodusModule.h"
//#include "vtkInteractionWidgetsModule.h"
//#include "vtkPVVTKExtensionsCoreModule.h"
//#include "vtkPVVTKExtensionsRenderingModule.h"

#if (__cplusplus >= 201103L) || ( defined(_MSC_VER) && _MSC_VER >= 1700 )
# define VTK_OVERRIDE override
# define VTK_FINAL final
#else
# define VTK_OVERRIDE
# define VTK_FINAL
#endif

#if (__cplusplus >= 201103L) || ( defined(_MSC_VER) && _MSC_VER >= 1800 )
# define VTK_DELETE_FUNCTION =delete
#else
# define VTK_DELETE_FUNCTION
#endif

#endif
