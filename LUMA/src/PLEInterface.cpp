
#include "appdefs_rt.h"

// This will only compile on Linux since PreCICE cannot be built on Windows
#ifdef HAVE_LIBPRECICE

#include "PreciceInterface.h"

template <typename T>
PreciceInterface<T>::PreciceInterface
(
	precice::SolverInterface& precice,
	std::string meshName,
	std::vector<std::string>& readData,
	std::vector<std::string>& writeData,
	InOutRepo<T> &mesh
)
:
precice_(precice),
meshName_(meshName),
readData_(readData),
writeData_(writeData)
{

	// Get the number of points of the original mesh. Just to check that the user doesn't change them. It shouldn't. Also it doesn't 
	// guarantee that the user is not sending another InOutData that happens to have the same number of points in the mesh. 
	numDataLocations_ = mesh.getinitElements();
 
    // Get meshID from precice
    meshID_ = precice_.getMeshID(meshName_);

	// Create preCICE data IDs for the written and read data. 
	for (int i = 0; i < readData.size(); i++)
		readDataIDs_.push_back(precice_.getDataID(readData_.at(i), meshID_));

	for (int i = 0; i < writeData.size(); i++)
		writeDataIDs_.push_back(precice_.getDataID(writeData_.at(i), meshID_));
   
    // Configure the mesh (set the data locations)
    configureMesh(mesh);

 }

template <typename T>
void PreciceInterface<T>::configureMesh(InOutRepo<T>& mesh)
{

	// Array of the indices of the mesh vertices.
    // Each vertex has one index, but three coordinates.
    vertexIDs_ = new int[numDataLocations_];
    
    // Pass the mesh vertices information to preCICE
	//precice_.setMeshVertices(meshID_, numDataLocations_, static_cast<double *>(const_cast<T *>(mesh.getCoordinates())), vertexIDs_);
   precice_.setMeshVertices(meshID_, numDataLocations_, const_cast<T *>(mesh.getCoordinates()), vertexIDs_);
}

template <typename T>
void PreciceInterface<T>::readCouplingData(InOutRepo<T> &data)
{
    // Are new data available or is the participant subcycling?
    if (precice_.isReadDataAvailable())
    {
        // Read all the variables
        for (int i = 0; i < readData_.size(); i++)
        {
			if ( (readData_.at(i).find("v") != std::string::npos) || (readData_.at(i).find("U") != std::string::npos) || (readData_.at(i).find("Re") != std::string::npos)
          || (readData_.at(i).find("u") != std::string::npos))
			{
        
				precice_.readBlockVectorData
				(
					readDataIDs_.at(i),
					numDataLocations_,
					vertexIDs_,
					data.getVectorData(readData_.at(i))
				);
        
       // std::cout << "Data ID " << readDataIDs_.at(i) << " num of elements " << numDataLocations_ << std::endl; //<< " vertex ID 3 " << vertexIDs_[3] << " p at 3 " << data->getInitP()[3] << std::endl;
       // for (int i = 0; i< numDataLocations_ ; i++)
         // std::cout << data->getInitVelocity()[3*i] << " " << data->getInitVelocity()[3*i + 1] << " " << data->getInitVelocity()[3*i + 2] << std::endl;
        
			}
            
            else if ( (readData_.at(i).find("p") != std::string::npos) || (readData_.at(i).find("P") != std::string::npos) )
            {
                precice_.readBlockScalarData
                (
                    readDataIDs_.at(i),
                    numDataLocations_,
                    vertexIDs_,
					data.getScalarData(readData_.at(i))
                );
            }

			else
			{
				std::cout << "Data to read not recognised. Pressure and epsilon data should contain p in its name, velocity data should contain v, u or U and Re stress data should contain Re." << std::endl;
				std::cout << "Input data name: " << readData_.at(i) << std::endl;
				std::cout << "Aborting simulation." << std::endl;
				exit(123);
			}

            // Activate the flag so that GASCANS knows there is new data to read
			data.bReadyToWriteToLBM = true;
            //data->bNewBC = true;
        }
    }
}

template <typename T>
void PreciceInterface<T>::writeCouplingData(InOutRepo<T> &data)
{
	// Write the in InOut to preCICE. This is done every time step, since preCICE can't communicate directly with GASCANS. 
	// Read all the variables
	for (int i = 0; i < writeData_.size(); i++)
	{
		if ((writeData_.at(i).find("v") != std::string::npos) || (writeData_.at(i).find("U") != std::string::npos) || (writeData_.at(i).find("Re") != std::string::npos)
        || (writeData_.at(i).find("u") != std::string::npos))
		{
			precice_.writeBlockVectorData
			(
				writeDataIDs_.at(i),
				numDataLocations_,
				vertexIDs_,
				data.getVectorData(writeData_.at(i))
			);
		}

		else if ((writeData_.at(i).find("p") != std::string::npos) || (writeData_.at(i).find("P") != std::string::npos))
		{
     // std::cout << "Data ID " << readDataIDs_.at(i) << " num of elements " << numDataLocations_ << " vertex ID 3 " << vertexIDs_[3] << " p at 3 " << data->getInitP()[3] << std::endl;
      //for (int i = 0; i< numDataLocations_ ; i++)
      //  std::cout << data->getInitP()[i] << std::endl;
         
			precice_.writeBlockScalarData
			(
				writeDataIDs_.at(i),
				numDataLocations_,
				vertexIDs_,
				data.getScalarData(writeData_.at(i))
			);
      
		}
		else
		{
			std::cout << "Data to write not recognised. Pressure and epsilon data should contain p in its name, and velocity data should contain v, u or U and Re stress data should contain Re." << std::endl;
			std::cout << "Input data name: " << writeData_.at(i) << std::endl;
			std::cout << "Aborting simulation." << std::endl;
			exit(23);
		}

		// Activate the flag so that LUMIS knows there is new data to read
		//data->bOutputBC = true;
		data.bReadyToReadFromLBM = true;
	}
}

template <typename T>
PreciceInterface<T>::~PreciceInterface()
{
    // Delete the vertexIDs_
    delete [] vertexIDs_;
}

// ************************************************************************** //

// Explicit instantiations for single and double precision
//template class PreciceInterface<float>; //It can't be constructed as a float. preCICE doesn't like it. I'll have to fix that. 
template class PreciceInterface<double>;

#endif