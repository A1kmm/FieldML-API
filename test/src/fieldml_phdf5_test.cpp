/* \file
 * $Id$
 * \author Alan Wu
 * \brief
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is FieldML
 *
 * The Initial Developer of the Original Code is Auckland Uniservices Ltd,
 * Auckland, New Zealand. Portions created by the Initial Developer are
 * Copyright (C) 2010 the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include "math.h"
#include "FieldmlIoApi.h"
#include "fieldml_api.h"
#include "hdf5.h"
#include <time.h>
#include <vector>

#define PI 3.14159265

using namespace std;

/*
 * Set some global variables
 *---------------------------------------------------------------------------*/
int mpi_rank, mpi_size;
double **highResCoordinates, **lowResCoordinates;
int size = 41, sizeY=3, numberOfTimes = 100;
int arraySize = size * size *size;
double maximumTime = 40.0;
int *arraySizeForEachNode, *numberOfRowsForEachNode;
int lowResSize = 11, lowResNumberOfNodes = 1;
int lowResArraySize = lowResSize * lowResSize * lowResSize;
int lowResElementSize = (lowResSize - 1) * (lowResSize - 1) * (lowResSize - 1);
int highResElementSize = (size - 1) * (size - 1) * (size - 1);
int randomOrdering = 0;

/*
 * Variables for domain decomposition for parallel-HDF5
 *-------------------------------------------------------------------------*/
int count[2], offset[2];
double *write_buf;

/* Find the size and offset of the array which will allow the potential to be
 * calculated
 */
int findUsefulArraySizeForNode(int ratioPerNode, int *offset, int nodeID)
{
	int requiredRows = 1;
	if (lowResNumberOfNodes != 1)
	{
		int upperRow = 0, lowerRow = 0;
		if (nodeID == 0)
		{
			lowerRow = 0;
			upperRow = numberOfRowsForEachNode[0] * ratioPerNode - 1;
		}
		else if (nodeID == lowResNumberOfNodes - 1)
		{
			int temp = -1;
			for (int i = 0; i < nodeID; i ++)
			{
				temp = temp + numberOfRowsForEachNode[i];
			}
			lowerRow = temp * ratioPerNode + 1;
			upperRow = size - 1;
		}
		else
		{
			int temp = -1;
			for (int i = 0; i < nodeID; i ++)
			{
				temp = temp + numberOfRowsForEachNode[i];
			}
			lowerRow = temp * ratioPerNode + 1;
			upperRow = (temp + numberOfRowsForEachNode[nodeID] + 1) * ratioPerNode - 1;
		}
		offset[0] = lowerRow * size *size;
		requiredRows = upperRow - lowerRow + 1;
	}
	else if (nodeID == 0)
	{
		requiredRows = size;
		offset[0] = 0;
	}
	return requiredRows;
}

/* Given the low res identifier and return an high res identifier for a data point*/
int mappedLowResIdenatifierToHighRes(int lowResIdentifier)
{
	return lowResCoordinates[0][lowResIdentifier] + lowResCoordinates[1][lowResIdentifier] * size + lowResCoordinates[2][lowResIdentifier] * size *size;
}

/* read the potential from the hdf5 file output earlier into an array */
double **readHDF5Potential(int ratioPerNode, int *highResOffset)
{
	bool testOk = true;

	printf("Testing HDF5 array read\n");

	int offsets[2];
	int sizes[2];
	sizes[1] = 1;
	double **potentialArray = 0;
	if (lowResNumberOfNodes > mpi_rank)
	{
		if (lowResNumberOfNodes > 1)
		{
			int requiredRows = findUsefulArraySizeForNode(ratioPerNode, &offsets[0], mpi_rank);
			sizes[0] = requiredRows * size * size;
			potentialArray = new double*[numberOfTimes];
			for( int i = 0; i < numberOfTimes ; i++)
			{
				potentialArray[i] = new double[sizes[0]];
			}
		}
		else
		{
			potentialArray = new double*[numberOfTimes];
			for( int i = 0; i < numberOfTimes ; i++)
			{
				potentialArray[i] = new double[arraySize];
			}
			sizes[0] = arraySize;
			offsets[0] = 0;
		}
	}

	FmlSessionHandle session = Fieldml_Create( "test", "test" );
	FmlObjectHandle resource = 0;
	if (randomOrdering == 0)
		resource = Fieldml_CreateHrefDataResource( session, "potential.resource2", "PHDF5", "potential.h5" );
	else
		resource = Fieldml_CreateHrefDataResource( session, "potential.resource2", "PHDF5", "randomPotential.h5" );
	FmlObjectHandle sourceD = Fieldml_CreateArrayDataSource( session, "potential.source2_double", resource, "potential", 2 );

	FmlObjectHandle reader = Fieldml_OpenReader( session, sourceD );

	if (lowResNumberOfNodes > mpi_rank)
	{
		printf("Reading: mpi_size %d, mpi_rank %d highresOffset %d ratio %d \n", lowResNumberOfNodes, mpi_rank,offsets[0], ratioPerNode);

		for( int i = 0; i < numberOfTimes ; i++)
		{
			offsets[1] = i;
			Fieldml_ReadDoubleSlab( reader, offsets, sizes, &(potentialArray[i][0]) );
		}
	}
	//for (int i = 0; i< arraySize; i++)
	//	printf("mpi_rank %d potentialArray %d %g\n", mpi_rank, i+offsets[0], potentialArray[0][i]);

	Fieldml_CloseReader( reader );

	Fieldml_Destroy( session );

	*highResOffset = offsets[0];

	if( testOk )
	{
		printf( "TestHdf5Read - ok\n" );
	}
	else
	{
		printf( "TestHdf5Read - failed\n" );
	}

	return potentialArray;
}

/* find if the data point of the given identifier is on a boundary */
int isNodeOnBoundary(int identifier)
{
	for (int i = 0; i < sizeY; i++)
	{
		double coordinates = lowResCoordinates[i][identifier];
		if (coordinates == 0.0 || coordinates == (double)(size - 1))
			return 1;
	}
	return 0;
}

/* this will get the average potential around the
 * area of a data point with the given identifier */
double findAveragePotentialLevel(double **potentialArray, int identifier,
	int time, double ratioPerNode)
{
	int lowerLimit = -(ratioPerNode - 1);
	int upperLimit = ratioPerNode - 1;
	int numberOfNodesTaken = 0;
	double totalPotential = 0.0;
	int localIdentifier = 0;
	for (int z = lowerLimit; z <= upperLimit; z++)
	{
		localIdentifier = identifier + (z * size * size);
		int current_z_id = localIdentifier;
		for (int y = lowerLimit; y <= upperLimit; y++)
		{
			localIdentifier = current_z_id + (y * size);
			int current_y_id = localIdentifier;
			for (int x = lowerLimit; x <= upperLimit; x++)
			{
				localIdentifier = current_y_id + x;
				totalPotential += potentialArray[time][localIdentifier];
				numberOfNodesTaken++;
			}
		}
	}
	return (totalPotential/(double)numberOfNodesTaken);
}

/* find the activation time for all the data in this node and put into the activationTimes Array */
void findActivationTimes(double **potentialArray,double *activationTimesArray,
	int localOffset, int highResOffset, double ratioPerNode)
{
	printf("findActivationTimes localOffset %d\n", localOffset);
	for (int i = 0; i < arraySizeForEachNode[mpi_rank]; i++)
	{
		activationTimesArray[i] = 100.0;
		int found = 0;
		int highResIdentifier = mappedLowResIdenatifierToHighRes(localOffset + i);

		for (int j = 0; j < numberOfTimes && !found; j++)
		{
			if (isNodeOnBoundary(localOffset + i))
			{
				if (potentialArray[j][highResIdentifier - highResOffset] > (double)0.8)
				{
					found = 1;
				}
			}
			else
			{
				if (findAveragePotentialLevel(potentialArray, highResIdentifier - highResOffset, j, ratioPerNode) > (double)0.8)
				{
					found = 1;
				}
			}
			if (found)
			{
				//printf("highResIdentifier %d, localIdentifier %d potential %d value %g\n", highResIdentifier, localOffset + i, j, potentialArray[j][highResIdentifier - highResOffset]);
				activationTimesArray[i] = (double)j * maximumTime / (numberOfTimes - 1);
			}
		}
	}
}

/* the following function will reorder the array
 * and make it easier to find the activation time */
double **reorderPotentialArray(int ratioPerNode, int highResOffset, double **potentialArray)
{
	bool testOk = true;

	double **reorderPotentialArray = 0;

	printf("Testing HDF5 array read\n");

	int offsets = 0;
	int sizes = 0;
	int *orderingArray = 0;
	if (lowResNumberOfNodes > mpi_rank)
	{
		orderingArray = new int[arraySize];
		sizes = arraySize;
		offsets = 0;
	}

	/* First read the ordering into an array */
	FmlSessionHandle session = Fieldml_Create( "test", "test" );

	FmlObjectHandle indexResource = Fieldml_CreateHrefDataResource( session, "index.resource", "PHDF5", "randomIndex.h5" );

	FmlObjectHandle indexSourceD = Fieldml_CreateArrayDataSource( session, "index.source_int", indexResource, "id", 1 );

	FmlObjectHandle reader = Fieldml_OpenReader( session, indexSourceD );

	if (lowResNumberOfNodes > mpi_rank)
	{
		Fieldml_ReadIntSlab( reader, &offsets, &sizes, orderingArray);
	}

	Fieldml_CloseReader( reader );

	Fieldml_Destroy( session );

	/* the following part will allow each node to communicate with each other
	 * and allow each node to reconstruct an array sufficient to get the activation time */
	if (lowResNumberOfNodes > mpi_rank)
	{
		int requiredRows = findUsefulArraySizeForNode(ratioPerNode, &offsets, mpi_rank);
		sizes = requiredRows * size * size;
		int localStart = offsets, localEnd = offsets + sizes - 1;
		reorderPotentialArray = new double*[numberOfTimes];
		for( int i = 0; i < numberOfTimes ; i++)
		{
			reorderPotentialArray[i] = new double[sizes];
		}
		vector<int> idToSend[lowResNumberOfNodes];
		vector<int> idToReceive[lowResNumberOfNodes];
	   int ierr;
		if ( lowResNumberOfNodes > 1)
		{
			double **sendbuff,**recvbuff;
			sendbuff = new double*[lowResNumberOfNodes];
			recvbuff = new double*[lowResNumberOfNodes];
		   MPI_Status status;
		   MPI_Request	send_request[lowResNumberOfNodes],recv_request[lowResNumberOfNodes];
			int currentID = 0;
			int buffersize = 0;
			for (int i = 0; i < lowResNumberOfNodes; i++)
			{
				requiredRows = findUsefulArraySizeForNode(ratioPerNode, &offsets, i);
				sizes = requiredRows * size * size;
				int remoteStart = offsets, remoteEnd = offsets +sizes - 1;
				idToSend[i].clear();
				idToReceive[i].clear();
				printf("mpi_rank %d, remote rank %d local start %d local end %d remote start %d remote end %d \n",
					mpi_rank, i, localStart, localEnd, remoteStart, remoteEnd);
				vector<int> localArrayID;
				localArrayID.clear();
				for (int j = localStart; j <= localEnd; j++)
				{
					currentID = orderingArray[j] - 1;
					if ( remoteStart <= currentID && currentID <= remoteEnd)
					{
						idToSend[i].push_back(currentID);
						localArrayID.push_back(j - localStart);
					}
				}
				buffersize = idToSend[i].size() * numberOfTimes;
				sendbuff[i] = new double[buffersize];
				vector<int>::iterator pos = idToSend[i].begin();
				int currentOrder = 0;
				/* construct an array to send value to other nodes */
				while (pos != idToSend[i].end())
				{
					for (int k = 0; k < numberOfTimes; k ++)
					{
						sendbuff[i][currentOrder * 100 + k] = potentialArray[k][localArrayID[currentOrder]];
					}
					//printf("mpi_rank %d to %d currentID %d local array %d value %g \n", mpi_rank, i, *pos, localArrayID[currentOrder], potentialArray[0][localArrayID[currentOrder]]);
					++currentOrder;
					++pos;
				}
				if (i != mpi_rank)
				{
					ierr=MPI_Isend(sendbuff[i],buffersize,MPI_DOUBLE,
						i,mpi_rank,MPI_COMM_WORLD,&send_request[i]);
					printf("mpi_rank %d, to send %d to remote rank %d  \n",
						mpi_rank, buffersize, i);
				}

				localArrayID.clear();
				for (int j = remoteStart; j <= remoteEnd; j++)
				{
					currentID = orderingArray[j] - 1;
					if ( localStart <= currentID && currentID <= localEnd)
					{
						idToReceive[i].push_back(currentID);
						localArrayID.push_back(j - remoteStart);
					}
				}
				buffersize = idToReceive[i].size() * numberOfTimes;
				recvbuff[i] = new double[buffersize];
				if (i != mpi_rank)
				{
					ierr=MPI_Irecv(recvbuff[i],buffersize,MPI_DOUBLE,
						i,i,MPI_COMM_WORLD,&recv_request[i]);
					printf("mpi_rank %d, to receive %d from remote rank %d  \n",
						mpi_rank, buffersize, i);
				}
				else
				{
					/* it is local, assign value here */
					vector<int>::iterator pos = idToReceive[i].begin();
					/* construct an array to send value to other nodes */
					currentOrder = 0;
					while (pos != idToReceive[i].end())
					{
						for (int k = 0; k < numberOfTimes; k ++)
						{
							reorderPotentialArray[k][*pos - localStart] = potentialArray[k][localArrayID[currentOrder]];
						}
				//		printf("mpi_rank %d to %d currentID %d local array %d value %g \n", mpi_rank, i, *pos, localArrayID[currentOrder], potentialArray[0][localArrayID[currentOrder]]);
						++currentOrder;
						++pos;
					}
				}
			}
			for (int i = 0; i < lowResNumberOfNodes; i++)
			{
				if (i != mpi_rank)
				{
					ierr=MPI_Wait(&send_request[i],&status);
					ierr=MPI_Wait(&recv_request[i],&status);
				}
				if (sendbuff[i])
				{
					delete[] sendbuff[i];
				}

				if (i != mpi_rank)
				{
					vector<int>::iterator pos = idToReceive[i].begin();
					/* construct an array to send value to other nodes */
					int currentOrder = 0;
					while (pos != idToReceive[i].end())
					{
						for (int k = 0; k < numberOfTimes; k ++)
						{
							reorderPotentialArray[k][*pos-localStart] = recvbuff[i][currentOrder * 100 + k];
						}
					//	printf("currentOrder %d mpi_rank %d from %d pos %d local start %d value %g\n", currentOrder, mpi_rank, i, *pos, localStart, recvbuff[i][currentOrder * 100 + 0]);
						++currentOrder;
						++pos;
					}
				}
				if (recvbuff[i])
				{
					delete[] recvbuff[i];
				}
			}
			if (sendbuff)
				delete[] sendbuff;
			if (recvbuff)
				delete[] recvbuff;
		}
		else
		{
			for (int j = localStart; j <= localEnd; j++)
			{
				int currentID = orderingArray[j] - 1;
				for (int k = 0; k < numberOfTimes; k++)
				{
					reorderPotentialArray[k][currentID] = potentialArray[k][j];
				}
			}
		}
	}

	if (orderingArray)
		delete[] orderingArray;

	if( testOk )
	{
		printf( "TestHdf5Read - ok\n" );
	}
	else
	{
		printf( "TestHdf5Read - failed\n" );
	}


	return reorderPotentialArray;
}

/* calculate the analytical result and write them out */
int writeLowResActivationTime(FmlSessionHandle session, FmlObjectHandle cType, FmlObjectHandle sourceD)
{
	bool testOk = true;
	int offsets;
	int sizes;
	int ratioPerNode = (size - 1) / (lowResSize - 1);
	int highResOffset = 0;
	double *activationTimes = NULL;

	/* read the potential from the hdf5 file output earlier into an array */
	double **potentialArray = readHDF5Potential(ratioPerNode, &highResOffset);

	/* if random ordering is enabled, the following part will reorder the array
	 * and make it easier to find the activation time */
	if (randomOrdering)
	{
		double **temp = NULL;
		temp = potentialArray;
		potentialArray = reorderPotentialArray(ratioPerNode, highResOffset, temp);
		if (temp)
		{
			for ( int i = 0; i < numberOfTimes; i++ )
			{
				delete[] temp[i];
			}
			delete[] temp;
		}
	}

	sizes = lowResArraySize;

	FmlObjectHandle writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, &sizes, 1 );

	if (lowResNumberOfNodes > mpi_rank)
	{
		sizes = arraySizeForEachNode[mpi_rank];
		offsets = 0;

		for (int i = 0; i < mpi_rank; i ++)
		{
			offsets = offsets + arraySizeForEachNode[i];
		}
		activationTimes = new double[arraySizeForEachNode[mpi_rank]];
		/* find the activation time for all the data in this node and put into the activationTimes Array */
		findActivationTimes(potentialArray, &activationTimes[0], offsets, highResOffset, ratioPerNode);
		Fieldml_WriteDoubleSlab( writer, &offsets, &sizes, activationTimes );
	}

	Fieldml_CloseWriter( writer );

	if (potentialArray)
	{

		for ( int i = 0; i < numberOfTimes; i++ )
		{
			delete[] potentialArray[i];
		}
		delete[] potentialArray;
	}
	if (activationTimes)
		delete[] activationTimes;

	return testOk;
}


/*===========================================================================
  WRITE_LOWRES_COORDINATES

  C++ subroutine that writes the coordinate data to hard disk using the
  FieldML HDF5 API.
  ===========================================================================*/

bool write_lowres_coordinates( FmlSessionHandle session, FmlObjectHandle cType, 
                               FmlObjectHandle sourceD ) {

     bool testOk;
     int count[2], offset[2], remainder;
     FmlObjectHandle writer;
     FmlIoErrorNumber status;

  /*
   * Open access to the HDF5 file containing the coordinate data
   *-------------------------------------------------------------------------*/
     count[1] = lowResArraySize;
     count[0] = sizeY;
     writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, count, 2 );

  /*
   * Determine how many coordinate data values are written by the local MPI
   * task and its corresponding starting/offset point.
   *-------------------------------------------------------------------------*/
     count[1] = lowResArraySize / mpi_size;
     count[0] = sizeY;
     offset[1] = count[0]*mpi_rank;
     offset[0] = 0;

     remainder = lowResArraySize % mpi_size;
     if ( remainder>0 ) {
        if ( mpi_rank<remainder ) { count[1]++; }
        else { offset[1] += remainder-1; }
     }

  /*
   * Write the coordinate data to hard disk
   *-------------------------------------------------------------------------*/
     status = Fieldml_WriteDoubleSlab( writer, offset, count, 
                                 &(lowResCoordinates[offset[0]][offset[1]]) );    
     if ( status==FML_IOERR_NO_ERROR ) { testOk = true; }
     else { testOk = false; }

     Fieldml_CloseWriter( writer ); 
     return testOk;

}


int findElementsRangeToExport(int numebrOfNodes, int nodeId, int elementSize,int *lowElemID, int *highElemID)
{
	int remainder = elementSize % numebrOfNodes;
	int lowSize =  (int)(elementSize / numebrOfNodes);
	int lowElemLocal = 1;
	int highElemLocal = 0;
	for (int i = 0; i < numebrOfNodes; i++)
	{
		highElemLocal += lowSize;
		if (remainder != 0)
		{
			highElemLocal++;
			remainder--;
		}
		if (nodeId == i)
		{
			*lowElemID = lowElemLocal;
			*highElemID = highElemLocal;
			return *highElemID - *lowElemID + 1;
		}
		lowElemLocal = highElemLocal + 1;
	}
	return 0;
}

int getFirstNodeNumberForElement(int elemNo, int lowResSize)
{
	int numberOfElemsPerRow = lowResSize - 1;
	int elemNoMinusOne = elemNo - 1;
	int numOnZ = (int) ((elemNoMinusOne) / ((numberOfElemsPerRow) * (numberOfElemsPerRow)));
	int remainingElem = (elemNoMinusOne % ((numberOfElemsPerRow) * (numberOfElemsPerRow)));
	int numOnY = (int) (remainingElem / (numberOfElemsPerRow));
	int numOnX = (remainingElem % (numberOfElemsPerRow));
	return 1 + numOnX + numOnY * lowResSize + numOnZ * lowResSize * lowResSize;

}

/* construct an array with data of connectivity */
int **constructConnectivityArray(int nodesPerElement, int *numberOfElements, int *firstElementNo, int numberCompNodes, int resSize, int elementSize)
{
	int firstElemID = 1, lastElemID = 1, sizeToExport= 1;
	int **connectivityArray = NULL;

	if (numberCompNodes > mpi_rank)
	{
		sizeToExport = findElementsRangeToExport(numberCompNodes, mpi_rank, elementSize, &firstElemID, &lastElemID);
		printf("Node No: %d Element array size %d firstElemID %d lastElemID %d \n",
			mpi_rank, sizeToExport, firstElemID, lastElemID);

		connectivityArray = new int*[nodesPerElement];
		for(int i = 0; i < nodesPerElement ; i++)
		{
			connectivityArray[i] = new int[sizeToExport];
		}
		int j = 0;
		for (int i = firstElemID; i <= lastElemID; i++)
		{
			connectivityArray[0][j] = getFirstNodeNumberForElement(i, resSize);
			connectivityArray[1][j] = connectivityArray[0][j] + 1;
			connectivityArray[2][j] = connectivityArray[0][j] + resSize;
			connectivityArray[3][j] = connectivityArray[2][j] + 1;
			connectivityArray[4][j] = connectivityArray[0][j] + resSize *  resSize;
			connectivityArray[5][j] = connectivityArray[4][j] + 1;
			connectivityArray[6][j] = connectivityArray[4][j] + resSize;
			connectivityArray[7][j] = connectivityArray[6][j] + 1;
			j++;
		}
		*numberOfElements = sizeToExport;
		*firstElementNo = firstElemID;
	}
	return connectivityArray;
}

/* write out the connectivity */
int writeConnectivity(FmlSessionHandle session, FmlObjectHandle emsembleType, FmlObjectHandle sourceD,
	int elementSize, int numberCompNodes, int resSize)
{
	bool testOk = true;
	int offsets[2];
	int sizes[2];
	int firstElementNo = 1;
	int nodesPerElement = pow(2.0, sizeY);
	sizes[0] = elementSize;
	sizes[1] = nodesPerElement;
	int numberOfElements = 1;
	int **connectivityArray = constructConnectivityArray(nodesPerElement, &numberOfElements, &firstElementNo, numberCompNodes, resSize, elementSize);

	FmlObjectHandle writer = Fieldml_OpenArrayWriter( session, sourceD, emsembleType, 0, sizes, 2 );

	if (numberCompNodes > mpi_rank)
	{
		sizes[0] = numberOfElements;
		sizes[1] = 1;
		offsets[0] = firstElementNo - 1;

		for( int i = 0; i < nodesPerElement ; i++)
		{
			offsets[1] = i;
			Fieldml_WriteIntSlab( writer, offsets, sizes, &(connectivityArray[i][0]) );
		}
	}

	Fieldml_CloseWriter( writer );

	if (connectivityArray)
	{
		for ( int i = 0; i < nodesPerElement; i++ )
		{
			delete[] connectivityArray[i];
		}
		delete[] connectivityArray;
	}

	return testOk;
}

/* Function to find the array size per computational node, it is done
 * differently from the high res one. This function find the number of rows
 * each node will write instead of the total number of data to write out */
void findLowResArraySizePerNode(int resSize, int numebrOfNodes)
{
	arraySizeForEachNode = new int[numebrOfNodes];
	numberOfRowsForEachNode = new int[numebrOfNodes];
	int remainder = resSize % numebrOfNodes;
	int lowSize =  (int)(resSize / numebrOfNodes);
	for (int i = 0; i < numebrOfNodes; i++)
	{
		numberOfRowsForEachNode[i] = lowSize;
		if (remainder != 0)
		{
			numberOfRowsForEachNode[i]++;
			remainder--;
		}
		arraySizeForEachNode[i] = numberOfRowsForEachNode[i] * resSize * resSize;
	}
}

/* Function to analyse the data in high resolution and write out
 * the solution in lower resolution */

/*===========================================================================
  LOWRES_HDF5_WRITE

  C++ subroutine that writes the coordinate data to hard disk using the
  FieldML HDF5 API.
  ===========================================================================*/

bool lowResHdf5Write( FmlSessionHandle session, FmlObjectHandle realType,
	              FmlObjectHandle sourceD, FmlObjectHandle potentialSourceD,
	              FmlObjectHandle emsembleType, FmlObjectHandle connectivityD ) {

     bool testOk = write_lowres_coordinates(session, realType, sourceD);
/*


	findLowResArraySizePerNode(lowResSize, lowResNumberOfNodes);

	if (lowResNumberOfNodes > mpi_rank)
	{
		printf("mpi_size %d, mpi_rank %d number of rows %d size %d \n", lowResNumberOfNodes, mpi_rank,
			numberOfRowsForEachNode[mpi_rank], arraySizeForEachNode[mpi_rank]);
	}

	if (arraySizeForEachNode[mpi_rank] > 0)
	{

		printf("mpi_rank %d low res HDF5 array write1\n", mpi_rank);

		printf("mpi_rank %d low res HDF5 array write2\n", mpi_rank);
		testOk = testOk && writeLowResActivationTime(session, realType, potentialSourceD);
		printf("mpi_rank %d low res HDF5 array write3\n", mpi_rank);
		testOk = testOk && writeConnectivity(session, emsembleType, connectivityD, lowResElementSize, lowResNumberOfNodes, lowResSize);
	}

	
	delete[] arraySizeForEachNode;

	delete[] numberOfRowsForEachNode;
*/
     return testOk;
}

/* Function to write out the solution in lower resolution, it will handle
 * arbitrary and regular ordering data differently */
bool writeLowResFieldMLSolution() {

     bool testOk;
     int importHandle;
     FmlSessionHandle session;
     FmlObjectHandle chart3dArgumentHandle, realType, shapeType, 
                     trilinearPointsArgumentHandle, trilinearLagrangeParametersHandle, 
                     trilinearLagrangeParameteArgumentrsHandle, trilinearInterpolatorHandle, 
                     coordinatesRC3DComponenentHandle, coordinatesRC3DCHandle,
                     myEnsmebleHandle, nodesArgumentHandle, cubeMeshHandle, 
                     meshChartHandle, meshEnsembleHandle;

     session = Fieldml_Create( "test", "test" );
     importHandle = Fieldml_AddImportSource( session, "../FieldML_Library_0.5.xml",
                                            "library" );
     chart3dArgumentHandle = Fieldml_AddImport( session, importHandle, 
                                                "chart.3d.argument", 
                                                "chart.3d.argument" );
     realType = Fieldml_AddImport( session, importHandle, "real.type", 
                                   "real.1d" );
     shapeType = Fieldml_AddImport( session, importHandle, "shape.unit.cube", 
                                   "shape.unit.cube" );
     trilinearPointsArgumentHandle = Fieldml_AddImport( session, importHandle,
   	                                     "trilinearLagrange.points.argument", 
                       "parameters.3d.unit.trilinearLagrange.component.argument" );
     trilinearLagrangeParametersHandle = Fieldml_AddImport( session, importHandle,
   	                                           "trilinearLagrange.parameters", 
                                           "parameters.3d.unit.trilinearLagrange" );
     trilinearLagrangeParameteArgumentrsHandle = Fieldml_AddImport( session, 
                                                    importHandle,
   	                                  "trilinearLagrange.parameters.argument", 
                                  "parameters.3d.unit.trilinearLagrange.argument" );
     trilinearInterpolatorHandle = Fieldml_AddImport( session, importHandle,
                                                 "trilinearLagrange.interpolator", 
                                         "interpolator.3d.unit.trilinearLagrange" );
     coordinatesRC3DComponenentHandle = Fieldml_AddImport( session, importHandle,
                                           "coordinates.rc.3d.component.argument", 
                                           "coordinates.rc.3d.component.argument" );
     coordinatesRC3DCHandle = Fieldml_AddImport( session, importHandle,
                                                 "coordinates.rc.3d", 
                                                 "coordinates.rc.3d" );

     myEnsmebleHandle = Fieldml_CreateEnsembleType( session, "cube.nodes" );
     Fieldml_SetEnsembleMembersRange( session, myEnsmebleHandle, 1, 
                                      lowResArraySize, 1 );

     nodesArgumentHandle = Fieldml_CreateArgumentEvaluator( session, 
                                  "cube.nodes.argument", myEnsmebleHandle );

     cubeMeshHandle = Fieldml_CreateMeshType( session, "cube.mesh" );
     meshEnsembleHandle = Fieldml_CreateMeshElementsType( session, 
                                        cubeMeshHandle, "elements" );
     Fieldml_SetEnsembleMembersRange( session, meshEnsembleHandle, 1, 
                                      lowResElementSize, 1 );
     meshChartHandle = Fieldml_CreateMeshChartType( session, cubeMeshHandle, 
                                                   "chart" );
     Fieldml_CreateContinuousTypeComponents( session, meshChartHandle, 
                                           "cube.mesh.chart.component", 3 );
     Fieldml_SetMeshShapes( session, cubeMeshHandle, shapeType );

     Fieldml_CreateArgumentEvaluator( session, "cube.mesh.argument", 
                                      cubeMeshHandle );

   FmlObjectHandle cubeDOFsNodeArgumentHandle = Fieldml_CreateArgumentEvaluator( session, "cube.dofs.node.argument", realType );
   Fieldml_AddArgument( session, cubeDOFsNodeArgumentHandle, nodesArgumentHandle );

	FmlObjectHandle connectivityResource = Fieldml_CreateHrefDataResource( session,
		"cube.component1.trilinearLagrange.connectivity.resource", "PHDF5", "lowResConnectivity.h5" );
	FmlObjectHandle sourceConnectivityD = Fieldml_CreateArrayDataSource( session,
		"cube.component1.trilinearLagrange.connectivity.resource.data", connectivityResource, "1", 2 );
	int rawSizes[2];
	rawSizes[0]= lowResElementSize;
	rawSizes[1]= (int)pow(2.0, sizeY);
	Fieldml_SetArrayDataSourceRawSizes( session, sourceConnectivityD, rawSizes );
	Fieldml_SetArrayDataSourceSizes( session, sourceConnectivityD, rawSizes );

	FmlObjectHandle connectivityParameterHandle = Fieldml_CreateParameterEvaluator( session,
		"cube.component1.trilinearLagrange.connectivity", myEnsmebleHandle );
	Fieldml_SetParameterDataDescription( session, connectivityParameterHandle, FML_DATA_DESCRIPTION_DENSE_ARRAY );
	Fieldml_SetDataSource( session, connectivityParameterHandle, sourceConnectivityD );
	Fieldml_AddDenseIndexEvaluator( session, connectivityParameterHandle,
		Fieldml_GetObjectByName( session, "cube.mesh.argument.elements"), FML_INVALID_HANDLE );
	Fieldml_AddDenseIndexEvaluator( session, connectivityParameterHandle, trilinearPointsArgumentHandle, FML_INVALID_HANDLE );

	FmlObjectHandle cubeTrilinearLagrangeParameters = Fieldml_CreateAggregateEvaluator( session,
		"cube.component1.trilinearLagrange.parameters", trilinearLagrangeParametersHandle );
	Fieldml_SetIndexEvaluator( session, cubeTrilinearLagrangeParameters, 1, trilinearPointsArgumentHandle );
	Fieldml_SetDefaultEvaluator( session, cubeTrilinearLagrangeParameters, cubeDOFsNodeArgumentHandle );
	Fieldml_SetBind(session, cubeTrilinearLagrangeParameters, nodesArgumentHandle, connectivityParameterHandle);

	FmlObjectHandle cubeReferenceHandle = Fieldml_CreateReferenceEvaluator( session,
		"cube.component1.trilinearLagrange", trilinearInterpolatorHandle );
	Fieldml_SetBind(session, cubeReferenceHandle, chart3dArgumentHandle,
		Fieldml_GetObjectByName( session, "cube.mesh.argument.chart"));
	Fieldml_SetBind(session, cubeReferenceHandle, trilinearLagrangeParameteArgumentrsHandle,
		cubeTrilinearLagrangeParameters);

	FmlObjectHandle cubeTemplatehandle  = Fieldml_CreatePiecewiseEvaluator( session, "cube.component1.template", realType );
	Fieldml_SetIndexEvaluator( session, cubeTemplatehandle, 1, Fieldml_GetObjectByName( session, "cube.mesh.argument.elements") );
	Fieldml_SetDefaultEvaluator( session, cubeTemplatehandle, cubeReferenceHandle );

	FmlObjectHandle nodeResource = Fieldml_CreateHrefDataResource( session,
		"cube.geometric.dofs.node.resource", "PHDF5", "lowrescoordinates.h5" );
	FmlObjectHandle nodeSourceD = Fieldml_CreateArrayDataSource( session,
		"cube.geometric.dofs.node.data", nodeResource, "1", 2 );
	rawSizes[0]= lowResArraySize;
	rawSizes[1]= 3;
	Fieldml_SetArrayDataSourceRawSizes( session, nodeSourceD, rawSizes );
	Fieldml_SetArrayDataSourceSizes( session, nodeSourceD, rawSizes );

	FmlObjectHandle cubeNodeParameterHandle = Fieldml_CreateParameterEvaluator( session,
		"cube.geometric.dofs.node", realType );
	Fieldml_SetParameterDataDescription( session, cubeNodeParameterHandle, FML_DATA_DESCRIPTION_DENSE_ARRAY );
	Fieldml_SetDataSource( session, cubeNodeParameterHandle, nodeSourceD );
	Fieldml_AddDenseIndexEvaluator( session, cubeNodeParameterHandle,
		nodesArgumentHandle, FML_INVALID_HANDLE );
	Fieldml_AddDenseIndexEvaluator( session, cubeNodeParameterHandle, coordinatesRC3DComponenentHandle, FML_INVALID_HANDLE );

	FmlObjectHandle cubeGiometricParameters = Fieldml_CreateAggregateEvaluator( session,
		"cube.geometric.parameters", coordinatesRC3DCHandle );
	Fieldml_SetIndexEvaluator( session, cubeGiometricParameters, 1, coordinatesRC3DComponenentHandle );
	Fieldml_SetDefaultEvaluator( session, cubeGiometricParameters, cubeTemplatehandle );
	Fieldml_SetBind(session, cubeGiometricParameters, cubeDOFsNodeArgumentHandle, cubeNodeParameterHandle);

	FmlObjectHandle potentialResource = Fieldml_CreateHrefDataResource( session, "cube.geometric.dofs.activationtime.resource",
		"PHDF5", "lowresActivationTime.h5" );
	FmlObjectHandle potentialSourceD = Fieldml_CreateArrayDataSource( session, "cube.geometric.dofs.activationtime.data",
		potentialResource, "1", 1 );
	Fieldml_SetArrayDataSourceRawSizes( session, potentialSourceD, &lowResArraySize );
	Fieldml_SetArrayDataSourceSizes( session, potentialSourceD, &lowResArraySize );

	FmlObjectHandle nodeActivationParameterHandle = Fieldml_CreateParameterEvaluator( session,
		"cube.activationTime.dofs", realType );
	Fieldml_SetParameterDataDescription( session, nodeActivationParameterHandle, FML_DATA_DESCRIPTION_DENSE_ARRAY );
	Fieldml_SetDataSource( session, nodeActivationParameterHandle, potentialSourceD );
	Fieldml_AddDenseIndexEvaluator( session, nodeActivationParameterHandle,
		nodesArgumentHandle, FML_INVALID_HANDLE );

	FmlObjectHandle activationReferenceHandle = Fieldml_CreateReferenceEvaluator( session,
		"cube.activationTime.source", cubeTemplatehandle );
	Fieldml_SetBind(session, activationReferenceHandle, Fieldml_GetObjectByName( session, "cube.dofs.node.argument"),
		nodeActivationParameterHandle);

     testOk = lowResHdf5Write( session, realType, nodeSourceD, 
                               potentialSourceD, myEnsmebleHandle, 
                               sourceConnectivityD );

     if ( mpi_rank==0 )
        Fieldml_WriteFile( session, "phdf5_test.xml" );

     Fieldml_Destroy( session );
     return testOk;
}

/* calculate the potential at a given coordinates and time */
double findPotential(double x, double y, double z, double time)
{
	double propagation_speed = maximumTime / ( numberOfTimes - 1 );
	double radius = (( size - 1) / 2 );
	double distanceX = x - radius;
	double distanceY = y - radius;
	double distanceZ = z - radius;
	double distance = sqrt(distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
	double activation_level = cos ( (M_PI / (double)2.0 * ( distance - (propagation_speed * time)) / radius ) );
	if (0.0 > activation_level )
		activation_level = 0.0;
	return activation_level;
}


/*===========================================================================
  WRITE_HIGHRES_POTENTIAL

  C++ subroutine that writes the potential data to hard disk using the
  FieldML HDF5 API.
  ===========================================================================*/

bool write_highres_potential( FmlSessionHandle session, FmlObjectHandle cType, 
                              FmlObjectHandle sourceD ) {

     bool testOk=true;
     int i, j, count[2], offset[2], remainder, nodeID;
     FmlObjectHandle writer;
     FmlIoErrorNumber status;
     double time, increment, *potentialArray;

  /*
   * Open access to the HDF5 file containing the potential data
   *-------------------------------------------------------------------------*/
     count[0] = arraySize;
     count[1] = numberOfTimes;
     writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, count, 2 );

  /*
   * Determine how many data values are written by the local MPI
   * task and its corresponding starting/offset point.
   *-------------------------------------------------------------------------*/
     count[0] = arraySize / mpi_size;
     count[1] = 1;
     offset[0] = count[0]*mpi_rank;

     remainder = arraySize % mpi_size;
     if ( remainder>0 ) {
        if ( mpi_rank<remainder ) { count[0]++; }
        else { offset[0] += remainder-1; }
     }

  /*
   * Construct the local potential data for each MPI task 
   *-------------------------------------------------------------------------*/
     potentialArray = new double[count[0]];
     increment = maximumTime / (double)( numberOfTimes - 1 );

     for ( i=0; i<numberOfTimes; i++ ) {
         time = increment * (double)i;
         for ( j=0; j<count[0]; j++ ) {
             nodeID = offset[0] + j;
             potentialArray[j] = 1.234;
//             potentialArray[j] = findPotential( highResCoordinates[0][nodeID], 
//                                                highResCoordinates[1][nodeID], 
//                                                highResCoordinates[2][nodeID], 
//                                                time );
         }
         offset[1] = i;
         status = Fieldml_WriteDoubleSlab( writer, offset, count, 
                                           potentialArray );
         if ( status!=FML_IOERR_NO_ERROR ) { 
            testOk = false;
            return testOk; 
         }
     }

     Fieldml_CloseWriter( writer ); 
     return testOk;

}

/* Routine to write the potential out using FieldML hdf5 APIs */
bool writeHighResPotential( FmlSessionHandle session, FmlObjectHandle cType, 
                            FmlObjectHandle sourceD ) {

	bool testOk = true;
	int offsets[2];
	int sizes[1];

	sizes[0] = arraySize;
	sizes[1] = numberOfTimes;

	FmlObjectHandle writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, sizes, 2 );

	sizes[0] = arraySizeForEachNode[mpi_rank];
	sizes[1] = 1;
	offsets[0] = 0;
	for (int i = 0; i < mpi_rank; i ++)
	{
		offsets[0] = offsets[0] + arraySizeForEachNode[i];
	}
	double *potentialArray = new double[arraySizeForEachNode[mpi_rank]];
	double increment = maximumTime / ( numberOfTimes - 1 );
	for( int i = 0; i < numberOfTimes ; i++)
	{
		double time = increment * i;
		for (int j = 0; j < arraySizeForEachNode[mpi_rank]; j++)
		{
			int nodeID = offsets[0] + j;
			potentialArray[j] = findPotential(highResCoordinates[0][nodeID], highResCoordinates[1][nodeID], highResCoordinates[2][nodeID], time );
		}
		offsets[1] = i;
		Fieldml_WriteDoubleSlab( writer, offsets, sizes, potentialArray );
	}

	Fieldml_CloseWriter( writer );

	return testOk;
}

/*===========================================================================
  HIGHRES_WRITE

  C++ subroutine that writes the coordinate data to hard disk using the
  FieldML HDF5 API.
  ===========================================================================*/

bool highres_write( FmlSessionHandle session, FmlObjectHandle cType, 
                    FmlObjectHandle sourceD, int num_slices ) {

     bool testOk;
     int count[2], offset[2], remainder;
     FmlObjectHandle writer;
     FmlIoErrorNumber status;

  /*
   * Open access to the HDF5 file containing the coordinate data
   *-------------------------------------------------------------------------*/
     count[1] = arraySize;
     count[0] = num_slices;
     writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, count, 2 );

  /*
   * Determine how many coordinate data values are written by the local MPI
   * task and its corresponding starting/offset point.
   *-------------------------------------------------------------------------*/
     count[1] = arraySize / mpi_size;
     count[0] = num_slices;
     offset[1] = count[1]*mpi_rank;
     offset[0] = 0;

     remainder = arraySize % mpi_size;
     if ( remainder>0 ) {
        if ( mpi_rank<remainder ) { count[1]++; }
        else { offset[1] += remainder-1; }
     }

  /*
   * Write the coordinate data to hard disk
   *-------------------------------------------------------------------------*/
     status = Fieldml_WriteDoubleSlab( writer, offset, count, 
                                 &(highResCoordinates[offset[0]][offset[1]]) );    
     if ( status==FML_IOERR_NO_ERROR ) { testOk = true; }
     else { testOk = false; }

     Fieldml_CloseWriter( writer ); 
     return testOk;

}


/**============================================================================
 ** highResHdf5Write
 ** 
 ** C subroutine that creates two HDF5 files that have data written to them in
 ** parallel.
 **============================================================================*/

bool highResHdf5Write() {

     bool testOk;
     FmlSessionHandle session;
     FmlObjectHandle cType, resource, sourceD, potentialResource, 
                     potentialSourceD;

     session = Fieldml_Create( "test", "test" );
     cType = Fieldml_CreateContinuousType( session, "test.scalar_real" );

 /*
  * Define a HDF5 file for coordinate data and write data to it. 
  *----------------------------------------------------------------------------*/
     resource = Fieldml_CreateHrefDataResource( session, "coordinates.resource",
                                               "PHDF5", "coordinates.h5" );
     sourceD = Fieldml_CreateArrayDataSource( session, "coordinates.source_double", 
                                              resource, "coordinates", 2 );
     testOk = highres_write( session, cType, sourceD, sizeY );

 /*
  * Define a HDF5 file for potential field data and write data to it. 
  *----------------------------------------------------------------------------*/
     potentialResource = Fieldml_CreateHrefDataResource( session, 
                                "potential.resource", "PHDF5", "potential.h5" );
     potentialSourceD = Fieldml_CreateArrayDataSource( session, 
                                "potential.source_double", potentialResource, 
                                "potential", 2 );
     testOk = testOk && highres_write( session, cType, potentialSourceD, 
                                       numberOfTimes );

     Fieldml_Destroy( session );
//     delete[] arraySizeForEachNode;
     return testOk;
}


/*===========================================================================
  SETUP_HIGHRES_COORDINATES 

  C++ subroutine that creates the locations of each coordinate point in the
  global domain. 
  ===========================================================================*/

void setup_highres_coordinates() {

     int j, x, y, z;
     double zValue, yValue, tmp = (double)size / ((double)size - 1.0);

     highResCoordinates = new double*[sizeY];
     for ( x=0; x<sizeY; x++ )
          highResCoordinates[x] = new double[arraySize]; 

     j = 0;
     for ( z=0; z<size ; z++ ) {
         zValue = tmp * (double)z;
         for ( y=0; y<size; y++ ) {
             yValue = tmp * (double)y;
             for ( x=0; x<size; x++ ) {
                 highResCoordinates[0][j] = tmp * (double)x;
                 highResCoordinates[1][j] = yValue;
                 highResCoordinates[2][j] = zValue;
                 j++;
             }
         }
     }
     return;
}

/*===========================================================================
  SETUP_LOWRES_COORDINATES 

  C++ subroutine that creates the locations of each coordinate point in the
  lower resolution global domain. 
  ===========================================================================*/
 
void setup_lowres_coordinates() {

     int j, x, y, z;
     double yValue, zValue, increment;

     lowResCoordinates = new double*[sizeY];
     for ( x=0; x<sizeY; x++ ) 
         lowResCoordinates[x] = new double[lowResArraySize];

     j = 0;
     increment = (double)(size - 1) / (double)(lowResSize - 1);
     for ( z=0; z<lowResSize; z++ ) { 
         zValue = increment * z;
         for ( y=0; y<lowResSize; y++ ) {
             yValue = increment * y;
             for ( x=0; x<lowResSize; x++ ) {
                 lowResCoordinates[0][j] = increment * x;
                 lowResCoordinates[1][j] = yValue;
                 lowResCoordinates[2][j] = zValue;
                 j++;
             }
         }
     }
     return;
}


/*
 *  This demo first generates coordinates and its associated potential
 * at different time steps then write them out in the hdf5 format using
 * FieldML phdf5 implementation.
 *
 * These hdf5 files is then read back into this demo using phdf5
 * and the data will be analyse and written out with a much lower resolution.
 *
 */
int main( int argc, char *argv[] ) {

    int i;
    bool flag;
    double inittime, totaltime, maxtime;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

    if ( mpi_rank==0 ) { printf( "\n   FIELDML PARALLEL IO TEST\n------------------------------\n\n" ); } 

    inittime = MPI_Wtime();

   /*
    * Setup the coordinate lists
    *------------------------------------------------------------------------*/
    setup_highres_coordinates();
    setup_lowres_coordinates();

    if ( mpi_size>1 ) { lowResNumberOfNodes = mpi_size-1; }

   /*
    * Perform a HDF5 write test 
    *------------------------------------------------------------------------*/
    flag = highResHdf5Write();
    if ( mpi_rank==0 ) {
       if ( flag ) 
          printf( "  High-resolution write... PASSED\n" );
       else 
          printf( "  High-resolution write... FAILED\n" );
    } 

   /*
    * The program will then read the data and analyse it 
    *------------------------------------------------------------------------*/
/*    flag = writeLowResFieldMLSolution();
    if ( mpi_rank==0 ) {
       if ( flag ) 
          printf( "  Low-resolution write...  PASSED\n" );
       else 
          printf( "  Low-resolution write...  FAILED\n" );
    }
*/
    totaltime = MPI_Wtime() - inittime;
    MPI_Reduce( &totaltime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( mpi_rank==0 ) { printf("\n  Time I/O was %5.2f sec\n\n", maxtime); };

   /*
    * Clean-up 
    *------------------------------------------------------------------------*/
    for ( i=0; i<sizeY; i++ ) {
        delete[] highResCoordinates[i];
        delete[] lowResCoordinates[i];
    }

    delete[] highResCoordinates;
    delete[] lowResCoordinates;
    MPI_Finalize();
    return 0;
}
