#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>

typedef struct trajectory
{
	int atomID, atomType, molType, ix, iy, iz;
	float x, y, z;
	int isEndGroup;
} TRAJECTORY;

typedef struct vector
{
	float x1, y1, z1;
	float x2, y2, z2;
	float xc, yc, zc;
} VECTOR;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct rdf
{
	float rlo, rhi, gofr;
} RDF;

typedef struct stats
{
	float average, standardDeviation;
} STATS;

typedef struct orderParameterBins
{
	float orderParameter, rlo, rhi, count;
} ORDERPARAMETER_BINS;

typedef struct distanceBins
{
	float rlo, rhi, count;
} DIST_BINS;

int countNAtoms (FILE *file_dump, int *nAtomEntries)
{
	int nAtoms, currentAtomID, nAtomsFixed;
	char lineString[2000];
	rewind (file_dump);

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, file_dump); }

	sscanf (lineString, "%d\n", &nAtoms);
	(*nAtomEntries) = nAtoms;
	rewind (file_dump);
	nAtomsFixed = nAtoms;

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > nAtoms) {
			nAtomsFixed = currentAtomID; }
	}

	return nAtomsFixed;
}

SIMULATION_BOUNDARY readDumpBoundary (FILE *file_dump, SIMULATION_BOUNDARY boundary)
{
	rewind (file_dump);
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}

	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.xlo, &boundary.xhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.ylo, &boundary.yhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.zlo, &boundary.zhi);
	rewind (file_dump);

	boundary.xLength = boundary.xhi - boundary.xlo;
	boundary.yLength = boundary.yhi - boundary.ylo;
	boundary.zLength = boundary.zhi - boundary.zlo;

	return boundary;
}

TRAJECTORY *initializeAtoms (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].molType = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].isEndGroup = 0;
	}

	return atoms;
}

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtomEntries, SIMULATION_BOUNDARY *boundary)
{
	char lineString[2000];
	int currentAtomID = 1;

	for (int i = 0; i < 5; ++i) {
		fgets (lineString, 2000, file_dump); }

	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).xlo, &(*boundary).xhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).ylo, &(*boundary).yhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).zlo, &(*boundary).zhi);
	fgets (lineString, 2000, file_dump);

	for (int i = 0; i < nAtomEntries; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);
		sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
		atoms[currentAtomID - 1].isEndGroup = 0;
	}

	return atoms;
}

DIST_BINS *initializeDistBins (DIST_BINS *coordDist, int dist_nBins, float dist_binWidth)
{
	for (int i = 0; i < dist_nBins; ++i)
	{
		if (i == 0) {
			coordDist[0].rlo = 0;
			coordDist[0].rhi = dist_binWidth; }
		else {
			coordDist[i].rlo = coordDist[i - 1].rhi;
			coordDist[i].rhi = coordDist[i].rlo + dist_binWidth; }

		coordDist[i].count = 0;
	}

	return coordDist;
}

float translatePeriodicDistance (float x1, float x2, float simBoxLength, float newR)
{
	if (fabs (x1 - x2) > (simBoxLength / (float)2))
	{
		if (x1 >= x2) {
			newR = x1 - simBoxLength; }
		else if (x1 < x2) {
			newR = x1 + simBoxLength; }

		return newR;
	}
	else
	{
		return x1;
	}
}

bool checkIfWithinBin (bool withinBin, TRAJECTORY *atoms, int i, int j, int nAtoms, SIMULATION_BOUNDARY boundary, DIST_BINS *coordDist, int k, float dist_cutoff)
{
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);
	float distance;
	float newX, newY, newZ;

	newX = translatePeriodicDistance (atoms[i].x, atoms[j].x, xLength, newX);
	newY = translatePeriodicDistance (atoms[i].y, atoms[j].y, yLength, newY);
	newZ = translatePeriodicDistance (atoms[i].z, atoms[j].z, zLength, newZ);

	distance = sqrt (
		(newX - atoms[j].x) * (newX - atoms[j].x) +
		(newY - atoms[j].y) * (newY - atoms[j].y) +
		(newZ - atoms[j].z) * (newZ - atoms[j].z)
		);

	if (distance < coordDist[k].rhi && distance >= coordDist[k].rlo)
	{
		withinBin = true;
		return withinBin;
	}
	else
	{
		withinBin = false;
	}

	return withinBin;
}

DIST_BINS *computeCoordination (DIST_BINS *coordDist, int *nCoordination, int dist_nBins, TRAJECTORY *atoms, int atomType1, int atomType2, float dist_cutoff, int nAtoms, SIMULATION_BOUNDARY boundary)
{
	bool withinBin;
	(*nCoordination) = 0;

	for (int k = 0; k < dist_nBins; ++k)
	{
		for (int i = 0; i < nAtoms; ++i)
		{
			if ((atoms[i].atomType == atomType1) || (atomType1 == -1))
			{
				for (int j = i + 1; j < nAtoms; ++j)
				{
					if ((atoms[j].atomType == atomType2) || (atomType2 == -1))
					{
						withinBin = checkIfWithinBin (withinBin, atoms, i, j, nAtoms, boundary, coordDist, k, dist_cutoff);

						if (withinBin) {
							(*nCoordination)++;
							coordDist[k].count++; }
					}
				}
			}
		}		
	}


	return coordDist;
}

void printCoordinationNumber (FILE *file_stats, int nCoordination)
{
	fprintf(file_stats, "%d\n", nCoordination);
}

void printCoordinationDistribution (FILE *file_dist, DIST_BINS *coordDist, int dist_nBins)
{
	for (int i = 0; i < dist_nBins; ++i)
	{
		fprintf(file_dist, "%f\t", coordDist[i].count);
	}

	fprintf(file_dist, "\n");
}

void printCoordinationDistributionHeader (FILE *file_dist, DIST_BINS *coordDist, int dist_nBins)
{
	fprintf(file_dist, "# rlo ==> ");
	
	for (int i = 0; i < dist_nBins; ++i)
	{
		fprintf(file_dist, "%f\t", coordDist[i].rlo);
	}

	fprintf(file_dist, "\n");
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		fprintf(stdout, "REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = program\n {~} argv[1] = input dump file name\n {~} argv[2] = atom type 1\n {~} argv[3] = atom type 2\n {~} argv[4] = cutoff distance for coordination\n {~} argv[5] = distance bin width\n {~} argv[6] = output stats file\n\n");
		fflush (stdout);
		exit (1);
	}

	FILE *file_dump, *file_dist, *file_stats;
	char *pipeString;
	pipeString = (char *) malloc (200 * sizeof (char));

	if (strstr (argv[1], ".xz")) {
		snprintf (pipeString, 200, "xzcat %s", argv[1]);
		file_dump = popen (pipeString, "r"); }
	else {
		file_dump = fopen (argv[1], "r"); }

	file_dist = fopen ("coordinationDistance.stats", "w");
	file_stats = fopen (argv[6], "w");

	int nAtomEntries, nAtoms = countNAtoms (file_dump, &nAtomEntries), atomType1 = atoi (argv[2]), atomType2 = atoi (argv[3]), file_status;
	SIMULATION_BOUNDARY boundary;
	boundary = readDumpBoundary (file_dump, boundary);

	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	printf("Number of atoms in the trajectory file: %d\n", nAtoms);

	atoms = initializeAtoms (atoms, nAtoms);

	rewind (file_dump);
	file_status = fgetc (file_dump);

	DIST_BINS *coordDist;
	float dist_cutoff = atof (argv[4]);
	float dist_binWidth = atof (argv[5]);
	int dist_nBins = ceil (dist_cutoff / dist_binWidth);
	coordDist = (DIST_BINS *) malloc (dist_nBins * sizeof (DIST_BINS));
	int nCoordination;
	int currentTimestep = 0;


	while (file_status != EOF)
	{
		nCoordination = 0;
		coordDist = initializeDistBins (coordDist, dist_nBins, dist_binWidth);
		fprintf(stdout, "Scanning timestep: %d                           \r", currentTimestep);
		fflush (stdout);

		if (currentTimestep == 0) {
			printCoordinationDistributionHeader (file_dist, coordDist, dist_nBins); }

		atoms = readTimestep (file_dump, atoms, nAtomEntries, &boundary);
		coordDist = computeCoordination (coordDist, &nCoordination, dist_nBins, atoms, atomType1, atomType2, dist_cutoff, nAtoms, boundary);

		printCoordinationNumber (file_stats, nCoordination);
		printCoordinationDistribution (file_dist, coordDist, dist_nBins);

		currentTimestep++;
		file_status = fgetc (file_dump);
	}

	fclose (file_dump);
	fclose (file_dist);
	fclose (file_stats);
	return 0;
}